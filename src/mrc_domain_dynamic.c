
#include "mrc_domain_private.h"
#include "mrc_params.h"
#include "mrc_ddc.h"
#include "mrc_io.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#define TAG_SCAN_OFF (1000)

static inline struct mrc_domain_multi *
mrc_domain_multi(struct mrc_domain *domain)
{
  return domain->obj.subctx;
}

// ======================================================================
// map
// 
// maps between global patch index (contiguous) and 1D SFC idx
// (potentially non-contiguous)
// if sfc_indices is NULL, the map will be the indentity

static void
map_create(struct mrc_domain *domain, int *sfc_indices, int nr_gpatches)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  if (!sfc_indices) {
    return;
  }
  multi->gp = malloc(sizeof(int) * multi->nr_global_patches);
  for (int i = 0; i < nr_gpatches; i++) {
    multi->gp[i] = sfc_indices[i];
  }

  //Create the bintree for performant searching
  int vals[nr_gpatches];
  for (int i = 0; i < nr_gpatches; i++) {
    vals[i] = i;
  }
  bintree_create_from_ordered_list(&multi->g_patches, sfc_indices, vals, nr_gpatches);
}

static void
map_destroy(struct mrc_domain *domain)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  if (!multi->gp) {
    return;
  }
  free(multi->gp);
  bintree_destroy(&multi->g_patches);
}

static int
map_sfc_idx_to_gpatch(struct mrc_domain *domain, int sfc_idx)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  if (!multi->gp) {
    return sfc_idx;
  }
  int retval;
  int rc = bintree_get(&multi->g_patches, sfc_idx, &retval);
  if (rc == 0) {
    return -1;
  }
  return retval;
}

static int
map_gpatch_to_sfc_idx(struct mrc_domain *domain, int gpatch)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  if (!multi->gp) {
    return gpatch;
  }
  return multi->gp[gpatch];
}

// ======================================================================

static void
sfc_idx_to_rank_patch(struct mrc_domain *domain, int sfc_idx,
		      int *rank, int *patch)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  int gpatch = map_sfc_idx_to_gpatch(domain, sfc_idx);
  if (gpatch < 0) {
    *rank = -1;
    *patch = -1;
    return;
  }
  
  // FIXME, this can be done much more efficiently using binary search...
  for (int i = 0; i < domain->size; i++) {
    if (gpatch < multi->gpatch_off_all[i+1]) {
      *rank = i;
      *patch = gpatch - multi->gpatch_off_all[i];
      break;
    }
  }
}

// ======================================================================

static void
mrc_domain_multi_view(struct mrc_domain *domain)
{
}

static void
mrc_domain_multi_get_global_patch_info(struct mrc_domain *domain, int gpatch,
				       struct mrc_patch_info *info)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  assert(gpatch < multi->nr_global_patches);
  info->global_patch = gpatch;
  int sfc_idx = map_gpatch_to_sfc_idx(domain, gpatch);
  sfc_idx_to_rank_patch(domain, sfc_idx, &info->rank, &info->patch);
  
  assert(info->rank >= 0);
  
  int p3[3];
  sfc_idx_to_idx3(&multi->sfc, sfc_idx, p3);
  for (int d = 0; d < 3; d++) {
    info->ldims[d] = multi->ldims[d][p3[d]];
    info->off[d] = multi->off[d][p3[d]];
    info->idx3[d] = p3[d];
  }
}

static void
mrc_domain_multi_get_local_patch_info(struct mrc_domain *domain, int patch,
				      struct mrc_patch_info *info)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  mrc_domain_multi_get_global_patch_info(domain, multi->gpatch_off + patch,
					 info);
}

static void
setup_gpatch_off_all(struct mrc_domain *domain)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  multi->gpatch_off_all = calloc(domain->size + 1, sizeof(*multi->gpatch_off_all));
  int nr_global_patches = multi->nr_global_patches;

  if (multi->nr_patches >= 0) {
    // prescribed mapping patch <-> proc
    MPI_Comm comm = mrc_domain_comm(domain);
    int *nr_patches_all = calloc(domain->size, sizeof(*nr_patches_all));
    MPI_Gather(&multi->nr_patches, 1, MPI_INT, nr_patches_all, 1, MPI_INT,
	       0, comm);
    MPI_Bcast(nr_patches_all, domain->size, MPI_INT, 0, comm);

    for (int i = 1; i <= domain->size; i++) {
      int nr_patches = nr_patches_all[i-1];
      multi->gpatch_off_all[i] = multi->gpatch_off_all[i-1] + nr_patches;
    }
    free(nr_patches_all);
  } else {
    // map patch <-> proc uniformly (roughly)
    int patches_per_proc = nr_global_patches / domain->size;
    int patches_per_proc_rmndr = nr_global_patches % domain->size;
    for (int i = 1; i <= domain->size; i++) {
      int nr_patches = patches_per_proc + ((i-1) < patches_per_proc_rmndr);
      multi->gpatch_off_all[i] = multi->gpatch_off_all[i-1] + nr_patches;
    }
  }
}

static void
mrc_domain_multi_setup_patches(struct mrc_domain *domain)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  int sfc_indices[multi->nr_global_patches];
  
  setup_gpatch_off_all(domain);

  multi->gpatch_off = multi->gpatch_off_all[domain->rank];
  multi->nr_patches = multi->gpatch_off_all[domain->rank+1] - multi->gpatch_off;

  multi->patches = malloc(sizeof(*multi->patches) * multi->nr_patches);
 
  int activerank = 0;
  int npatches = 0;
  
  //TODO Find a smarter way than iterating over all possible patches
  int npt = multi->np[0] * multi->np[1] * multi->np[2];
  for(int i = 0; i < npt; i++) {
    int idx[3];
    sfc_idx_to_idx3(&multi->sfc, i, idx);
    if(bitfield3d_isset(&multi->activepatches, idx)) {
      //Calculate rank
      if (npatches >= multi->gpatch_off_all[activerank+1]) {
	activerank++;
      }
      
      //Register the patch
      sfc_indices[npatches] = i;

      if (activerank == domain->rank) { // Create the patch on multi processor
	//Setup patches[lpatch]
	int lpatch = npatches - multi->gpatch_off;
	for(int d = 0; d < 3; d++) {
	  multi->patches[lpatch].ldims[d] = multi->ldims[d][idx[d]];
	  multi->patches[lpatch].off[d] = multi->off[d][idx[d]];
	}
      }
      npatches++;
    }
  }
  
  map_create(domain, sfc_indices, multi->nr_global_patches);
}

static void
mrc_domain_multi_setup(struct mrc_domain *domain)
{
  assert(!domain->is_setup);
  domain->is_setup = true;

  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  MPI_Comm comm = mrc_domain_comm(domain);
  MPI_Comm_rank(comm, &domain->rank);
  MPI_Comm_size(comm, &domain->size);
  
  //Copy the activepatch-list
  bitfield3d_copy(&multi->activepatches, multi->p_activepatches);
  
  multi->nr_global_patches = bitfield3d_count_bits_set(&multi->activepatches);
  
  int *np = multi->np;
  for (int d = 0; d < 3; d++) {
    int ldims[3], rmndr[3];
    ldims[d] = multi->gdims[d] / np[d];
    rmndr[d] = multi->gdims[d] % np[d];

    multi->ldims[d] = calloc(np[d], sizeof(*multi->ldims[d]));
    multi->off[d] = calloc(np[d], sizeof(*multi->off[d]));
    for (int i = 0; i < np[d]; i++) {
      multi->ldims[d][i] = ldims[d] + (i < rmndr[d]);
      if (i > 0) {
	multi->off[d][i] = multi->off[d][i-1] + multi->ldims[d][i-1];
      }
    }
  }
  
  //Create list of patches
  sfc_setup(&multi->sfc, multi->np);
  mrc_domain_multi_setup_patches(domain);
}

static void
mrc_domain_multi_destroy(struct mrc_domain *domain)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  for (int d = 0; d < 3; d++) {
    free(multi->ldims[d]);
    free(multi->off[d]);
  }
  free(multi->gpatch_off_all);
  free(multi->patches);
  bitfield3d_destroy(&multi->activepatches);
  map_destroy(domain);
}

static struct mrc_patch *
mrc_domain_multi_get_patches(struct mrc_domain *domain, int *nr_patches)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);
  if (nr_patches) {
    *nr_patches = multi->nr_patches;
  }
  return multi->patches;
}

static void
mrc_domain_multi_get_global_dims(struct mrc_domain *domain, int *dims)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  for (int d = 0; d < 3; d++) {
    dims[d] = multi->gdims[d];
  }
}

static void
mrc_domain_multi_get_nr_procs(struct mrc_domain *domain, int *nr_procs)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  for (int d = 0; d < 3; d++) {
    nr_procs[d] = multi->np[d];
  }
}

static void
mrc_domain_multi_get_bc(struct mrc_domain *domain, int *bc)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  for (int d = 0; d < 3; d++) {
    bc[d] = multi->bc[d];
  }
}

static void
mrc_domain_multi_get_nr_global_patches(struct mrc_domain *domain, int *nr_global_patches)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  *nr_global_patches = multi->nr_global_patches;
}

static void
mrc_domain_multi_get_idx3_patch_info(struct mrc_domain *domain, int idx[3],
				     struct mrc_patch_info *info)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);
  //Check if the patch is active
  if (!bitfield3d_isset(&multi->activepatches, idx)) {
    info->rank = -1;
    info->patch = -1;
    info->global_patch = -1;
    for(int d = 0; d < 3; d++) {
      info->ldims[d] = multi->ldims[d][idx[d]];
      info->off[d] = multi->off[d][idx[d]];
      info->idx3[d] = idx[d];
    }
    return;
  }

  int sfc_idx = sfc_idx3_to_idx(&multi->sfc, idx);
  int gpatch = map_sfc_idx_to_gpatch(domain, sfc_idx);
  mrc_domain_multi_get_global_patch_info(domain, gpatch, info);
}

static void
mrc_domain_multi_write(struct mrc_domain *domain, struct mrc_io *io)
{
  int nr_global_patches;
  mrc_domain_multi_get_nr_global_patches(domain, &nr_global_patches);
  mrc_io_write_attr_int(io, mrc_domain_name(domain), "nr_global_patches",
			nr_global_patches);
  
  //Iterate over all global patches
  for (int gp = 0; gp < nr_global_patches; gp++) {
    struct mrc_patch_info info;
    mrc_domain_multi_get_global_patch_info(domain, gp, &info);
    char path[strlen(mrc_domain_name(domain)) + 10];
    sprintf(path, "%s/p%d", mrc_domain_name(domain), gp);
    mrc_io_write_attr_int3(io, path, "ldims", info.ldims);
    mrc_io_write_attr_int3(io, path, "off", info.off);
    mrc_io_write_attr_int3(io, path, "idx3", info.idx3);
    int sfc_idx = map_gpatch_to_sfc_idx(domain, gp);
    mrc_io_write_attr_int(io, path, "sfc_idx", sfc_idx);
  }
}

static void
mrc_domain_multi_plot(struct mrc_domain *domain)
{
}

static struct mrc_ddc *
mrc_domain_multi_create_ddc(struct mrc_domain *domain)
{
  struct mrc_ddc *ddc = mrc_ddc_create(domain->obj.comm);
  mrc_ddc_set_type(ddc, "multi");
  mrc_ddc_set_domain(ddc, domain);
  return ddc;
}

static struct mrc_param_select bc_descr[] = {
  { .val = BC_NONE       , .str = "none"     },
  { .val = BC_PERIODIC   , .str = "periodic" },
  {},
};

static struct mrc_param_select curve_descr[] = {
  { .val = CURVE_BYDIM   , .str = "bydim"    },
  { .val = CURVE_MORTON  , .str = "morton"   },
  { .val = CURVE_HILBERT , .str = "hilbert"  },
  {},
};

#define VAR(x) (void *)offsetof(struct mrc_domain_multi, x)
static struct param mrc_domain_multi_params_descr[] = {
  { "m"               , VAR(gdims)           , PARAM_INT3(32, 32, 32) },
  { "np"              , VAR(np)              , PARAM_INT3(1, 1, 1)    },
  { "bcx"             , VAR(bc[0])           , PARAM_SELECT(BC_NONE,
							    bc_descr) },
  { "bcy"             , VAR(bc[1])           , PARAM_SELECT(BC_NONE,
							    bc_descr) },
  { "bcz"             , VAR(bc[2])           , PARAM_SELECT(BC_NONE,
							    bc_descr) },
  { "curve_type"      , VAR(sfc.curve_type)  , PARAM_SELECT(CURVE_BYDIM,
							    curve_descr) },
  { "nr_patches"      , VAR(nr_patches)      , PARAM_INT(-1) },
  { "activepatches"   , VAR(p_activepatches) , PARAM_PTR(NULL) },
  {},
};
#undef VAR

struct mrc_domain_ops mrc_domain_dynamic_ops = {
  .name                  = "dynamic",
  .size                  = sizeof(struct mrc_domain_multi),
  .param_descr           = mrc_domain_multi_params_descr,
  .setup                 = mrc_domain_multi_setup,
  .view                  = mrc_domain_multi_view,
  .write                 = mrc_domain_multi_write,
  .destroy               = mrc_domain_multi_destroy,
  .get_patches           = mrc_domain_multi_get_patches,
  .get_global_dims       = mrc_domain_multi_get_global_dims,
  .get_nr_procs          = mrc_domain_multi_get_nr_procs,
  .get_bc                = mrc_domain_multi_get_bc,
  .get_nr_global_patches = mrc_domain_multi_get_nr_global_patches,
  .get_global_patch_info = mrc_domain_multi_get_global_patch_info,
  .get_local_patch_info  = mrc_domain_multi_get_local_patch_info,
  .get_idx3_patch_info   = mrc_domain_multi_get_idx3_patch_info,
  .plot                  = mrc_domain_multi_plot,
  .create_ddc            = mrc_domain_multi_create_ddc,
};
