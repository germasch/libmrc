
#include <mrc_ts.h>
#include <mrc_fld.h>
#include <mrc_domain.h>
#include <mrc_params.h>
#include <mrc_bits.h>

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

// ======================================================================

#define BND 2

enum {
  U,
  NR_FLDS,
};

struct kdv {
  struct mrc_obj obj;

  struct mrc_domain *domain;
};

MRC_CLASS_DECLARE(kdv, struct kdv);

// ======================================================================

// FIXME BND
#define CRDX(ix) (MRC_CRDX(crds, (ix)+BND))

static void
_kdv_create(struct kdv *kdv)
{
  kdv->domain = mrc_domain_create(kdv_comm(kdv));
  kdv_add_child(kdv, (struct mrc_obj *) kdv->domain);
  // set defaults
  mrc_domain_set_param_int3(kdv->domain, "m", (int [3]) { 160, 1, 1 });
  struct mrc_crds *crds = mrc_domain_get_crds(kdv->domain);
  mrc_crds_set_param_int(crds, "sw", BND);
  mrc_crds_set_param_float3(crds, "l", (float[3]) { -8., 0., 0. });
  mrc_crds_set_param_float3(crds, "h", (float[3]) {  8., 0., 0. });
}

static struct mrc_fld *
kdv_get_fld(struct kdv *kdv, int nr_comps, const char *name)
{
  struct mrc_fld *x = mrc_domain_f1_create(kdv->domain);
  mrc_fld_set_name(x, name);
  mrc_fld_set_sw(x, BND);
  mrc_fld_set_nr_comps(x, nr_comps);
  mrc_fld_setup(x);
  return x;
}

static void
kdv_fill_ghosts(struct kdv *kdv, struct mrc_fld *x, int m_x)
{
  int mx = mrc_fld_dims(x)[0];
  MRC_F1(x, m_x , -2  ) = MRC_F1(x, m_x , mx-2);
  MRC_F1(x, m_x , -1  ) = MRC_F1(x, m_x , mx-1);
  MRC_F1(x, m_x , mx  ) = MRC_F1(x, m_x , 0);
  MRC_F1(x, m_x , mx+1) = MRC_F1(x, m_x , 1);
}

#define Dx(x, m_x, ix)							\
  ((MRC_F1(x, m_x, ix+1) - MRC_F1(x, m_x, ix-1)) / (CRDX(ix+1) - CRDX(ix-1)))

// assumes uniform coordinates!
#define Dxxx(x, m_x, ix)						\
  ((MRC_F1(x, m_x, ix+2) - 2.*MRC_F1(x, m_x, ix+1) + 2.*MRC_F1(x, m_x, ix-1) - MRC_F1(x, m_x, ix-2)) / (2.*powf(CRDX(ix+1) - CRDX(ix), 3.)))

static void
kdv_rhsf(void *ctx, struct mrc_fld *rhs, float time, struct mrc_fld *x)
{
  struct kdv *kdv = ctx;
  struct mrc_crds *crds = mrc_domain_get_crds(kdv->domain);

  kdv_fill_ghosts(kdv, x, U);

  mrc_f1_foreach(x, ix, 0, 0) {
    MRC_F1(rhs, U, ix) =  - Dxxx(x, U, ix) + 6. * MRC_F1(x, U, ix) * Dx(x, U, ix);
  } mrc_f1_foreach_end;

  mrc_obj_print_class_info(CLASS_INFO_VERB_DIFF);
}

// ======================================================================

static struct mrc_obj_method kdv_methods[] = {
  MRC_OBJ_METHOD("rhsf", kdv_rhsf),
  {},
};

struct mrc_class_kdv mrc_class_kdv = {
  .name             = "kdv",
  .size             = sizeof(struct kdv),
  .create           = _kdv_create,
  .methods          = kdv_methods,
};

// ======================================================================

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  mrc_class_mrc_fld.watch = true;
  libmrc_params_init(argc, argv);

  struct kdv *kdv = kdv_create(MPI_COMM_WORLD);
  kdv_set_from_options(kdv);
  kdv_setup(kdv);
  kdv_view(kdv);

  // i.c.
  struct mrc_crds *crds = mrc_domain_get_crds(kdv->domain);
  struct mrc_fld *x = kdv_get_fld(kdv, NR_FLDS, "x");
  mrc_fld_set_comp_name(x, U, "u");

  // setup initial equilibrium and perturbation
  mrc_f1_foreach(x, ix, 0, 0) {
    //    MRC_F1(x, U, ix) = sin(2.*M_PI * CRDX(ix));
    MRC_F1(x, U, ix) = -12. * 1./sqr(cosh(CRDX(ix))); // 3 solitons
  } mrc_f1_foreach_end;

  // run time integration
  struct mrc_ts *ts = mrc_ts_create_std(MPI_COMM_WORLD, NULL, NULL);
  mrc_ts_set_context(ts, kdv_to_mrc_obj(kdv));
  mrc_ts_set_solution(ts, mrc_fld_to_mrc_obj(x));
  mrc_ts_set_from_options(ts);
  mrc_ts_setup(ts);
  mrc_ts_solve(ts);
  mrc_ts_view(ts);
  mrc_ts_destroy(ts);

  mrc_fld_destroy(x);

  kdv_destroy(kdv);

  libmrc_finalize(true, CLASS_INFO_VERB_ACTIVE);
  MPI_Finalize();
  return 0;
}
