
#include "ggcm_mhd_diag_item_private.h"

#include "ggcm_mhd_diag.h"

#include <stdio.h>
#include <assert.h>

#define ggcm_mhd_diag_item_ops(item) ((struct ggcm_mhd_diag_item_ops *) (item)->obj.ops)

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item_run

void
ggcm_mhd_diag_item_run(struct ggcm_mhd_diag_item *item, struct mrc_io *io,
		       struct mrc_f3 *f, int diag_type, float plane)
{
  struct ggcm_mhd_diag_item_ops *ops = ggcm_mhd_diag_item_ops(item);
  assert(ops && ops->run);
  ops->run(item, io, f, diag_type, plane);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_init

static void
ggcm_mhd_diag_init()
{
#if 0
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_v);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_rr);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_pp);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_b);
#endif

  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_rr1);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_uu1);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_rv1);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_b1);

#if 0
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_e);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_j);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_xtra);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_zmask);
#endif

  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_divb);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item description

#define VAR(x) (void *)offsetof(struct ggcm_mhd_diag_item, x)
static struct param ggcm_mhd_diag_item_descr[] = {
  { "diag"            , VAR(diag)            , PARAM_OBJ(ggcm_mhd_diag) },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item class

struct mrc_class_ggcm_mhd_diag_item mrc_class_ggcm_mhd_diag_item = {
  .name             = "ggcm_mhd_diag_item",
  .size             = sizeof(struct ggcm_mhd_diag_item),
  .param_descr      = ggcm_mhd_diag_item_descr,
  .init             = ggcm_mhd_diag_init,
};
