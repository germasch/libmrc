// revised
#include <ggcm_mhd_private.h>
#include <ggcm_mhd_step.h>
#include <ggcm_mhd_ic_private.h>
#include <ggcm_mhd_crds_private.h>
#include <ggcm_mhd_bnd.h>
#include <ggcm_mhd_diag.h>

#include <mrc_fld_as_double.h>
#include <mrc_domain.h>

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

// ======================================================================
// ggcm_mhd_ic subclass "alfven"

struct ggcm_mhd_ic_alfven
{
  // params
  double rr0;
  double pp0;
  double B0;    // background field
  double theta; // background field angle
  double amp;   // perturbation amplitude
  int wave;     // wave type

  // state
  double kx;
  double ky;
  double kz;
};

#define ggcm_mhd_ic_alfven(ic) mrc_to_subobj(ic, struct ggcm_mhd_ic_alfven)

enum
{
  OPT_WAVE_ALFVEN,
  OPT_WAVE_FAST,
  OPT_WAVE_SLOW
};

// ----------------------------------------------------------------------
// ggcm_mhd_ic_alfven_setup

static void ggcm_mhd_ic_alfven_setup(struct ggcm_mhd_ic* ic)
{
  struct ggcm_mhd_ic_alfven* sub = ggcm_mhd_ic_alfven(ic);
  struct mrc_crds* crds = mrc_domain_get_crds(ic->mhd->domain);
  const double *lo = mrc_crds_lo(crds), *hi = mrc_crds_hi(crds);

  ggcm_mhd_ic_setup_super(ic);

  sub->kx = 0.; // 2. * M_PI / (hi[0] - lo[0]);
  sub->ky = 0.;
  sub->kz = 2. * M_PI / (hi[2] - lo[2]);
}

// ----------------------------------------------------------------------
// pert_B

static void pert_B(struct ggcm_mhd_ic* ic, double v[3], double om, double b[3])
{
  struct ggcm_mhd_ic_alfven* sub = ggcm_mhd_ic_alfven(ic);

  double k = sub->kz;
  double B0 = sub->B0;
  double theta = 2. * M_PI / 360. * sub->theta;

  b[0] = (k * B0 * cos(theta) * v[0] + k * v[2] * B0 * sin(theta)) / om;
  b[1] = (k * B0 * cos(theta) * v[1]) / om;
  b[2] = (k * B0 * cos(theta) * v[2] + k * v[2] * B0 * cos(theta)) / om;
}

// ----------------------------------------------------------------------
// pert_rr

static void pert_rr(struct ggcm_mhd_ic* ic, double v[3], double om, double* rr)
{
  struct ggcm_mhd_ic_alfven* sub = ggcm_mhd_ic_alfven(ic);

  double k = sub->kz;
  double rr0 = sub->rr0;

  *rr = rr0 * k * v[2] / om;
}

// ----------------------------------------------------------------------
// pert_pp

static void pert_pp(struct ggcm_mhd_ic* ic, double v[3], double om, double* pp)
{
  struct ggcm_mhd* mhd = ic->mhd;
  struct ggcm_mhd_ic_alfven* sub = ggcm_mhd_ic_alfven(ic);

  double k = sub->kz;
  double pp0 = sub->pp0;

  *pp = mhd->par.gamm * pp0 * k * v[2] / om;
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_alfven_primitive

static double ggcm_mhd_ic_alfven_primitive(struct ggcm_mhd_ic* ic, int m,
                                           double crd[3])
{
  struct ggcm_mhd* mhd = ic->mhd;
  struct ggcm_mhd_ic_alfven* sub = ggcm_mhd_ic_alfven(ic);

  double amp = sub->amp;
  double k = sub->kz;
  double rr0 = sub->rr0, pp0 = sub->pp0, B0 = sub->B0;
  double theta = 2. * M_PI / 360. * sub->theta;
  double zz = crd[2];

  double vA = sqrt(sqr(B0) / rr0);
  double vS = sqrt(mhd->par.gamm * pp0 / rr0);
  double vA2 = sqr(vA), vS2 = sqr(vS);
  double vp = sqrt(
    .5 * (vA2 + vS2 + sqrt(sqr(vA2 + vS2) - 4. * vA2 * vS2 * sqr(cos(theta)))));

  double om;
  double v[3];
  switch (sub->wave) {
    case OPT_WAVE_ALFVEN:
      om = k * vA * cos(theta);
      v[0] = 0.;
      v[1] = amp;
      v[2] = 0.;
      break;
    case OPT_WAVE_FAST:
      om = k * vp;
      double vx = amp * (sqr(vp) - sqr(cos(theta) * sqr(vS)));
      double vz = amp * (cos(theta) * sin(theta) * sqr(vS));
      v[0] = cos(theta) * vx + sin(theta) * vz;
      v[1] = 0.;
      v[2] = -sin(theta) * vx + cos(theta) * vz;
      break;
    default: assert(0);
  }

  double b[3], rr1, pp1;
  pert_B(ic, v, om, b);
  pert_rr(ic, v, om, &rr1);
  pert_pp(ic, v, om, &pp1);

  switch (m) {
    case RR: return rr1 * sin(k * zz) + rr0;
    case PP: return pp1 * sin(k * zz) + pp0;

    case VX: return v[0] * sin(k * zz);
    case VY: return v[1] * sin(k * zz);
    case VZ: return v[2] * sin(k * zz);

    case BX: return b[0] * sin(k * zz) + B0 * sin(theta);
    case BY: return b[1] * sin(k * zz);
    case BZ: return b[2] * sin(k * zz) + B0 * cos(theta);

    default: return 0.;
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_alfven_descr

static struct mrc_param_select opt_wave_descr[] _mrc_unused = {
  {.val = OPT_WAVE_ALFVEN, .str = "alfven"},
  {.val = OPT_WAVE_FAST, .str = "fast"},
  {.val = OPT_WAVE_SLOW, .str = "slow"},
  {},
};

#define VAR(x) (void*)offsetof(struct ggcm_mhd_ic_alfven, x)
static struct param ggcm_mhd_ic_alfven_descr[] = {
  {"rr0", VAR(rr0), PARAM_DOUBLE(1.)},
  {"pp0", VAR(pp0), PARAM_DOUBLE(1.)},
  {"B0", VAR(B0), PARAM_DOUBLE(1.)},
  {"amp", VAR(amp), PARAM_DOUBLE(0.001)},
  {"theta", VAR(theta), PARAM_DOUBLE(0.)},
  {"wave", VAR(wave), PARAM_SELECT(OPT_WAVE_ALFVEN, opt_wave_descr)},
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_ic_alfven_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_alfven_ops = {
  .name = "alfven",
  .size = sizeof(struct ggcm_mhd_ic_alfven),
  .param_descr = ggcm_mhd_ic_alfven_descr,
  .setup = ggcm_mhd_ic_alfven_setup,
  .primitive = ggcm_mhd_ic_alfven_primitive,
};

// ======================================================================
// ggcm_mhd subclass "ot"

// ----------------------------------------------------------------------
// ggcm_mhd_ot_create

static void ggcm_mhd_ot_create(struct ggcm_mhd* mhd)
{
  ggcm_mhd_default_box(mhd);

  // default mesh size
  mrc_domain_set_param_int3(mhd->domain, "m", (int[3]){1, 1, 64});

  // default domain size
  struct mrc_crds* crds = mrc_domain_get_crds(mhd->domain);
  mrc_crds_set_type(crds, "uniform");
  mrc_crds_set_param_int(crds, "sw", SW_2); // 'stencil width'
  mrc_crds_set_param_double3(crds, "l", (double[3]){0.0, 0.0, 0.0});
  mrc_crds_set_param_double3(crds, "h", (double[3]){0.01, 0.01, 1.0});
}

// ----------------------------------------------------------------------
// ggcm_mhd_ot_ops

static struct ggcm_mhd_ops ggcm_mhd_ot_ops = {
  .name = "ot",
  .create = ggcm_mhd_ot_create,
};

//================================================================
// main

extern struct ggcm_mhd_diag_ops ggcm_mhd_diag_c_ops;

int main(int argc, char** argv)
{
  mrc_class_register_subclass(&mrc_class_ggcm_mhd, &ggcm_mhd_ot_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag, &ggcm_mhd_diag_c_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_alfven_ops);

  return ggcm_mhd_main(&argc, &argv);
}
