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

  // state
  double kx;
  double ky;
  double kz;
};

#define ggcm_mhd_ic_alfven(ic) mrc_to_subobj(ic, struct ggcm_mhd_ic_alfven)

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
// aw_omega
// Linear MHD waves (no dispersion)
// ksq: k dot k
// k_par: k dot B / |B|???

static double aw_omega(double ksq, double kz, double B0, double rr0, double pp0)
{
  double Bsq;
  Bsq = sqr(B0); // sqr(B0[0]) + sqr(B0[1]) + sqr(B0[2]);
  double vA = sqrt(Bsq / rr0);
  double vs = sqrt(5. * pp0 / (3. * rr0));
  // double q = (sqr(vA * kz));
  double a = ((sqr(vA)) + (sqr(vs)));
  double q =
    (ksq * a) - (sqrt(sqr(ksq * a) - sqr(2. * sqrt(ksq) * kz * vs * vA)));
  return sqrt(q / 2);
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_alfven_primitive

static double ggcm_mhd_ic_alfven_primitive(struct ggcm_mhd_ic* ic, int m,
                                           double crd[3])
{
  struct ggcm_mhd_ic_alfven* sub = ggcm_mhd_ic_alfven(ic);

  double amp = sub->amp;
  double rr0 = sub->rr0, pp0 = sub->pp0;
  double B0 = sub->B0, kx = sub->kx, ky = sub->ky, kz = sub->kz;
  double theta = 2. * M_PI / 360. * sub->theta;
  double xx = crd[0], yy = crd[1], zz = crd[2];

  double vA = sqrt(sqr(B0) / rr0);
  double rr1 = 0.;
  double ksq = sqr(kx) + sqr(ky) + sqr(kz);
  double k = sqrt(ksq);
  double om = k * vA * cos(theta);

  double vy1 = amp * sin((kx * xx) + (ky * yy) + (kz * zz));
  double by1 = k * B0 * cos(theta) / om * vy1;
  switch (m) {
    case RR: return rr0 + (rr1)*sin((kx * xx) + (ky * yy) + (kz * zz));
    case PP: return pp0;

    case VX: return 0.;
    case VY: return vy1;
    case VZ: return 0.;

    case BX: return B0 * sin(theta);
    case BY: return by1;
    case BZ: return B0 * cos(theta);

    default: return 0.;
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_alfven_descr

#define VAR(x) (void*)offsetof(struct ggcm_mhd_ic_alfven, x)
static struct param ggcm_mhd_ic_alfven_descr[] = {
  {"rr0", VAR(rr0), PARAM_DOUBLE(1.)},
  {"pp0", VAR(pp0), PARAM_DOUBLE(1.)},
  {"B0", VAR(B0), PARAM_DOUBLE(1.)},
  {"amp", VAR(amp), PARAM_DOUBLE(0.01)},
  {"theta", VAR(theta), PARAM_DOUBLE(0.)},
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
