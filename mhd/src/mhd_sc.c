
// ----------------------------------------------------------------------
// zmaskn

static void __unused
zmaskn(struct ggcm_mhd *mhd, struct mrc_fld *x)
{
  float va02i = 1.f / sqr(mhd->par.speedlimit / mhd->par.vvnorm);
  float eps   = 1e-15f;

  mrc_fld_foreach(x, ix,iy,iz, 1, 1) {
    mrc_fld_data_t bb = (sqr(.5f * (BX(x, ix,iy,iz) + BX(x, ix-1,iy,iz))) +
			 sqr(.5f * (BY(x, ix,iy,iz) + BY(x, ix,iy-1,iz))) +
			 sqr(.5f * (BZ(x, ix,iy,iz) + BZ(x, ix,iy,iz-1))));
    float rrm = fmaxf(eps, bb * va02i);
    F3(x, _ZMASK, ix,iy,iz) = F3(x, _YMASK, ix,iy,iz) *
      fminf(1.f, RR(x, ix,iy,iz) / rrm);
  } mrc_fld_foreach_end;
}

// ----------------------------------------------------------------------
// newstep_sc

static mrc_fld_data_t __unused
newstep_sc(struct ggcm_mhd *mhd, struct mrc_fld *x)
{
  float *fd1x = ggcm_mhd_crds_get_crd(mhd->crds, 0, FD1);
  float *fd1y = ggcm_mhd_crds_get_crd(mhd->crds, 1, FD1);
  float *fd1z = ggcm_mhd_crds_get_crd(mhd->crds, 2, FD1);

  mrc_fld_data_t splim2 = sqr(mhd->par.speedlimit / mhd->par.vvnorm);
  mrc_fld_data_t gamm   = mhd->par.gamm;
  mrc_fld_data_t d_i    = mhd->par.d_i;

  mrc_fld_data_t thx    = mhd->par.thx;
  mrc_fld_data_t eps    = 1e-9f;
  mrc_fld_data_t dt     = 1e10f;

  mrc_fld_data_t two_pi_d_i = 2. * M_PI * d_i;
  bool have_hall = d_i > 0.f;

  mrc_fld_foreach(x, ix, iy, iz, 0, 0) {
    mrc_fld_data_t hh = mrc_fld_max(mrc_fld_max(fd1x[ix], fd1y[iy]), fd1z[iz]);
    mrc_fld_data_t rri = 1.f / mrc_fld_abs(RR(x, ix,iy,iz)); // FIXME abs necessary?
    mrc_fld_data_t bb = (sqr(.5f * (BX(x, ix,iy,iz) + BX(x, ix-1,iy,iz))) + 
			 sqr(.5f * (BY(x, ix,iy,iz) + BY(x, ix,iy-1,iz))) +
			 sqr(.5f * (BZ(x, ix,iy,iz) + BZ(x, ix,iy,iz-1))));
    if (have_hall) {
      bb *= 1 + sqr(two_pi_d_i * hh);
    }      

    mrc_fld_data_t vv1 = mrc_fld_min(bb * rri, splim2);
    mrc_fld_data_t rrvv = (sqr(RVX(x, ix,iy,iz)) + 
			   sqr(RVY(x, ix,iy,iz)) +
			   sqr(RVZ(x, ix,iy,iz)));
    mrc_fld_data_t pp = (gamm - 1.f) * (UU(x, ix,iy,iz) - .5f * rrvv * rri);
    mrc_fld_data_t vv2 = gamm * pp * rri;
    mrc_fld_data_t vv = mrc_fld_sqrt(vv1 + vv2) + mrc_fld_sqrt(rrvv) * rri;
    vv = mrc_fld_max(eps, vv);

    mrc_fld_data_t tt = thx / mrc_fld_max(eps, hh*vv*F3(x, _ZMASK, ix,iy,iz));
    dt = mrc_fld_min(dt, tt);
  } mrc_fld_foreach_end;

  mrc_fld_data_t dtn;
  MPI_Allreduce(&dt, &dtn, 1, MPI_MRC_FLD_DATA_T, MPI_MIN, ggcm_mhd_comm(mhd));

  return dtn;
}

// ----------------------------------------------------------------------
// newstep_sc_ggcm
//
// like newstep_c, but doesn't need any of the prep work
// (no primvar, primbb, ymask, zmask)

static mrc_fld_data_t __unused
newstep_sc_ggcm(struct ggcm_mhd *mhd, struct mrc_fld *x)
{
  static int PR;
  if (!PR) {
    PR = prof_register("newstep_sc_ggcm", 1., 0, 0);
  }
  prof_start(PR);

  float *fd1x = ggcm_mhd_crds_get_crd(mhd->crds, 0, FD1);
  float *fd1y = ggcm_mhd_crds_get_crd(mhd->crds, 1, FD1);
  float *fd1z = ggcm_mhd_crds_get_crd(mhd->crds, 2, FD1);
  float *fx2x = ggcm_mhd_crds_get_crd(mhd->crds, 0, FX2);
  float *fx2y = ggcm_mhd_crds_get_crd(mhd->crds, 1, FX2);
  float *fx2z = ggcm_mhd_crds_get_crd(mhd->crds, 2, FX2);

  mrc_fld_data_t splim2 = sqr(mhd->par.speedlimit / mhd->par.vvnorm);
  mrc_fld_data_t isphere2 = sqr(mhd->par.isphere);
  mrc_fld_data_t gamm   = mhd->par.gamm;
  mrc_fld_data_t d_i    = mhd->par.d_i;
  mrc_fld_data_t thx    = mhd->par.thx;
  mrc_fld_data_t eps    = 1e-9f;
  mrc_fld_data_t dt     = 1e10f;
  mrc_fld_data_t va02i  = 1.f / sqr(mhd->par.speedlimit / mhd->par.vvnorm);
  mrc_fld_data_t epsz   = 1e-15f;
  mrc_fld_data_t s      = gamm - 1.f;

  mrc_fld_data_t two_pi_d_i = 2. * M_PI * d_i;
  bool have_hall = d_i > 0.f;
  mrc_fld_foreach(x, ix, iy, iz, 0, 0) {
    mrc_fld_data_t hh = fmaxf(fmaxf(fd1x[ix], fd1y[iy]), fd1z[iz]);
    mrc_fld_data_t rri = 1.f / fabsf(RR(x, ix,iy,iz)); // FIXME abs necessary?
    mrc_fld_data_t bb = 
      sqr(.5f * (BX(x, ix,iy,iz) + BX(x, ix-1,iy,iz))) + 
      sqr(.5f * (BY(x, ix,iy,iz) + BY(x, ix,iy-1,iz))) +
      sqr(.5f * (BZ(x, ix,iy,iz) + BZ(x, ix,iy,iz-1)));
    if (have_hall) {
      bb *= 1 + sqr(two_pi_d_i * hh);
    }
    mrc_fld_data_t vv1 = fminf(bb * rri, splim2);
    
    mrc_fld_data_t rv2 = 
      sqr(RVX(x, ix,iy,iz)) + sqr(RVY(x, ix,iy,iz)) + sqr(RVZ(x, ix,iy,iz));
    mrc_fld_data_t rvv = rri * rv2;
    mrc_fld_data_t pp = s * (UU(x, ix,iy,iz) - .5f * rvv);
    mrc_fld_data_t vv2 = gamm * fmaxf(0.f, pp) * rri;
    mrc_fld_data_t vv3 = rri * sqrtf(rv2);
    mrc_fld_data_t vv = sqrtf(vv1 + vv2) + vv3;
    vv = fmaxf(eps, vv);
    
    mrc_fld_data_t ymask = 1.f;
    if (fx2x[ix] + fx2y[iy] + fx2z[iz] < isphere2)
      ymask = 0.f;
    
    mrc_fld_data_t rrm = fmaxf(epsz, bb * va02i);
    mrc_fld_data_t zmask = ymask * fminf(1.f, RR(x, ix,iy,iz) / rrm);
    
    mrc_fld_data_t tt = thx / fmaxf(eps, hh*vv*zmask);
    dt = fminf(dt, tt);
  } mrc_fld_foreach_end;

  mrc_fld_data_t dtn;
  MPI_Allreduce(&dt, &dtn, 1, MPI_MRC_FLD_DATA_T, MPI_MIN, ggcm_mhd_comm(mhd));

  prof_stop(PR);

  return dtn;
}

