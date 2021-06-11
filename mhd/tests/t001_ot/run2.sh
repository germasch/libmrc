#!/bin/sh 

mpirun -n 1 ../mhd_ot \
    --mrc_crds_lx 0. --mrc_crds_hx 1.0 \
    --mrc_crds_lz 0. --mrc_crds_hz 1.0 \
    --mrc_crds_ly -0.01 --mrc_crds_hy 0.01 \
    \
    --mrc_domain_mx 1 --mrc_domain_my 1 --mrc_domain_mz 128 \
    --mrc_domain_npx 1 --mrc_domain_npy 1 \
    \
    --mrc_ts_output_every_time .01  \
    --ggcm_mhd_diag_fields rr1:uu1:rv1:j:b1:divb:rr:pp:v:b \
    \
    --ggcm_mhd_ic_B0 2. \
    --ggcm_mhd_ic_theta 90. \
    --ggcm_mhd_ic_wave fast \
    --ggcm_mhd_ic_amp 1e-4 \
    \
    --mrc_ts_max_time 4.0 \
    --ggcm_mhd_magdiffu const \
    --ggcm_mhd_step_type mhdcc_double \
    --ggcm_mhd_step_limiter minmod \
    --xggcm_mhd_step_debug_dump \
    --ggcm_mhd_step_legacy_dt_handling false \
    \
    --timelo 1000. \
    \
    2>&1 | tee log

#./plot.py
