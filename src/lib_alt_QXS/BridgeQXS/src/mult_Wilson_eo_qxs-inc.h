/*!
      @file    mult_Wilson_eo_qxs-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#ifndef MULT_WILSON_EO_QXS_INCLUDED
#define MULT_WILSON_EO_QXS_INCLUDED

#include "mult_common_th-inc.h"

//#define IMPLE  0   // 0: original, 1: Nitadori-san's
#define IMPLE    1 // 0: original, 1: Nitadori-san's
//#define IMPLE  2   // setting of predicate/index only once

//====================================================================
void BridgeQXS::mult_wilson_eo_bulk_dirac(real_t *v2, real_t *up,
                                          real_t *v1, real_t *xp,
                                          real_t kappa, int *bc,
                                          int *Nsize, int *do_comm,
                                          int *Leo, const int ieo,
                                          const int iflag)
{
  int Nx2v = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];

  int Nst2v = Nx2v * Ny * Nz * Nt;
  int Nst2  = Nst2v * VLEN;

  int Nxy2  = Nx2v * Ny;
  int Nxyz2 = Nx2v * Ny * Nz;

  svbool_t pg1e_xp, pg2e_xp, pg3e_xp, pg1e_xm, pg2e_xm, pg3e_xm;
  svbool_t pg1o_xp, pg2o_xp, pg3o_xp, pg1o_xm, pg2o_xm, pg3o_xm;
  svbool_t pg1_yp, pg2_yp, pg1_ym, pg2_ym;

#if IMPLE == 0
  // original
  set_predicate_xp_eo(pg1e_xp, pg2e_xp, pg3e_xp, 0);
  set_predicate_xp_eo(pg1o_xp, pg2o_xp, pg3o_xp, 1);
  set_predicate_xm_eo(pg1e_xm, pg2e_xm, pg3e_xm, 0);
  set_predicate_xm_eo(pg1o_xm, pg2o_xm, pg3o_xm, 1);
  set_predicate_yp(pg1_yp, pg2_yp);
  set_predicate_ym(pg1_ym, pg2_ym);
#else
  // nitadori-san
  svuint_t idx1e_xp, idx1o_xp, idx1e_xm, idx1o_xm;
  svuint_t idx1_yp, idx1_ym;
  set_idx_predicate_xp_eo(pg1e_xp, idx1e_xp, 0);
  set_idx_predicate_xp_eo(pg1o_xp, idx1o_xp, 1);
  set_idx_predicate_xm_eo(pg1e_xm, idx1e_xm, 0);
  set_idx_predicate_xm_eo(pg1o_xm, idx1o_xm, 1);
  set_idx_predicate_yp(pg1_yp, idx1_yp);
  set_idx_predicate_ym(pg1_ym, idx1_ym);
#endif

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nst2v);

  for (int site = is; site < ns; ++site) {
    int ix   = site % Nx2v;
    int iyzt = site / Nx2v;
    int iy   = iyzt % Ny;
    int izt  = site / Nxy2;
    int iz   = izt % Nz;
    int it   = izt / Nz;
    int ixy  = ix + Nx2v * iy;
    int ixyz = ixy + Nxy2 * iz;
    int jeo  = (ieo + Leo[VLENY * iyzt]) % 2;

    Vsimd_t v2v[NVCD];
    clear_vec(v2v, NVCD);

    real_t zL[VLEN * NVCD], uL[VLEN * NDF];

    if ((ix < Nx2v - 1) || (do_comm[0] == 0)) {
      int nei = ix + 1 + Nx2v * iyzt;
      if (ix == Nx2v - 1) nei = 0 + Nx2v * iyzt;
      real_t *u = &up[NDF * Nst2 * (ieo + 2 * 0)];

      if (jeo == 0) {
#if IMPLE == 0
        // original
        mult_wilson_eo_xpb(pg1e_xp, pg2e_xp, pg3e_xp,
                           v2v, &u[VLEN * NDF * site],
                           &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
#else
        // nitadori-san
#if IMPLE < 2
        set_idx_predicate_xp_eo(pg1e_xp, idx1e_xp, 0);
#endif
        mult_wilson_eo_xpb(pg1e_xp, idx1e_xp,
                           v2v, &u[VLEN * NDF * site],
                           &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
#endif
      } else {
#if IMPLE == 0
        // original
        mult_wilson_eo_xpb(pg1o_xp, pg2o_xp, pg3o_xp,
                           v2v, &u[VLEN * NDF * site],
                           &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
#else
        // nitadori-san
#if IMPLE < 2
        set_idx_predicate_xp_eo(pg1o_xp, idx1o_xp, 1);
#endif
        mult_wilson_eo_xpb(pg1o_xp, idx1o_xp,
                           v2v, &u[VLEN * NDF * site],
                           &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
#endif
      }
    }

    if ((ix > 0) || (do_comm[0] == 0)) {
      int nei = ix - 1 + Nx2v * iyzt;
      if (ix == 0) nei = Nx2v - 1 + Nx2v * iyzt;
      real_t *u = &up[NDF * Nst2 * (1 - ieo + 2 * 0)];
      if (jeo == 0) {
#if IMPLE == 0
        // original
        mult_wilson_eo_xmb(pg1e_xm, pg2e_xm, pg3e_xm, v2v,
                           &u[VLEN * NDF * site], &u[VLEN * NDF * nei],
                           &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
#else
        // nitadori-san
#if IMPLE < 2
        set_idx_predicate_xm_eo(pg1e_xm, idx1e_xm, 0);
#endif
        mult_wilson_eo_xmb(pg1e_xm, idx1e_xm, v2v,
                           &u[VLEN * NDF * site], &u[VLEN * NDF * nei],
                           &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
#endif
      } else {
#if IMPLE == 0
        // original
        mult_wilson_eo_xmb(pg1o_xm, pg2o_xm, pg3o_xm, v2v,
                           &u[VLEN * NDF * site], &u[VLEN * NDF * nei],
                           &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
#else
        // nitadori-san
#if IMPLE < 2
        set_idx_predicate_xm_eo(pg1o_xm, idx1o_xm, 1);
#endif
        mult_wilson_eo_xmb(pg1o_xm, idx1o_xm, v2v,
                           &u[VLEN * NDF * site], &u[VLEN * NDF * nei],
                           &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
#endif
      }
    }

    if ((iy < Ny - 1) || (do_comm[1] == 0)) {
      int    iy2 = (iy + 1) % Ny;
      int    nei = ix + Nx2v * (iy2 + Ny * izt);
      real_t *u  = &up[NDF * Nst2 * (ieo + 2 * 1)];
#if IMPLE == 0
      mult_wilson_ypb(pg1_yp, pg2_yp, v2v,
                      &u[VLEN * NDF * site],
                      &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
#else
#if IMPLE < 2
      set_idx_predicate_yp(pg1_yp, idx1_yp);
#endif
      mult_wilson_ypb(pg1_yp, idx1_yp, v2v,
                      &u[VLEN * NDF * site],
                      &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
#endif
    }

    if ((iy > 0) || (do_comm[1] == 0)) {
      int    iy2 = (iy - 1 + Ny) % Ny;
      int    nei = ix + Nx2v * (iy2 + Ny * izt);
      real_t *u  = &up[NDF * Nst2 * (1 - ieo + 2 * 1)];
#if IMPLE == 0
      mult_wilson_ymb(pg1_ym, pg2_ym, v2v,
                      &u[VLEN * NDF * site], &u[VLEN * NDF * nei],
                      &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
#else
#if IMPLE < 2
      set_idx_predicate_ym(pg1_ym, idx1_ym);
#endif
      mult_wilson_ymb(pg1_ym, idx1_ym, v2v,
                      &u[VLEN * NDF * site], &u[VLEN * NDF * nei],
                      &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
#endif
    }

    if ((iz < Nz - 1) || (do_comm[2] == 0)) {
      int    iz2 = (iz + 1) % Nz;
      int    nei = ixy + Nxy2 * (iz2 + Nz * it);
      real_t *u  = &up[NDF * Nst2 * (ieo + 2 * 2)];
      mult_wilson_zpb(v2v, &u[VLEN * NDF * site], &v1[VLEN * NVCD * nei]);
    }

    if ((iz > 0) || (do_comm[2] == 0)) {
      int    iz2 = (iz - 1 + Nz) % Nz;
      int    nei = ixy + Nxy2 * (iz2 + Nz * it);
      real_t *u  = &up[NDF * Nst2 * (1 - ieo + 2 * 2)];
      mult_wilson_zmb(v2v, &u[VLEN * NDF * nei], &v1[VLEN * NVCD * nei]);
    }

    if ((it < Nt - 1) || (do_comm[3] == 0)) {
      int    it2 = (it + 1) % Nt;
      int    nei = ixyz + Nxyz2 * it2;
      real_t *u  = &up[NDF * Nst2 * (ieo + 2 * 3)];
      mult_wilson_tpb_dirac(v2v, &u[VLEN * NDF * site], &v1[VLEN * NVCD * nei]);
    }

    if ((it > 0) || (do_comm[3] == 0)) {
      int    it2 = (it - 1 + Nt) % Nt;
      int    nei = ixyz + Nxyz2 * it2;
      real_t *u  = &up[NDF * Nst2 * (1 - ieo + 2 * 3)];
      mult_wilson_tmb_dirac(v2v, &u[VLEN * NDF * nei], &v1[VLEN * NVCD * nei]);
    }

    svbool_t pg = set_predicate();
    if (iflag == 0) {
      real_t *vv2 = &v2[VLEN * NVCD * site];
      for (int i = 0; i < NVCD; ++i) {
        svreal_t v2t;
        load_vec(pg, v2t, &v2v[i].v[0]);
        scal_vec(pg, v2t, -kappa);
        save_vec(pg, &vv2[VLEN * i], v2t);
      }
    } else {
      mult_wilson_aypx_save(&v2[VLEN * NVCD * site],
                            kappa, v2v, &xp[VLEN * NVCD * site]);
    }
  }
}


//====================================================================
void BridgeQXS::mult_wilson_eo_1_dirac(real_t *buf_xp, real_t *buf_xm,
                                       real_t *buf_yp, real_t *buf_ym,
                                       real_t *buf_zp, real_t *buf_zm,
                                       real_t *buf_tp, real_t *buf_tm,
                                       real_t *up, real_t *v1, int *bc,
                                       int *Nsize, int *do_comm, int *Leo,
                                       const int ieo, const int iflag)
{
  int Nx2v = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];

  int Nst2v = Nx2v * Ny * Nz * Nt;
  int Nst2  = Nst2v * VLEN;

  int Nxy2  = Nx2v * Ny;
  int Nxyz2 = Nx2v * Ny * Nz;

  svbool_t pg1e_xp, pg2e_xp, pg3e_xp, pg1e_xm, pg2e_xm, pg3e_xm;
  svbool_t pg1o_xp, pg2o_xp, pg3o_xp, pg1o_xm, pg2o_xm, pg3o_xm;
  set_predicate_xp_eo(pg1e_xp, pg2e_xp, pg3e_xp, 0);
  set_predicate_xp_eo(pg1o_xp, pg2o_xp, pg3o_xp, 1);
  set_predicate_xm_eo(pg1e_xm, pg2e_xm, pg3e_xm, 0);
  set_predicate_xm_eo(pg1o_xm, pg2o_xm, pg3o_xm, 1);
  svbool_t pg1_yp, pg2_yp, pg1_ym, pg2_ym;
  set_predicate_yp(pg1_yp, pg2_yp);
  set_predicate_ym(pg1_ym, pg2_ym);
#if IMPLE == 0
  svint_t svidx_xp, svidx_xm;
#else
  svuint_t svidx_xp, svidx_xm;
#endif
  set_index_xp_eo(svidx_xp);
  set_index_xm_eo(svidx_xm);

  int Nskipx = (VLENY + 1) / 2;


  if (do_comm[0] == 1) {
    int idir = 0;

    int Nyzt = Ny * Nz * Nt;

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, Nyzt);

    for (int iyzt = is; iyzt < ns; ++iyzt) {
      int jeo = (ieo + Leo[VLENY * iyzt]) % 2;
      {
        int ix     = 0;
        int site   = ix + Nx2v * iyzt;
        int ibf_up = Nskipx * NVC * ND2 * iyzt;
        if (jeo == 0) {
#if VLENY > 1
#if IMPLE == 0
          set_index_xm_eo(svidx_xm);
          mult_wilson_eo_xp1(pg2o_xm, svidx_xm,
                             &buf_xp[ibf_up], &v1[VLEN * NVCD * site]);
#else
          mult_wilson_eo_xp1(pg2o_xm,
                             &buf_xp[ibf_up], &v1[VLEN * NVCD * site]);
#endif
#endif
        } else {
          if (VLENY == 1) ibf_up = Nskipx * NVC * ND2 * (iyzt / 2);
#if IMPLE == 0
          set_index_xm_eo(svidx_xm);
          mult_wilson_eo_xp1(pg2e_xm, svidx_xm,
                             &buf_xp[ibf_up], &v1[VLEN * NVCD * site]);
#else
          mult_wilson_eo_xp1(pg2e_xm,
                             &buf_xp[ibf_up], &v1[VLEN * NVCD * site]);
#endif
        }
      }

      {
        int    ix     = Nx2v - 1;
        int    site   = ix + Nx2v * iyzt;
        int    ibf_dn = Nskipx * NVC * ND2 * iyzt;
        real_t *u     = &up[NDF * Nst2 * (1 - ieo + 2 * idir)];
        if (jeo == 0) {
          if (VLENY == 1) ibf_dn = Nskipx * NVC * ND2 * (iyzt / 2);
#if IMPLE == 0
          set_index_xp_eo(svidx_xp);
          mult_wilson_eo_xm1(pg2o_xp, svidx_xp, &buf_xm[ibf_dn],
                             &u[VLEN * NDF * site], &v1[VLEN * NVCD * site]);
#else
          mult_wilson_eo_xm1(pg2o_xp, &buf_xm[ibf_dn],
                             &u[VLEN * NDF * site], &v1[VLEN * NVCD * site]);
#endif
        } else {
#if VLENY > 1
#if IMPLE == 0
          set_index_xp_eo(svidx_xp);
          mult_wilson_eo_xm1(pg2e_xp, svidx_xp, &buf_xm[ibf_dn],
                             &u[VLEN * NDF * site], &v1[VLEN * NVCD * site]);
#else
          mult_wilson_eo_xm1(pg2e_xp, &buf_xm[ibf_dn],
                             &u[VLEN * NDF * site], &v1[VLEN * NVCD * site]);
#endif
#endif
        }
      }
    }
  }

  if (do_comm[1] > 0) {
    int idir = 1;

    int Nxzt = Nx2v * Nz * Nt;

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, Nxzt);

    for (int ixzt = is; ixzt < ns; ++ixzt) {
      int ix  = ixzt % Nx2v;
      int izt = ixzt / Nx2v;
      {
        int iy   = 0;
        int site = ix + Nx2v * (iy + Ny * izt);
        int ibf  = VLENX * NVC * ND2 * (ix + Nx2v * izt);
        mult_wilson_yp1(pg2_ym,
                        &buf_yp[ibf], &v1[VLEN * NVCD * site]);
      }
      {
        int    iy   = Ny - 1;
        int    site = ix + Nx2v * (iy + Ny * izt);
        real_t *u   = &up[NDF * Nst2 * ((1 - ieo) + 2 * idir)];
        int    ibf  = VLENX * NVC * ND2 * (ix + Nx2v * izt);
        mult_wilson_ym1(pg2_yp,
                        &buf_ym[ibf], &u[VLEN * NDF * site],
                        &v1[VLEN * NVCD * site]);
      }
    }
  }

  if (do_comm[2] > 0) {
    int idir  = 2;
    int Nxyt2 = Nxy2 * Nt;

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, Nxyt2);

    for (int ixyt = is; ixyt < ns; ++ixyt) {
      int ixy = ixyt % Nxy2;
      int it  = ixyt / Nxy2;
      {
        int iz   = 0;
        int site = ixy + Nxy2 * (iz + Nz * it);
        int ibf  = VLEN * NVC * ND2 * (ixy + Nxy2 * it);
        mult_wilson_zp1(&buf_zp[ibf], &v1[VLEN * NVCD * site]);
      }
      {
        int    iz   = Nz - 1;
        int    site = ixy + Nxy2 * (iz + Nz * it);
        int    ibf  = VLEN * NVC * ND2 * (ixy + Nxy2 * it);
        real_t *u   = &up[NDF * Nst2 * (1 - ieo + 2 * idir)];
        mult_wilson_zm1(&buf_zm[ibf], &u[VLEN * NDF * site],
                        &v1[VLEN * NVCD * site]);
      }
    }
  }

  if (do_comm[3] > 0) {
    int idir = 3;

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, Nxyz2);

    for (int ixyz = is; ixyz < ns; ++ixyz) {
      {
        int it   = 0;
        int site = ixyz + Nxyz2 * it;
        mult_wilson_tp1_dirac(&buf_tp[VLEN * NVC * ND2 * ixyz],
                              &v1[VLEN * NVCD * site]);
      }
      {
        int    it   = Nt - 1;
        int    site = ixyz + Nxyz2 * it;
        real_t *u   = &up[NDF * Nst2 * (1 - ieo + 2 * idir)];
        mult_wilson_tm1_dirac(&buf_tm[VLEN * NVC * ND2 * ixyz],
                              &u[VLEN * NDF * site], &v1[VLEN * NVCD * site]);
      }
    }
  }
}


//====================================================================
void BridgeQXS::mult_wilson_eo_2_dirac(real_t *v2, real_t *up, real_t *v1,
                                       real_t *xp,
                                       real_t *buf_xp, real_t *buf_xm,
                                       real_t *buf_yp, real_t *buf_ym,
                                       real_t *buf_zp, real_t *buf_zm,
                                       real_t *buf_tp, real_t *buf_tm,
                                       real_t kappa, int *bc,
                                       int *Nsize, int *do_comm, int *Leo,
                                       const int ieo, const int iflag)
{
  int Nx2v = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];

  int Nst2v = Nx2v * Ny * Nz * Nt;
  int Nst2  = Nst2v * VLEN;

  int Nxy2  = Nx2v * Ny;
  int Nxyz2 = Nx2v * Ny * Nz;

  svbool_t pg1e_xp, pg2e_xp, pg3e_xp, pg1e_xm, pg2e_xm, pg3e_xm;
  svbool_t pg1o_xp, pg2o_xp, pg3o_xp, pg1o_xm, pg2o_xm, pg3o_xm;
  set_predicate_xp_eo(pg1e_xp, pg2e_xp, pg3e_xp, 0);
  set_predicate_xp_eo(pg1o_xp, pg2o_xp, pg3o_xp, 1);
  set_predicate_xm_eo(pg1e_xm, pg2e_xm, pg3e_xm, 0);
  set_predicate_xm_eo(pg1o_xm, pg2o_xm, pg3o_xm, 1);
  svbool_t pg1_yp, pg2_yp, pg1_ym, pg2_ym;
  set_predicate_yp(pg1_yp, pg2_yp);
  set_predicate_ym(pg1_ym, pg2_ym);
#if IMPLE == 0
  svint_t svidx_xp, svidx_xm;
#else
  svuint_t svidx_xp, svidx_xm;
#endif
  set_index_xp_eo(svidx_xp);
  set_index_xm_eo(svidx_xm);
  svbool_t pg = set_predicate();

  real_t kappa_eo = kappa;
  if (iflag == 0) {
    kappa_eo = -kappa;
  }

  int Nskipx = (VLENY + 1) / 2;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nst2v);

  for (int site = is; site < ns; ++site) {
    int ix   = site % Nx2v;
    int iyzt = site / Nx2v;
    int iy   = iyzt % Ny;
    int izt  = site / Nxy2;
    int iz   = izt % Nz;
    int it   = izt / Nz;
    int ixy  = ix + Nx2v * iy;
    int ixyz = ixy + Nxy2 * iz;
    int jeo  = (ieo + Leo[VLENY * iyzt]) % 2;

    Vsimd_t v2v[NVCD];
    clear_vec(v2v, NVCD);

    real_t zL[VLEN * NVCD], uL[VLEN * NDF];

    if ((ix == Nx2v - 1) && (do_comm[0] > 0)) {
      int ibf_up = Nskipx * NVC * ND2 * iyzt;
      if (VLENY == 1) ibf_up = Nskipx * NVC * ND2 * (iyzt / 2);
      real_t *u = &up[NDF * Nst2 * (ieo + 2 * 0)];
      if (jeo == 0) {
        set_index_xp_eo(svidx_xp);
        mult_wilson_eo_xp2(pg1e_xp, pg2e_xp, pg3e_xp, svidx_xp,
                           v2v, &u[VLEN * NDF * site],
                           &v1[VLEN * NVCD * site], &buf_xp[ibf_up]);
      } else {
        set_index_xp_eo(svidx_xp);
        mult_wilson_eo_xp2(pg1o_xp, pg2o_xp, pg3o_xp, svidx_xp,
                           v2v, &u[VLEN * NDF * site],
                           &v1[VLEN * NVCD * site], &buf_xp[ibf_up]);
      }
    }

    if ((ix == 0) && (do_comm[0] > 0)) {
      int ibf_dn = Nskipx * NVC * ND2 * iyzt;
      if (VLENY == 1) ibf_dn = Nskipx * NVC * ND2 * (iyzt / 2);
      real_t *u = &up[NDF * Nst2 * (1 - ieo + 2 * 0)];
      if (jeo == 0) {
        set_index_xm_eo(svidx_xm);
        mult_wilson_eo_xm2(pg1e_xm, pg2e_xm, pg3e_xm, svidx_xm,
                           v2v, &u[VLEN * NDF * site],
                           &v1[VLEN * NVCD * site], &buf_xm[ibf_dn]);
      } else {
        set_index_xm_eo(svidx_xm);
        mult_wilson_eo_xm2(pg1o_xm, pg2o_xm, pg3o_xm, svidx_xm,
                           v2v, &u[VLEN * NDF * site],
                           &v1[VLEN * NVCD * site], &buf_xm[ibf_dn]);
      }
    }

    if ((iy == Ny - 1) && (do_comm[1] > 0)) {
      int    ibf = VLENX * NVC * ND2 * (ix + Nx2v * izt);
      real_t *u  = &up[NDF * Nst2 * (ieo + 2 * 1)];
      mult_wilson_yp2(pg1_yp, pg2_yp,
                      v2v, &u[VLEN * NDF * site],
                      &v1[VLEN * NVCD * site], &buf_yp[ibf]);
    }

    if ((iy == 0) && (do_comm[1] > 0)) {
      int    ibf = VLENX * NVC * ND2 * (ix + Nx2v * izt);
      real_t *u  = &up[NDF * Nst2 * (1 - ieo + 2 * 1)];
      mult_wilson_ym2(pg1_ym, pg2_ym,
                      v2v, &u[VLEN * NDF * site],
                      &v1[VLEN * NVCD * site], &buf_ym[ibf]);
    }

    if ((iz == Nz - 1) && (do_comm[2] > 0)) {
      int    ibf = VLEN * NVC * ND2 * (ixy + Nxy2 * it);
      real_t *u  = &up[NDF * Nst2 * (ieo + 2 * 2)];
      mult_wilson_zp2(v2v, &u[VLEN * NDF * site], &buf_zp[ibf]);
    }

    if ((iz == 0) && (do_comm[2] > 0)) {
      int ibf = VLEN * NVC * ND2 * (ixy + Nxy2 * it);
      mult_wilson_zm2(v2v, &buf_zm[ibf]);
    }

    if ((it == Nt - 1) && (do_comm[3] > 0)) {
      real_t *u = &up[NDF * Nst2 * (ieo + 2 * 3)];
      mult_wilson_tp2_dirac(v2v, &u[VLEN * NDF * site],
                            &buf_tp[VLEN * NVC * ND2 * ixyz]);
    }

    if ((it == 0) && (do_comm[3] > 0)) {
      mult_wilson_tm2_dirac(v2v, &buf_tm[VLEN * NVC * ND2 * ixyz]);
    }


    real_t *ww2 = &v2v[0].v[0];
    real_t *vv2 = &v2[VLEN * NVCD * site];
    for (int i = 0; i < NVCD; i += 8) {
      svreal_t wt1, wt2, wt3, wt4;
      svreal_t vt1, vt2, vt3, vt4;
      svreal_t wt5, wt6, wt7, wt8;
      svreal_t vt5, vt6, vt7, vt8;
      load_vec(pg, vt1, &vv2[VLEN * (i)]);
      load_vec(pg, vt2, &vv2[VLEN * (i + 1)]);
      load_vec(pg, vt3, &vv2[VLEN * (i + 2)]);
      load_vec(pg, vt4, &vv2[VLEN * (i + 3)]);

      load_vec(pg, wt1, &ww2[VLEN * (i)]);
      load_vec(pg, wt2, &ww2[VLEN * (i + 1)]);
      load_vec(pg, wt3, &ww2[VLEN * (i + 2)]);
      load_vec(pg, wt4, &ww2[VLEN * (i + 3)]);

      axpy_vec(pg, vt1, kappa_eo, wt1);
      axpy_vec(pg, vt2, kappa_eo, wt2);
      axpy_vec(pg, vt3, kappa_eo, wt3);

      load_vec(pg, vt5, &vv2[VLEN * (i + 4)]);
      load_vec(pg, vt6, &vv2[VLEN * (i + 5)]);
      load_vec(pg, vt7, &vv2[VLEN * (i + 6)]);
      load_vec(pg, vt8, &vv2[VLEN * (i + 7)]);

      load_vec(pg, wt5, &ww2[VLEN * (i + 4)]);
      load_vec(pg, wt6, &ww2[VLEN * (i + 5)]);
      load_vec(pg, wt7, &ww2[VLEN * (i + 6)]);
      load_vec(pg, wt8, &ww2[VLEN * (i + 7)]);

      axpy_vec(pg, vt4, kappa_eo, wt4);
      axpy_vec(pg, vt5, kappa_eo, wt5);
      axpy_vec(pg, vt6, kappa_eo, wt6);
      axpy_vec(pg, vt7, kappa_eo, wt7);
      axpy_vec(pg, vt8, kappa_eo, wt8);

      save_vec(pg, &vv2[VLEN * (i)], vt1);
      save_vec(pg, &vv2[VLEN * (i + 1)], vt2);
      save_vec(pg, &vv2[VLEN * (i + 2)], vt3);
      save_vec(pg, &vv2[VLEN * (i + 3)], vt4);
      save_vec(pg, &vv2[VLEN * (i + 4)], vt5);
      save_vec(pg, &vv2[VLEN * (i + 5)], vt6);
      save_vec(pg, &vv2[VLEN * (i + 6)], vt7);
      save_vec(pg, &vv2[VLEN * (i + 7)], vt8);
    }
  }
}


#endif
//============================================================END=====
