/*!
      @file    mult_Clover_qxs-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#ifndef MULT_CLOVER_QXS_INCLUDED
#define MULT_CLOVER_QXS_INCLUDED

#include "mult_common_th-inc.h"

//====================================================================
void BridgeQXS::mult_clover_bulk_dirac(real_t *v2, real_t *up,
                                       real_t *ct, real_t *v1,
                                       real_t kappa, int *bc,
                                       int *Nsize, int *do_comm)
{
  int Nxv  = Nsize[0];
  int Nyv  = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Nyv * Nz * Nt;
  int Nst  = Nstv * VLEN;

  svbool_t pg1_xp, pg2_xp, pg1_xm, pg2_xm;
  svbool_t pg1_yp, pg2_yp, pg1_ym, pg2_ym;
  set_predicate_xp(pg1_xp, pg2_xp);
  set_predicate_xm(pg1_xm, pg2_xm);
  set_predicate_yp(pg1_yp, pg2_yp);
  set_predicate_ym(pg1_ym, pg2_ym);

  int Nxy  = Nxv * Nyv;
  int Nxyz = Nxv * Nyv * Nz;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstv);

  for (int site = is; site < ns; ++site) {
    int ix   = site % Nxv;
    int iyzt = site / Nxv;
    int iy   = iyzt % Nyv;
    int izt  = site / Nxy;
    int iz   = izt % Nz;
    int it   = izt / Nz;
    int ixy  = ix + Nxv * iy;
    int ixyz = ixy + Nxy * iz;

    Vsimd_t v2v[NVCD];
    clear_vec(v2v, NVCD);

    real_t zL[VLEN * NVCD];
    real_t uL[VLEN * NDF];

    if (ix < Nxv - 1) {
      real_t *u  = &up[NDF * Nst * 0];
      int    nei = ix + 1 + Nxv * iyzt;
      mult_wilson_xpb(pg1_xp, pg2_xp, v2v, &u[VLEN * NDF * site],
                      &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
    } else if (do_comm[0] == 0) {  // ix = Nxv-1
      real_t *u  = &up[NDF * Nst * 0];
      int    nei = 0 + Nxv * iyzt;
      mult_wilson_xpb(pg1_xp, pg2_xp, v2v, &u[VLEN * NDF * site],
                      &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
    }

    if (ix > 0) {
      real_t *u  = &up[NDF * Nst * 0];
      int    nei = ix - 1 + Nxv * iyzt;
      mult_wilson_xmb(pg1_xm, pg2_xm, v2v,
                      &u[VLEN * NDF * site], &u[VLEN * NDF * nei],
                      &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
    } else if (do_comm[0] == 0) {   // ix = 0
      real_t *u  = &up[NDF * Nst * 0];
      int    nei = Nxv - 1 + Nxv * iyzt;
      mult_wilson_xmb(pg1_xm, pg2_xm, v2v,
                      &u[VLEN * NDF * site], &u[VLEN * NDF * nei],
                      &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
    }

    if (iy < Nyv - 1) {
      int    iy2 = (iy + 1) % Nyv;
      int    nei = ix + Nxv * (iy2 + Nyv * izt);
      real_t *u  = &up[NDF * Nst * 1];
      mult_wilson_ypb(pg1_yp, pg2_yp, v2v,
                      &u[VLEN * NDF * site],
                      &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
    } else if (do_comm[1] == 0) {  // iy = Nyv-1
      int    iy2 = (iy + 1) % Nyv;
      int    nei = ix + Nxv * (iy2 + Nyv * izt);
      real_t *u  = &up[NDF * Nst * 1];
      mult_wilson_ypb(pg1_yp, pg2_yp, v2v,
                      &u[VLEN * NDF * site],
                      &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
    }

    if (iy > 0) {
      int    iy2 = (iy - 1 + Nyv) % Nyv;
      int    nei = ix + Nxv * (iy2 + Nyv * izt);
      real_t *u  = &up[NDF * Nst * 1];
      mult_wilson_ymb(pg1_ym, pg2_ym, v2v,
                      &u[VLEN * NDF * site], &u[VLEN * NDF * nei],
                      &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
    } else if (do_comm[1] == 0) {  // iy = 0
      int    iy2 = (iy - 1 + Nyv) % Nyv;
      int    nei = ix + Nxv * (iy2 + Nyv * izt);
      real_t *u  = &up[NDF * Nst * 1];
      mult_wilson_ymb(pg1_ym, pg2_ym, v2v,
                      &u[VLEN * NDF * site], &u[VLEN * NDF * nei],
                      &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
    }

    if ((iz < Nz - 1) || (do_comm[2] == 0)) {
      int    iz2 = (iz + 1) % Nz;
      int    nei = ixy + Nxy * (iz2 + Nz * it);
      real_t *u  = &up[NDF * Nst * 2];
      mult_wilson_zpb(v2v, &u[VLEN * NDF * site], &v1[VLEN * NVCD * nei]);
    }

    if ((iz > 0) || (do_comm[2] == 0)) {
      int    iz2 = (iz - 1 + Nz) % Nz;
      int    nei = ixy + Nxy * (iz2 + Nz * it);
      real_t *u  = &up[NDF * Nst * 2];
      mult_wilson_zmb(v2v, &u[VLEN * NDF * nei], &v1[VLEN * NVCD * nei]);
    }

    if ((it < Nt - 1) || (do_comm[3] == 0)) {
      int    it2 = (it + 1) % Nt;
      int    nei = ixyz + Nxyz * it2;
      real_t *u  = &up[NDF * Nst * 3];
      mult_wilson_tpb_dirac(v2v, &u[VLEN * NDF * site],
                            &v1[VLEN * NVCD * nei]);
    }

    if ((it > 0) || (do_comm[3] == 0)) {
      int    it2 = (it - 1 + Nt) % Nt;
      int    nei = ixyz + Nxyz * it2;
      real_t *u  = &up[NDF * Nst * 3];
      mult_wilson_tmb_dirac(v2v, &u[VLEN * NDF * nei],
                            &v1[VLEN * NVCD * nei]);
    }

    mult_clover_csw_aypx(&v2[VLEN * NVCD * site], -kappa, &v2v[0],
                         &ct[VLEN * NDF * 2 * ND * site], &v1[VLEN * NVCD * site]);
  }
}


//====================================================================
void BridgeQXS::mult_clover_bulk_dirac_chrot(
  real_t *v2, real_t *up,
  real_t *ct, real_t *v1,
  real_t kappa, int *bc,
  int *Nsize, int *do_comm)
{
  int Nxv  = Nsize[0];
  int Nyv  = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Nyv * Nz * Nt;
  int Nst  = Nstv * VLEN;

  svbool_t pg1_xp, pg2_xp, pg1_xm, pg2_xm;
  svbool_t pg1_yp, pg2_yp, pg1_ym, pg2_ym;
  set_predicate_xp(pg1_xp, pg2_xp);
  set_predicate_xm(pg1_xm, pg2_xm);
  set_predicate_yp(pg1_yp, pg2_yp);
  set_predicate_ym(pg1_ym, pg2_ym);

  int Nxy  = Nxv * Nyv;
  int Nxyz = Nxv * Nyv * Nz;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstv);

  for (int site = is; site < ns; ++site) {
    int ix   = site % Nxv;
    int iyzt = site / Nxv;
    int iy   = iyzt % Nyv;
    int izt  = site / Nxy;
    int iz   = izt % Nz;
    int it   = izt / Nz;
    int ixy  = ix + Nxv * iy;
    int ixyz = ixy + Nxy * iz;

    Vsimd_t v2v[NVCD];
    clear_vec(v2v, NVCD);


    if (ix < Nxv - 1) {
      real_t *u  = &up[NDF * Nst * 0];
      int    nei = ix + 1 + Nxv * iyzt;
      mult_wilson_xpb(pg1_xp, pg2_xp, &v2v[0].v[0], &u[VLEN * NDF * site],
                      &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
    } else if (do_comm[0] == 0) {  // ix = Nxv-1
      real_t *u  = &up[NDF * Nst * 0];
      int    nei = 0 + Nxv * iyzt;
      mult_wilson_xpb(pg1_xp, pg2_xp, v2v, &u[VLEN * NDF * site],
                      &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
    }

    if (ix > 0) {
      real_t *u  = &up[NDF * Nst * 0];
      int    nei = ix - 1 + Nxv * iyzt;
      mult_wilson_xmb(pg1_xm, pg2_xm, v2v,
                      &u[VLEN * NDF * site], &u[VLEN * NDF * nei],
                      &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
    } else if (do_comm[0] == 0) {   // ix = 0
      real_t *u  = &up[NDF * Nst * 0];
      int    nei = Nxv - 1 + Nxv * iyzt;
      mult_wilson_xmb(pg1_xm, pg2_xm, v2v,
                      &u[VLEN * NDF * site], &u[VLEN * NDF * nei],
                      &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
    }

    if (iy < Nyv - 1) {
      int    iy2 = (iy + 1) % Nyv;
      int    nei = ix + Nxv * (iy2 + Nyv * izt);
      real_t *u  = &up[NDF * Nst * 1];
      mult_wilson_ypb(pg1_yp, pg2_yp, v2v,
                      &u[VLEN * NDF * site],
                      &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
    } else if (do_comm[1] == 0) {  // iy = Nyv-1
      int    iy2 = (iy + 1) % Nyv;
      int    nei = ix + Nxv * (iy2 + Nyv * izt);
      real_t *u  = &up[NDF * Nst * 1];
      mult_wilson_ypb(pg1_yp, pg2_yp, v2v,
                      &u[VLEN * NDF * site],
                      &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
    }

    if (iy > 0) {
      int    iy2 = (iy - 1 + Nyv) % Nyv;
      int    nei = ix + Nxv * (iy2 + Nyv * izt);
      real_t *u  = &up[NDF * Nst * 1];
      mult_wilson_ymb(pg1_ym, pg2_ym, v2v,
                      &u[VLEN * NDF * site], &u[VLEN * NDF * nei],
                      &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
    } else if (do_comm[1] == 0) {  // iy = 0
      int    iy2 = (iy - 1 + Nyv) % Nyv;
      int    nei = ix + Nxv * (iy2 + Nyv * izt);
      real_t *u  = &up[NDF * Nst * 1];
      mult_wilson_ymb(pg1_ym, pg2_ym, v2v,
                      &u[VLEN * NDF * site], &u[VLEN * NDF * nei],
                      &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
    }

    if ((iz < Nz - 1) || (do_comm[2] == 0)) {
      int    iz2 = (iz + 1) % Nz;
      int    nei = ixy + Nxy * (iz2 + Nz * it);
      real_t *u  = &up[NDF * Nst * 2];
      mult_wilson_zpb(v2v, &u[VLEN * NDF * site], &v1[VLEN * NVCD * nei]);
    }

    if ((iz > 0) || (do_comm[2] == 0)) {
      int    iz2 = (iz - 1 + Nz) % Nz;
      int    nei = ixy + Nxy * (iz2 + Nz * it);
      real_t *u  = &up[NDF * Nst * 2];
      mult_wilson_zmb(v2v, &u[VLEN * NDF * nei], &v1[VLEN * NVCD * nei]);
    }

    if ((it < Nt - 1) || (do_comm[3] == 0)) {
      int    it2 = (it + 1) % Nt;
      int    nei = ixyz + Nxyz * it2;
      real_t *u  = &up[NDF * Nst * 3];
      mult_wilson_tpb_dirac(v2v, &u[VLEN * NDF * site],
                            &v1[VLEN * NVCD * nei]);
    }

    if ((it > 0) || (do_comm[3] == 0)) {
      int    it2 = (it - 1 + Nt) % Nt;
      int    nei = ixyz + Nxyz * it2;
      real_t *u  = &up[NDF * Nst * 3];
      mult_wilson_tmb_dirac(v2v, &u[VLEN * NDF * nei],
                            &v1[VLEN * NVCD * nei]);
    }

    mult_clover_csw_aypx_chrot(
      &v2[VLEN * NVCD * site], -kappa, &v2v[0].v[0],
      &ct[VLEN * NDF * 2 * ND * site], &v1[VLEN * NVCD * site]);
  }
}


//====================================================================
void BridgeQXS::mult_clover_cswinv_dirac_chrot(real_t *v2,
                                               real_t *ct, real_t *v1,
                                               int *Nsize)
{
  int Nxv  = Nsize[0];
  int Nyv  = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Nyv * Nz * Nt;
  int Nst  = Nstv * VLEN;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstv);

  for (int site = is; site < ns; ++site) {
    mult_cswinv_chrot(&v2[VLEN * NVCD * site],
                      &ct[VLEN * NDF * ND * site], &v1[VLEN * NVCD * site]);
  }
}


//============================================================END=====
#endif
