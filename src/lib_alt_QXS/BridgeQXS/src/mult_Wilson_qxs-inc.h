/*!
      @file    mult_Wilson_qxs-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#ifndef MULT_WILSON_QXS_INCLUDED
#define MULT_WILSON_QXS_INCLUDED

#include "mult_common_th-inc.h"
//#include "mult_Wilson_qxs_parts-inc.h"

//====================================================================
void BridgeQXS::mult_wilson_bulk_dirac(real_t *v2, real_t *up, real_t *v1,
                                       real_t kappa, int *bc, int *Nsize, int *do_comm)
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

  Vsimd_t v2v[NVCD];

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

    mult_wilson_aypx_save(&v2[VLEN * NVCD * site],
                          -kappa, v2v, &v1[VLEN * NVCD * site]);
  }
}


//====================================================================
void BridgeQXS::mult_wilson_1_dirac(
  real_t *buf_xp, real_t *buf_xm,
  real_t *buf_yp, real_t *buf_ym,
  real_t *buf_zp, real_t *buf_zm,
  real_t *buf_tp, real_t *buf_tm,
  real_t *up, real_t *v1,
  int *bc, int *Nsize, int *do_comm)
{
  int Nxv  = Nsize[0];
  int Nyv  = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Nyv * Nz * Nt;
  int Nst  = Nstv * VLEN;

  int Nxy  = Nxv * Nyv;
  int Nxyz = Nxv * Nyv * Nz;

  svbool_t pg1_xp, pg2_xp, pg1_xm, pg2_xm;
  svbool_t pg1_yp, pg2_yp, pg1_ym, pg2_ym;
  set_predicate_xp(pg1_xp, pg2_xp);
  set_predicate_xm(pg1_xm, pg2_xm);
  set_predicate_yp(pg1_yp, pg2_yp);
  set_predicate_ym(pg1_ym, pg2_ym);
  svint_t svidx_xp, svidx_xm;
  set_index_xp(svidx_xp);
  set_index_xm(svidx_xm);

  if (do_comm[0] > 0) {
    int    idir = 0;
    real_t *u   = &up[NDF * Nst * idir];

    int Nyzt = Nyv * Nz * Nt;

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, Nyzt);

    for (int iyzt = is; iyzt < ns; ++iyzt) {
      {
        int ix   = 0;
        int site = ix + Nxv * iyzt;
        int ibf  = VLENY * NVC * ND2 * iyzt;
        set_index_xm(svidx_xm);
        mult_wilson_xp1(pg2_xm, svidx_xm,
                        &buf_xp[ibf], &v1[VLEN * NVCD * site]);
      }
      {
        int ix   = Nxv - 1;
        int site = ix + Nxv * iyzt;
        int ibf  = VLENY * NVC * ND2 * iyzt;
        set_index_xp(svidx_xp);
        mult_wilson_xm1(pg2_xp, svidx_xp,
                        &buf_xm[ibf], &u[VLEN * NDF * site],
                        &v1[VLEN * NVCD * site]);
      }
    }
  }

  if (do_comm[1] > 0) {
    int    idir = 1;
    real_t *u   = &up[NDF * Nst * idir];

    int Nxzt = Nxv * Nz * Nt;

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, Nxzt);

    for (int ixzt = is; ixzt < ns; ++ixzt) {
      int ix  = ixzt % Nxv;
      int izt = ixzt / Nxv;
      {
        int iy   = 0;
        int site = ix + Nxv * (iy + Nyv * izt);
        int ibf  = VLENX * NVC * ND2 * ixzt;
        mult_wilson_yp1(pg2_ym,
                        &buf_yp[ibf], &v1[VLEN * NVCD * site]);
      }
      {
        int iy   = Nyv - 1;
        int site = ix + Nxv * (iy + Nyv * izt);
        int ibf  = VLENX * NVC * ND2 * ixzt;
        mult_wilson_ym1(pg2_yp,
                        &buf_ym[ibf], &u[VLEN * NDF * site],
                        &v1[VLEN * NVCD * site]);
      }
    }
  }

  if (do_comm[2] > 0) {
    int    idir = 2;
    real_t *u   = &up[NDF * Nst * idir];

    int Nxyt = Nxv * Nyv * Nt;

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, Nxyt);

    for (int ixyt = is; ixyt < ns; ++ixyt) {
      int ixy = ixyt % Nxy;
      int it  = ixyt / Nxy;
      {
        int iz   = 0;
        int site = ixy + Nxy * (iz + Nz * it);
        int ibf  = VLEN * NVC * ND2 * (ixy + Nxy * it);
        mult_wilson_zp1(&buf_zp[ibf], &v1[VLEN * NVCD * site]);
      }
      {
        int iz   = Nz - 1;
        int site = ixy + Nxy * (iz + Nz * it);
        int ibf  = VLEN * NVC * ND2 * (ixy + Nxy * it);
        mult_wilson_zm1(&buf_zm[ibf], &u[VLEN * NDF * site],
                        &v1[VLEN * NVCD * site]);
      }
    }
  }

  if (do_comm[3] > 0) {
    int    idir = 3;
    real_t *u   = &up[NDF * Nst * idir];

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, Nxyz);
    {
      int it = 0;
      for (int ixyz = is; ixyz < ns; ++ixyz) {
        int site = ixyz + Nxyz * it;
        mult_wilson_tp1_dirac(&buf_tp[VLEN * NVC * ND2 * ixyz],
                              &v1[VLEN * NVCD * site]);
      }
    }
    {
      int it = Nt - 1;
      for (int ixyz = is; ixyz < ns; ++ixyz) {
        int site = ixyz + Nxyz * it;
        mult_wilson_tm1_dirac(&buf_tm[VLEN * NVC * ND2 * ixyz],
                              &u[VLEN * NDF * site], &v1[VLEN * NVCD * site]);
      }
    }
  }
}


//====================================================================
void BridgeQXS::mult_wilson_2_dirac(real_t *v2, real_t *up, real_t *v1,
                                    real_t *buf_xp, real_t *buf_xm,
                                    real_t *buf_yp, real_t *buf_ym,
                                    real_t *buf_zp, real_t *buf_zm,
                                    real_t *buf_tp, real_t *buf_tm,
                                    real_t kappa, int *bc, int *Nsize, int *do_comm)
{
  int Nxv  = Nsize[0];
  int Nyv  = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Nyv * Nz * Nt;
  int Nst  = Nstv * VLEN;

  int Nxy  = Nxv * Nyv;
  int Nxyz = Nxv * Nyv * Nz;

  svbool_t pg1_xp, pg2_xp, pg1_xm, pg2_xm;
  svbool_t pg1_yp, pg2_yp, pg1_ym, pg2_ym;
  set_predicate_xp(pg1_xp, pg2_xp);
  set_predicate_xm(pg1_xm, pg2_xm);
  set_predicate_yp(pg1_yp, pg2_yp);
  set_predicate_ym(pg1_ym, pg2_ym);
  svint_t svidx_xp, svidx_xm;
  set_index_xp(svidx_xp);
  set_index_xm(svidx_xm);

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
    int opr_any = 0;

    if ((ix == Nxv - 1) && (do_comm[0] > 0)) {
      real_t *u  = &up[NDF * Nst * 0];
      int    ibf = VLENY * NVC * ND2 * iyzt;
      set_index_xp(svidx_xp);
      mult_wilson_xp2(pg1_xp, pg2_xp, svidx_xp,
                      v2v, &u[VLEN * NDF * site],
                      &v1[VLEN * NVCD * site], &buf_xp[ibf]);
      ++opr_any;
    }

    if ((ix == 0) && (do_comm[0] > 0)) {
      real_t *u  = &up[NDF * Nst * 0];
      int    ibf = VLENY * NVC * ND2 * iyzt;
      set_index_xm(svidx_xm);
      mult_wilson_xm2(pg1_xm, pg2_xm, svidx_xm,
                      v2v, &u[VLEN * NDF * site],
                      &v1[VLEN * NVCD * site], &buf_xm[ibf]);
      ++opr_any;
    }

    if ((iy == Nyv - 1) && (do_comm[1] > 0)) {
      real_t *u  = &up[NDF * Nst * 1];
      int    ibf = VLENX * NVC * ND2 * (ix + Nxv * izt);
      mult_wilson_yp2(pg1_yp, pg2_yp,
                      v2v, &u[VLEN * NDF * site],
                      &v1[VLEN * NVCD * site], &buf_yp[ibf]);
      ++opr_any;
    }

    if ((iy == 0) && (do_comm[1] > 0)) {
      real_t *u  = &up[NDF * Nst * 1];
      int    ibf = VLENX * NVC * ND2 * (ix + Nxv * izt);
      mult_wilson_ym2(pg1_ym, pg2_ym,
                      v2v, &u[VLEN * NDF * site],
                      &v1[VLEN * NVCD * site], &buf_ym[ibf]);
      ++opr_any;
    }

    if ((iz == Nz - 1) && (do_comm[2] > 0)) {
      int    ibf = VLEN * NVC * ND2 * (ixy + Nxy * it);
      real_t *u  = &up[NDF * Nst * 2];
      mult_wilson_zp2(v2v, &u[VLEN * NDF * site], &buf_zp[ibf]);
      ++opr_any;
    }

    if ((iz == 0) && (do_comm[2] > 0)) {
      int ibf = VLEN * NVC * ND2 * (ixy + Nxy * it);
      mult_wilson_zm2(v2v, &buf_zm[ibf]);
      ++opr_any;
    }

    if ((it == Nt - 1) && (do_comm[3] > 0)) {
      real_t *u = &up[NDF * Nst * 3];
      mult_wilson_tp2_dirac(v2v, &u[VLEN * NDF * site],
                            &buf_tp[VLEN * NVC * ND2 * ixyz]);
      ++opr_any;
    }

    if ((it == 0) && (do_comm[3] > 0)) {
      mult_wilson_tm2_dirac(v2v, &buf_tm[VLEN * NVC * ND2 * ixyz]);
      ++opr_any;
    }

    if (opr_any > 0) {
      mult_wilson_aypx_save(&v2[VLEN * NVCD * site],
                            -kappa, v2v, &v2[VLEN * NVCD * site]);
    }
  }
}


void BridgeQXS::mult_wilson_gm5_dirac(real_t *v2, real_t *v1,
                                      int *Nsize)
{
  int Nxv  = Nsize[0];
  int Nyv  = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Nyv * Nz * Nt;

  svbool_t pg = set_predicate();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstv);

  for (int site = is; site < ns; ++site) {
    real_t *vv1 = v1 + VLEN * 2 * ND * NC * site;
    real_t *vv2 = v2 + VLEN * 2 * ND * NC * site;
    for (int ic = 0; ic < NC; ++ic) {
      mult_gm5_dirac_vec(pg, &vv2[VLEN * 2 * ND * ic],
                         &vv1[VLEN * 2 * ND * ic]);
    }
  }
}


#endif
//============================================================END=====
