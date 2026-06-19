/*!
      @file    mult_Staggered_eo_qxs-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#ifndef MULT_STAGGERED_EO_QXS_INCLUDED
#define MULT_STAGGERED_EO_QXS_INCLUDED

#include "mult_common_th-inc.h"

//====================================================================
void BridgeQXS::mult_staggered_eo_bulk(real_t *v2, real_t *up,
                                       real_t *v1, real_t *xp,
                                       real_t mq, int jd,
                                       int *Nsize, int *do_comm,
                                       int *Leo, int ieo, int iflag)
{
  int Nx2v = Nsize[0];
  int Nyv  = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];

  int    Nst2v = Nx2v * Nyv * Nz * Nt;
  int    Nst2  = Nst2v * VLEN;
  real_t fac   = 0.5 * real_t(jd);

  svbool_t pg = set_predicate();
  svbool_t pg1e_xp, pg2e_xp, pg3e_xp, pg1e_xm, pg2e_xm, pg3e_xm;
  svbool_t pg1o_xp, pg2o_xp, pg3o_xp, pg1o_xm, pg2o_xm, pg3o_xm;
  svbool_t pg1_yp, pg2_yp, pg1_ym, pg2_ym;
  set_predicate_xp_eo(pg1e_xp, pg2e_xp, pg3e_xp, 0);
  set_predicate_xp_eo(pg1o_xp, pg2o_xp, pg3o_xp, 1);
  set_predicate_xm_eo(pg1e_xm, pg2e_xm, pg3e_xm, 0);
  set_predicate_xm_eo(pg1o_xm, pg2o_xm, pg3o_xm, 1);
  set_predicate_yp(pg1_yp, pg2_yp);
  set_predicate_ym(pg1_ym, pg2_ym);

  int Nxy  = Nx2v * Nyv;
  int Nxyz = Nx2v * Nyv * Nz;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nst2v);

  for (int site = is; site < ns; ++site) {
    int ix   = site % Nx2v;
    int iyzt = site / Nx2v;
    int iy   = iyzt % Nyv;
    int izt  = site / Nxy;
    int iz   = izt % Nz;
    int it   = izt / Nz;
    int ixy  = ix + Nx2v * iy;
    int ixyz = ixy + Nxy * iz;
    int jeo  = (ieo + Leo[VLENY * iyzt]) % 2;

    svreal_t vt0, vt1, vt2, vt3, vt4, vt5;
    clear_vec(pg, vt0);
    clear_vec(pg, vt1);
    clear_vec(pg, vt2);
    clear_vec(pg, vt3);
    clear_vec(pg, vt4);
    clear_vec(pg, vt5);

    if ((ix < Nx2v - 1) || (do_comm[0] == 0)) {
      int nei = ix + 1 + Nx2v * iyzt;
      if (ix == Nx2v - 1) nei = 0 + Nx2v * iyzt;
      real_t *u = &up[NDF * Nst2 * (ieo + 2 * 0)];
      if (jeo == 0) {
        mult_staggered_eo_xp(pg, pg1e_xp, pg2e_xp, pg3e_xp,
                             vt0, vt1, vt2, vt3, vt4, vt5,
                             &u[VLEN * NDF * site],
                             &v1[VLEN * NVC * site], &v1[VLEN * NVC * nei]);
      } else {
        mult_staggered_eo_xp(pg, pg1o_xp, pg2o_xp, pg3o_xp,
                             vt0, vt1, vt2, vt3, vt4, vt5,
                             &u[VLEN * NDF * site],
                             &v1[VLEN * NVC * site], &v1[VLEN * NVC * nei]);
      }
    }

    if ((ix > 0) || (do_comm[0] == 0)) {
      real_t *u  = &up[NDF * Nst2 * (1 - ieo + 2 * 0)];
      int    nei = ix - 1 + Nx2v * iyzt;
      if (ix == 0) nei = Nx2v - 1 + Nx2v * iyzt;

      svreal_t wt0, wt1, wt2, wt3, wt4, wt5;
      if (jeo == 0) {
        mult_staggered_eo_xm(pg, pg1e_xm, pg2e_xm, pg3e_xm,
                             vt0, vt1, vt2, vt3, vt4, vt5,
                             &u[VLEN * NDF * site], &u[VLEN * NDF * nei],
                             &v1[VLEN * NVC * site], &v1[VLEN * NVC * nei]);
      } else {
        mult_staggered_eo_xm(pg, pg1o_xm, pg2o_xm, pg3o_xm,
                             vt0, vt1, vt2, vt3, vt4, vt5,
                             &u[VLEN * NDF * site], &u[VLEN * NDF * nei],
                             &v1[VLEN * NVC * site], &v1[VLEN * NVC * nei]);
      }
    }

    if ((iy < Nyv - 1) || (do_comm[1] == 0)) {
      int    iy2 = (iy + 1) % Nyv;
      int    nei = ix + Nx2v * (iy2 + Nyv * izt);
      real_t *u  = &up[NDF * Nst2 * (ieo + 2 * 1)];
      mult_staggered_yp(pg, pg1_yp, pg2_yp,
                        vt0, vt1, vt2, vt3, vt4, vt5,
                        &u[VLEN * NDF * site],
                        &v1[VLEN * NVC * site], &v1[VLEN * NVC * nei]);
    }

    if ((iy > 0) || (do_comm[1] == 0)) {
      int    iy2 = (iy - 1 + Nyv) % Nyv;
      int    nei = ix + Nx2v * (iy2 + Nyv * izt);
      real_t *u  = &up[NDF * Nst2 * (1 - ieo + 2 * 1)];
      mult_staggered_ym(pg, pg1_ym, pg2_ym,
                        vt0, vt1, vt2, vt3, vt4, vt5,
                        &u[VLEN * NDF * site], &u[VLEN * NDF * nei],
                        &v1[VLEN * NVC * site], &v1[VLEN * NVC * nei]);
    }

    if ((iz < Nz - 1) || (do_comm[2] == 0)) {
      int iz2 = (iz + 1) % Nz;
      int nei = ixy + Nxy * (iz2 + Nz * it);
      mult_staggered_up(pg, vt0, vt1, vt2, vt3, vt4, vt5,
                        &up[VLEN * NDF * (site + Nst2v * (ieo + 2 * 2))],
                        &v1[VLEN * NVC * nei]);
    }

    if ((iz > 0) || (do_comm[2] == 0)) {
      int iz2 = (iz - 1 + Nz) % Nz;
      int nei = ixy + Nxy * (iz2 + Nz * it);
      mult_staggered_dn(pg, vt0, vt1, vt2, vt3, vt4, vt5,
                        &up[VLEN * NDF * (nei + Nst2v * (1 - ieo + 2 * 2))],
                        &v1[VLEN * NVC * nei]);
    }

    if ((it < Nt - 1) || (do_comm[3] == 0)) {
      int it2 = (it + 1) % Nt;
      int nei = ixyz + Nxyz * it2;
      mult_staggered_up(pg, vt0, vt1, vt2, vt3, vt4, vt5,
                        &up[VLEN * NDF * (site + Nst2v * (ieo + 2 * 3))],
                        &v1[VLEN * NVC * nei]);
    }

    if ((it > 0) || (do_comm[3] == 0)) {
      int it2 = (it - 1 + Nt) % Nt;
      int nei = ixyz + Nxyz * it2;
      mult_staggered_dn(pg, vt0, vt1, vt2, vt3, vt4, vt5,
                        &up[VLEN * NDF * (nei + Nst2v * (1 - ieo + 2 * 3))],
                        &v1[VLEN * NVC * nei]);
    }

    scal_vec(pg, vt0, fac);
    scal_vec(pg, vt1, fac);
    scal_vec(pg, vt2, fac);
    scal_vec(pg, vt3, fac);
    scal_vec(pg, vt4, fac);
    scal_vec(pg, vt5, fac);

    if (iflag == 0) {
      save_vec(pg, &v2[VLEN * (0 + NVC * site)], vt0);
      save_vec(pg, &v2[VLEN * (1 + NVC * site)], vt1);
      save_vec(pg, &v2[VLEN * (2 + NVC * site)], vt2);
      save_vec(pg, &v2[VLEN * (3 + NVC * site)], vt3);
      save_vec(pg, &v2[VLEN * (4 + NVC * site)], vt4);
      save_vec(pg, &v2[VLEN * (5 + NVC * site)], vt5);
    } else {
      svreal_t wt0, wt1, wt2, wt3, wt4, wt5;
      load_vec(pg, wt0, &xp[VLEN * (0 + NVC * site)]);
      load_vec(pg, wt1, &xp[VLEN * (1 + NVC * site)]);
      load_vec(pg, wt2, &xp[VLEN * (2 + NVC * site)]);
      load_vec(pg, wt3, &xp[VLEN * (3 + NVC * site)]);
      load_vec(pg, wt4, &xp[VLEN * (4 + NVC * site)]);
      load_vec(pg, wt5, &xp[VLEN * (5 + NVC * site)]);

      real_t fac2 = -1.0 / (mq * mq);
      axpy_vec(pg, wt0, fac2, vt0);
      axpy_vec(pg, wt1, fac2, vt1);
      axpy_vec(pg, wt2, fac2, vt2);
      axpy_vec(pg, wt3, fac2, vt3);
      axpy_vec(pg, wt4, fac2, vt4);
      axpy_vec(pg, wt5, fac2, vt5);

      save_vec(pg, &v2[VLEN * (0 + NVC * site)], wt0);
      save_vec(pg, &v2[VLEN * (1 + NVC * site)], wt1);
      save_vec(pg, &v2[VLEN * (2 + NVC * site)], wt2);
      save_vec(pg, &v2[VLEN * (3 + NVC * site)], wt3);
      save_vec(pg, &v2[VLEN * (4 + NVC * site)], wt4);
      save_vec(pg, &v2[VLEN * (5 + NVC * site)], wt5);
    }
  }
}


//====================================================================
void BridgeQXS::mult_staggered_eo_1(real_t *buf_xp, real_t *buf_xm,
                                    real_t *buf_yp, real_t *buf_ym,
                                    real_t *buf_zp, real_t *buf_zm,
                                    real_t *buf_tp, real_t *buf_tm,
                                    real_t *up, real_t *v1,
                                    int *Nsize, int *do_comm,
                                    int *Leo, int ieo)
{
  int Nx2v  = Nsize[0];
  int Nyv   = Nsize[1];
  int Nz    = Nsize[2];
  int Nt    = Nsize[3];
  int Nst2v = Nx2v * Nyv * Nz * Nt;
  int Nst2  = Nst2v * VLEN;

  int Nxy  = Nx2v * Nyv;
  int Nxyz = Nx2v * Nyv * Nz;

  svbool_t pg = set_predicate();
  svbool_t pg1e_xp, pg2e_xp, pg3e_xp, pg1e_xm, pg2e_xm, pg3e_xm;
  svbool_t pg1o_xp, pg2o_xp, pg3o_xp, pg1o_xm, pg2o_xm, pg3o_xm;
  set_predicate_xp_eo(pg1e_xp, pg2e_xp, pg3e_xp, 0);
  set_predicate_xp_eo(pg1o_xp, pg2o_xp, pg3o_xp, 1);
  set_predicate_xm_eo(pg1e_xm, pg2e_xm, pg3e_xm, 0);
  set_predicate_xm_eo(pg1o_xm, pg2o_xm, pg3o_xm, 1);
  svbool_t pg1_yp, pg2_yp, pg1_ym, pg2_ym;
  set_predicate_yp(pg1_yp, pg2_yp);
  set_predicate_ym(pg1_ym, pg2_ym);
  svint_t svidx_xp, svidx_xm;
  set_index_xp_eo(svidx_xp);
  set_index_xm_eo(svidx_xm);

  int Nskipx = (VLENY + 1) / 2;

  if (do_comm[0] > 0) {
    int    Nyzt = Nyv * Nz * Nt;
    real_t *u   = &up[NDF * Nst2 * (1 - ieo + 2 * 0)];

    int ith, nth, isx, nsx;
    set_threadtask(ith, nth, isx, nsx, Nyzt);

    for (int iyzt = isx; iyzt < nsx; ++iyzt) {
      {
        int    ix   = 0;
        int    site = ix + Nx2v * iyzt;
        real_t *buf = &buf_xp[Nskipx * NVC * iyzt];
        int    jeo  = (ieo + Leo[VLENY * iyzt]) % 2;
        if (jeo == 0) {
#if VLENY > 1
          set_index_xm_eo(svidx_xm);
          for (int ivc = 0; ivc < NVC; ++ivc) {
            svreal_t wt;
            load_vec(pg2o_xm, wt, &v1[VLEN * (ivc + NVC * site)]);
            svst1_scatter_index(pg2o_xm, &buf[Nskipx * ivc], svidx_xm, wt);
          }
#endif
        } else {
#if VLENY == 1
          buf = &buf_xp[Nskipx * NVC * (iyzt / 2)];
#endif
          set_index_xm_eo(svidx_xm);
          for (int ivc = 0; ivc < NVC; ++ivc) {
            svreal_t wt;
            load_vec(pg2e_xm, wt, &v1[VLEN * (ivc + NVC * site)]);
            svst1_scatter_index(pg2e_xm, &buf[Nskipx * ivc], svidx_xm, wt);
          }
        }
      }
      {
        int    ix   = Nx2v - 1;
        int    site = ix + Nx2v * iyzt;
        real_t *buf = &buf_xm[Nskipx * NVC * iyzt];
        int    jeo  = (ieo + Leo[VLENY * iyzt]) % 2;
        if (jeo == 0) {
#if VLENY == 1
          buf = &buf_xm[Nskipx * NVC * (iyzt / 2)];
#endif
          set_index_xp_eo(svidx_xp);
          mult_staggered_eo_xm1(pg2o_xp, svidx_xp, buf,
                                &u[VLEN * NDF * site], &v1[VLEN * NVC * site]);
        } else {
#if VLENY > 1
          set_index_xp_eo(svidx_xp);
          mult_staggered_eo_xm1(pg2e_xp, svidx_xp, buf,
                                &u[VLEN * NDF * site], &v1[VLEN * NVC * site]);
#endif
        }
      }
    }
  }

  if (do_comm[1] > 0) {
    int    Nxzt = Nx2v * Nz * Nt;
    real_t *u   = &up[NDF * Nst2 * (1 - ieo + 2 * 1)];

    int ith, nth, isy, nsy;
    set_threadtask(ith, nth, isy, nsy, Nxzt);

    for (int ixzt = isy; ixzt < nsy; ++ixzt) {
      int ix  = ixzt % Nx2v;
      int izt = ixzt / Nx2v;
      {
        int    iy   = 0;
        int    site = ix + Nx2v * (iy + Nyv * izt);
        real_t *buf = &buf_yp[VLENX * NVC * ixzt];
        for (int ic = 0; ic < NC; ++ic) {
          svreal_t wtr, wti;
          load_vec(pg2_ym, wtr, &v1[VLEN * (2 * ic + NVC * site)]);
          load_vec(pg2_ym, wti, &v1[VLEN * (2 * ic + 1 + NVC * site)]);
          //svst1(pg2_ym, &buf[VLENX*(2*ic)  ], wtr);
          //svst1(pg2_ym, &buf[VLENX*(2*ic+1)], wti);
          save_vec(pg2_ym, &buf[VLENX * (2 * ic)], wtr);
          save_vec(pg2_ym, &buf[VLENX * (2 * ic + 1)], wti);
        }
      }
      {
        int      iy = Nyv - 1;
        int      site = ix + Nx2v * (iy + Nyv * izt);
        real_t   *buf = &buf_ym[VLENX * NVC * ixzt];
        svreal_t wt0, wt1, wt2, wt3, wt4, wt5;
        load_vec(pg2_yp, wt0, &v1[VLEN * (0 + NVC * site)]);
        load_vec(pg2_yp, wt1, &v1[VLEN * (1 + NVC * site)]);
        load_vec(pg2_yp, wt2, &v1[VLEN * (2 + NVC * site)]);
        load_vec(pg2_yp, wt3, &v1[VLEN * (3 + NVC * site)]);
        load_vec(pg2_yp, wt4, &v1[VLEN * (4 + NVC * site)]);
        load_vec(pg2_yp, wt5, &v1[VLEN * (5 + NVC * site)]);

        int offset = -VLENX * (VLENY - 1);
        for (int ic = 0; ic < NC; ++ic) {
          svreal_t ut0, ut1, ut2, ut3, ut4, ut5;
          svreal_t xtr, xti;
          load_udag(pg2_yp, ut0, ut1, ut2, ut3, ut4, ut5,
                    &u[VLEN * (NVC * ic + NDF * site)]);
          mult_udv(pg2_yp, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
                   wt0, wt1, wt2, wt3, wt4, wt5);
          //svst1(pg2_yp, &buf[offset + VLENX*(2*ic)  ], xtr);
          //svst1(pg2_yp, &buf[offset + VLENX*(2*ic+1)], xti);
          save_vec(pg2_yp, &buf[offset + VLENX * (2 * ic)], xtr);
          save_vec(pg2_yp, &buf[offset + VLENX * (2 * ic + 1)], xti);
        }
      }
    }
  }

  if (do_comm[2] > 0) {
    int    Nxyt = Nx2v * Nyv * Nt;
    real_t *u   = &up[NDF * Nst2 * (1 - ieo + 2 * 2)];

    int ith, nth, isz, nsz;
    set_threadtask(ith, nth, isz, nsz, Nxyt);

    for (int ixyt = isz; ixyt < nsz; ++ixyt) {
      int ixy = ixyt % Nxy;
      int it  = ixyt / Nxy;
      {
        int    iz   = 0;
        int    site = ixy + Nxy * (iz + Nz * it);
        real_t *buf = &buf_zp[VLEN * NVC * ixyt];
        for (int ic = 0; ic < NC; ++ic) {
          svreal_t wtr, wti;
          load_vec(pg, wtr, &v1[VLEN * (2 * ic + NVC * site)]);
          load_vec(pg, wti, &v1[VLEN * (2 * ic + 1 + NVC * site)]);
          //svst1(pg, &buf[VLEN * (2*ic)  ], wtr);
          //svst1(pg, &buf[VLEN * (2*ic+1)], wti);
          save_vec(pg, &buf[VLEN * (2 * ic)], wtr);
          save_vec(pg, &buf[VLEN * (2 * ic + 1)], wti);
        }
      }
      {
        int iz   = Nz - 1;
        int site = ixy + Nxy * (iz + Nz * it);
        mult_staggered_dn1(pg, &buf_zm[VLEN * NVC * ixyt],
                           &u[VLEN * NDF * site], &v1[VLEN * NVC * site]);
      }
    }
  }

  if (do_comm[3] > 0) {
    int    Nxyz = Nx2v * Nyv * Nz;
    real_t *u   = &up[NDF * Nst2 * (1 - ieo + 2 * 3)];

    int ith, nth, ist, nst;
    set_threadtask(ith, nth, ist, nst, Nxyz);

    for (int ixyz = ist; ixyz < nst; ++ixyz) {
      {
        int    it   = 0;
        int    site = ixyz + Nxyz * it;
        real_t *buf = &buf_tp[VLEN * NVC * ixyz];
        for (int ic = 0; ic < NC; ++ic) {
          svreal_t wtr, wti;
          load_vec(pg, wtr, &v1[VLEN * (2 * ic + NVC * site)]);
          load_vec(pg, wti, &v1[VLEN * (2 * ic + 1 + NVC * site)]);
          //svst1(pg, &buf[VLEN * (2*ic)  ], wtr);
          //svst1(pg, &buf[VLEN * (2*ic+1)], wti);
          save_vec(pg, &buf[VLEN * (2 * ic)], wtr);
          save_vec(pg, &buf[VLEN * (2 * ic + 1)], wti);
        }
      }
      {
        int it   = Nt - 1;
        int site = ixyz + Nxyz * it;
        mult_staggered_dn1(pg, &buf_tm[VLEN * NVC * ixyz],
                           &u[VLEN * NDF * site], &v1[VLEN * NVC * site]);
      }
    }
  }
}


//====================================================================
void BridgeQXS::mult_staggered_eo_2(real_t *v2, real_t *up, real_t *v1,
                                    real_t *buf_xp, real_t *buf_xm,
                                    real_t *buf_yp, real_t *buf_ym,
                                    real_t *buf_zp, real_t *buf_zm,
                                    real_t *buf_tp, real_t *buf_tm,
                                    real_t mq, int *Nsize, int *do_comm,
                                    int *Leo, int ieo, int iflag)
{
  int Nx2v  = Nsize[0];
  int Nyv   = Nsize[1];
  int Nz    = Nsize[2];
  int Nt    = Nsize[3];
  int Nst2v = Nx2v * Nyv * Nz * Nt;
  int Nst2  = Nst2v * VLEN;

  real_t fac = 0.5;
  if (iflag != 0) fac = -0.5 / (mq * mq);

  int Nxy  = Nx2v * Nyv;
  int Nxyz = Nx2v * Nyv * Nz;

  svbool_t pg = set_predicate();
  svbool_t pg1e_xp, pg2e_xp, pg3e_xp, pg1e_xm, pg2e_xm, pg3e_xm;
  svbool_t pg1o_xp, pg2o_xp, pg3o_xp, pg1o_xm, pg2o_xm, pg3o_xm;
  set_predicate_xp_eo(pg1e_xp, pg2e_xp, pg3e_xp, 0);
  set_predicate_xp_eo(pg1o_xp, pg2o_xp, pg3o_xp, 1);
  set_predicate_xm_eo(pg1e_xm, pg2e_xm, pg3e_xm, 0);
  set_predicate_xm_eo(pg1o_xm, pg2o_xm, pg3o_xm, 1);
  svbool_t pg1_yp, pg2_yp, pg1_ym, pg2_ym;
  set_predicate_yp(pg1_yp, pg2_yp);
  set_predicate_ym(pg1_ym, pg2_ym);
  svint_t svidx_xp, svidx_xm;
  set_index_xp_eo(svidx_xp);
  set_index_xm_eo(svidx_xm);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nst2v);

  for (int site = is; site < ns; ++site) {
    int ix   = site % Nx2v;
    int iyzt = site / Nx2v;
    int iy   = iyzt % Nyv;
    int izt  = site / Nxy;
    int iz   = izt % Nz;
    int it   = izt / Nz;
    int ixy  = ix + Nx2v * iy;
    int ixyz = ixy + Nxy * iz;
    int jeo  = (ieo + Leo[VLENY * iyzt]) % 2;

    int Nskipx = (VLENY + 1) / 2;

    svreal_t vt0, vt1, vt2, vt3, vt4, vt5;
    clear_vec(pg, vt0);
    clear_vec(pg, vt1);
    clear_vec(pg, vt2);
    clear_vec(pg, vt3);
    clear_vec(pg, vt4);
    clear_vec(pg, vt5);

    int opr_any = 0;

    if ((ix == Nx2v - 1) && (do_comm[0] > 0)) {
      real_t *u   = &up[NDF * Nst2 * (ieo + 2 * 0)];
      real_t *buf = &buf_xp[Nskipx * NVC * iyzt];
#if VLENY == 1
      buf = &buf_xp[Nskipx * NVC * (iyzt / 2)];
#endif
      svreal_t wt0, wt1, wt2, wt3, wt4, wt5;
      set_index_xp_eo(svidx_xp);
      if (jeo == 0) {
        mult_staggered_eo_xp2(pg1e_xp, pg2e_xp, pg3e_xp, svidx_xp,
                              vt0, vt1, vt2, vt3, vt4, vt5,
                              &u[VLEN * NDF * site],
                              &v1[VLEN * NVC * site], buf);
      } else {
        mult_staggered_eo_xp2(pg1o_xp, pg2o_xp, pg3o_xp, svidx_xp,
                              vt0, vt1, vt2, vt3, vt4, vt5,
                              &u[VLEN * NDF * site],
                              &v1[VLEN * NVC * site], buf);
      }
      ++opr_any;
    }

    if ((ix == 0) && (do_comm[0] > 0)) {
      real_t *u   = &up[NDF * Nst2 * (1 - ieo + 2 * 0)];
      real_t *buf = &buf_xm[Nskipx * NVC * iyzt];
#if VLENY == 1
      buf = &buf_xm[Nskipx * NVC * (iyzt / 2)];
#endif
      set_index_xm_eo(svidx_xm);
      if (jeo == 0) {
        mult_staggered_eo_xm2(pg1e_xm, pg2e_xm, pg3e_xm, svidx_xm,
                              vt0, vt1, vt2, vt3, vt4, vt5,
                              &u[VLEN * NDF * site],
                              &v1[VLEN * NVC * site], buf);
      } else {
        mult_staggered_eo_xm2(pg1o_xm, pg2o_xm, pg3o_xm, svidx_xm,
                              vt0, vt1, vt2, vt3, vt4, vt5,
                              &u[VLEN * NDF * site],
                              &v1[VLEN * NVC * site], buf);
      }
      ++opr_any;
    }

    if ((iy == Nyv - 1) && (do_comm[1] > 0)) {
      real_t *u   = &up[NDF * Nst2 * (ieo + 2 * 1)];
      int    ixzt = ix + Nx2v * izt;
      real_t *buf = &buf_yp[VLENX * NVC * ixzt];

      svreal_t wt0, wt1, wt2, wt3, wt4, wt5;
      load_vec(pg1_yp, wt0, &v1[VLEN * (0 + NVC * site) + VLENX]);
      load_vec(pg1_yp, wt1, &v1[VLEN * (1 + NVC * site) + VLENX]);
      load_vec(pg1_yp, wt2, &v1[VLEN * (2 + NVC * site) + VLENX]);
      load_vec(pg1_yp, wt3, &v1[VLEN * (3 + NVC * site) + VLENX]);
      load_vec(pg1_yp, wt4, &v1[VLEN * (4 + NVC * site) + VLENX]);
      load_vec(pg1_yp, wt5, &v1[VLEN * (5 + NVC * site) + VLENX]);

      int offset = -VLENX * (VLENY - 1);
      load_add(pg2_yp, wt0, &buf[offset + VLENX * 0]);
      load_add(pg2_yp, wt1, &buf[offset + VLENX * 1]);
      load_add(pg2_yp, wt2, &buf[offset + VLENX * 2]);
      load_add(pg2_yp, wt3, &buf[offset + VLENX * 3]);
      load_add(pg2_yp, wt4, &buf[offset + VLENX * 4]);
      load_add(pg2_yp, wt5, &buf[offset + VLENX * 5]);

      svreal_t ut0, ut1, ut2, ut3, ut4, ut5;
      svreal_t xtr, xti;
      load_u(pg, ut0, ut1, ut2, ut3, ut4, ut5,
             &u[VLEN * (0 + NDF * site)]);
      mult_uv(pg, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
              wt0, wt1, wt2, wt3, wt4, wt5);
      add_vec(pg, vt0, xtr);
      add_vec(pg, vt1, xti);
      load_u(pg, ut0, ut1, ut2, ut3, ut4, ut5,
             &u[VLEN * (2 + NDF * site)]);
      mult_uv(pg, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
              wt0, wt1, wt2, wt3, wt4, wt5);
      add_vec(pg, vt2, xtr);
      add_vec(pg, vt3, xti);
      load_u(pg, ut0, ut1, ut2, ut3, ut4, ut5,
             &u[VLEN * (4 + NDF * site)]);
      mult_uv(pg, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
              wt0, wt1, wt2, wt3, wt4, wt5);
      add_vec(pg, vt4, xtr);
      add_vec(pg, vt5, xti);

      ++opr_any;
    }

    if ((iy == 0) && (do_comm[1] > 0)) {
      real_t *u   = &up[NDF * Nst2 * (1 - ieo + 2 * 1)];
      int    ixzt = ix + Nx2v * izt;
      real_t *buf = &buf_ym[VLENX * NVC * ixzt];

      svreal_t wt0, wt1, wt2, wt3, wt4, wt5;
      load_vec(pg1_ym, wt0, &v1[VLEN * (0 + NVC * site) - VLENX]);
      load_vec(pg1_ym, wt1, &v1[VLEN * (1 + NVC * site) - VLENX]);
      load_vec(pg1_ym, wt2, &v1[VLEN * (2 + NVC * site) - VLENX]);
      load_vec(pg1_ym, wt3, &v1[VLEN * (3 + NVC * site) - VLENX]);
      load_vec(pg1_ym, wt4, &v1[VLEN * (4 + NVC * site) - VLENX]);
      load_vec(pg1_ym, wt5, &v1[VLEN * (5 + NVC * site) - VLENX]);

      svreal_t ut0, ut1, ut2, ut3, ut4, ut5;
      svreal_t xtr, xti;
      load_udag(pg1_ym, ut0, ut1, ut2, ut3, ut4, ut5,
                &u[VLEN * (0 + NDF * site) - VLENX]);
      mult_udv(pg1_ym, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
               wt0, wt1, wt2, wt3, wt4, wt5);
      load_add(pg2_ym, xtr, &buf[VLENX * 0]);
      load_add(pg2_ym, xti, &buf[VLENX * 1]);
      sub_vec(pg, vt0, xtr);
      sub_vec(pg, vt1, xti);

      load_udag(pg1_ym, ut0, ut1, ut2, ut3, ut4, ut5,
                &u[VLEN * (NVC + NDF * site) - VLENX]);
      mult_udv(pg1_ym, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
               wt0, wt1, wt2, wt3, wt4, wt5);
      load_add(pg2_ym, xtr, &buf[VLENX * 2]);
      load_add(pg2_ym, xti, &buf[VLENX * 3]);
      sub_vec(pg, vt2, xtr);
      sub_vec(pg, vt3, xti);

      load_udag(pg1_ym, ut0, ut1, ut2, ut3, ut4, ut5,
                &u[VLEN * (2 * NVC + NDF * site) - VLENX]);
      mult_udv(pg1_ym, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
               wt0, wt1, wt2, wt3, wt4, wt5);
      load_add(pg2_ym, xtr, &buf[VLENX * 4]);
      load_add(pg2_ym, xti, &buf[VLENX * 5]);
      sub_vec(pg, vt4, xtr);
      sub_vec(pg, vt5, xti);

      ++opr_any;
    }

    if ((iz == Nz - 1) && (do_comm[2] > 0)) {
      int ixyt = ixy + Nxy * it;
      mult_staggered_up(pg, vt0, vt1, vt2, vt3, vt4, vt5,
                        &up[VLEN * NDF * (site + Nst2v * (ieo + 2 * 2))],
                        &buf_zp[VLEN * NVC * ixyt]);
      ++opr_any;
    }

    if ((iz == 0) && (do_comm[2] > 0)) {
      int    ixyt = ixy + Nxy * it;
      real_t *buf = &buf_zm[VLEN * NVC * ixyt];

      svreal_t wt0, wt1, wt2, wt3, wt4, wt5;
      load_vec(pg, wt0, &buf[VLEN * 0]);
      load_vec(pg, wt1, &buf[VLEN * 1]);
      load_vec(pg, wt2, &buf[VLEN * 2]);
      load_vec(pg, wt3, &buf[VLEN * 3]);
      load_vec(pg, wt4, &buf[VLEN * 4]);
      load_vec(pg, wt5, &buf[VLEN * 5]);

      sub_vec(pg, vt0, wt0);
      sub_vec(pg, vt1, wt1);
      sub_vec(pg, vt2, wt2);
      sub_vec(pg, vt3, wt3);
      sub_vec(pg, vt4, wt4);
      sub_vec(pg, vt5, wt5);

      ++opr_any;
    }

    if ((it == Nt - 1) && (do_comm[3] > 0)) {
      mult_staggered_up(pg, vt0, vt1, vt2, vt3, vt4, vt5,
                        &up[VLEN * NDF * (site + Nst2v * (ieo + 2 * 3))],
                        &buf_tp[VLEN * NVC * ixyz]);
      ++opr_any;
    }

    if ((it == 0) && (do_comm[3] > 0)) {
      real_t *buf = &buf_tm[VLEN * NVC * ixyz];

      svreal_t wt0, wt1, wt2, wt3, wt4, wt5;
      load_vec(pg, wt0, &buf[VLEN * 0]);
      load_vec(pg, wt1, &buf[VLEN * 1]);
      load_vec(pg, wt2, &buf[VLEN * 2]);
      load_vec(pg, wt3, &buf[VLEN * 3]);
      load_vec(pg, wt4, &buf[VLEN * 4]);
      load_vec(pg, wt5, &buf[VLEN * 5]);

      sub_vec(pg, vt0, wt0);
      sub_vec(pg, vt1, wt1);
      sub_vec(pg, vt2, wt2);
      sub_vec(pg, vt3, wt3);
      sub_vec(pg, vt4, wt4);
      sub_vec(pg, vt5, wt5);

      ++opr_any;
    }

    if (opr_any > 0) {
      svreal_t wt0, wt1, wt2, wt3, wt4, wt5;
      load_vec(pg, wt0, &v2[VLEN * (0 + NVC * site)]);
      load_vec(pg, wt1, &v2[VLEN * (1 + NVC * site)]);
      load_vec(pg, wt2, &v2[VLEN * (2 + NVC * site)]);
      load_vec(pg, wt3, &v2[VLEN * (3 + NVC * site)]);
      load_vec(pg, wt4, &v2[VLEN * (4 + NVC * site)]);
      load_vec(pg, wt5, &v2[VLEN * (5 + NVC * site)]);

      scal_vec(pg, vt0, fac);
      scal_vec(pg, vt1, fac);
      scal_vec(pg, vt2, fac);
      scal_vec(pg, vt3, fac);
      scal_vec(pg, vt4, fac);
      scal_vec(pg, vt5, fac);

      add_vec(pg, vt0, wt0);
      add_vec(pg, vt1, wt1);
      add_vec(pg, vt2, wt2);
      add_vec(pg, vt3, wt3);
      add_vec(pg, vt4, wt4);
      add_vec(pg, vt5, wt5);

      //svst1(pg, &v2[VLEN * (0 + NVC * site)], vt0);
      //svst1(pg, &v2[VLEN * (1 + NVC * site)], vt1);
      //svst1(pg, &v2[VLEN * (2 + NVC * site)], vt2);
      //svst1(pg, &v2[VLEN * (3 + NVC * site)], vt3);
      //svst1(pg, &v2[VLEN * (4 + NVC * site)], vt4);
      //svst1(pg, &v2[VLEN * (5 + NVC * site)], vt5);

      save_vec(pg, &v2[VLEN * (0 + NVC * site)], vt0);
      save_vec(pg, &v2[VLEN * (1 + NVC * site)], vt1);
      save_vec(pg, &v2[VLEN * (2 + NVC * site)], vt2);
      save_vec(pg, &v2[VLEN * (3 + NVC * site)], vt3);
      save_vec(pg, &v2[VLEN * (4 + NVC * site)], vt4);
      save_vec(pg, &v2[VLEN * (5 + NVC * site)], vt5);
    }
  }
}


#endif
//============================================================END=====
