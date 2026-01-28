/*!
      @file    mult_Staggered_qxs-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#ifndef MULT_STAGGERED_QXS_INCLUDED
#define MULT_STAGGERED_QXS_INCLUDED

#include "mult_common_th-inc.h"

//====================================================================
void BridgeQXS::mult_staggered_phase(real_t *v, real_t *ph,
                                     int *Nsize, int Nin)
{
  int Nxv  = Nsize[0];
  int Nyv  = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Nyv * Nz * Nt;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstv);

  svbool_t pg = set_predicate();

  for (int site = is; site < ns; ++site) {
    svreal_t vph;
    load_vec(pg, vph, &ph[VLEN * site]);

    for (int in = 0; in < Nin; ++in) {
      svreal_t vt;
      load_vec(pg, vt, &v[VLEN * (in + Nin * site)]);
      scal_vec(pg, vt, vph);
      save_vec(pg, &v[VLEN * (in + Nin * site)], vt);
    }
  }
}


//====================================================================
void BridgeQXS::mult_staggered_clear(real_t *v, int *Nsize, int Nin)
{
  int Nxv  = Nsize[0];
  int Nyv  = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Nyv * Nz * Nt;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstv);

  svbool_t pg = set_predicate();
  svreal_t vz;
  clear_vec(pg, vz);

  for (int site = is; site < ns; ++site) {
    for (int in = 0; in < Nin; ++in) {
      //svst1(pg, &v[VLEN * (in + Nin * site)], vz);
      save_vec(pg, &v[VLEN * (in + Nin * site)], vz);
    }
  }
}


//====================================================================
void BridgeQXS::mult_staggered_axpby(real_t b, real_t *v,
                                     real_t a, real_t *w,
                                     int *Nsize, int Nin)
{
  int Nxv  = Nsize[0];
  int Nyv  = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Nyv * Nz * Nt;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstv);

  svbool_t pg = set_predicate();

  for (int site = is; site < ns; ++site) {
    for (int in = 0; in < Nin; ++in) {
      svreal_t vt, wt;
      load_vec(pg, wt, &w[VLEN * (in + Nin * site)]);
      load_vec(pg, vt, &v[VLEN * (in + Nin * site)]);
      scal_vec(pg, vt, b);
      axpy_vec(pg, vt, a, wt);
      //svst1(pg, &v[VLEN * (in + Nin * site)], vt);
      save_vec(pg, &v[VLEN * (in + Nin * site)], vt);
    }
  }
}


//====================================================================
void BridgeQXS::mult_staggered_mult_Gn(real_t *v, real_t *u,
                                       real_t *w, int *Nsize)
{
  int Nxv  = Nsize[0];
  int Nyv  = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Nyv * Nz * Nt;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstv);

  svbool_t pg = set_predicate();

  for (int site = is; site < ns; ++site) {
    svreal_t wt0, wt1, wt2, wt3, wt4, wt5;
    load_vec(pg, wt0, &w[VLEN * (0 + NVC * site)]);
    load_vec(pg, wt1, &w[VLEN * (1 + NVC * site)]);
    load_vec(pg, wt2, &w[VLEN * (2 + NVC * site)]);
    load_vec(pg, wt3, &w[VLEN * (3 + NVC * site)]);
    load_vec(pg, wt4, &w[VLEN * (4 + NVC * site)]);
    load_vec(pg, wt5, &w[VLEN * (5 + NVC * site)]);

    for (int ic = 0; ic < NC; ++ic) {
      svreal_t ut0, ut1, ut2, ut3, ut4, ut5;
      load_u(pg, ut0, ut1, ut2, ut3, ut4, ut5,
             &u[VLEN * (2 * ic + NDF * site)]);

      svreal_t xtr, xti;
      mult_uv(pg, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
              wt0, wt1, wt2, wt3, wt4, wt5);

      //svst1(pg, &v[VLEN * (2*ic   + NVC * site)], xtr);
      //svst1(pg, &v[VLEN * (2*ic+1 + NVC * site)], xti);
      save_vec(pg, &v[VLEN * (2 * ic + NVC * site)], xtr);
      save_vec(pg, &v[VLEN * (2 * ic + 1 + NVC * site)], xti);
    }
  }
}


//====================================================================
void BridgeQXS::mult_staggered_mult_Gd(real_t *v, real_t *u, real_t *w,
                                       int *Nsize)
{
  int Nxv  = Nsize[0];
  int Nyv  = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Nyv * Nz * Nt;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstv);

  svbool_t pg = set_predicate();

  for (int site = is; site < ns; ++site) {
    svreal_t wt0, wt1, wt2, wt3, wt4, wt5;
    load_vec(pg, wt0, &w[VLEN * (0 + NVC * site)]);
    load_vec(pg, wt1, &w[VLEN * (1 + NVC * site)]);
    load_vec(pg, wt2, &w[VLEN * (2 + NVC * site)]);
    load_vec(pg, wt3, &w[VLEN * (3 + NVC * site)]);
    load_vec(pg, wt4, &w[VLEN * (4 + NVC * site)]);
    load_vec(pg, wt5, &w[VLEN * (5 + NVC * site)]);

    for (int ic = 0; ic < NC; ++ic) {
      svreal_t ut0, ut1, ut2, ut3, ut4, ut5;
      load_udag(pg, ut0, ut1, ut2, ut3, ut4, ut5,
                &u[VLEN * (NVC * ic + NDF * site)]);

      svreal_t xtr, xti;
      mult_udv(pg, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
               wt0, wt1, wt2, wt3, wt4, wt5);

      // svst1(pg, &v[VLEN * (2*ic   + NVC * site)], xtr);
      // svst1(pg, &v[VLEN * (2*ic+1 + NVC * site)], xti);
      save_vec(pg, &v[VLEN * (2 * ic + NVC * site)], xtr);
      save_vec(pg, &v[VLEN * (2 * ic + 1 + NVC * site)], xti);
    }
  }
}


//====================================================================
void BridgeQXS::mult_staggered_bulk(real_t *v2, real_t *up, real_t *v1,
                                    real_t mq, int jd,
                                    int *Nsize, int *do_comm)
{
  int    Nxv  = Nsize[0];
  int    Nyv  = Nsize[1];
  int    Nz   = Nsize[2];
  int    Nt   = Nsize[3];
  int    Nstv = Nxv * Nyv * Nz * Nt;
  int    Nst  = Nstv * VLEN;
  real_t fac  = 0.5 * real_t(jd);

  svbool_t pg = set_predicate();
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

    svreal_t vt0, vt1, vt2, vt3, vt4, vt5;
    clear_vec(pg, vt0);
    clear_vec(pg, vt1);
    clear_vec(pg, vt2);
    clear_vec(pg, vt3);
    clear_vec(pg, vt4);
    clear_vec(pg, vt5);

    if ((ix < Nxv - 1) || (do_comm[0] == 0)) {
      real_t *u  = &up[NDF * Nst * 0];
      int    nei = ix + 1 + Nxv * iyzt;
      if (ix == Nxv - 1) nei = 0 + Nxv * iyzt;
      mult_staggered_xp(pg, pg1_xp, pg2_xp,
                        vt0, vt1, vt2, vt3, vt4, vt5,
                        &u[VLEN * NDF * site],
                        &v1[VLEN * NVC * site], &v1[VLEN * NVC * nei]);
    }

    if ((ix > 0) || (do_comm[0] == 0)) {
      real_t *u  = &up[NDF * Nst * 0];
      int    nei = ix - 1 + Nxv * iyzt;
      if (ix == 0) nei = Nxv - 1 + Nxv * iyzt;
      mult_staggered_xm(pg, pg1_xm, pg2_xm,
                        vt0, vt1, vt2, vt3, vt4, vt5,
                        &u[VLEN * NDF * site], &u[VLEN * NDF * nei],
                        &v1[VLEN * NVC * site], &v1[VLEN * NVC * nei]);
    }

    if ((iy < Nyv - 1) || (do_comm[1] == 0)) {
      int    iy2 = (iy + 1) % Nyv;
      int    nei = ix + Nxv * (iy2 + Nyv * izt);
      real_t *u  = &up[NDF * Nst * 1];
      mult_staggered_yp(pg, pg1_yp, pg2_yp,
                        vt0, vt1, vt2, vt3, vt4, vt5,
                        &u[VLEN * NDF * site],
                        &v1[VLEN * NVC * site], &v1[VLEN * NVC * nei]);
    }

    if ((iy > 0) || (do_comm[1] == 0)) {
      int    iy2 = (iy - 1 + Nyv) % Nyv;
      int    nei = ix + Nxv * (iy2 + Nyv * izt);
      real_t *u  = &up[NDF * Nst * 1];
      mult_staggered_ym(pg, pg1_ym, pg2_ym,
                        vt0, vt1, vt2, vt3, vt4, vt5,
                        &u[VLEN * NDF * site], &u[VLEN * NDF * nei],
                        &v1[VLEN * NVC * site], &v1[VLEN * NVC * nei]);
    }

    if ((iz < Nz - 1) || (do_comm[2] == 0)) {
      int iz2 = (iz + 1) % Nz;
      int nei = ixy + Nxy * (iz2 + Nz * it);
      mult_staggered_up(pg, vt0, vt1, vt2, vt3, vt4, vt5,
                        &up[VLEN * NDF * (site + Nstv * 2)],
                        &v1[VLEN * NVC * nei]);
    }

    if ((iz > 0) || (do_comm[2] == 0)) {
      int iz2 = (iz - 1 + Nz) % Nz;
      int nei = ixy + Nxy * (iz2 + Nz * it);
      mult_staggered_dn(pg, vt0, vt1, vt2, vt3, vt4, vt5,
                        &up[VLEN * NDF * (nei + Nstv * 2)],
                        &v1[VLEN * NVC * nei]);
    }

    if ((it < Nt - 1) || (do_comm[3] == 0)) {
      int it2 = (it + 1) % Nt;
      int nei = ixyz + Nxyz * it2;
      mult_staggered_up(pg, vt0, vt1, vt2, vt3, vt4, vt5,
                        &up[VLEN * NDF * (site + Nstv * 3)],
                        &v1[VLEN * NVC * nei]);
    }

    if ((it > 0) || (do_comm[3] == 0)) {
      int it2 = (it - 1 + Nt) % Nt;
      int nei = ixyz + Nxyz * it2;
      mult_staggered_dn(pg, vt0, vt1, vt2, vt3, vt4, vt5,
                        &up[VLEN * NDF * (nei + Nstv * 3)],
                        &v1[VLEN * NVC * nei]);
    }

    svreal_t wt0, wt1, wt2, wt3, wt4, wt5;
    load_vec(pg, wt0, &v1[VLEN * (0 + NVC * site)]);
    load_vec(pg, wt1, &v1[VLEN * (1 + NVC * site)]);
    load_vec(pg, wt2, &v1[VLEN * (2 + NVC * site)]);
    load_vec(pg, wt3, &v1[VLEN * (3 + NVC * site)]);
    load_vec(pg, wt4, &v1[VLEN * (4 + NVC * site)]);
    load_vec(pg, wt5, &v1[VLEN * (5 + NVC * site)]);

    scal_vec(pg, vt0, fac);
    scal_vec(pg, vt1, fac);
    scal_vec(pg, vt2, fac);
    scal_vec(pg, vt3, fac);
    scal_vec(pg, vt4, fac);
    scal_vec(pg, vt5, fac);

    axpy_vec(pg, vt0, mq, wt0);
    axpy_vec(pg, vt1, mq, wt1);
    axpy_vec(pg, vt2, mq, wt2);
    axpy_vec(pg, vt3, mq, wt3);
    axpy_vec(pg, vt4, mq, wt4);
    axpy_vec(pg, vt5, mq, wt5);

    save_vec(pg, &v2[VLEN * (0 + NVC * site)], vt0);
    save_vec(pg, &v2[VLEN * (1 + NVC * site)], vt1);
    save_vec(pg, &v2[VLEN * (2 + NVC * site)], vt2);
    save_vec(pg, &v2[VLEN * (3 + NVC * site)], vt3);
    save_vec(pg, &v2[VLEN * (4 + NVC * site)], vt4);
    save_vec(pg, &v2[VLEN * (5 + NVC * site)], vt5);
  }
}


//====================================================================
void BridgeQXS::mult_staggered_1(real_t *buf_xp, real_t *buf_xm,
                                 real_t *buf_yp, real_t *buf_ym,
                                 real_t *buf_zp, real_t *buf_zm,
                                 real_t *buf_tp, real_t *buf_tm,
                                 real_t *up, real_t *v1,
                                 int *Nsize, int *do_comm)
{
  int Nxv  = Nsize[0];
  int Nyv  = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Nyv * Nz * Nt;
  int Nst  = Nstv * VLEN;

  int Nxy  = Nxv * Nyv;
  int Nxyz = Nxv * Nyv * Nz;

  svbool_t pg = set_predicate();
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
    int ith, nth, isx, nsx;
    set_threadtask(ith, nth, isx, nsx, Nyzt);

    for (int iyzt = isx; iyzt < nsx; ++iyzt) {
      {
        int    ix   = 0;
        int    site = ix + Nxv * iyzt;
        real_t *buf = &buf_xp[VLENY * NVC * iyzt];

        set_index_xm(svidx_xm);
        for (int ivc = 0; ivc < NVC; ++ivc) {
          svreal_t wt;
          load_vec(pg2_xm, wt, &v1[VLEN * (ivc + NVC * site)]);
          //svst1_scatter_index(pg2_xm, &buf[VLENY*ivc], svidx_xm, wt);
          save_vec_scatter(pg2_xm, &buf[VLENY * ivc], wt, svidx_xm);
        }
      }
      {
        int    ix   = Nxv - 1;
        int    site = ix + Nxv * iyzt;
        real_t *buf = &buf_xm[VLENY * NVC * iyzt];

        svreal_t wt0, wt1, wt2, wt3, wt4, wt5;
        load_vec(pg2_xp, wt0, &v1[VLEN * (0 + NVC * site)]);
        load_vec(pg2_xp, wt1, &v1[VLEN * (1 + NVC * site)]);
        load_vec(pg2_xp, wt2, &v1[VLEN * (2 + NVC * site)]);
        load_vec(pg2_xp, wt3, &v1[VLEN * (3 + NVC * site)]);
        load_vec(pg2_xp, wt4, &v1[VLEN * (4 + NVC * site)]);
        load_vec(pg2_xp, wt5, &v1[VLEN * (5 + NVC * site)]);

        set_index_xp(svidx_xp);
        for (int ic = 0; ic < NC; ++ic) {
          svreal_t ut0, ut1, ut2, ut3, ut4, ut5;
          svreal_t xtr, xti;
          load_udag(pg2_xp, ut0, ut1, ut2, ut3, ut4, ut5,
                    &u[VLEN * (NVC * ic + NDF * site)]);
          mult_udv(pg2_xp, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
                   wt0, wt1, wt2, wt3, wt4, wt5);
          //svst1_scatter_index(pg2_xp, &buf[VLENY*(2*ic)  ], svidx_xp, xtr);
          //svst1_scatter_index(pg2_xp, &buf[VLENY*(2*ic+1)], svidx_xp, xti);
          save_vec_scatter(pg2_xp, &buf[VLENY * (2 * ic)], xtr, svidx_xp);
          save_vec_scatter(pg2_xp, &buf[VLENY * (2 * ic + 1)], xti, svidx_xp);
        }
      }
    }
  }

  if (do_comm[1] > 0) {
    int    idir = 1;
    real_t *u   = &up[NDF * Nst * idir];

    int Nxzt = Nxv * Nz * Nt;
    int ith, nth, isy, nsy;
    set_threadtask(ith, nth, isy, nsy, Nxzt);

    for (int ixzt = isy; ixzt < nsy; ++ixzt) {
      int ix  = ixzt % Nxv;
      int izt = ixzt / Nxv;
      {
        int    iy   = 0;
        int    site = ix + Nxv * (iy + Nyv * izt);
        real_t *buf = &buf_yp[VLENX * NVC * ixzt];
        for (int ivc = 0; ivc < NVC; ++ivc) {
          svreal_t wt;
          load_vec(pg2_ym, wt, &v1[VLEN * (ivc + NVC * site)]);
          //svst1(pg2_ym, &buf[VLENX*ivc], wt);
          save_vec(pg2_ym, &buf[VLENX * ivc], wt);
        }
      }
      {
        int    iy   = Nyv - 1;
        int    site = ix + Nxv * (iy + Nyv * izt);
        real_t *buf = &buf_ym[VLENX * NVC * ixzt];

        svreal_t wt0, wt1, wt2, wt3, wt4, wt5;
        load_vec(pg2_yp, wt0, &v1[VLEN * (0 + NVC * site)]);
        load_vec(pg2_yp, wt1, &v1[VLEN * (1 + NVC * site)]);
        load_vec(pg2_yp, wt2, &v1[VLEN * (2 + NVC * site)]);
        load_vec(pg2_yp, wt3, &v1[VLEN * (3 + NVC * site)]);
        load_vec(pg2_yp, wt4, &v1[VLEN * (4 + NVC * site)]);
        load_vec(pg2_yp, wt5, &v1[VLEN * (5 + NVC * site)]);

        for (int ic = 0; ic < NC; ++ic) {
          svreal_t ut0, ut1, ut2, ut3, ut4, ut5;
          svreal_t xtr, xti;
          load_udag(pg2_yp, ut0, ut1, ut2, ut3, ut4, ut5,
                    &u[VLEN * (NVC * ic + NDF * site)]);
          mult_udv(pg2_yp, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
                   wt0, wt1, wt2, wt3, wt4, wt5);
          int offset = -VLENX * (VLENY - 1);
          //svst1(pg2_yp, &buf[offset + VLENX*(2*ic)  ], xtr);
          //svst1(pg2_yp, &buf[offset + VLENX*(2*ic+1)], xti);
          save_vec(pg2_yp, &buf[offset + VLENX * (2 * ic)], xtr);
          save_vec(pg2_yp, &buf[offset + VLENX * (2 * ic + 1)], xti);
        }
      }
    }
  }

  if (do_comm[2] > 0) {
    int idir = 2;

    int Nxyt = Nxv * Nyv * Nt;
    int ith, nth, isz, nsz;
    set_threadtask(ith, nth, isz, nsz, Nxyt);

    for (int ixyt = isz; ixyt < nsz; ++ixyt) {
      int ixy = ixyt % Nxy;
      int it  = ixyt / Nxy;
      {
        int    iz   = 0;
        int    site = ixy + Nxy * (iz + Nz * it);
        real_t *buf = &buf_zp[VLEN * NVC * ixyt];
        for (int ivc = 0; ivc < NVC; ++ivc) {
          svreal_t wt;
          load_vec(pg, wt, &v1[VLEN * (ivc + NVC * site)]);
          //svst1(pg, &buf[VLEN * ivc], wt);
          save_vec(pg, &buf[VLEN * ivc], wt);
        }
      }
      {
        int iz   = Nz - 1;
        int site = ixy + Nxy * (iz + Nz * it);
        mult_staggered_dn1(pg, &buf_zm[VLEN * NVC * ixyt],
                           &up[VLEN * NDF * (site + Nstv * 2)],
                           &v1[VLEN * NVC * site]);
      }
    }
  }

  if (do_comm[3] > 0) {
    int idir = 3;

    int ith, nth, ist, nst;
    set_threadtask(ith, nth, ist, nst, Nxyz);

    for (int ixyz = ist; ixyz < nst; ++ixyz) {
      {
        int    it   = 0;
        int    site = ixyz + Nxyz * it;
        real_t *buf = &buf_tp[VLEN * NVC * ixyz];
        for (int ivc = 0; ivc < NVC; ++ivc) {
          svreal_t wt;
          load_vec(pg, wt, &v1[VLEN * (ivc + NVC * site)]);
          //svst1(pg, &buf[VLEN * ivc], wt);
          save_vec(pg, &buf[VLEN * ivc], wt);
        }
      }
      {
        int it   = Nt - 1;
        int site = ixyz + Nxyz * it;
        mult_staggered_dn1(pg, &buf_tm[VLEN * NVC * ixyz],
                           &up[VLEN * NDF * (site + Nstv * 3)],
                           &v1[VLEN * NVC * site]);
      }
    }
  }
}


//====================================================================
void BridgeQXS::mult_staggered_2(real_t *v2, real_t *up, real_t *v1,
                                 real_t *buf_xp, real_t *buf_xm,
                                 real_t *buf_yp, real_t *buf_ym,
                                 real_t *buf_zp, real_t *buf_zm,
                                 real_t *buf_tp, real_t *buf_tm,
                                 real_t qm, int jd,
                                 int *Nsize, int *do_comm)
{
  int Nxv  = Nsize[0];
  int Nyv  = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Nyv * Nz * Nt;
  int Nst  = Nstv * VLEN;

  int Nxy  = Nxv * Nyv;
  int Nxyz = Nxv * Nyv * Nz;

  real_t fac = 0.5 * real_t(jd);

  svbool_t pg = set_predicate();
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

    svreal_t vt0, vt1, vt2, vt3, vt4, vt5;
    clear_vec(pg, vt0);
    clear_vec(pg, vt1);
    clear_vec(pg, vt2);
    clear_vec(pg, vt3);
    clear_vec(pg, vt4);
    clear_vec(pg, vt5);

    int opr_any = 0;

    if ((ix == Nxv - 1) && (do_comm[0] > 0)) {
      real_t *u   = &up[NDF * Nst * 0];
      real_t *buf = &buf_xp[VLENY * NVC * iyzt];

      set_index_xp(svidx_xp);
      svreal_t wt0, wt1, wt2, wt3, wt4, wt5;
      load_vec(pg1_xp, wt0, &v1[VLEN * (0 + NVC * site) + 1]);
      load_vec(pg1_xp, wt1, &v1[VLEN * (1 + NVC * site) + 1]);
      load_vec(pg1_xp, wt2, &v1[VLEN * (2 + NVC * site) + 1]);
      load_vec(pg1_xp, wt3, &v1[VLEN * (3 + NVC * site) + 1]);
      load_vec(pg1_xp, wt4, &v1[VLEN * (4 + NVC * site) + 1]);
      load_vec(pg1_xp, wt5, &v1[VLEN * (5 + NVC * site) + 1]);

      load_add_gather(pg2_xp, wt0, &buf[VLENY * 0], svidx_xp);
      load_add_gather(pg2_xp, wt1, &buf[VLENY * 1], svidx_xp);
      load_add_gather(pg2_xp, wt2, &buf[VLENY * 2], svidx_xp);
      load_add_gather(pg2_xp, wt3, &buf[VLENY * 3], svidx_xp);
      load_add_gather(pg2_xp, wt4, &buf[VLENY * 4], svidx_xp);
      load_add_gather(pg2_xp, wt5, &buf[VLENY * 5], svidx_xp);

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

    if ((ix == 0) && (do_comm[0] > 0)) {
      real_t *u   = &up[NDF * Nst * 0];
      real_t *buf = &buf_xm[VLENY * NVC * iyzt];

      svreal_t wt0, wt1, wt2, wt3, wt4, wt5;
      load_vec(pg1_xm, wt0, &v1[VLEN * (0 + NVC * site) - 1]);
      load_vec(pg1_xm, wt1, &v1[VLEN * (1 + NVC * site) - 1]);
      load_vec(pg1_xm, wt2, &v1[VLEN * (2 + NVC * site) - 1]);
      load_vec(pg1_xm, wt3, &v1[VLEN * (3 + NVC * site) - 1]);
      load_vec(pg1_xm, wt4, &v1[VLEN * (4 + NVC * site) - 1]);
      load_vec(pg1_xm, wt5, &v1[VLEN * (5 + NVC * site) - 1]);

      set_index_xm(svidx_xm);
      svreal_t ut0, ut1, ut2, ut3, ut4, ut5;
      svreal_t xtr, xti;
      load_udag(pg1_xm, ut0, ut1, ut2, ut3, ut4, ut5,
                &u[VLEN * (0 + NDF * site) - 1]);
      mult_udv(pg1_xm, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
               wt0, wt1, wt2, wt3, wt4, wt5);
      load_add_gather(pg2_xm, xtr, &buf[VLENY * 0], svidx_xm);
      load_add_gather(pg2_xm, xti, &buf[VLENY * 1], svidx_xm);
      sub_vec(pg, vt0, xtr);
      sub_vec(pg, vt1, xti);

      load_udag(pg1_xm, ut0, ut1, ut2, ut3, ut4, ut5,
                &u[VLEN * (NVC + NDF * site) - 1]);
      mult_udv(pg1_xm, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
               wt0, wt1, wt2, wt3, wt4, wt5);
      load_add_gather(pg2_xm, xtr, &buf[VLENY * 2], svidx_xm);
      load_add_gather(pg2_xm, xti, &buf[VLENY * 3], svidx_xm);
      sub_vec(pg, vt2, xtr);
      sub_vec(pg, vt3, xti);

      load_udag(pg1_xm, ut0, ut1, ut2, ut3, ut4, ut5,
                &u[VLEN * (2 * NVC + NDF * site) - 1]);
      mult_udv(pg1_xm, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
               wt0, wt1, wt2, wt3, wt4, wt5);
      load_add_gather(pg2_xm, xtr, &buf[VLENY * 4], svidx_xm);
      load_add_gather(pg2_xm, xti, &buf[VLENY * 5], svidx_xm);
      sub_vec(pg, vt4, xtr);
      sub_vec(pg, vt5, xti);

      ++opr_any;
    }

    if ((iy == Nyv - 1) && (do_comm[1] > 0)) {
      real_t *u   = &up[NDF * Nst * 1];
      int    ixzt = ix + Nxv * izt;
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
      real_t *u   = &up[NDF * Nst * 1];
      int    ixzt = ix + Nxv * izt;
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
                        &up[VLEN * NDF * (site + Nstv * 2)],
                        &buf_zp[VLEN * NVC * ixyt]);
      ++opr_any;
    }

    if ((iz == 0) && (do_comm[2] > 0)) {
      int    ixyt = ixy + Nxy * it;
      real_t *buf = &buf_zm[VLEN * NVC * ixyt];
      real_t *u   = &up[NDF * Nst * 2];

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
                        &up[VLEN * NDF * (site + Nstv * 3)],
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
