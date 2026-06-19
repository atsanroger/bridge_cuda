/*!
      @file    mult_Doainwall_5din_qxs-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2581 $
*/

#ifndef MULT_DOMAINWALL_5DIN_QXS_INCLUDED
#define MULT_DOMAINWALL_5DIN_QXS_INCLUDED

#include "mult_common_th-inc.h"

//====================================================================
void BridgeQXS::mult_domainwall_5din_5dir_dirac(
  real_t *vp, real_t *yp, real_t *wp,
  real_t mq, real_t M0, int Ns, int *bc,
  real_t *b, real_t *c,
  int *Nsize, int *do_comm)
{
  int Nxv  = Nsize[0];
  int Nyv  = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Nyv * Nz * Nt;
  int Nst  = Nstv * VLEN;

  int Nin4 = VLEN * NVCD;
  int Nin5 = Nin4 * Ns;

  int ith, nth, site0, site1;
  set_threadtask(ith, nth, site0, site1, Nstv);

  for (int site = site0; site < site1; ++site) {
    real_t *vp2 = &vp[Nin5 * site];
    real_t *wp2 = &wp[Nin5 * site];
    real_t *yp2 = &yp[Nin5 * site];

    for (int is = 0; is < Ns; ++is) {
      svbool_t pg = set_predicate();

      for (int ic = 0; ic < NC; ++ic) {
        svreal_t vt1r, vt1i, vt2r, vt2i, vt3r, vt3i, vt4r, vt4i;

        int    is_up = (is + 1) % Ns;
        real_t Fup   = 0.5;
        if (is == Ns - 1) Fup = -0.5 * mq;
        set_aPm5_dirac_vec(vt1r, vt1i, vt2r, vt2i,
                           vt3r, vt3i, vt4r, vt4i,
                           Fup, wp2, is_up, ic);

        int    is_dn = (is - 1 + Ns) % Ns;
        real_t Fdn   = 0.5;
        if (is == 0) Fdn = -0.5 * mq;
        add_aPp5_dirac_vec(vt1r, vt1i, vt2r, vt2i,
                           vt3r, vt3i, vt4r, vt4i,
                           Fdn, wp2, is_dn, ic);

        svreal_t xt;
        real_t   FF1    = b[is] * (4.0 - M0) + 1.0;
        real_t   FF2    = c[is] * (4.0 - M0) - 1.0;
        int      offset = 2 * ND * ic + NVCD * is;
        dw_5dir_axpy(pg, vp2, yp2, wp2, FF1, FF2, b[is], c[is], vt1r, 0 + offset);
        dw_5dir_axpy(pg, vp2, yp2, wp2, FF1, FF2, b[is], c[is], vt1i, 1 + offset);
        dw_5dir_axpy(pg, vp2, yp2, wp2, FF1, FF2, b[is], c[is], vt2r, 2 + offset);
        dw_5dir_axpy(pg, vp2, yp2, wp2, FF1, FF2, b[is], c[is], vt2i, 3 + offset);
        dw_5dir_axpy(pg, vp2, yp2, wp2, FF1, FF2, b[is], c[is], vt3r, 4 + offset);
        dw_5dir_axpy(pg, vp2, yp2, wp2, FF1, FF2, b[is], c[is], vt3i, 5 + offset);
        dw_5dir_axpy(pg, vp2, yp2, wp2, FF1, FF2, b[is], c[is], vt4r, 6 + offset);
        dw_5dir_axpy(pg, vp2, yp2, wp2, FF1, FF2, b[is], c[is], vt4i, 7 + offset);
      }
    }
  }
}


//====================================================================
void BridgeQXS::accum_mult_domainwall_5din_5dir_dirac(
  real_t *vp, real_t *yp, real_t *wp,
  real_t mq, real_t M0, int Ns, int *bc,
  real_t *b, real_t *c,
  int *Nsize, int *do_comm)
{
  int Nxv  = Nsize[0];
  int Nyv  = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Nyv * Nz * Nt;
  int Nst  = Nstv * VLEN;

  int Nin4 = VLEN * NVCD;
  int Nin5 = Nin4 * Ns;

  int ith, nth, site0, site1;
  set_threadtask(ith, nth, site0, site1, Nstv);

  for (int site = site0; site < site1; ++site) {
    real_t *vp2 = &vp[Nin5 * site];
    real_t *wp2 = &wp[Nin5 * site];
    real_t *yp2 = &yp[Nin5 * site];

    for (int is = 0; is < Ns; ++is) {
      svbool_t pg = set_predicate();

      for (int ic = 0; ic < NC; ++ic) {
        svreal_t vt1r, vt1i, vt2r, vt2i, vt3r, vt3i, vt4r, vt4i;

        int    is_up = (is + 1) % Ns;
        real_t Fup   = 0.5;
        if (is == Ns - 1) Fup = -0.5 * mq;
        set_aPm5_dirac_vec(vt1r, vt1i, vt2r, vt2i,
                           vt3r, vt3i, vt4r, vt4i,
                           Fup, wp2, is_up, ic);

        int    is_dn = (is - 1 + Ns) % Ns;
        real_t Fdn   = 0.5;
        if (is == 0) Fdn = -0.5 * mq;
        add_aPp5_dirac_vec(vt1r, vt1i, vt2r, vt2i,
                           vt3r, vt3i, vt4r, vt4i,
                           Fdn, wp2, is_dn, ic);

        svreal_t xt;
        real_t   FF1    = b[is] * (4.0 - M0) + 1.0;
        real_t   FF2    = c[is] * (4.0 - M0) - 1.0;
        int      offset = 2 * ND * ic + NVCD * is;
        accum_dw_5dir_axpy(pg, vp2, yp2, wp2, FF1, FF2, b[is], c[is], vt1r, 0 + offset);
        accum_dw_5dir_axpy(pg, vp2, yp2, wp2, FF1, FF2, b[is], c[is], vt1i, 1 + offset);
        accum_dw_5dir_axpy(pg, vp2, yp2, wp2, FF1, FF2, b[is], c[is], vt2r, 2 + offset);
        accum_dw_5dir_axpy(pg, vp2, yp2, wp2, FF1, FF2, b[is], c[is], vt2i, 3 + offset);
        accum_dw_5dir_axpy(pg, vp2, yp2, wp2, FF1, FF2, b[is], c[is], vt3r, 4 + offset);
        accum_dw_5dir_axpy(pg, vp2, yp2, wp2, FF1, FF2, b[is], c[is], vt3i, 5 + offset);
        accum_dw_5dir_axpy(pg, vp2, yp2, wp2, FF1, FF2, b[is], c[is], vt4r, 6 + offset);
        accum_dw_5dir_axpy(pg, vp2, yp2, wp2, FF1, FF2, b[is], c[is], vt4i, 7 + offset);
      }
    }
  }
}


//====================================================================
void BridgeQXS::mult_domainwall_5din_5dirdag_dirac(
  real_t *vp, real_t *yp, real_t *wp,
  real_t mq, real_t M0, int Ns, int *bc,
  real_t *b, real_t *c,
  int *Nsize, int *do_comm)
{
  int Nxv  = Nsize[0];
  int Nyv  = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Nyv * Nz * Nt;
  int Nst  = Nstv * VLEN;

  int Nin4 = VLEN * NVCD;
  int Nin5 = Nin4 * Ns;

  int ith, nth, site0, site1;
  set_threadtask(ith, nth, site0, site1, Nstv);

  for (int site = site0; site < site1; ++site) {
    real_t *vp2 = &vp[Nin5 * site];
    real_t *wp2 = &wp[Nin5 * site];
    real_t *yp2 = &yp[Nin5 * site];

    Vsimd_t vL[NVCD], yL[NVCD], xL[NVCD];

    real_t FF1, FF2;

    for (int is = 0; is < Ns; ++is) {
      svbool_t pg = set_predicate();

      for (int ic = 0; ic < NC; ++ic) {
        svreal_t vt1r, vt1i, vt2r, vt2i, vt3r, vt3i, vt4r, vt4i;
        FF1 = b[is] * (4.0 - M0) + 1.0;
        real_t a1    = -0.5 * b[is];
        int    index = 2 * ND * ic + NVCD * is;
        // dw_5dir_dag(pg, vt1r, vt1i, vt2r, vt2i, vt3r, vt3i, vt4r, vt4i,
        //             wp2, yp2, FF1, a1, index);
        // note: origianl operation multiply gm5 here.
        dw_5dir_dag(pg, vt3r, vt3i, vt4r, vt4i, vt1r, vt1i, vt2r, vt2i,
                    wp2, yp2, -FF1, -a1, index);

        int is_up = (is + 1) % Ns;
        FF2 = c[is_up] * (4.0 - M0) - 1.0;
        real_t   aup = -0.5 * c[is_up];
        int      index_up = 2 * ND * ic + NVCD * is_up;
        svreal_t xt1r, xt1i, xt2r, xt2i, xt3r, xt3i, xt4r, xt4i;
        dw_5dir_dag(pg, xt1r, xt1i, xt2r, xt2i, xt3r, xt3i, xt4r, xt4i,
                    wp2, yp2, FF2, aup, index_up);

        real_t Fup = 0.5;
        if (is == Ns - 1) Fup = -0.5 * mq;
        add_aPp5_dirac_vec(pg, vt1r, vt1i, vt2r, vt2i, vt3r, vt3i, vt4r, vt4i,
                           Fup, xt1r, xt1i, xt2r, xt2i, xt3r, xt3i, xt4r, xt4i);

        int is_dn = (is - 1 + Ns) % Ns;
        FF2 = c[is_dn] * (4.0 - M0) - 1.0;
        real_t adn      = -0.5 * c[is_dn];
        int    index_dn = 2 * ND * ic + NVCD * is_dn;
        dw_5dir_dag(pg, xt1r, xt1i, xt2r, xt2i, xt3r, xt3i, xt4r, xt4i,
                    wp2, yp2, FF2, adn, index_dn);

        real_t Fdn = 0.5;
        if (is == 0) Fdn = -0.5 * mq;
        add_aPm5_dirac_vec(pg, vt1r, vt1i, vt2r, vt2i, vt3r, vt3i, vt4r, vt4i,
                           -Fdn, xt1r, xt1i, xt2r, xt2i, xt3r, xt3i, xt4r, xt4i);

        int offset = 2 * ND * ic + NVCD * is;
        save_vec(pg, &vp2[VLEN * (0 + offset)], vt1r);
        save_vec(pg, &vp2[VLEN * (1 + offset)], vt1i);
        save_vec(pg, &vp2[VLEN * (2 + offset)], vt2r);
        save_vec(pg, &vp2[VLEN * (3 + offset)], vt2i);
        save_vec(pg, &vp2[VLEN * (4 + offset)], vt3r);
        save_vec(pg, &vp2[VLEN * (5 + offset)], vt3i);
        save_vec(pg, &vp2[VLEN * (6 + offset)], vt4r);
        save_vec(pg, &vp2[VLEN * (7 + offset)], vt4i);
      }
    }
  }
}


//====================================================================
void BridgeQXS::mult_domainwall_5din_mult_gm5_dirac(real_t *vp, real_t *wp,
                                                    int Ns, int *Nsize)
{
  int Nxv  = Nsize[0];
  int Nyv  = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Nyv * Nz * Nt;

  int Nin4 = VLEN * NVCD;
  int Nin5 = Nin4 * Ns;

  int ith, nth, site0, site1;
  set_threadtask(ith, nth, site0, site1, Nstv);

  for (int site = site0; site < site1; ++site) {
    real_t *vp2 = &vp[Nin5 * site];
    real_t *wp2 = &wp[Nin5 * site];

    svbool_t pg = set_predicate();
    for (int is = 0; is < Ns; ++is) {
      for (int ic = 0; ic < NC; ++ic) {
        svreal_t vt1r, vt1i, vt2r, vt2i, vt3r, vt3i, vt4r, vt4i;
        int      index = 2 * ND * ic + NVCD * is;
        load_vec(pg, vt1r, &wp2[VLEN * (4 + index)]);
        load_vec(pg, vt1i, &wp2[VLEN * (5 + index)]);
        load_vec(pg, vt2r, &wp2[VLEN * (6 + index)]);
        load_vec(pg, vt2i, &wp2[VLEN * (7 + index)]);
        load_vec(pg, vt3r, &wp2[VLEN * (0 + index)]);
        load_vec(pg, vt3i, &wp2[VLEN * (1 + index)]);
        load_vec(pg, vt4r, &wp2[VLEN * (2 + index)]);
        load_vec(pg, vt4i, &wp2[VLEN * (3 + index)]);

        flip_sign(pg, vt1r);  // vt1r := -vt1r
        flip_sign(pg, vt1i);
        flip_sign(pg, vt2r);
        flip_sign(pg, vt2i);
        flip_sign(pg, vt3r);
        flip_sign(pg, vt3i);
        flip_sign(pg, vt4r);
        flip_sign(pg, vt4i);
        save_vec(pg, &vp2[VLEN * (0 + index)], vt1r);
        save_vec(pg, &vp2[VLEN * (1 + index)], vt1i);
        save_vec(pg, &vp2[VLEN * (2 + index)], vt2r);
        save_vec(pg, &vp2[VLEN * (3 + index)], vt2i);
        save_vec(pg, &vp2[VLEN * (4 + index)], vt3r);
        save_vec(pg, &vp2[VLEN * (5 + index)], vt3i);
        save_vec(pg, &vp2[VLEN * (6 + index)], vt4r);
        save_vec(pg, &vp2[VLEN * (7 + index)], vt4i);
      }
    }
  }
}


//====================================================================
void BridgeQXS::mult_domainwall_5din_clear(real_t *vp, int Ns, int *Nsize)
{
  int Nxv  = Nsize[0];
  int Nyv  = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Nyv * Nz * Nt;

  int Nin4 = VLEN * NVCD;
  int Nin5 = Nin4 * Ns;

  int ith, nth, site0, site1;
  set_threadtask(ith, nth, site0, site1, Nstv);

  svbool_t pg = set_predicate();
  svreal_t vzero;
  clear_vec(pg, vzero);

  for (int site = site0; site < site1; ++site) {
    real_t *vp2 = &vp[Nin5 * site];
    for (int is = 0; is < Ns; ++is) {
      for (int in = 0; in < NVCD; ++in) {
        save_vec(pg, &vp2[Nin4*is + VLEN*in], vzero);
      }
    }
  }

  /*
  Vsimd_t yL[NVCD];
  clear_vec(yL, NVCD);

  for (int site = site0; site < site1; ++site) {
    real_t *vp2 = &vp[Nin5 * site];

    for (int is = 0; is < Ns; ++is) {

        save_vec(pg, &vp2[VLEN * (0 + index)], vt1r);

      save_vec(&vp2[Nin4 * is], yL, NVCD);
    }
  }
  */
}


//====================================================================
void BridgeQXS::mult_domainwall_5din_hopb_dirac(
  real_t *vp, real_t *up, real_t *wp,
  real_t mq, real_t M0, int Ns, int *bc,
  real_t *b, real_t *c,
  int *Nsize, int *do_comm)
{
  int Nxv  = Nsize[0];
  int Nyv  = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Nyv * Nz * Nt;
  int Nst  = Nstv * VLEN;

  int Nin4 = VLEN * NVCD;
  int Nin5 = Nin4 * Ns;

  int NvU = NDF * Nst;

  int Nxy  = Nxv * Nyv;
  int Nxyz = Nxv * Nyv * Nz;

  svbool_t pg1_xp, pg2_xp, pg1_xm, pg2_xm;
  svbool_t pg1_yp, pg2_yp, pg1_ym, pg2_ym;
  set_predicate_xp(pg1_xp, pg2_xp);
  set_predicate_xm(pg1_xm, pg2_xm);
  set_predicate_yp(pg1_yp, pg2_yp);
  set_predicate_ym(pg1_ym, pg2_ym);

  int ith, nth, site0, site1;
  set_threadtask(ith, nth, site0, site1, Nstv);

  real_t bufL[Nin5];

  for (int site = site0; site < site1; ++site) {
    int ix   = site % Nxv;
    int iyzt = site / Nxv;
    int ixy  = site % Nxy;
    int iy   = iyzt % Nyv;
    int izt  = site / Nxy;
    int iz   = izt % Nz;
    int it   = izt / Nz;
    int ixyz = site % Nxyz;

    int idir, nei;

    real_t *wp2 = &wp[Nin5 * site];
    real_t *vp2 = &vp[Nin5 * site];

    real_t  z4[VLEN * NVCD];
    Vsimd_t vL[NVCD * Ns];
    Vsimd_t uL[NDF];

    load_vec(vL, vp2, NVCD * Ns);

    idir = 0;

    if ((ix < Nxv - 1) || (do_comm[idir] == 0)) {
      int ix2 = (ix + 1) % Nxv;
      nei = ix2 + Nxv * iyzt;
      real_t *wpn = &wp[Nin5 * nei];
      real_t *u   = &up[VLEN * NDF * site + NvU * idir];
      for (int is = 0; is < Ns; ++is) {
        mult_wilson_xpb(pg1_xp, pg2_xp, &vL[NVCD * is].v[0], u,
                        &wp2[Nin4 * is], &wpn[Nin4 * is]);
      }
    }

    if ((ix > 0) || (do_comm[idir] == 0)) {
      int ix2 = (ix - 1 + Nxv) % Nxv;
      nei = ix2 + Nxv * iyzt;
      real_t *wpn = &wp[Nin5 * nei];
      real_t *ux  = &up[VLEN * NDF * site + NvU * idir];
      real_t *un  = &up[VLEN * NDF * nei + NvU * idir];
      for (int is = 0; is < Ns; ++is) {
        mult_wilson_xmb(pg1_xm, pg2_xm, &vL[NVCD * is], ux, un,
                        &wp2[Nin4 * is], &wpn[Nin4 * is]);
      }
    }

    idir = 1;

    if ((iy < Nyv - 1) || (do_comm[idir] == 0)) {
      int iy2 = (iy + 1) % Nyv;
      nei = ix + Nxv * (iy2 + Nyv * izt);
      real_t *wpn = &wp[Nin5 * nei];
      real_t *u   = &up[VLEN * NDF * site + NvU * idir];
      for (int is = 0; is < Ns; ++is) {
        mult_wilson_ypb(pg1_yp, pg2_yp, &vL[NVCD * is].v[0], u,
                        &wp2[Nin4 * is], &wpn[Nin4 * is]);
      }
    }

    if ((iy > 0) || (do_comm[idir] == 0)) {
      int iy2 = (iy - 1 + Nyv) % Nyv;
      nei = ix + Nxv * (iy2 + Nyv * izt);
      real_t *wpn = &wp[Nin5 * nei];
      real_t *ux  = &up[VLEN * NDF * site + NvU * idir];
      real_t *un  = &up[VLEN * NDF * nei + NvU * idir];
      for (int is = 0; is < Ns; ++is) {
        mult_wilson_ymb(pg1_ym, pg2_ym, &vL[NVCD * is].v[0], ux, un,
                        &wp2[Nin4 * is], &wpn[Nin4 * is]);
      }
    }

    idir = 2;

    if ((iz < Nz - 1) || (do_comm[idir] == 0)) {
      int iz2 = (iz + 1) % Nz;
      nei = ixy + Nxy * (iz2 + Nz * it);
      real_t *wpn = &wp[Nin5 * nei];
      real_t *u   = &up[VLEN * NDF * site + NvU * idir];
      for (int is = 0; is < Ns; ++is) {
        mult_wilson_zpb(&vL[NVCD * is].v[0], u, &wpn[Nin4 * is]);
      }
    }

    if ((iz > 0) || (do_comm[idir] == 0)) {
      int iz2 = (iz - 1 + Nz) % Nz;
      nei = ixy + Nxy * (iz2 + Nz * it);
      real_t *wpn = &wp[Nin5 * nei];
      real_t *u   = &up[VLEN * NDF * nei + NvU * idir];
      for (int is = 0; is < Ns; ++is) {
        mult_wilson_zmb(&vL[NVCD * is].v[0], u, &wpn[Nin4 * is]);
      }
    }

    idir = 3;

    if ((it < Nt - 1) || (do_comm[idir] == 0)) {
      int it2 = (it + 1) % Nt;
      nei = ixyz + Nxyz * it2;
      real_t *wpn = &wp[Nin5 * nei];
      real_t *u   = &up[VLEN * NDF * site + NvU * idir];
      for (int is = 0; is < Ns; ++is) {
        mult_wilson_tpb_dirac(&vL[NVCD * is].v[0], u, &wpn[Nin4 * is]);
      }
    }

    if ((it > 0) || (do_comm[idir] == 0)) {
      int it2 = (it - 1 + Nt) % Nt;
      nei = ixyz + Nxyz * it2;
      real_t *wpn = &wp[Nin5 * nei];
      real_t *u   = &up[VLEN * NDF * nei + NvU * idir];
      for (int is = 0; is < Ns; ++is) {
        mult_wilson_tmb_dirac(&vL[NVCD * is].v[0], u, &wpn[Nin4 * is]);
      }
    }

    save_vec(vp2, vL, NVCD * Ns);
  }
}


//====================================================================
void BridgeQXS::mult_domainwall_5din_hop1_dirac(
  real_t *buf1_xp, real_t *buf1_xm,
  real_t *buf1_yp, real_t *buf1_ym,
  real_t *buf1_zp, real_t *buf1_zm,
  real_t *buf1_tp, real_t *buf1_tm,
  real_t *up, real_t *wp,
  real_t mq, real_t M0, int Ns, int *bc,
  int *Nsize, int *do_comm)
{
  int Nxv  = Nsize[0];
  int Nyv  = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Nyv * Nz * Nt;
  int Nst  = Nstv * VLEN;
  int NvU  = NDF * Nst;

  int Nin4  = VLEN * NVCD;
  int Nin5  = Nin4 * Ns;
  int Nin4H = VLEN * NVC * ND2;
  int Nin5H = Nin4H * Ns;

  svbool_t pg1_xp, pg2_xp, pg1_xm, pg2_xm;
  svbool_t pg1_yp, pg2_yp, pg1_ym, pg2_ym;
  set_predicate_xp(pg1_xp, pg2_xp);
  set_predicate_xm(pg1_xm, pg2_xm);
  set_predicate_yp(pg1_yp, pg2_yp);
  set_predicate_ym(pg1_ym, pg2_ym);
  svint_t svidx_xp, svidx_xm;
  set_index_xp(svidx_xp);
  set_index_xm(svidx_xm);

  if (do_comm[0] == 1) {
    int idir = 0;

    int Nyzt = Nyv * Nz * Nt;

    int ith, nth, site0, site1;
    set_threadtask(ith, nth, site0, site1, Nyzt);

    for (int iyzt = site0; iyzt < site1; ++iyzt) {
      {
        int    ix   = 0;
        int    site = ix + Nxv * iyzt;
        real_t *wp2 = &wp[Nin5 * site];
        int    ibf  = VLENY * NVC * ND2 * Ns * iyzt;
        set_index_xm(svidx_xm);
        for (int is = 0; is < Ns; ++is) {
          mult_wilson_xp1(pg2_xm, svidx_xm,
                          &buf1_xp[ibf + VLENY * NVC * ND2 * is], &wp2[Nin4 * is]);
        }
      }

      {
        int    ix   = Nxv - 1;
        int    site = ix + Nxv * iyzt;
        real_t *wp2 = &wp[Nin5 * site];
        int    ibf  = VLENY * NVC * ND2 * Ns * iyzt;
        real_t *u2  = &up[VLEN * NDF * site + NvU * idir];
        set_index_xp(svidx_xp);
        for (int is = 0; is < Ns; ++is) {
          mult_wilson_xm1(pg2_xp, svidx_xp,
                          &buf1_xm[ibf + VLENY * NVC * ND2 * is], u2, &wp2[Nin4 * is]);
        }
      }
    }
  } // do_comm[0]

  if (do_comm[1] == 1) {
    int idir = 1;
    int Nxzt = Nxv * Nz * Nt;

    int ith, nth, site0, site1;
    set_threadtask(ith, nth, site0, site1, Nxzt);

    for (int ixzt = site0; ixzt < site1; ++ixzt) {
      int ix  = ixzt % Nxv;
      int izt = ixzt / Nxv;

      {
        int    iy   = 0;
        int    site = ix + Nxv * (iy + Nyv * izt);
        real_t *wp2 = &wp[Nin5 * site];

        int ibf = VLENX * NVC * ND2 * Ns * (ix + Nxv * izt);
        for (int is = 0; is < Ns; ++is) {
          mult_wilson_yp1(pg2_ym,
                          &buf1_yp[ibf + VLENX * NVC * ND2 * is], &wp2[Nin4 * is]);
        }
      }

      {
        int    iy   = Nyv - 1;
        int    site = ix + Nxv * (iy + Nyv * izt);
        real_t *wp2 = &wp[Nin5 * site];

        int    ibf = VLENX * NVC * ND2 * Ns * (ix + Nxv * izt);
        real_t *u2 = &up[VLEN * NDF * site + NvU * idir];
        for (int is = 0; is < Ns; ++is) {
          mult_wilson_ym1(pg2_yp,
                          &buf1_ym[ibf + VLENX * NVC * ND2 * is], u2, &wp2[Nin4 * is]);
        }
      }
    }
  } // do_comm[1]


  if (do_comm[2] == 1) {
    int idir = 2;
    int Nxy  = Nxv * Nyv;
    int Nxyt = Nxv * Nyv * Nt;

    int ith, nth, site0, site1;
    set_threadtask(ith, nth, site0, site1, Nxyt);

    for (int ixyt = site0; ixyt < site1; ++ixyt) {
      int ixy = ixyt % Nxy;
      int it  = ixyt / Nxy;

      {
        int    iz   = 0;
        int    site = ixy + Nxy * (iz + Nz * it);
        real_t *wp2 = &wp[Nin5 * site];
        int    ibf  = Nin5H * (ixy + Nxy * it);
        for (int is = 0; is < Ns; ++is) {
          mult_wilson_zp1(&buf1_zp[ibf + Nin4H * is], &wp2[Nin4 * is]);
        }
      }
      {
        int    iz   = Nz - 1;
        int    site = ixy + Nxy * (iz + Nz * it);
        real_t *wp2 = &wp[Nin5 * site];
        int    ibf  = Nin5H * (ixy + Nxy * it);
        real_t *u2  = &up[VLEN * NDF * site + NvU * idir];
        for (int is = 0; is < Ns; ++is) {
          mult_wilson_zm1(&buf1_zm[ibf + Nin4H * is], u2, &wp2[Nin4 * is]);
        }
      }
    }
  }

  if (do_comm[3] == 1) {
    int idir = 3;
    int Nxyz = Nxv * Nyv * Nz;

    int ith, nth, site0, site1;
    set_threadtask(ith, nth, site0, site1, Nxyz);

    for (int ixyz = site0; ixyz < site1; ++ixyz) {
      {
        int    it   = 0;
        int    site = ixyz + Nxyz * it;
        real_t *wp2 = &wp[Nin5 * site];
        int    ibf  = Nin5H * ixyz;
        for (int is = 0; is < Ns; ++is) {
          mult_wilson_tp1_dirac(&buf1_tp[ibf + Nin4H * is], &wp2[Nin4 * is]);
        }
      }

      {
        int    it   = Nt - 1;
        int    site = ixyz + Nxyz * it;
        real_t *wp2 = &wp[Nin5 * site];
        int    ibf  = Nin5H * ixyz;
        real_t *u2  = &up[VLEN * NDF * site + NvU * idir];
        for (int is = 0; is < Ns; ++is) {
          mult_wilson_tm1_dirac(&buf1_tm[ibf + Nin4H * is], u2, &wp2[Nin4 * is]);
        }
      }
    }
  }
}


//====================================================================
void BridgeQXS::mult_domainwall_5din_hop2_dirac(
  real_t *vp, real_t *up, real_t *wp,
  real_t *buf2_xp, real_t *buf2_xm,
  real_t *buf2_yp, real_t *buf2_ym,
  real_t *buf2_zp, real_t *buf2_zm,
  real_t *buf2_tp, real_t *buf2_tm,
  real_t mq, real_t M0, int Ns, int *bc,
  int *Nsize, int *do_comm)


{
  int Nxv  = Nsize[0];
  int Nyv  = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Nyv * Nz * Nt;
  int Nst  = Nstv * VLEN;

  int Nin4  = VLEN * NVCD;
  int Nin5  = Nin4 * Ns;
  int Nin4H = VLEN * NVC * ND2;
  int Nin5H = Nin4H * Ns;
  int NvU   = NDF * Nst;

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

  int ith, nth, site0, site1;
  set_threadtask(ith, nth, site0, site1, Nstv);

  real_t bufL[Nin5];

  for (int site = site0; site < site1; ++site) {
    int ix   = site % Nxv;
    int iyzt = site / Nxv;
    int ixy  = site % Nxy;
    int iy   = iyzt % Nyv;
    int izt  = site / Nxy;
    int iz   = izt % Nz;
    int it   = izt / Nz;
    int ixyz = site % Nxyz;

    int idir, nei;

    real_t *wp2 = &wp[Nin5 * site];
    real_t *vp2 = &vp[Nin5 * site];

    real_t z4[VLEN * NVCD];

    Vsimd_t vL[NVCD * Ns];

    clear_vec(vL, NVCD * Ns);

    int opr_any = 0;

    idir = 0;
    if (do_comm[idir] == 1) {
      if (ix == Nxv - 1) {
        real_t *u2 = &up[VLEN * NDF * site + NvU * idir];
        int    ibf = VLENY * NVC * ND2 * Ns * iyzt;
        for (int is = 0; is < Ns; ++is) {
          set_index_xp(svidx_xp);
          mult_wilson_xp2(pg1_xp, pg2_xp, svidx_xp,
                          &vL[NVCD * is].v[0], u2,
                          &wp2[Nin4 * is], &buf2_xp[ibf + VLENY * NVC * ND2 * is]);
        }
        ++opr_any;
      }

      if (ix == 0) {
        real_t *u2 = &up[VLEN * NDF * site + NvU * idir];
        int    ibf = VLENY * NVC * ND2 * Ns * iyzt;
        set_index_xm(svidx_xm);
        for (int is = 0; is < Ns; ++is) {
          mult_wilson_xm2(pg1_xm, pg2_xm, svidx_xm,
                          &vL[NVCD * is].v[0], u2,
                          &wp2[Nin4 * is], &buf2_xm[ibf + VLENY * NVC * ND2 * is]);
        }
        ++opr_any;
      }
    }

    idir = 1;
    if (do_comm[idir] == 1) {
      if (iy == Nyv - 1) {
        int    ibf = VLENX * NVC * ND2 * Ns * (ix + Nxv * izt);
        real_t *u2 = &up[VLEN * NDF * site + NvU * idir];
        for (int is = 0; is < Ns; ++is) {
          mult_wilson_yp2(pg1_yp, pg2_yp,
                          &vL[NVCD * is].v[0], u2,
                          &wp2[Nin4 * is], &buf2_yp[ibf + VLENX * NVC * ND2 * is]);
        }
        ++opr_any;
      }

      if (iy == 0) {
        int    ibf = VLENX * NVC * ND2 * Ns * (ix + Nxv * izt);
        real_t *u2 = &up[VLEN * NDF * site + NvU * idir];
        for (int is = 0; is < Ns; ++is) {
          mult_wilson_ym2(pg1_ym, pg2_ym,
                          &vL[NVCD * is].v[0], u2,
                          &wp2[Nin4 * is], &buf2_ym[ibf + VLENX * NVC * ND2 * is]);
        }
        ++opr_any;
      }
    }

    idir = 2;
    if (do_comm[idir] == 1) {
      if (iz == Nz - 1) {
        int    ibf = Nin5H * (ixy + Nxy * it);
        real_t *u2 = &up[VLEN * NDF * site + NvU * idir];
        for (int is = 0; is < Ns; ++is) {
          mult_wilson_zp2(&vL[NVCD * is].v[0], u2, &buf2_zp[ibf + Nin4H * is]);
        }
        ++opr_any;
      }

      if (iz == 0) {
        int ibf = Nin5H * (ixy + Nxy * it);
        for (int is = 0; is < Ns; ++is) {
          mult_wilson_zm2(&vL[NVCD * is].v[0], &buf2_zm[ibf + Nin4H * is]);
        }
        ++opr_any;
      }
    }

    idir = 3;
    if (do_comm[idir] == 1) {
      if (it == Nt - 1) {
        int    ibf = Nin5H * ixyz;
        real_t *u2 = &up[VLEN * NDF * site + NvU * idir];
        for (int is = 0; is < Ns; ++is) {
          mult_wilson_tp2_dirac(&vL[NVCD * is].v[0], u2, &buf2_tp[ibf + Nin4H * is]);
        }
        ++opr_any;
      }

      if (it == 0) {
        int ibf = Nin5H * ixyz;
        for (int is = 0; is < Ns; ++is) {
          mult_wilson_tm2_dirac(&vL[NVCD * is].v[0], &buf2_tm[ibf + Nin4H * is]);
        }
        ++opr_any;
      }
    }

    // if(opr_any > 0) save_vec(vp2, vL, NVCD * Ns);

    if (opr_any > 0) {
      svbool_t pg = set_predicate();
      for (int ivs = 0; ivs < NVCD * Ns; ivs += 4) {
        svreal_t vt1, vt2, vt3, vt4;
        svreal_t yt1, yt2, yt3, yt4;
        load_vec(pg, yt1, &vp2[VLEN * (ivs)]);
        load_vec(pg, yt2, &vp2[VLEN * (ivs + 1)]);
        load_vec(pg, yt3, &vp2[VLEN * (ivs + 2)]);
        load_vec(pg, yt4, &vp2[VLEN * (ivs + 3)]);
        load_vec(pg, vt1, &(vL[ivs].v[0]));
        load_vec(pg, vt2, &(vL[ivs + 1].v[0]));
        load_vec(pg, vt3, &(vL[ivs + 2].v[0]));
        load_vec(pg, vt4, &(vL[ivs + 3].v[0]));

        add_vec(pg, yt1, vt1);
        add_vec(pg, yt2, vt2);
        add_vec(pg, yt3, vt3);
        add_vec(pg, yt4, vt4);
        save_vec(pg, &vp2[VLEN * (ivs)], yt1);
        save_vec(pg, &vp2[VLEN * (ivs + 1)], yt2);
        save_vec(pg, &vp2[VLEN * (ivs + 2)], yt3);
        save_vec(pg, &vp2[VLEN * (ivs + 3)], yt4);
      }
    }
  }
}


//====================================================================
void BridgeQXS::mult_domainwall_5din_L_inv_dirac(
  real_t *__restrict vp,
  real_t *__restrict wp,
  int Ns, int *Nsize,
  real_t *e, real_t *dpinv, real_t *dm)
{
  BridgeQXS::mult_domainwall_5din_eo_L_inv_dirac(vp, wp, Ns, Nsize,
                                                 e, dpinv, dm);
}


//====================================================================
void BridgeQXS::mult_domainwall_5din_U_inv_dirac(
  real_t *__restrict vp,
  real_t *__restrict wp,
  int Ns, int *Nsize,
  real_t *f, real_t *dpinv, real_t *dm)
{
  BridgeQXS::mult_domainwall_5din_eo_U_inv_dirac(vp, wp, Ns, Nsize,
                                                 f, dpinv, dm);
}


//====================================================================
void BridgeQXS::mult_domainwall_5din_Ldag_inv_dirac(
  real_t *__restrict vp,
  real_t *__restrict wp,
  int Ns, int *Nsize,
  real_t *e, real_t *dpinv, real_t *dm)
{
  BridgeQXS::mult_domainwall_5din_eo_Ldag_inv_dirac(vp, wp, Ns, Nsize,
                                                    e, dpinv, dm);
}


//====================================================================
void BridgeQXS::mult_domainwall_5din_Udag_inv_dirac(
  real_t *__restrict vp,
  real_t *__restrict wp,
  int Ns, int *Nsize,
  real_t *f, real_t *dpinv, real_t *dm)
{
  BridgeQXS::mult_domainwall_5din_eo_Udag_inv_dirac(vp, wp, Ns, Nsize,
                                                    f, dpinv, dm);
}


#endif
//============================================================END=====
