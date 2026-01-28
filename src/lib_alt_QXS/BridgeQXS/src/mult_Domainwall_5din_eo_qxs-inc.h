/*!
      @file    mult_Doainwall_5din_eo_qxs-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2492 $
*/

/*
 Copyright
   Bridge++ project and RIKEN (2022)

 Licence: GPL
    see README.txt and LICENSE for more details
*/

#ifndef MULT_DOMAINWALL_5DIN_EO_QXS_INCLUDED
#define MULT_DOMAINWALL_5DIN_EO_QXS_INCLUDED

#include "mult_common_th-inc.h"
#include "prefetch.h"
//====================================================================
void BridgeQXS::mult_domainwall_5din_eo_5dir_dirac(
  real_t *__restrict yp,
  real_t *__restrict wp,
  real_t mq, real_t M0, int Ns, int *bc,
  real_t *__restrict b,
  real_t *__restrict c,
  int *Nsize, int *do_comm)
{
  int Nx2v  = Nsize[0];
  int Ny    = Nsize[1];
  int Nz    = Nsize[2];
  int Nt    = Nsize[3];
  int Nst2v = Nx2v * Ny * Nz * Nt;

  int Nin4 = VLEN * NVCD;
  int Nin5 = Nin4 * Ns;

  int ith, nth, site0, site1;
  set_threadtask(ith, nth, site0, site1, Nst2v);

  for (int site = site0; site < site1; ++site) {
    real_t *wp2 = &wp[Nin5 * site];
    real_t *yp2 = &yp[Nin5 * site];

    svbool_t pg = set_predicate();
    for (int is = 0; is < Ns; ++is) {
      for (int ic = 0; ic < NC; ++ic) {
        svreal_t vt1r, vt1i, vt2r, vt2i, vt3r, vt3i, vt4r, vt4i;
        int      offset = 2 * ND * ic + NVCD * is;
        real_t   factor = -0.5;
        real_t   bb     = b[is] * factor;

        load_vec(pg, vt1r, &wp2[VLEN * (offset + 0)]);
        load_vec(pg, vt1i, &wp2[VLEN * (offset + 1)]);
        load_vec(pg, vt2r, &wp2[VLEN * (offset + 2)]);
        load_vec(pg, vt2i, &wp2[VLEN * (offset + 3)]);
        load_vec(pg, vt3r, &wp2[VLEN * (offset + 4)]);
        load_vec(pg, vt3i, &wp2[VLEN * (offset + 5)]);
        load_vec(pg, vt4r, &wp2[VLEN * (offset + 6)]);
        load_vec(pg, vt4i, &wp2[VLEN * (offset + 7)]);
        scal_vec(pg, vt1r, bb);
        scal_vec(pg, vt1i, bb);
        scal_vec(pg, vt2r, bb);
        scal_vec(pg, vt2i, bb);
        scal_vec(pg, vt3r, bb);
        scal_vec(pg, vt3i, bb);
        scal_vec(pg, vt4r, bb);
        scal_vec(pg, vt4i, bb);

        int    is_up = (is + 1) % Ns;
        real_t Fup   = 0.5 * c[is] * factor;
        if (is == Ns - 1) Fup *= -mq;
        add_aPm5_dirac_vec(vt1r, vt1i, vt2r, vt2i,
                           vt3r, vt3i, vt4r, vt4i,
                           Fup, wp2, is_up, ic);

        int    is_dn = (is - 1 + Ns) % Ns;
        real_t Fdn   = 0.5 * c[is] * factor;
        if (is == 0) Fdn *= -mq;
        add_aPp5_dirac_vec(vt1r, vt1i, vt2r, vt2i,
                           vt3r, vt3i, vt4r, vt4i,
                           Fdn, wp2, is_dn, ic);
        save_vec(pg, &yp2[VLEN * (offset + 0)], vt1r);
        save_vec(pg, &yp2[VLEN * (offset + 1)], vt1i);
        save_vec(pg, &yp2[VLEN * (offset + 2)], vt2r);
        save_vec(pg, &yp2[VLEN * (offset + 3)], vt2i);
        save_vec(pg, &yp2[VLEN * (offset + 4)], vt3r);
        save_vec(pg, &yp2[VLEN * (offset + 5)], vt3i);
        save_vec(pg, &yp2[VLEN * (offset + 6)], vt4r);
        save_vec(pg, &yp2[VLEN * (offset + 7)], vt4i);
      }
    } // is
  }   // site
}


//====================================================================
void BridgeQXS::mult_domainwall_5din_eo_mult_gm5_dirac(
  real_t *__restrict vp,
  real_t *__restrict wp,
  int Ns, int *Nsize)
{
  int Nxv  = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Ny * Nz * Nt;

  int Nin4 = VLEN * NVCD;
  int Nin5 = Nin4 * Ns;

  int ith, nth, site0, site1;
  set_threadtask(ith, nth, site0, site1, Nstv);

  for (int site = site0; site < site1; ++site) {
    real_t   *vp2 = &vp[Nin5 * site];
    real_t   *wp2 = &wp[Nin5 * site];
    svbool_t pg   = set_predicate();
    for (int is = 0; is < Ns; ++is) {
      for (int ic = 0; ic < NC; ++ic) {
        svreal_t vt1r, vt1i, vt2r, vt2i;
        svreal_t vt3r, vt3i, vt4r, vt4i;
        int      offset = 2 * ND * ic + NVCD * is;

        load_mult_gm5_dirac_vec(pg, vt1r, vt1i, vt2r, vt2i,
                                vt3r, vt3i, vt4r, vt4i, &wp2[VLEN * offset]);
        save_vec(pg, &vp2[VLEN * (offset + 0)], vt1r);
        save_vec(pg, &vp2[VLEN * (offset + 1)], vt1i);
        save_vec(pg, &vp2[VLEN * (offset + 2)], vt2r);
        save_vec(pg, &vp2[VLEN * (offset + 3)], vt2i);
        save_vec(pg, &vp2[VLEN * (offset + 4)], vt3r);
        save_vec(pg, &vp2[VLEN * (offset + 5)], vt3i);
        save_vec(pg, &vp2[VLEN * (offset + 6)], vt4r);
        save_vec(pg, &vp2[VLEN * (offset + 7)], vt4i);
      }
    }
  }
}


//====================================================================
void BridgeQXS::mult_domainwall_5din_eo_clear(real_t *vp, int Ns, int *Nsize)
{
  int Nxv  = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Ny * Nz * Nt;

  int Nin4 = VLEN * NVCD;
  int Nin5 = Nin4 * Ns;

  int ith, nth, site0, site1;
  set_threadtask(ith, nth, site0, site1, Nstv);

  svreal_t y;
  svbool_t pg;
  pg = set_predicate();
  clear_vec(pg, y);
  for (int site = site0; site < site1; ++site) {
    real_t *vp2 = &vp[Nin5 * site];

    for (int is = 0; is < Ns; ++is) {
      for (int i = 0; i < NVCD; ++i) {
        save_vec(pg, &vp2[VLEN * (is * NVCD + i)], y);
      }
    }
  }
}


//====================================================================
void BridgeQXS::mult_domainwall_5din_eo_5dirdag_dirac(
  real_t *__restrict vp,
  real_t *__restrict yp,
  real_t mq, real_t M0, int Ns, int *bc,
  real_t *b, real_t *c,
  int *Nsize, int *do_comm)
{
  int Nx2v  = Nsize[0];
  int Ny    = Nsize[1];
  int Nz    = Nsize[2];
  int Nt    = Nsize[3];
  int Nst2v = Nx2v * Ny * Nz * Nt;
  int Nst2  = Nst2v * VLEN;

  int Nin4 = VLEN * NVCD;
  int Nin5 = Nin4 * Ns;

  int ith, nth, site0, site1;
  set_threadtask(ith, nth, site0, site1, Nst2v);

  for (int site = site0; site < site1; ++site) {
    real_t *vp2 = &vp[Nin5 * site];
    real_t *yp2 = &yp[Nin5 * site];

    svbool_t pg = set_predicate();
    for (int is = 0; is < Ns; ++is) {
      for (int ic = 0; ic < NC; ++ic) {
        svreal_t vt1r, vt1i, vt2r, vt2i, vt3r, vt3i, vt4r, vt4i;
        int      offset = 2 * ND * ic + NVCD * is;
        real_t   bb     = 0.5 * b[is];

        load_vec(pg, vt3r, &yp2[VLEN * (offset + 0)]);
        load_vec(pg, vt3i, &yp2[VLEN * (offset + 1)]);
        load_vec(pg, vt4r, &yp2[VLEN * (offset + 2)]);
        load_vec(pg, vt4i, &yp2[VLEN * (offset + 3)]);
        load_vec(pg, vt1r, &yp2[VLEN * (offset + 4)]);
        load_vec(pg, vt1i, &yp2[VLEN * (offset + 5)]);
        load_vec(pg, vt2r, &yp2[VLEN * (offset + 6)]);
        load_vec(pg, vt2i, &yp2[VLEN * (offset + 7)]);
        scal_vec(pg, vt1r, bb);
        scal_vec(pg, vt1i, bb);
        scal_vec(pg, vt2r, bb);
        scal_vec(pg, vt2i, bb);
        scal_vec(pg, vt3r, bb);
        scal_vec(pg, vt3i, bb);
        scal_vec(pg, vt4r, bb);
        scal_vec(pg, vt4i, bb);

        int    is_up = (is + 1) % Ns;
        real_t Fup   = 0.5 * (-0.5) * c[is_up];
        if (is == Ns - 1) Fup *= -mq;
        add_aPp5_dirac_vec(vt1r, vt1i, vt2r, vt2i,
                           vt3r, vt3i, vt4r, vt4i,
                           Fup, yp2, is_up, ic);

        int    is_dn = (is - 1 + Ns) % Ns;
        real_t Fdn   = -0.5 * (-0.5) * c[is_dn];
        if (is == 0) Fdn *= -mq;
        add_aPm5_dirac_vec(vt1r, vt1i, vt2r, vt2i,
                           vt3r, vt3i, vt4r, vt4i,
                           Fdn, yp2, is_dn, ic);
        save_vec(pg, &vp2[VLEN * (offset + 0)], vt1r);
        save_vec(pg, &vp2[VLEN * (offset + 1)], vt1i);
        save_vec(pg, &vp2[VLEN * (offset + 2)], vt2r);
        save_vec(pg, &vp2[VLEN * (offset + 3)], vt2i);
        save_vec(pg, &vp2[VLEN * (offset + 4)], vt3r);
        save_vec(pg, &vp2[VLEN * (offset + 5)], vt3i);
        save_vec(pg, &vp2[VLEN * (offset + 6)], vt4r);
        save_vec(pg, &vp2[VLEN * (offset + 7)], vt4i);
      } // ic
    }   //is
  }     //site
}


//====================================================================
void BridgeQXS::mult_domainwall_5din_eo_hopb_dirac(
  real_t *vp, real_t *up, real_t *wp,
  real_t mq, real_t M0, int Ns, int *bc,
  real_t *b, real_t *c,
  int *Leo, int *Nsize, int *do_comm,
  const int ieo)
{
  int Nx2v  = Nsize[0];
  int Ny    = Nsize[1];
  int Nz    = Nsize[2];
  int Nt    = Nsize[3];
  int Nst2v = Nx2v * Ny * Nz * Nt;
  int Nst2  = Nst2v * VLEN;

  int Nin4 = VLEN * NVCD;
  int Nin5 = Nin4 * Ns;

  int NvU2 = NDF * Nst2;

  int Nxy2  = Nx2v * Ny;
  int Nxyz2 = Nx2v * Ny * Nz;

  svbool_t pg1e_xp, pg2e_xp, pg3e_xp, pg1e_xm, pg2e_xm, pg3e_xm;
  svbool_t pg1o_xp, pg2o_xp, pg3o_xp, pg1o_xm, pg2o_xm, pg3o_xm;
  svbool_t pg1_yp, pg2_yp, pg1_ym, pg2_ym;
  set_predicate_xp_eo(pg1e_xp, pg2e_xp, pg3e_xp, 0);
  set_predicate_xp_eo(pg1o_xp, pg2o_xp, pg3o_xp, 1);
  set_predicate_xm_eo(pg1e_xm, pg2e_xm, pg3e_xm, 0);
  set_predicate_xm_eo(pg1o_xm, pg2o_xm, pg3o_xm, 1);
  set_predicate_yp(pg1_yp, pg2_yp);
  set_predicate_ym(pg1_ym, pg2_ym);

  int ith, nth, site0, site1;
  set_threadtask(ith, nth, site0, site1, Nst2v);

  real_t bufL[Nin5];

  for (int site = site0; site < site1; ++site) {
    int ix   = site % Nx2v;
    int iyzt = site / Nx2v;
    int iy   = iyzt % Ny;
    int izt  = site / Nxy2;
    int iz   = izt % Nz;
    int it   = izt / Nz;
    int ixy  = ix + Nx2v * iy;
    int ixyz = ixy + Nxy2 * iz;
    int jeo  = (ieo + Leo[VLENY * iyzt]) % 2;

    // index for prefetch
    int ix_p1   = (site + 1) % Nx2v;
    int iyzt_p1 = (site + 1) / Nx2v;
    int iy_p1   = iyzt_p1 % Ny;
    int izt_p1  = (site + 1) / Nxy2;
    int iz_p1   = izt_p1 % Nz;
    int it_p1   = izt_p1 / Nz;
    int ixy_p1  = ix_p1 + Nx2v * iy_p1;
    int ixyz_p1 = ixy_p1 + Nxy2 * iz_p1;
    int jeo_p1  = (ieo + Leo[VLENY * iyzt_p1]) % 2;

    int idir;

    real_t *wp2 = &wp[Nin5 * site];
    real_t *vp2 = &vp[Nin5 * site];

    real_t  z4[VLEN * NVCD];
    Vsimd_t vL[NVCD];

    int nei_p1;

    idir = 0;

    __prefetch_load_hop_u_l2(up, ieo + 2 * idir, site + 1);

    if ((ix < Nx2v - 1) || (do_comm[0] == 0)) {
      int nei = ix + 1 + Nx2v * iyzt;

      nei_p1 = ix - 1 + Nx2v * iyzt;
      if (ix == Nx2v - 1) nei = 0 + Nx2v * iyzt;

      real_t *u2  = &up[VLEN * NDF * site + NvU2 * (ieo + 2 * idir)];
      real_t *wpn = &wp[Nin5 * nei];

      if (jeo == 0) {
        for (int is = 0; is < Ns; ++is) {
          __prefetch_load_hop_vec_l2(wp, site + 1, is);
          __prefetch_write_hop_vec_l2(vp, site + 1, is);
          __prefetch_load_hop_vec_l1(wp, nei_p1, is);

          mult_wilson_eo_xpb(pg1e_xp, pg2e_xp, pg3e_xp,
                             &vp2[Nin4 * is], u2, &wp2[Nin4 * is], &wpn[Nin4 * is]);
        }
      } else {
        for (int is = 0; is < Ns; ++is) {
          __prefetch_load_hop_vec_l2(wp, site + 1, is);
          __prefetch_write_hop_vec_l2(vp, site + 1, is);
          __prefetch_load_hop_vec_l1(wp, nei_p1, is);

          mult_wilson_eo_xpb(pg1o_xp, pg2o_xp, pg3o_xp,
                             &vp2[Nin4 * is], u2, &wp2[Nin4 * is], &wpn[Nin4 * is]);
        }
      }
    }


    __prefetch_load_hop_u_l2(up, 1 - ieo + 2 * idir, site + 1);

    if ((ix > 0) || (do_comm[0] == 0)) {
      int nei = ix - 1 + Nx2v * iyzt;

      int iy2_p1 = (iy + 1) % Ny;
      nei_p1 = ix + Nx2v * (iy2_p1 + Ny * izt);

      if (ix == 0) nei = Nx2v - 1 + Nx2v * iyzt;

      real_t *u2  = &up[VLEN * NDF * site + NvU2 * (1 - ieo + 2 * idir)];
      real_t *un  = &up[VLEN * NDF * nei + NvU2 * (1 - ieo + 2 * idir)];
      real_t *wpn = &wp[Nin5 * nei];
      if (jeo == 0) {
        for (int is = 0; is < Ns; ++is) {
          __prefetch_load_hop_vec_l2(wp, site + 1, is);
          __prefetch_write_hop_vec_l2(vp, site + 1, is);
          __prefetch_load_hop_vec_l1(wp, nei_p1, is);

          mult_wilson_eo_xmb(pg1e_xm, pg2e_xm, pg3e_xm, &vp2[Nin4 * is],
                             u2, un, &wp2[Nin4 * is], &wpn[Nin4 * is]);
        }
      } else {
        for (int is = 0; is < Ns; ++is) {
          __prefetch_load_hop_vec_l2(wp, site + 1, is);
          __prefetch_write_hop_vec_l2(vp, site + 1, is);
          __prefetch_load_hop_vec_l1(wp, nei_p1, is);

          mult_wilson_eo_xmb(pg1o_xm, pg2o_xm, pg3o_xm, &vp2[Nin4 * is],
                             u2, un, &wp2[Nin4 * is], &wpn[Nin4 * is]);
        }
      }
    }


    idir = 1;

    __prefetch_load_hop_u_l2(up, ieo + 2 * idir, site + 1);

    if ((iy < Ny - 1) || (do_comm[idir] == 0)) {
      int iy2 = (iy + 1) % Ny;
      int nei = ix + Nx2v * (iy2 + Ny * izt);

      int iy2_p1 = (iy - 1 + Ny) % Ny;
      nei_p1 = ix + Nx2v * (iy2_p1 + Ny * izt);


      real_t *wpn = &wp[Nin5 * nei];
      real_t *u2  = &up[VLEN * NDF * site + NvU2 * (ieo + 2 * idir)];
      for (int is = 0; is < Ns; ++is) {
        __prefetch_load_hop_vec_l1(wp, nei_p1, is);

        mult_wilson_ypb(pg1_yp, pg2_yp,
                        &vp2[Nin4 * is], u2, &wp2[Nin4 * is], &wpn[Nin4 * is]);
      }
    }

    __prefetch_load_hop_u_l2(up, 1 - ieo + 2 * idir, site + 1);

    if ((iy > 0) || (do_comm[idir] == 0)) {
      int iy2 = (iy - 1 + Ny) % Ny;
      int nei = ix + Nx2v * (iy2 + Ny * izt);

      int iz2_p1 = (iz + 1) % Nz;
      nei_p1 = ixy + Nxy2 * (iz2_p1 + Nz * it);

      real_t *wpn = &wp[Nin5 * nei];
      real_t *u2  = &up[VLEN * NDF * site + NvU2 * (1 - ieo + 2 * idir)];
      real_t *un  = &up[VLEN * NDF * nei + NvU2 * (1 - ieo + 2 * idir)];
      for (int is = 0; is < Ns; ++is) {
        __prefetch_load_hop_vec_l1(wp, nei_p1, is);
        mult_wilson_ymb(pg1_ym, pg2_ym,
                        &vp2[Nin4 * is], u2, un, &wp2[Nin4 * is], &wpn[Nin4 * is]);
      }
    }

    idir = 2;

    __prefetch_load_hop_u_l2(up, ieo + 2 * idir, site + 1);

    if ((iz < Nz - 1) || (do_comm[idir] == 0)) {
      int iz2 = (iz + 1) % Nz;
      int nei = ixy + Nxy2 * (iz2 + Nz * it);

      int iz2_p1 = (iz - 1 + Nz) % Nz;
      nei_p1 = ixy + Nxy2 * (iz2_p1 + Nz * it);

      real_t *wpn = &wp[Nin5 * nei];
      real_t *u2  = &up[VLEN * NDF * site + NvU2 * (ieo + 2 * idir)];
      for (int is = 0; is < Ns; ++is) {
        __prefetch_load_hop_vec_l1(wp, nei_p1, is);

        mult_wilson_zpb(&vp2[Nin4 * is], u2, &wpn[Nin4 * is]);
      }
    }

    __prefetch_load_hop_u_l2(up, 1 - ieo + 2 * idir, site + 1);

    if ((iz > 0) || (do_comm[idir] == 0)) {
      int iz2 = (iz - 1 + Nz) % Nz;
      int nei = ixy + Nxy2 * (iz2 + Nz * it);

      int it2_p1 = (it + 1) % Nt;
      nei_p1 = ixyz + Nxyz2 * it2_p1;

      real_t *wpn = &wp[Nin5 * nei];
      real_t *u2  = &up[VLEN * NDF * nei + NvU2 * (1 - ieo + 2 * idir)];
      for (int is = 0; is < Ns; ++is) {
        __prefetch_load_hop_vec_l1(wp, nei_p1, is);

        mult_wilson_zmb(&vp2[Nin4 * is], u2, &wpn[Nin4 * is]);
      }
    }

    idir = 3;

    __prefetch_load_hop_u_l2(up, ieo + 2 * idir, site + 1);

    if ((it < Nt - 1) || (do_comm[idir] == 0)) {
      int it2 = (it + 1) % Nt;
      int nei = ixyz + Nxyz2 * it2;

      int it2_p1 = (it - 1 + Nt) % Nt;
      nei_p1 = ixyz + Nxyz2 * it2_p1;

      real_t *wpn = &wp[Nin5 * nei];
      real_t *u2  = &up[VLEN * NDF * site + NvU2 * (ieo + 2 * idir)];
      for (int is = 0; is < Ns; ++is) {
        __prefetch_load_hop_vec_l1(wp, nei_p1, is);

        mult_wilson_tpb_dirac(&vp2[Nin4 * is], u2, &wpn[Nin4 * is]);
      }
    }

    __prefetch_load_hop_u_l2(up, 1 - ieo + 2 * idir, site + 1);

    if ((it > 0) || (do_comm[idir] == 0)) {
      int it2 = (it - 1 + Nt) % Nt;
      int nei = ixyz + Nxyz2 * it2;

      nei_p1 = ix_p1 + 1 + Nx2v * iyzt_p1;

      real_t *wpn = &wp[Nin5 * nei];
      real_t *u2  = &up[VLEN * NDF * nei + NvU2 * (1 - ieo + 2 * idir)];
      for (int is = 0; is < Ns; ++is) {
        __prefetch_load_hop_vec_l1(wp, nei_p1, is);

        mult_wilson_tmb_dirac(&vp2[Nin4 * is], u2, &wpn[Nin4 * is]);
      }
    }
  }
}


//====================================================================
void BridgeQXS::mult_domainwall_5din_eo_hop1_dirac(
  real_t *buf1_xp, real_t *buf1_xm,
  real_t *buf1_yp, real_t *buf1_ym,
  real_t *buf1_zp, real_t *buf1_zm,
  real_t *buf1_tp, real_t *buf1_tm,
  real_t *up, real_t *wp,
  real_t mq, real_t M0, int Ns, int *bc,
  int *Leo, int *Nsize, int *do_comm,
  const int ieo)
{
  int Nx2v  = Nsize[0];
  int Ny    = Nsize[1];
  int Nz    = Nsize[2];
  int Nt    = Nsize[3];
  int Nst2v = Nx2v * Ny * Nz * Nt;
  int Nst2  = Nst2v * VLEN;

  int Nin4  = VLEN * NVCD;
  int Nin5  = Nin4 * Ns;
  int Nin4H = VLEN * NVC * ND2;
  int Nin5H = Nin4H * Ns;
  int NvU2  = NDF * Nst2;

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
  svint_t svidx_xp, svidx_xm;
  set_index_xp_eo(svidx_xp);
  set_index_xm_eo(svidx_xm);

  int idir = 0;
  if (do_comm[idir] == 1) {
    int Nyzt = Ny * Nz * Nt;

    int ith, nth, site0, site1;
    set_threadtask(ith, nth, site0, site1, Nyzt);

    for (int iyzt = site0; iyzt < site1; ++iyzt) {
      int iy = iyzt % Ny;
      int iz = (iyzt / Ny) % Nz;
      int it = iyzt / (Ny * Nz);

      int jeo = (ieo + Leo[VLENY * iyzt]) % 2;

      int Nskipx = (VLENY + 1) / 2;

      // index for prefetch
#if VLENY > 1
      int ibf_up_xp1   = Nskipx * NVC * ND2 * Ns * (iyzt + 1);
      int ibf_dn_xp1   = Nskipx * NVC * ND2 * Ns * (iyzt + 1);
      int ibf_dn_xp1_2 = Nskipx * NVC * ND2 * Ns * iyzt;
      int ibf_up_xp1_2 = Nskipx * NVC * ND2 * Ns * (iyzt + 1);
#else
      int ibf_up_xp1   = Nskipx * NVC * ND2 * Ns * ((iyzt + 1) / 2);
      int ibf_dn_xp1   = Nskipx * NVC * ND2 * Ns * ((iyzt + 1) / 2);
      int ibf_dn_xp1_2 = Nskipx * NVC * ND2 * Ns * (iyzt / 2);
      int ibf_up_xp1_2 = Nskipx * NVC * ND2 * Ns * ((iyzt + 1) / 2);
#endif

      int i_p1 = 1;
      int i_p2 = 2;

#if VLENY > 1
      int ibf_up = Nskipx * NVC * ND2 * Ns * iyzt;
      int ibf_dn = Nskipx * NVC * ND2 * Ns * iyzt;
#else
      int ibf_up = Nskipx * NVC * ND2 * Ns * (iyzt / 2);
      int ibf_dn = Nskipx * NVC * ND2 * Ns * (iyzt / 2);
#endif

      {
        int    ix   = 0;
        int    site = ix + Nx2v * iyzt;
        real_t *wp2 = &wp[Nin5 * site];

        // index for prefetch
        int site_xp1   = ix + Nx2v * (iyzt + 1);
        int site_xp1_2 = Nx2v - 1 + Nx2v * (iyzt + 1);

        set_index_xm_eo(svidx_xm);
        if (jeo == 0) {
#if VLENY > 1
          for (int is = 0; is < Ns; ++is) {
            __prefetch_load_hop_vec_l2(wp, site_xp1, is);
            __prefetch_write_hop1_buf_x_l2(buf1_xp, ibf_up_xp1, is,
                                           Nskipx * NVC * ND2);

            if ((!((it == 0) || (it == Nt - 1))) &&
                (!((iz == 0) || (iz == Nz - 1))) &&
                (!((iy == 0) || (iy == Ny - 1)))) {
              __prefetch_load_hop_vec_l1(wp, site_xp1_2, is);
              __prefetch_write_hop1_buf_x_l1(buf1_xm, ibf_dn_xp1_2, is,
                                             Nskipx * NVC * ND2);
            }

            mult_wilson_eo_xp1(pg2o_xm, svidx_xm,
                               &buf1_xp[ibf_up + Nskipx * NVC * ND2 * is],
                               &wp2[Nin4 * is]);
          }
#endif
        } else {
          for (int is = 0; is < Ns; ++is) {
            __prefetch_load_hop_vec_l2(wp, site_xp1, is);
            __prefetch_write_hop1_buf_x_l2(buf1_xp, ibf_up_xp1, is,
                                           Nskipx * NVC * ND2);

            if ((!((it == 0) || (it == Nt - 1))) &&
                (!((iz == 0) || (iz == Nz - 1))) &&
                (!((iy == 0) || (iy == Ny - 1)))) {
              __prefetch_load_hop_vec_l1(wp, site_xp1_2, is);
              __prefetch_write_hop1_buf_x_l1(buf1_xm, ibf_dn_xp1_2, is,
                                             Nskipx * NVC * ND2);
            }

            mult_wilson_eo_xp1(pg2e_xm, svidx_xm,
                               &buf1_xp[ibf_up + Nskipx * NVC * ND2 * is],
                               &wp2[Nin4 * is]);
          }
        }
      }

      {
        int    ix   = Nx2v - 1;
        int    site = ix + Nx2v * iyzt;
        real_t *wp2 = &wp[Nin5 * site];
        real_t *u2  = &up[VLEN * NDF * site + NvU2 * (1 - ieo + 2 * idir)];

        int site_xp1   = ix + Nx2v * (iyzt + 1);
        int site_xp1_2 = 0 + Nx2v * (iyzt + 1);

        __prefetch_load_hop_u_l2(up, 1 - ieo + 2 * 0, site_xp1);
        // __prefetch_load_hop_u_l1(up,1-ieo+2*0,site_xp1);

        set_index_xp_eo(svidx_xp);
        if (jeo == 0) {
          for (int is = 0; is < Ns; ++is) {
            __prefetch_load_hop_vec_l2(wp, site_xp1, is);
            __prefetch_write_hop1_buf_x_l2(buf1_xm, ibf_dn_xp1, is,
                                           Nskipx * NVC * ND2);

            if ((!((it == 0) || (it == Nt - 1))) &&
                (!((iz == 0) || (iz == Nz - 1))) &&
                (!((iy == 0) || (iy == Ny - 1)))) {
              __prefetch_load_hop_vec_l1(wp, site_xp1_2, is);
              __prefetch_write_hop1_buf_x_l1(buf1_xp, ibf_up_xp1_2, is,
                                             Nskipx * NVC * ND2);
            }

            mult_wilson_eo_xm1(pg2o_xp, svidx_xp,
                               &buf1_xm[ibf_dn + Nskipx * NVC * ND2 * is],
                               u2, &wp2[Nin4 * is]);
          }
        } else {
#if VLENY > 1
          for (int is = 0; is < Ns; ++is) {
            __prefetch_load_hop_vec_l2(wp, site_xp1, is);
            __prefetch_write_hop1_buf_x_l2(buf1_xm, ibf_dn_xp1, is,
                                           Nskipx * NVC * ND2);

            if ((!((it == 0) || (it == Nt - 1))) &&
                (!((iz == 0) || (iz == Nz - 1))) &&
                (!((iy == 0) || (iy == Ny - 1)))) {
              __prefetch_load_hop_vec_l1(wp, site_xp1_2, is);
              __prefetch_write_hop1_buf_x_l1(buf1_xp, ibf_up_xp1_2, is,
                                             Nskipx * NVC * ND2);
            }

            mult_wilson_eo_xm1(pg2e_xp, svidx_xp,
                               &buf1_xm[ibf_dn + Nskipx * NVC * ND2 * is],
                               u2, &wp2[Nin4 * is]);
          }
#endif
        }
      }
    } // iyzt loop
  }

  idir = 1;
  if (do_comm[idir] == 1) {
    int Nxzt2 = Nx2v * Nz * Nt;

    int ith, nth, site0, site1;
    set_threadtask(ith, nth, site0, site1, Nxzt2);

    for (int ixzt = site0; ixzt < site1; ++ixzt) {
      int ix  = ixzt % Nx2v;
      int izt = ixzt / Nx2v;
      int iz  = izt % Nz;
      int it  = izt / Nz;

      // index for prefetch
      int i_p1 = 1;
      int i_p2 = 2;

      int ibf_yp1 = VLENX * NVC * ND2 * Ns *
                    ((ix + i_p1) % Nx2v + Nx2v * (izt + (ix + i_p1) / Nx2v));
      int ibf_yp2 = VLENX * NVC * ND2 * Ns *
                    ((ix + i_p2) % Nx2v + Nx2v * (izt + (ix + i_p2) / Nx2v));

      {
        int    iy   = 0;
        int    site = ix + Nx2v * (iy + Ny * izt);
        int    ibf  = VLENX * NVC * ND2 * Ns * ixzt;
        real_t *wp2 = &wp[Nin5 * site];

        int iy_p1    = ((iy + (ix + i_p1) / Nx2v) % Ny) * (Ny - 1);
        int izt_p1   = izt + (iy + (ix + i_p1) / Nx2v) / Ny;
        int site_yp1 = (ix + i_p1) % Nx2v + Nx2v * (izt_p1 * Ny + iy_p1);
        int iy_p2    = ((iy + (ix + i_p2) / Nx2v) % Ny) * (Ny - 1);
        int izt_p2   = izt + (iy + (ix + i_p2) / Nx2v) / Ny;
        int site_yp2 = (ix + i_p2) % Nx2v + Nx2v * (izt_p2 * Ny + iy_p2);

        for (int is = 0; is < Ns; ++is) {
          __prefetch_load_hop_vec_l2(wp, site_yp2, is);
          __prefetch_write_hop1_buf_y_l2(buf1_yp, ibf_yp2, is, VLENX * NVC * ND2);

          if ((!((it == 0) || (it == Nt - 1))) &&
              (!((iz == 0) || (iz == Nz - 1)))) {
            __prefetch_load_hop_vec_l1(wp, site_yp1, is);
            __prefetch_write_hop1_buf_y_l1(buf1_yp, ibf_yp1, is,
                                           VLENX * NVC * ND2);
          }
          mult_wilson_yp1(pg2_ym,
                          &buf1_yp[ibf + VLENX * NVC * ND2 * is], &wp2[Nin4 * is]);
        }
      }
      {
        int    iy   = Ny - 1;
        int    site = ix + Nx2v * (iy + Ny * izt);
        int    ibf  = VLENX * NVC * ND2 * Ns * ixzt;
        real_t *wp2 = &wp[Nin5 * site];
        real_t *u2  = &up[VLEN * NDF * site + NvU2 * (1 - ieo + 2 * idir)];

        int iy_p1    = ((iy + (ix + i_p1) / Nx2v) % Ny) * (Ny - 1);
        int izt_p1   = izt + (iy + (ix + i_p1) / Nx2v) / Ny;
        int site_yp1 = (ix + i_p1) % Nx2v + Nx2v * (izt_p1 * Ny + iy_p1);
        int iy_p2    = ((iy + (ix + i_p2) / Nx2v) % Ny) * (Ny - 1);
        int izt_p2   = izt + (iy + (ix + i_p2) / Nx2v) / Ny;
        int site_yp2 = (ix + i_p2) % Nx2v + Nx2v * (izt_p2 * Ny + iy_p2);

        __prefetch_load_hop_u_l2(up, 1 - ieo + 2 * idir, site_yp2);

        for (int is = 0; is < Ns; ++is) {
          __prefetch_load_hop_vec_l2(wp, site_yp2, is);
          __prefetch_write_hop1_buf_y_l2(buf1_ym, ibf_yp2, is, VLENX * NVC * ND2);

          if ((!((it == 0) || (it == Nt - 1))) &&
              (!((iz == 0) || (iz == Nz - 1)))) {
            __prefetch_load_hop_vec_l1(wp, site_yp1, is);
            __prefetch_write_hop1_buf_y_l1(buf1_ym, ibf_yp1, is,
                                           VLENX * NVC * ND2);
          }

          mult_wilson_ym1(pg2_yp,
                          &buf1_ym[ibf + VLENX * NVC * ND2 * is], u2, &wp2[Nin4 * is]);
        }
      }
    }
  }

  idir = 2;
  if (do_comm[idir] == 1) {
    int Nxyt2 = Nxy2 * Nt;

    int ith, nth, site0, site1;
    set_threadtask(ith, nth, site0, site1, Nxyt2);

    for (int ixyt = site0; ixyt < site1; ++ixyt) {
      int ixy = ixyt % Nxy2;
      int it  = ixyt / Nxy2;

      int i_p1    = 1;
      int i_p2    = 2;
      int ibf_zp1 = Nin5H * ((ixy + i_p1) % Nxy2 + Nxy2 * (it + (ixy + i_p1) / Nxy2));
      int ibf_zp2 = Nin5H * ((ixy + i_p2) % Nxy2 + Nxy2 * (it + (ixy + i_p2) / Nxy2));

      int ibf = Nin5H * (ixy + Nxy2 * it);

      {
        int    iz   = 0;
        int    site = ixy + Nxy2 * (iz + Nz * it);
        real_t *wp2 = &wp[Nin5 * site];

        int iz_p1    = ((iz + (ixy + i_p1) / Nxy2) % Nz) * (Nz - 1);
        int it_p1    = it + (iz + (ixy + i_p1) / Nxy2) / Nz;
        int site_zp1 = (ixy + i_p1) % Nxy2 + Nxy2 * (it_p1 * Nz + iz_p1);
        int iz_p2    = ((iz + (ixy + 2) / Nxy2) % Nz) * (Nz - 1);
        int it_p2    = it + (iz + (ixy + i_p2) / Nxy2) / Nz;
        int site_zp2 = (ixy + 2) % Nxy2 + Nxy2 * (it_p2 * Nz + iz_p2);

        for (int is = 0; is < Ns; ++is) {
          __prefetch_load_hop_vec_l2(wp, site_zp2, is);
          __prefetch_write_hop1_buf_zt_l2(buf1_zp, ibf_zp2, is, Nin4H);

          if (!((it == 0) || (it == Nt - 1))) {
            __prefetch_load_hop_vec_l1(wp, site_zp1, is);
            __prefetch_write_hop1_buf_zt_l1(buf1_zp, ibf_zp1, is, Nin4H);
          }

          mult_wilson_zp1(&buf1_zp[ibf + Nin4H * is], &wp2[Nin4 * is]);
        }
      }

      {
        int    iz   = Nz - 1;
        int    site = ixy + Nxy2 * (iz + Nz * it);
        real_t *wp2 = &wp[Nin5 * site];
        real_t *u2  = &up[VLEN * NDF * site + NvU2 * (1 - ieo + 2 * idir)];

        int iz_p1    = ((iz + (ixy + i_p1) / Nxy2) % Nz) * (Nz - 1);
        int it_p1    = it + (iz + (ixy + i_p1) / Nxy2) / Nz;
        int site_zp1 = (ixy + i_p1) % Nxy2 + Nxy2 * (it_p1 * Nz + iz_p1);
        int iz_p2    = ((iz + (ixy + 2) / Nxy2) % Nz) * (Nz - 1);
        int it_p2    = it + (iz + (ixy + i_p2) / Nxy2) / Nz;
        int site_zp2 = (ixy + 2) % Nxy2 + Nxy2 * (it_p2 * Nz + iz_p2);

        __prefetch_load_hop_u_l2(up, 1 - ieo + 2 * idir, site_zp2);

        for (int is = 0; is < Ns; ++is) {
          __prefetch_load_hop_vec_l2(wp, site_zp2, is);
          __prefetch_write_hop1_buf_zt_l2(buf1_zm, ibf_zp2, is, Nin4H);

          if (!((it == 0) || (it == Nt - 1))) {
            __prefetch_load_hop_vec_l1(wp, site_zp1, is);
            __prefetch_write_hop1_buf_zt_l1(buf1_zm, ibf_zp1, is, Nin4H);
          }

          mult_wilson_zm1(&buf1_zm[ibf + Nin4H * is], u2, &wp2[Nin4 * is]);
        }
      }
    }
  }

  idir = 3;
  if (do_comm[idir] == 1) {
    int ith, nth, site0, site1;
    set_threadtask(ith, nth, site0, site1, Nxyz2);

    for (int ixyz = site0; ixyz < site1; ++ixyz) {
      int ibf = Nin5H * ixyz;

      int i_p1 = 1;
      int i_p2 = 2;

      // index for prefetch
      int ibf_tp1 = Nin5H * (ixyz + i_p1);
      int ibf_tp2 = Nin5H * (ixyz + i_p2);

      {
        int    it   = 0;
        int    site = ixyz + Nxyz2 * it;
        real_t *wp2 = &wp[Nin5 * site];

        // index for prefetch
        int site_tp1 = ixyz + i_p1 + it * Nxy2 * Nz;
        int site_tp2 = ixyz + i_p2 + it * Nxy2 * Nz;

        for (int is = 0; is < Ns; ++is) {
          __prefetch_load_hop_vec_l2(wp, site_tp2, is);
          __prefetch_write_hop1_buf_zt_l2(buf1_tp, ibf_tp2, is, Nin4H);

          __prefetch_load_hop_vec_l1(wp, site_tp1, is);
          __prefetch_write_hop1_buf_zt_l1(buf1_tp, ibf_tp1, is, Nin4H);

          mult_wilson_tp1_dirac(&buf1_tp[ibf + Nin4H * is], &wp2[Nin4 * is]);
        }
      }

      {
        int    it   = Nt - 1;
        int    site = ixyz + Nxyz2 * it;
        real_t *wp2 = &wp[Nin5 * site];

        // index for prefetch
        int site_tp1 = ixyz + i_p1 + it * Nxy2 * Nz;
        int site_tp2 = ixyz + i_p2 + it * Nxy2 * Nz;

        __prefetch_load_hop_u_l2(up, 1 - ieo + 2 * idir, site_tp2);

        real_t *u2 = &up[VLEN * NDF * site + NvU2 * (1 - ieo + 2 * idir)];
        for (int is = 0; is < Ns; ++is) {
          __prefetch_load_hop_vec_l2(wp, site_tp2, is);
          __prefetch_write_hop1_buf_zt_l2(buf1_tm, ibf_tp2, is, Nin4H);

          __prefetch_load_hop_vec_l1(wp, site_tp1, is);
          __prefetch_write_hop1_buf_zt_l1(buf1_tm, ibf_tp1, is, Nin4H);

          mult_wilson_tm1_dirac(&buf1_tm[ibf + Nin4H * is], u2, &wp2[Nin4 * is]);
        }
      }
    }
  }
}


//====================================================================
void BridgeQXS::mult_domainwall_5din_eo_hop2_dirac(
  real_t *vp, real_t *up, real_t *wp,
  real_t *buf2_xp, real_t *buf2_xm,
  real_t *buf2_yp, real_t *buf2_ym,
  real_t *buf2_zp, real_t *buf2_zm,
  real_t *buf2_tp, real_t *buf2_tm,
  real_t mq, real_t M0, int Ns, int *bc,
  int *Leo, int *Nsize, int *do_comm,
  const int ieo)
{
  int Nx2v  = Nsize[0];
  int Ny    = Nsize[1];
  int Nz    = Nsize[2];
  int Nt    = Nsize[3];
  int Nst2v = Nx2v * Ny * Nz * Nt;
  int Nst2  = Nst2v * VLEN;

  int Nin4  = VLEN * NVCD;
  int Nin5  = Nin4 * Ns;
  int Nin4H = VLEN * NVC * ND2;
  int Nin5H = Nin4H * Ns;
  int NvU2  = NDF * Nst2;

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
  svint_t svidx_xp, svidx_xm;
  set_index_xp_eo(svidx_xp);
  set_index_xm_eo(svidx_xm);

  int Nskipx = (VLENY + 1) / 2;

  int ith, nth, site0, site1;
  set_threadtask(ith, nth, site0, site1, Nst2v);

  real_t bufL[Nin5];

  /* #pragma loop noprefetch */
  /* #pragma fj loop prefetch_sequential soft */
  for (int site = site0; site < site1; ++site) {
    int ix   = site % Nx2v;
    int iyzt = site / Nx2v;
    int iy   = iyzt % Ny;
    int izt  = site / Nxy2;
    int iz   = izt % Nz;
    int it   = izt / Nz;
    int ixy  = ix + Nx2v * iy;
    int ixyz = ixy + Nxy2 * iz;
    int jeo  = (ieo + Leo[VLENY * iyzt]) % 2;

    // index for prefetch
    int ibf_up_xp1, ibf_dn_xp1, site_xp1, site_xp1_2, ibf_up_xp1_2, ibf_dn_xp1_2;

    ibf_up_xp1 = Nskipx * NVC * ND2 * Ns * (iyzt + 1);
    if (VLENY == 1) {
      ibf_up_xp1 = Nskipx * NVC * ND2 * Ns * ((iyzt + 1) / 2);
    }
    ibf_dn_xp1 = Nskipx * NVC * ND2 * Ns * (iyzt + 1);
    if (VLENY == 1) {
      ibf_dn_xp1 = Nskipx * NVC * ND2 * Ns * ((iyzt + 1) / 2);
    }
    site_xp1 = (iyzt + 1) * Nx2v + ix;

    if (ix == 0) {
      site_xp1_2 = iyzt * Nx2v + Nx2v - 1;

      ibf_dn_xp1_2 = Nskipx * NVC * ND2 * Ns * iyzt;
      if (VLENY == 1) {
        ibf_dn_xp1_2 = Nskipx * NVC * ND2 * Ns * (iyzt / 2);
      }
    } else if (ix == Nx2v - 1) {
      site_xp1_2 = (iyzt + 1) * Nx2v + 0;

      ibf_up_xp1_2 = Nskipx * NVC * ND2 * Ns * (iyzt + 1);
      if (VLENY == 1) {
        ibf_up_xp1_2 = Nskipx * NVC * ND2 * Ns * ((iyzt + 1) / 2);
      }
    }

    int i_p1 = 1;
    int i_p2 = 2;

    int ibf_yp1  = VLENX * NVC * ND2 * Ns * ((ix + i_p1) % Nx2v + Nx2v * (izt + (ix + i_p1) / Nx2v));
    int iy_p1    = ((iy + (ix + i_p1) / Nx2v) % Ny) * (Ny - 1);
    int izt_p1   = izt + (iy + (ix + i_p1) / Nx2v) / Ny;
    int site_yp1 = (ix + i_p1) % Nx2v + Nx2v * (izt_p1 * Ny + iy_p1);

    int ibf_yp2  = VLENX * NVC * ND2 * Ns * ((ix + i_p2) % Nx2v + Nx2v * (izt + (ix + i_p2) / Nx2v));
    int iy_p2    = ((iy + (ix + i_p2) / Nx2v) % Ny) * (Ny - 1);
    int izt_p2   = izt + (iy + (ix + i_p2) / Nx2v) / Ny;
    int site_yp2 = (ix + i_p2) % Nx2v + Nx2v * (izt_p2 * Ny + iy_p2);

    int ibf_zp1  = Nin5H * ((ixy + i_p1) % Nxy2 + Nxy2 * (it + (ixy + i_p1) / Nxy2));
    int iz_p1    = ((iz + (ixy + i_p1) / Nxy2) % Nz) * (Nz - 1);
    int it_p1    = it + (iz + (ixy + i_p1) / Nxy2) / Nz;
    int site_zp1 = (ixy + i_p1) % Nxy2 + Nxy2 * (it_p1 * Nz + iz_p1);

    int ibf_zp2  = Nin5H * ((ixy + i_p2) % Nxy2 + Nxy2 * (it + (ixy + i_p2) / Nxy2));
    int iz_p2    = ((iz + (ixy + 2) / Nxy2) % Nz) * (Nz - 1);
    int it_p2    = it + (iz + (ixy + i_p2) / Nxy2) / Nz;
    int site_zp2 = (ixy + 2) % Nxy2 + Nxy2 * (it_p2 * Nz + iz_p2);

    int ibf_tp1  = Nin5H * (ixyz + i_p1);
    int site_tp1 = ixyz + i_p1 + it * Nxy2 * Nz;

    int ibf_tp2  = Nin5H * (ixyz + i_p2);
    int site_tp2 = ixyz + i_p2 + it * Nxy2 * Nz;
    //

    int idir, nei;

    real_t *wp2 = &wp[Nin5 * site];
    real_t *vp2 = &vp[Nin5 * site];

    idir = 0;
    if (do_comm[idir] == 1) {
      if (ix == Nx2v - 1) {
        int ibf_up = Nskipx * NVC * ND2 * Ns * iyzt;
        if (VLENY == 1) {
          ibf_up = Nskipx * NVC * ND2 * Ns * (iyzt / 2);
        }


        __prefetch_load_hop_u_l2(up, ieo + 2 * 0, site_xp1);
        real_t *u = &up[NDF * Nst2 * (ieo + 2 * 0)];
        set_index_xp_eo(svidx_xp);
        if (jeo == 0) {
          for (int is = 0; is < Ns; ++is) {
            __prefetch_load_hop_vec_l2(wp, site_xp1, is);
            __prefetch_write_hop_vec_l2(vp, site_xp1, is);
            __prefetch_load_hop2_buf_x_l2(buf2_xp, ibf_up_xp1, is, Nskipx * NVC * ND2);


            if ((!((it == 0) || (it == Nt - 1))) &&
                (!((iz == 0) || (iz == Nz - 1))) &&
                (!((iy == 0) || (iy == Ny - 1)))) {
              __prefetch_load_hop_vec_l1(wp, site_xp1_2, is);
              __prefetch_write_hop_vec_l1(vp, site_xp1_2, is);
              __prefetch_load_hop2_buf_x_l1(buf2_xm, ibf_dn_xp1_2, is, Nskipx * NVC * ND2);
            }

            mult_wilson_eo_xp2(pg1e_xp, pg2e_xp, pg3e_xp, svidx_xp,
                               &vp2[Nin4 * is], &u[VLEN * NDF * site], &wp2[Nin4 * is],
                               &buf2_xp[ibf_up + Nskipx * NVC * ND2 * is]);
          }
        } else {
          for (int is = 0; is < Ns; ++is) {
            __prefetch_load_hop_vec_l2(wp, site_xp1, is);
            __prefetch_write_hop_vec_l2(vp, site_xp1, is);
            __prefetch_load_hop2_buf_x_l2(buf2_xp, ibf_up_xp1, is, Nskipx * NVC * ND2);

            if ((!((it == 0) || (it == Nt - 1))) &&
                (!((iz == 0) || (iz == Nz - 1))) &&
                (!((iy == 0) || (iy == Ny - 1)))) {
              __prefetch_load_hop_vec_l1(wp, site_xp1_2, is);
              __prefetch_write_hop_vec_l1(vp, site_xp1_2, is);
              __prefetch_load_hop2_buf_x_l1(buf2_xm, ibf_dn_xp1_2, is, Nskipx * NVC * ND2);
            }

            mult_wilson_eo_xp2(pg1o_xp, pg2o_xp, pg3o_xp, svidx_xp,
                               &vp2[Nin4 * is], &u[VLEN * NDF * site], &wp2[Nin4 * is],
                               &buf2_xp[ibf_up + Nskipx * NVC * ND2 * is]);
          }
        }
      }

      if (ix == 0) {
        int ibf_dn = Nskipx * NVC * ND2 * Ns * iyzt;
        if (VLENY == 1) {
          ibf_dn = Nskipx * NVC * ND2 * Ns * (iyzt / 2);
        }


        __prefetch_load_hop_u_l2(up, 1 - ieo + 2 * 0, site_xp1);
        real_t *u = &up[NDF * Nst2 * (1 - ieo + 2 * 0)];

        set_index_xm_eo(svidx_xm);
        if (jeo == 0) {
          for (int is = 0; is < Ns; ++is) {
            __prefetch_load_hop_vec_l2(wp, site_xp1, is);
            __prefetch_write_hop_vec_l2(vp, site_xp1, is);
            __prefetch_load_hop2_buf_x_l2(buf2_xm, ibf_dn_xp1, is, Nskipx * NVC * ND2);

            if ((!((it == 0) || (it == Nt - 1))) &&
                (!((iz == 0) || (iz == Nz - 1))) &&
                (!((iy == 0) || (iy == Ny - 1)))) {
              __prefetch_load_hop_vec_l1(wp, site_xp1_2, is);
              __prefetch_write_hop_vec_l1(vp, site_xp1_2, is);
              __prefetch_write_hop1_buf_x_l1(buf2_xp, ibf_up_xp1_2, is, Nskipx * NVC * ND2)
            }

            mult_wilson_eo_xm2(pg1e_xm, pg2e_xm, pg3e_xm, svidx_xm,
                               &vp2[Nin4 * is], &u[VLEN * NDF * site], &wp2[Nin4 * is],
                               &buf2_xm[ibf_dn + Nskipx * NVC * ND2 * is]);
          }
        } else {
          for (int is = 0; is < Ns; ++is) {
            __prefetch_load_hop_vec_l2(wp, site_xp1, is);
            __prefetch_write_hop_vec_l2(vp, site_xp1, is);
            __prefetch_load_hop2_buf_x_l2(buf2_xm, ibf_dn_xp1, is, Nskipx * NVC * ND2);

            if ((!((it == 0) || (it == Nt - 1))) &&
                (!((iz == 0) || (iz == Nz - 1))) &&
                (!((iy == 0) || (iy == Ny - 1)))) {
              __prefetch_load_hop_vec_l1(wp, site_xp1_2, is);
              __prefetch_write_hop_vec_l1(vp, site_xp1_2, is);
              __prefetch_write_hop1_buf_x_l1(buf2_xp, ibf_up_xp1_2, is, Nskipx * NVC * ND2);
            }

            mult_wilson_eo_xm2(pg1o_xm, pg2o_xm, pg3o_xm, svidx_xm,
                               &vp2[Nin4 * is], &u[VLEN * NDF * site], &wp2[Nin4 * is],
                               &buf2_xm[ibf_dn + Nskipx * NVC * ND2 * is]);
          }
        }
      }
    }

    idir = 1;
    if (do_comm[idir] == 1) {
      if (iy == Ny - 1) {
        int ibf = VLENX * NVC * ND2 * Ns * (ix + Nx2v * izt);


        __prefetch_load_hop_u_l2(up, ieo + 2 * 1, site_yp2);
        real_t *u = &up[NDF * Nst2 * (ieo + 2 * 1)];
        for (int is = 0; is < Ns; ++is) {
          __prefetch_load_hop_vec_l2(wp, site_yp2, is);
          __prefetch_write_hop_vec_l2(vp, site_yp2, is);
          __prefetch_load_hop2_buf_y_l2(buf2_yp, ibf_yp2, is, VLENX * NVC * ND2);

          if ((!((it == 0) || (it == Nt - 1))) &&
              (!((iz == 0) || (iz == Nz - 1)))) {
            __prefetch_load_hop_vec_l1(wp, site_yp1, is);
            __prefetch_write_hop_vec_l1(vp, site_yp1, is);
            __prefetch_write_hop1_buf_y_l1(buf2_yp, ibf_yp1, is, VLENX * NVC * ND2);
          }

          mult_wilson_yp2(pg1_yp, pg2_yp,
                          &vp2[Nin4 * is], &u[VLEN * NDF * site],
                          &wp2[Nin4 * is], &buf2_yp[ibf + VLENX * NVC * ND2 * is]);
        }
      }

      if (iy == 0) {
        int ibf = VLENX * NVC * ND2 * Ns * (ix + Nx2v * izt);

        __prefetch_load_hop_u_l2(up, 1 - ieo + 2 * 1, site_yp2);
        real_t *u = &up[NDF * Nst2 * (1 - ieo + 2 * 1)];
        for (int is = 0; is < Ns; ++is) {
          __prefetch_load_hop_vec_l2(wp, site_yp2, is);
          __prefetch_write_hop_vec_l2(vp, site_yp2, is);
          __prefetch_load_hop2_buf_y_l2(buf2_ym, ibf_yp2, is, VLENX * NVC * ND2);

          if ((!((it == 0) || (it == Nt - 1))) &&
              (!((iz == 0) || (iz == Nz - 1)))) {
            __prefetch_load_hop_vec_l1(wp, site_yp1, is);
            __prefetch_write_hop_vec_l1(vp, site_yp1, is);

            __prefetch_write_hop1_buf_y_l1(buf2_ym, ibf_yp1, is, VLENX * NVC * ND2);
          }

          mult_wilson_ym2(pg1_ym, pg2_ym,
                          &vp2[Nin4 * is], &u[VLEN * NDF * site],
                          &wp2[Nin4 * is], &buf2_ym[ibf + VLENX * NVC * ND2 * is]);
        }
      }
    }

    idir = 2;
    if (do_comm[idir] == 1) {
      if (iz == Nz - 1) {
        int ibf = Nin5H * (ixy + Nxy2 * it);

        __prefetch_load_hop_u_l2(up, ieo + 2 * idir, site_zp2);
        real_t *u2 = &up[VLEN * NDF * site + NvU2 * (ieo + 2 * idir)];
        for (int is = 0; is < Ns; ++is) {
          __prefetch_write_hop_vec_l2(vp, site_zp2, is);
          __prefetch_load_hop2_buf_zt_l2(buf2_zp, ibf_zp2, is, Nin4H);

          if (!((it == 0) || (it == Nt - 1))) {
            __prefetch_write_hop_vec_l1(vp, site_zp1, is);
            __prefetch_write_hop1_buf_zt_l1(buf2_zp, ibf_zp1, is, Nin4H);
          }

          mult_wilson_zp2(&vp2[Nin4 * is], u2, &buf2_zp[ibf + Nin4H * is]);
        }
      }

      if (iz == 0) {
        int ibf = Nin5H * (ixy + Nxy2 * it);

        for (int is = 0; is < Ns; ++is) {
          __prefetch_write_hop_vec_l2(vp, site_zp2, is);
          __prefetch_load_hop2_buf_zt_l2(buf2_zm, ibf_zp2, is, Nin4H);

          if (!((it == 0) || (it == Nt - 1))) {
            __prefetch_write_hop_vec_l1(vp, site_zp1, is);
            __prefetch_write_hop1_buf_zt_l1(buf2_zm, ibf_zp1, is, Nin4H);
          }

          mult_wilson_zm2(&vp2[Nin4 * is], &buf2_zm[ibf + Nin4H * is]);
        }
      }
    }

    idir = 3;
    if (do_comm[idir] == 1) {
      if (it == Nt - 1) {
        int ibf = Nin5H * ixyz;

        __prefetch_load_hop_u_l2(up, ieo + 2 * idir, site_tp2);

        real_t *u2 = &up[VLEN * NDF * site + NvU2 * (ieo + 2 * idir)];
        for (int is = 0; is < Ns; ++is) {
          __prefetch_write_hop_vec_l2(vp, site_tp2, is);
          __prefetch_load_hop2_buf_zt_l2(buf2_tp, ibf_tp2, is, Nin4H);

          __prefetch_write_hop_vec_l1(vp, site_tp1, is);
          __prefetch_load_hop2_buf_zt_l1(buf2_tp, ibf_tp1, is, Nin4H);

          mult_wilson_tp2_dirac(&vp2[Nin4 * is], u2, &buf2_tp[ibf + Nin4H * is]);
        }
      }

      if (it == 0) {
        int ibf = Nin5H * ixyz;


        for (int is = 0; is < Ns; ++is) {
          __prefetch_write_hop_vec_l2(vp, site_tp2, is);
          __prefetch_load_hop2_buf_zt_l2(buf2_tp, ibf_tp2, is, Nin4H);

          __prefetch_write_hop_vec_l1(vp, site_tp1, is);
          __prefetch_load_hop2_buf_zt_l1(buf2_tp, ibf_tp1, is, Nin4H);

          mult_wilson_tm2_dirac(&vp2[Nin4 * is], &buf2_tm[ibf + Nin4H * is]);
        }
      }
    }
  }
}


//====================================================================
void BridgeQXS::mult_domainwall_5din_eo_bulk_dirac(
  real_t *vp, real_t *up, real_t *wp,
  real_t *yp,
  real_t mq, real_t M0, int Ns, int *bc,
  real_t *b, real_t *c,
  int *Leo, int *Nsize, int *do_comm,
  const int ieo)
{
  mult_domainwall_5din_eo_clear(vp, Ns, Nsize);

  mult_domainwall_5din_eo_5dir_dirac(yp, wp, mq, M0, Ns, bc, b, c,
                                     Nsize, do_comm);

  mult_domainwall_5din_eo_hopb_dirac(vp, up, yp, mq, M0, Ns, bc, b, c,
                                     Leo, Nsize, do_comm, ieo);
}


//====================================================================
void BridgeQXS::mult_domainwall_5din_eo_L_inv_dirac(
  real_t *__restrict vp,
  real_t *__restrict wp,
  int Ns, int *Nsize,
  real_t *e, real_t *dpinv, real_t *dm)
{
  int Nx2v  = Nsize[0];
  int Ny    = Nsize[1];
  int Nz    = Nsize[2];
  int Nt    = Nsize[3];
  int Nst2v = Nx2v * Ny * Nz * Nt;
  int Nst2  = Nst2v * VLEN;

  int Nin4 = VLEN * NVCD;
  int Nin5 = Nin4 * Ns;


  int ith, nth, site0, site1;
  set_threadtask(ith, nth, site0, site1, Nst2v);

  for (int site = site0; site < site1; ++site) {
    real_t   *vp2 = &vp[Nin5 * site];
    real_t   *wp2 = &wp[Nin5 * site];
    svbool_t pg   = set_predicate();

    for (int ic = 0; ic < NC; ++ic) {
      svreal_t x1r, x1i, x2r, x2i, x3r, x3i, x4r, x4i;
      svreal_t v1r, v1i, v2r, v2i, v3r, v3i, v4r, v4i;
      svreal_t y1r, y1i, y2r, y2i, y3r, y3i, y4r, y4i;
      int      offset0 = 2 * ND * ic;

      __prefetch_load_luinv_l1(wp, offset0);
      __prefetch_write_luinv_l1(vp, offset0);

      load_vec(pg, v1r, &wp2[VLEN * (offset0 + 0)]);
      load_vec(pg, v1i, &wp2[VLEN * (offset0 + 1)]);
      load_vec(pg, v2r, &wp2[VLEN * (offset0 + 2)]);
      load_vec(pg, v2i, &wp2[VLEN * (offset0 + 3)]);
      load_vec(pg, v3r, &wp2[VLEN * (offset0 + 4)]);
      load_vec(pg, v3i, &wp2[VLEN * (offset0 + 5)]);
      load_vec(pg, v4r, &wp2[VLEN * (offset0 + 6)]);
      load_vec(pg, v4i, &wp2[VLEN * (offset0 + 7)]);
      save_vec(pg, &vp2[VLEN * (offset0 + 0)], v1r);
      save_vec(pg, &vp2[VLEN * (offset0 + 1)], v1i);
      save_vec(pg, &vp2[VLEN * (offset0 + 2)], v2r);
      save_vec(pg, &vp2[VLEN * (offset0 + 3)], v2i);
      save_vec(pg, &vp2[VLEN * (offset0 + 4)], v3r);
      save_vec(pg, &vp2[VLEN * (offset0 + 5)], v3i);
      save_vec(pg, &vp2[VLEN * (offset0 + 6)], v4r);
      save_vec(pg, &vp2[VLEN * (offset0 + 7)], v4i);
      y1r = v1r;
      scal_vec(pg, y1r, e[0]);
      y1i = v1i;
      scal_vec(pg, y1i, e[0]);
      y2r = v2r;
      scal_vec(pg, y2r, e[0]);
      y2i = v2i;
      scal_vec(pg, y2i, e[0]);
      y3r = v3r;
      scal_vec(pg, y3r, e[0]);
      y3i = v3i;
      scal_vec(pg, y3i, e[0]);
      y4r = v4r;
      scal_vec(pg, y4r, e[0]);
      y4i = v4i;
      scal_vec(pg, y4i, e[0]);

      for (int is = 1; is < Ns - 1; ++is) {
        x1r = v1r;
        x1i = v1i;
        x2r = v2r;
        x2i = v2i;
        x3r = v3r;
        x3i = v3i;
        x4r = v4r;
        x4i = v4i;

        int offset = 2 * ND * ic + NVCD * is;


        __prefetch_load_luinv_l1(wp, offset);
        __prefetch_write_luinv_l1(vp, offset);

        load_vec(pg, v1r, &wp2[VLEN * (offset + 0)]);
        load_vec(pg, v1i, &wp2[VLEN * (offset + 1)]);
        load_vec(pg, v2r, &wp2[VLEN * (offset + 2)]);
        load_vec(pg, v2i, &wp2[VLEN * (offset + 3)]);
        load_vec(pg, v3r, &wp2[VLEN * (offset + 4)]);
        load_vec(pg, v3i, &wp2[VLEN * (offset + 5)]);
        load_vec(pg, v4r, &wp2[VLEN * (offset + 6)]);
        load_vec(pg, v4i, &wp2[VLEN * (offset + 7)]);

        real_t a = real_t(0.5) * dm[is] * dpinv[is - 1];

        add_aPp5_dirac_vec(pg,
                           v1r, v1i, v2r, v2i, v3r, v3i, v4r, v4i,
                           a,
                           x1r, x1i, x2r, x2i, x3r, x3i, x4r, x4i);
        save_vec(pg, &vp2[VLEN * (offset + 0)], v1r);
        save_vec(pg, &vp2[VLEN * (offset + 1)], v1i);
        save_vec(pg, &vp2[VLEN * (offset + 2)], v2r);
        save_vec(pg, &vp2[VLEN * (offset + 3)], v2i);
        save_vec(pg, &vp2[VLEN * (offset + 4)], v3r);
        save_vec(pg, &vp2[VLEN * (offset + 5)], v3i);
        save_vec(pg, &vp2[VLEN * (offset + 6)], v4r);
        save_vec(pg, &vp2[VLEN * (offset + 7)], v4i);
        axpy_vec(pg, y1r, e[is], v1r);
        axpy_vec(pg, y1i, e[is], v1i);
        axpy_vec(pg, y2r, e[is], v2r);
        axpy_vec(pg, y2i, e[is], v2i);
        axpy_vec(pg, y3r, e[is], v3r);
        axpy_vec(pg, y3i, e[is], v3i);
        axpy_vec(pg, y4r, e[is], v4r);
        axpy_vec(pg, y4i, e[is], v4i);
      }
      int is = Ns - 1;
      x1r = v1r;
      x1i = v1i;
      x2r = v2r;
      x2i = v2i;
      x3r = v3r;
      x3i = v3i;
      x4r = v4r;
      x4i = v4i;

      int offset = 2 * ND * ic + NVCD * is;
      load_vec(pg, v1r, &wp2[VLEN * (offset + 0)]);
      load_vec(pg, v1i, &wp2[VLEN * (offset + 1)]);
      load_vec(pg, v2r, &wp2[VLEN * (offset + 2)]);
      load_vec(pg, v2i, &wp2[VLEN * (offset + 3)]);
      load_vec(pg, v3r, &wp2[VLEN * (offset + 4)]);
      load_vec(pg, v3i, &wp2[VLEN * (offset + 5)]);
      load_vec(pg, v4r, &wp2[VLEN * (offset + 6)]);
      load_vec(pg, v4i, &wp2[VLEN * (offset + 7)]);
      real_t a = real_t(0.5) * dm[is] * dpinv[is - 1];
      add_aPp5_dirac_vec(pg,
                         v1r, v1i, v2r, v2i, v3r, v3i, v4r, v4i,
                         a,
                         x1r, x1i, x2r, x2i, x3r, x3i, x4r, x4i);

      add_aPm5_dirac_vec(pg,
                         v1r, v1i, v2r, v2i, v3r, v3i, v4r, v4i,
                         real_t(-0.5),
                         y1r, y1i, y2r, y2i, y3r, y3i, y4r, y4i);
      save_vec(pg, &vp2[VLEN * (offset + 0)], v1r);
      save_vec(pg, &vp2[VLEN * (offset + 1)], v1i);
      save_vec(pg, &vp2[VLEN * (offset + 2)], v2r);
      save_vec(pg, &vp2[VLEN * (offset + 3)], v2i);
      save_vec(pg, &vp2[VLEN * (offset + 4)], v3r);
      save_vec(pg, &vp2[VLEN * (offset + 5)], v3i);
      save_vec(pg, &vp2[VLEN * (offset + 6)], v4r);
      save_vec(pg, &vp2[VLEN * (offset + 7)], v4i);
    } //ic
  }   //site

#pragma omp barrier
}


//====================================================================
void BridgeQXS::mult_domainwall_5din_eo_U_inv_dirac(
  real_t *__restrict vp,
  real_t *__restrict wp,
  int Ns, int *Nsize,
  real_t *f, real_t *dpinv, real_t *dm)
{
  int Nx2v  = Nsize[0];
  int Ny    = Nsize[1];
  int Nz    = Nsize[2];
  int Nt    = Nsize[3];
  int Nst2v = Nx2v * Ny * Nz * Nt;
  int Nst2  = Nst2v * VLEN;

  int Nin4 = VLEN * NVCD;
  int Nin5 = Nin4 * Ns;


  int ith, nth, site0, site1;
  set_threadtask(ith, nth, site0, site1, Nst2v);

  for (int site = site0; site < site1; ++site) {
    real_t   *vp2 = &vp[Nin5 * site];
    real_t   *wp2 = &wp[Nin5 * site];
    svbool_t pg   = set_predicate();

    for (int ic = 0; ic < NC; ++ic) {
      svreal_t x1r, x1i, x2r, x2i, x3r, x3i, x4r, x4i;
      svreal_t v1r, v1i, v2r, v2i, v3r, v3i, v4r, v4i;
      svreal_t y1r, y1i, y2r, y2i, y3r, y3i, y4r, y4i;

      int offset0 = 2 * ND * ic + NVCD * (Ns - 1);

      __prefetch_load_luinv_l1(wp, offset0);
      __prefetch_write_luinv_l1(vp, offset0);

      load_vec(pg, v1r, &wp2[VLEN * (offset0 + 0)]);
      load_vec(pg, v1i, &wp2[VLEN * (offset0 + 1)]);
      load_vec(pg, v2r, &wp2[VLEN * (offset0 + 2)]);
      load_vec(pg, v2i, &wp2[VLEN * (offset0 + 3)]);
      load_vec(pg, v3r, &wp2[VLEN * (offset0 + 4)]);
      load_vec(pg, v3i, &wp2[VLEN * (offset0 + 5)]);
      load_vec(pg, v4r, &wp2[VLEN * (offset0 + 6)]);
      load_vec(pg, v4i, &wp2[VLEN * (offset0 + 7)]);

      real_t a = dpinv[Ns - 1];

      scal_vec(pg, v1r, a);
      scal_vec(pg, v1i, a);
      scal_vec(pg, v2r, a);
      scal_vec(pg, v2i, a);
      scal_vec(pg, v3r, a);
      scal_vec(pg, v3i, a);
      scal_vec(pg, v4r, a);
      scal_vec(pg, v4i, a);

      save_vec(pg, &vp2[VLEN * (offset0 + 0)], v1r);
      save_vec(pg, &vp2[VLEN * (offset0 + 1)], v1i);
      save_vec(pg, &vp2[VLEN * (offset0 + 2)], v2r);
      save_vec(pg, &vp2[VLEN * (offset0 + 3)], v2i);
      save_vec(pg, &vp2[VLEN * (offset0 + 4)], v3r);
      save_vec(pg, &vp2[VLEN * (offset0 + 5)], v3i);
      save_vec(pg, &vp2[VLEN * (offset0 + 6)], v4r);
      save_vec(pg, &vp2[VLEN * (offset0 + 7)], v4i);

      set_aPp5_dirac_vec(pg,
                         y1r, y1i, y2r, y2i, y3r, y3i, y4r, y4i,
                         real_t(0.5),
                         v1r, v1i, v2r, v2i, v3r, v3i, v4r, v4i);
      for (int is = Ns - 2; is >= 0; --is) {
        x1r = v1r;
        x1i = v1i;
        x2r = v2r;
        x2i = v2i;
        x3r = v3r;
        x3i = v3i;
        x4r = v4r;
        x4i = v4i;

        int offset = 2 * ND * ic + NVCD * is;

        __prefetch_load_luinv_l1(wp, offset);
        __prefetch_write_luinv_l1(vp, offset);

        load_vec(pg, v1r, &wp2[VLEN * (offset + 0)]);
        load_vec(pg, v1i, &wp2[VLEN * (offset + 1)]);
        load_vec(pg, v2r, &wp2[VLEN * (offset + 2)]);
        load_vec(pg, v2i, &wp2[VLEN * (offset + 3)]);
        load_vec(pg, v3r, &wp2[VLEN * (offset + 4)]);
        load_vec(pg, v3i, &wp2[VLEN * (offset + 5)]);
        load_vec(pg, v4r, &wp2[VLEN * (offset + 6)]);
        load_vec(pg, v4i, &wp2[VLEN * (offset + 7)]);

        real_t a = real_t(0.5) * dm[is];

        add_aPm5_dirac_vec(pg,
                           v1r, v1i, v2r, v2i, v3r, v3i, v4r, v4i,
                           a,
                           x1r, x1i, x2r, x2i, x3r, x3i, x4r, x4i);

        axpy_vec(pg, v1r, -f[is], y1r);
        axpy_vec(pg, v1i, -f[is], y1i);
        axpy_vec(pg, v2r, -f[is], y2r);
        axpy_vec(pg, v2i, -f[is], y2i);
        axpy_vec(pg, v3r, -f[is], y3r);
        axpy_vec(pg, v3i, -f[is], y3i);
        axpy_vec(pg, v4r, -f[is], y4r);
        axpy_vec(pg, v4i, -f[is], y4i);

        real_t aa = dpinv[is];
        scal_vec(pg, v1r, aa);
        scal_vec(pg, v1i, aa);
        scal_vec(pg, v2r, aa);
        scal_vec(pg, v2i, aa);
        scal_vec(pg, v3r, aa);
        scal_vec(pg, v3i, aa);
        scal_vec(pg, v4r, aa);
        scal_vec(pg, v4i, aa);

        save_vec(pg, &vp2[VLEN * (offset + 0)], v1r);
        save_vec(pg, &vp2[VLEN * (offset + 1)], v1i);
        save_vec(pg, &vp2[VLEN * (offset + 2)], v2r);
        save_vec(pg, &vp2[VLEN * (offset + 3)], v2i);
        save_vec(pg, &vp2[VLEN * (offset + 4)], v3r);
        save_vec(pg, &vp2[VLEN * (offset + 5)], v3i);
        save_vec(pg, &vp2[VLEN * (offset + 6)], v4r);
        save_vec(pg, &vp2[VLEN * (offset + 7)], v4i);
      }
    } //ic
  }   //site
#pragma omp barrier
}


//====================================================================
void BridgeQXS::mult_domainwall_5din_eo_Udag_inv_dirac(
  real_t *__restrict vp,
  real_t *__restrict wp,
  int Ns, int *Nsize,
  real_t *f, real_t *dpinv, real_t *dm)
{
  int Nx2v  = Nsize[0];
  int Ny    = Nsize[1];
  int Nz    = Nsize[2];
  int Nt    = Nsize[3];
  int Nst2v = Nx2v * Ny * Nz * Nt;
  int Nst2  = Nst2v * VLEN;

  int Nin4 = VLEN * NVCD;
  int Nin5 = Nin4 * Ns;


  int ith, nth, site0, site1;
  set_threadtask(ith, nth, site0, site1, Nst2v);

  for (int site = site0; site < site1; ++site) {
    real_t   *vp2 = &vp[Nin5 * site];
    real_t   *wp2 = &wp[Nin5 * site];
    svbool_t pg   = set_predicate();

    for (int ic = 0; ic < NC; ++ic) {
      svreal_t x1r, x1i, x2r, x2i, x3r, x3i, x4r, x4i;
      svreal_t v1r, v1i, v2r, v2i, v3r, v3i, v4r, v4i;
      svreal_t y1r, y1i, y2r, y2i, y3r, y3i, y4r, y4i;
      int      offset0 = 2 * ND * ic;

      __prefetch_load_luinv_l1(wp, offset0);
      __prefetch_write_luinv_l1(vp, offset0);

      load_vec(pg, v1r, &wp2[VLEN * (offset0 + 0)]);
      load_vec(pg, v1i, &wp2[VLEN * (offset0 + 1)]);
      load_vec(pg, v2r, &wp2[VLEN * (offset0 + 2)]);
      load_vec(pg, v2i, &wp2[VLEN * (offset0 + 3)]);
      load_vec(pg, v3r, &wp2[VLEN * (offset0 + 4)]);
      load_vec(pg, v3i, &wp2[VLEN * (offset0 + 5)]);
      load_vec(pg, v4r, &wp2[VLEN * (offset0 + 6)]);
      load_vec(pg, v4i, &wp2[VLEN * (offset0 + 7)]);

      real_t a = dpinv[0];
      scal_vec(pg, v1r, a);
      scal_vec(pg, v1i, a);
      scal_vec(pg, v2r, a);
      scal_vec(pg, v2i, a);
      scal_vec(pg, v3r, a);
      scal_vec(pg, v3i, a);
      scal_vec(pg, v4r, a);
      scal_vec(pg, v4i, a);

      save_vec(pg, &vp2[VLEN * (offset0 + 0)], v1r);
      save_vec(pg, &vp2[VLEN * (offset0 + 1)], v1i);
      save_vec(pg, &vp2[VLEN * (offset0 + 2)], v2r);
      save_vec(pg, &vp2[VLEN * (offset0 + 3)], v2i);
      save_vec(pg, &vp2[VLEN * (offset0 + 4)], v3r);
      save_vec(pg, &vp2[VLEN * (offset0 + 5)], v3i);
      save_vec(pg, &vp2[VLEN * (offset0 + 6)], v4r);
      save_vec(pg, &vp2[VLEN * (offset0 + 7)], v4i);
      y1r = v1r;
      scal_vec(pg, y1r, f[0]);
      y1i = v1i;
      scal_vec(pg, y1i, f[0]);
      y2r = v2r;
      scal_vec(pg, y2r, f[0]);
      y2i = v2i;
      scal_vec(pg, y2i, f[0]);
      y3r = v3r;
      scal_vec(pg, y3r, f[0]);
      y3i = v3i;
      scal_vec(pg, y3i, f[0]);
      y4r = v4r;
      scal_vec(pg, y4r, f[0]);
      y4i = v4i;
      scal_vec(pg, y4i, f[0]);

      for (int is = 1; is < Ns - 1; ++is) {
        x1r = v1r;
        x1i = v1i;
        x2r = v2r;
        x2i = v2i;
        x3r = v3r;
        x3i = v3i;
        x4r = v4r;
        x4i = v4i;

        int offset = 2 * ND * ic + NVCD * is;

        __prefetch_load_luinv_l1(wp, offset);
        __prefetch_write_luinv_l1(vp, offset);

        load_vec(pg, v1r, &wp2[VLEN * (offset + 0)]);
        load_vec(pg, v1i, &wp2[VLEN * (offset + 1)]);
        load_vec(pg, v2r, &wp2[VLEN * (offset + 2)]);
        load_vec(pg, v2i, &wp2[VLEN * (offset + 3)]);
        load_vec(pg, v3r, &wp2[VLEN * (offset + 4)]);
        load_vec(pg, v3i, &wp2[VLEN * (offset + 5)]);
        load_vec(pg, v4r, &wp2[VLEN * (offset + 6)]);
        load_vec(pg, v4i, &wp2[VLEN * (offset + 7)]);

        real_t a = real_t(0.5) * dm[is - 1];

        add_aPm5_dirac_vec(pg,
                           v1r, v1i, v2r, v2i, v3r, v3i, v4r, v4i,
                           a,
                           x1r, x1i, x2r, x2i, x3r, x3i, x4r, x4i);

        real_t aa = dpinv[is];
        scal_vec(pg, v1r, aa);
        scal_vec(pg, v1i, aa);
        scal_vec(pg, v2r, aa);
        scal_vec(pg, v2i, aa);
        scal_vec(pg, v3r, aa);
        scal_vec(pg, v3i, aa);
        scal_vec(pg, v4r, aa);
        scal_vec(pg, v4i, aa);

        save_vec(pg, &vp2[VLEN * (offset + 0)], v1r);
        save_vec(pg, &vp2[VLEN * (offset + 1)], v1i);
        save_vec(pg, &vp2[VLEN * (offset + 2)], v2r);
        save_vec(pg, &vp2[VLEN * (offset + 3)], v2i);
        save_vec(pg, &vp2[VLEN * (offset + 4)], v3r);
        save_vec(pg, &vp2[VLEN * (offset + 5)], v3i);
        save_vec(pg, &vp2[VLEN * (offset + 6)], v4r);
        save_vec(pg, &vp2[VLEN * (offset + 7)], v4i);
        axpy_vec(pg, y1r, f[is], v1r);
        axpy_vec(pg, y1i, f[is], v1i);
        axpy_vec(pg, y2r, f[is], v2r);
        axpy_vec(pg, y2i, f[is], v2i);
        axpy_vec(pg, y3r, f[is], v3r);
        axpy_vec(pg, y3i, f[is], v3i);
        axpy_vec(pg, y4r, f[is], v4r);
        axpy_vec(pg, y4i, f[is], v4i);
      }
      int is = Ns - 1;
      x1r = v1r;
      x1i = v1i;
      x2r = v2r;
      x2i = v2i;
      x3r = v3r;
      x3i = v3i;
      x4r = v4r;
      x4i = v4i;

      int offset = 2 * ND * ic + NVCD * is;
      load_vec(pg, v1r, &wp2[VLEN * (offset + 0)]);
      load_vec(pg, v1i, &wp2[VLEN * (offset + 1)]);
      load_vec(pg, v2r, &wp2[VLEN * (offset + 2)]);
      load_vec(pg, v2i, &wp2[VLEN * (offset + 3)]);
      load_vec(pg, v3r, &wp2[VLEN * (offset + 4)]);
      load_vec(pg, v3i, &wp2[VLEN * (offset + 5)]);
      load_vec(pg, v4r, &wp2[VLEN * (offset + 6)]);
      load_vec(pg, v4i, &wp2[VLEN * (offset + 7)]);
      a = real_t(0.5) * dm[is - 1];
      add_aPm5_dirac_vec(pg,
                         v1r, v1i, v2r, v2i, v3r, v3i, v4r, v4i,
                         a,
                         x1r, x1i, x2r, x2i, x3r, x3i, x4r, x4i);
      add_aPp5_dirac_vec(pg,
                         v1r, v1i, v2r, v2i, v3r, v3i, v4r, v4i,
                         real_t(-0.5),
                         y1r, y1i, y2r, y2i, y3r, y3i, y4r, y4i);
      real_t aa = dpinv[is];
      scal_vec(pg, v1r, aa);
      scal_vec(pg, v1i, aa);
      scal_vec(pg, v2r, aa);
      scal_vec(pg, v2i, aa);
      scal_vec(pg, v3r, aa);
      scal_vec(pg, v3i, aa);
      scal_vec(pg, v4r, aa);
      scal_vec(pg, v4i, aa);

      save_vec(pg, &vp2[VLEN * (offset + 0)], v1r);
      save_vec(pg, &vp2[VLEN * (offset + 1)], v1i);
      save_vec(pg, &vp2[VLEN * (offset + 2)], v2r);
      save_vec(pg, &vp2[VLEN * (offset + 3)], v2i);
      save_vec(pg, &vp2[VLEN * (offset + 4)], v3r);
      save_vec(pg, &vp2[VLEN * (offset + 5)], v3i);
      save_vec(pg, &vp2[VLEN * (offset + 6)], v4r);
      save_vec(pg, &vp2[VLEN * (offset + 7)], v4i);
    } //ic
  }   //site
#pragma omp barrier
}


//====================================================================
void BridgeQXS::mult_domainwall_5din_eo_Ldag_inv_dirac(
  real_t *__restrict vp,
  real_t *__restrict wp,
  int Ns, int *Nsize,
  real_t *e, real_t *dpinv, real_t *dm)
{
  int Nx2v  = Nsize[0];
  int Ny    = Nsize[1];
  int Nz    = Nsize[2];
  int Nt    = Nsize[3];
  int Nst2v = Nx2v * Ny * Nz * Nt;
  int Nst2  = Nst2v * VLEN;

  int Nin4 = VLEN * NVCD;
  int Nin5 = Nin4 * Ns;


  int ith, nth, site0, site1;
  set_threadtask(ith, nth, site0, site1, Nst2v);

  for (int site = site0; site < site1; ++site) {
    real_t   *vp2 = &vp[Nin5 * site];
    real_t   *wp2 = &wp[Nin5 * site];
    svbool_t pg   = set_predicate();

    for (int ic = 0; ic < NC; ++ic) {
      svreal_t x1r, x1i, x2r, x2i, x3r, x3i, x4r, x4i;
      svreal_t v1r, v1i, v2r, v2i, v3r, v3i, v4r, v4i;
      svreal_t y1r, y1i, y2r, y2i, y3r, y3i, y4r, y4i;

      int offset0 = 2 * ND * ic + NVCD * (Ns - 1);

      __prefetch_load_luinv_l1(wp, offset0);
      __prefetch_write_luinv_l1(vp, offset0);

      load_vec(pg, v1r, &wp2[VLEN * (offset0 + 0)]);
      load_vec(pg, v1i, &wp2[VLEN * (offset0 + 1)]);
      load_vec(pg, v2r, &wp2[VLEN * (offset0 + 2)]);
      load_vec(pg, v2i, &wp2[VLEN * (offset0 + 3)]);
      load_vec(pg, v3r, &wp2[VLEN * (offset0 + 4)]);
      load_vec(pg, v3i, &wp2[VLEN * (offset0 + 5)]);
      load_vec(pg, v4r, &wp2[VLEN * (offset0 + 6)]);
      load_vec(pg, v4i, &wp2[VLEN * (offset0 + 7)]);

      save_vec(pg, &vp2[VLEN * (offset0 + 0)], v1r);
      save_vec(pg, &vp2[VLEN * (offset0 + 1)], v1i);
      save_vec(pg, &vp2[VLEN * (offset0 + 2)], v2r);
      save_vec(pg, &vp2[VLEN * (offset0 + 3)], v2i);
      save_vec(pg, &vp2[VLEN * (offset0 + 4)], v3r);
      save_vec(pg, &vp2[VLEN * (offset0 + 5)], v3i);
      save_vec(pg, &vp2[VLEN * (offset0 + 6)], v4r);
      save_vec(pg, &vp2[VLEN * (offset0 + 7)], v4i);

      set_aPm5_dirac_vec(pg,
                         y1r, y1i, y2r, y2i, y3r, y3i, y4r, y4i,
                         real_t(0.5),
                         v1r, v1i, v2r, v2i, v3r, v3i, v4r, v4i);

      for (int is = Ns - 2; is >= 0; --is) {
        x1r = v1r;
        x1i = v1i;
        x2r = v2r;
        x2i = v2i;
        x3r = v3r;
        x3i = v3i;
        x4r = v4r;
        x4i = v4i;

        int offset = 2 * ND * ic + NVCD * is;

        __prefetch_load_luinv_l1(wp, offset);
        __prefetch_write_luinv_l1(vp, offset);

        load_vec(pg, v1r, &wp2[VLEN * (offset + 0)]);
        load_vec(pg, v1i, &wp2[VLEN * (offset + 1)]);
        load_vec(pg, v2r, &wp2[VLEN * (offset + 2)]);
        load_vec(pg, v2i, &wp2[VLEN * (offset + 3)]);
        load_vec(pg, v3r, &wp2[VLEN * (offset + 4)]);
        load_vec(pg, v3i, &wp2[VLEN * (offset + 5)]);
        load_vec(pg, v4r, &wp2[VLEN * (offset + 6)]);
        load_vec(pg, v4i, &wp2[VLEN * (offset + 7)]);

        real_t a = real_t(0.5) * dm[is + 1] * dpinv[is];

        add_aPp5_dirac_vec(pg,
                           v1r, v1i, v2r, v2i, v3r, v3i, v4r, v4i,
                           a,
                           x1r, x1i, x2r, x2i, x3r, x3i, x4r, x4i);

        axpy_vec(pg, v1r, -e[is], y1r);
        axpy_vec(pg, v1i, -e[is], y1i);
        axpy_vec(pg, v2r, -e[is], y2r);
        axpy_vec(pg, v2i, -e[is], y2i);
        axpy_vec(pg, v3r, -e[is], y3r);
        axpy_vec(pg, v3i, -e[is], y3i);
        axpy_vec(pg, v4r, -e[is], y4r);
        axpy_vec(pg, v4i, -e[is], y4i);

        save_vec(pg, &vp2[VLEN * (offset + 0)], v1r);
        save_vec(pg, &vp2[VLEN * (offset + 1)], v1i);
        save_vec(pg, &vp2[VLEN * (offset + 2)], v2r);
        save_vec(pg, &vp2[VLEN * (offset + 3)], v2i);
        save_vec(pg, &vp2[VLEN * (offset + 4)], v3r);
        save_vec(pg, &vp2[VLEN * (offset + 5)], v3i);
        save_vec(pg, &vp2[VLEN * (offset + 6)], v4r);
        save_vec(pg, &vp2[VLEN * (offset + 7)], v4i);
      }
    } //ic
  }   //site
#pragma omp barrier
}


#endif
//============================================================END=====
