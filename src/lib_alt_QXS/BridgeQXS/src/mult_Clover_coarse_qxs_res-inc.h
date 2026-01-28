/*!
      @file    mult_coarse_qxs-inc.h
      @brief
      @author  Issaku kanamori (kanamori)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2023-03-20 10:52:44 +0900 (Mon, 20 Mar 2023) $
      @version $LastChangedRevision: 2499 $
*/

#ifndef MULT_COARSE_QXS_INCLUDED
#define MULT_COARSE_QXS_INCLUDED

#include "mult_common_th-inc.h"


#define RUN_DIAG
#define RUN_HOP_XP
#define RUN_HOP_XM

#define RUN_HOP_YP
#define RUN_HOP_YM

#define RUN_HOP_ZP
#define RUN_HOP_ZM

#define RUN_HOP_TP
#define RUN_HOP_TM


namespace BridgeQXS {
  //====================================================================
  void mult_coarse_1(real_t *buf1_xp, real_t *buf1_xm,
                     real_t *buf1_yp, real_t *buf1_ym,
                     real_t *buf1_zp, real_t *buf1_zm,
                     real_t *buf1_tp, real_t *buf1_tm,
                     real_t *u0, real_t *v1, const int *Nsize,
                     int ncol, const int *do_comm)
  {
    int ith, nth, is, ns;
    int Nstv = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];
    int Nxv  = Nsize[0];
    int Nyv  = Nsize[1];
    int Nz   = Nsize[2];
    int Nt   = Nsize[3];
    int Nc   = ncol;
    int Nvc  = 2 * ncol; // 2 for complex
    int Nc2  = ncol * ncol;
    int Ndf  = 2 * Nc2;  // 2 for complex

#ifdef USE_QXS_ACLE
    svbool_t pg1_xp, pg2_xp, pg1_xm, pg2_xm;
    svbool_t pg1_yp, pg2_yp, pg1_ym, pg2_ym;
    set_predicate_xp(pg1_xp, pg2_xp);
    set_predicate_xm(pg1_xm, pg2_xm);
    set_predicate_yp(pg1_yp, pg2_yp);
    set_predicate_ym(pg1_ym, pg2_ym);
    svint_t svidx_xp, svidx_xm, svidx_yp, svidx_ym;
    set_index_xp(svidx_xp);
    set_index_xm(svidx_xm);
    set_index_yp(svidx_yp);
    set_index_ym(svidx_ym);
#endif

    int taskx      = (do_comm[0] > 0) ? (Nyv * Nz * Nt) : 0;
    int tasky      = (do_comm[1] > 0) ? (Nxv * Nz * Nt) : 0;
    int taskz      = (do_comm[2] > 0) ? (Nxv * Nyv * Nt) : 0;
    int taskt      = (do_comm[3] > 0) ? (Nxv * Nyv * Nz) : 0;
    int task_total = taskx + tasky + taskz + taskt;
    set_threadtask(ith, nth, is, ns, task_total);

    int isx = is;
    int nsx = (ns > taskx) ? taskx : ns;
    is -= taskx;
    ns -= taskx;
    int isy = (is < 0) ? 0 : is;
    int nsy = (ns > tasky) ? tasky : ns;
    is -= tasky;
    ns -= tasky;
    int isz = (is < 0) ? 0 : is;
    int nsz = (ns > taskz) ? taskz : ns;
    is -= taskz;
    ns -= taskz;
    int ist = (is < 0) ? 0 : is;
    int nst = (ns < 0) ? 0 : ns;

    for (int sitex = isx; sitex < nsx; ++sitex) {
      int    iyzt = sitex;
      int    ibf  = VLENY * Nvc * iyzt;
      int    idir = 0;
      real_t *u   = u0 + VLEN * Ndf * Nstv * idir;
      {
        int ix   = 0;
        int site = ix + Nxv * iyzt;
#ifdef USE_QXS_ACLE
        mult_coarse_xp1(pg2_xm, svidx_xm,
                        &buf1_xp[ibf], &v1[VLEN * Nvc * site], Nc);
#else
        mult_coarse_xp1(&buf1_xp[ibf], &v1[VLEN * Nvc * site], Nc);
#endif
      }
      {
        int ix   = Nxv - 1;
        int site = ix + Nxv * iyzt;
#ifdef USE_QXS_ACLE
        mult_coarse_xm1(pg2_xp, svidx_xp,
                        &buf1_xm[ibf], &u[VLEN * Ndf * site],
                        &v1[VLEN * Nvc * site], Nc);
#else
        mult_coarse_xm1(&buf1_xm[ibf], &u[VLEN * Ndf * site],
                        &v1[VLEN * Nvc * site], Nc);
#endif
      }
    }   // sitex

    for (int sitey = isy; sitey < nsy; sitey++) {
      int    ixzt = sitey;
      int    ix   = sitey % Nxv;
      int    izt  = sitey / Nxv;
      int    ibf  = VLENX * Nvc * ixzt;
      int    idir = 1;
      real_t *u   = u0 + VLEN * Ndf * Nstv * idir;
      {
        int iy   = 0;
        int site = ix + Nxv * iy + Nxv * Nyv * izt;
#ifdef USE_QXS_ACLE
        mult_coarse_yp1(pg2_ym, svidx_ym,
                        &buf1_yp[ibf], &v1[VLEN * Nvc * site], Nc);
#else
        mult_coarse_yp1(&buf1_yp[ibf], &v1[VLEN * Nvc * site], Nc);
#endif
      }
      {
        int iy   = Nyv - 1;
        int site = ix + Nxv * iy + Nxv * Nyv * izt;
#ifdef USE_QXS_ACLE
        mult_coarse_ym1(pg2_yp, svidx_yp,
                        &buf1_ym[ibf], &u[VLEN * Ndf * site],
                        &v1[VLEN * Nvc * site], Nc);
#else
        mult_coarse_ym1(&buf1_ym[ibf], &u[VLEN * Ndf * site],
                        &v1[VLEN * Nvc * site], Nc);
#endif
      }
    }   // sitey

    for (int sitez = isz; sitez < nsz; sitez++) {
      int    ixyt = sitez;
      int    ixy  = sitez % (Nxv * Nyv);
      int    it   = sitez / (Nxv * Nyv);
      int    idir = 2;
      real_t *u   = u0 + VLEN * Ndf * Nstv * idir;
      {
        int iz   = 0;
        int site = ixy + Nxv * Nyv * (iz + Nz * it);
        mult_coarse_zp1(&buf1_zp[VLEN * Nvc * ixyt], &v1[VLEN * Nvc * site], Nc);
      }
      {
        int iz   = Nz - 1;
        int site = ixy + Nxv * Nyv * (iz + Nz * it);
        mult_coarse_zm1(&buf1_zm[VLEN * Nvc * ixyt],
                        &u[VLEN * Ndf * site], &v1[VLEN * Nvc * site], Nc);
      }
    }   // sitez

    for (int sitet = ist; sitet < nst; sitet++) {
      int    ixyz = sitet;
      int    idir = 3;
      real_t *u   = u0 + VLEN * Ndf * Nstv * idir;
      {
        int it   = 0;
        int site = ixyz + Nxv * Nyv * Nz * it;
        mult_coarse_tp1(&buf1_tp[VLEN * Nvc * ixyz], &v1[VLEN * Nvc * site], Nc);
      }
      {
        int it   = Nt - 1;
        int site = ixyz + Nxv * Nyv * Nz * it;
        mult_coarse_tm1(&buf1_tm[VLEN * Nvc * ixyz],
                        &u[VLEN * Ndf * site], &v1[VLEN * Nvc * site], Nc);
      }
    }   // sitet
  }


//====================================================================
  void mult_coarse_b(real_t *v2,
                     real_t *u0, real_t *c0,
                     real_t *v1,
                     const int *Nsize, int ncol,
                     const int *do_comm, real_t *work)
  {
    int ith, nth, is, ns;
    int Nstv = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];
    int Nxv  = Nsize[0];
    int Nyv  = Nsize[1];
    int Nz   = Nsize[2];
    int Nt   = Nsize[3];
    int Nc   = ncol;
    int Nvc  = 2 * ncol; // 2 for complex
    int Nc2  = ncol * ncol;
    int Ndf  = 2 * Nc2;  // 2 for complex

#ifdef USE_QXS_ACLE
    svbool_t pg1_xp, pg2_xp, pg1_xm, pg2_xm;
    svbool_t pg1_yp, pg2_yp, pg1_ym, pg2_ym;
    set_predicate_xp(pg1_xp, pg2_xp);
    set_predicate_xm(pg1_xm, pg2_xm);
    set_predicate_yp(pg1_yp, pg2_yp);
    set_predicate_ym(pg1_ym, pg2_ym);
#endif


    int nv  = VLEN * Nvc;
    int nv2 = VLEN * Ndf;
    set_threadtask(ith, nth, is, ns, Nstv);

    for (int site = is; site < ns; ++site) {
      real_t *out = &v2[nv * site];

      // clover term
#ifdef RUN_DIAG
      set_mult_u(out, &v1[nv * site],
                 &c0[nv2 * site], Nc);
#else
      for (int i = 0; i < nv; i++) {
        out[i] = 0.0;
      }
#endif
      int ix   = site % Nxv;
      int iyzt = site / Nxv;
      {   // mult_xpb, mult_xmb
        int    idir = 0;
        real_t *u   = u0 + VLEN * Ndf * Nstv * idir;

#ifdef RUN_HOP_XP
        if ((ix < Nxv - 1) || (do_comm[0] == 0)) {
          int nei = (ix + 1) + Nxv * iyzt;
          if (ix == Nxv - 1) nei = 0 + Nxv * iyzt;
#ifdef USE_QXS_ACLE
          mult_coarse_xpb(pg1_xp, pg2_xp, out,
                          &u[nv2 * site],
                          &v1[nv * site], &v1[nv * nei], Nc, work);
#else
          mult_coarse_xpb(out,
                          &u[nv2 * site],
                          &v1[nv * site], &v1[nv * nei], Nc, work);
#endif
        }
#endif

#ifdef RUN_HOP_XM
        if ((ix > 0) || (do_comm[0] == 0)) {
          int ix2 = (ix - 1 + Nxv) % Nxv;
          int nei = ix2 + Nxv * iyzt;
#ifdef USE_QXS_ACLE
          mult_coarse_xmb(pg1_xm, pg2_xm, out,
                          &u[nv2 * site], &u[nv2 * nei],
                          &v1[nv * site], &v1[nv * nei],
                          Nc, work);
#else
          mult_coarse_xmb(out,
                          &u[nv2 * site], &u[nv2 * nei],
                          &v1[nv * site], &v1[nv * nei],
                          Nc, work);
#endif
        }
#endif
      }   // mult_xpb, mult_xmb, done

      int iy  = iyzt % Nyv;
      int izt = iyzt / Nyv;
      {   // mult_ypb, mult_ymb
        int    idir = 1;
        real_t *u   = u0 + VLEN * Ndf * Nstv * idir;
#ifdef RUN_HOP_YP
        if ((iy < Nyv - 1) || (do_comm[1] == 0)) {
          int iy2 = (iy + 1) % Nyv;
          int nei = ix + Nxv * (iy2 + Nyv * izt);
#ifdef USE_QXS_ACLE
          mult_coarse_ypb(pg1_yp, pg2_yp, out,
                          &u[nv2 * site],
                          &v1[nv * site], &v1[nv * nei],
                          Nc, work);
#else
          mult_coarse_ypb(out,
                          &u[nv2 * site],
                          &v1[nv * site], &v1[nv * nei],
                          Nc, work);
#endif
        }
#endif
#ifdef RUN_HOP_YM
        if ((iy != 0) || (do_comm[idir] == 0)) {
          int iy2 = (iy - 1 + Nyv) % Nyv;
          int nei = ix + Nxv * (iy2 + Nyv * izt);
#ifdef USE_QXS_ACLE
          mult_coarse_ymb(pg1_ym, pg2_ym, out,
                          &u[nv2 * site], &u[nv2 * nei],
                          &v1[nv * site], &v1[nv * nei],
                          Nc, work);
#else
          mult_coarse_ymb(out,
                          &u[nv2 * site], &u[nv2 * nei],
                          &v1[nv * site], &v1[nv * nei],
                          Nc, work);
#endif
        }
#endif
      }   // mult_ypb, mult_ymb, done

      int ixy  = ix + Nxv * iy;
      int iz   = izt % Nz;
      int it   = izt / Nz;
      int Nxyv = Nxv * Nyv;
      {   // mult_zpb, mult_zmb
        int    idir = 2;
        real_t *u   = u0 + VLEN * Ndf * Nstv * idir;

#ifdef RUN_HOP_ZP
        if ((iz != Nz - 1) || (do_comm[2] == 0)) {
          int iz2 = (iz + 1) % Nz;
          int nei = ixy + Nxyv * (iz2 + Nz * it);
          mult_coarse_zpb(out,
                          &u[nv2 * site], &v1[nv * nei], Nc);
        }
#endif
#ifdef RUN_HOP_ZM
        if ((iz > 0) || (do_comm[2] == 0)) {
          int iz2 = (iz - 1 + Nz) % Nz;
          int nei = ixy + Nxyv * (iz2 + Nz * it);
          mult_coarse_zmb(out,
                          &u[nv2 * nei], &v1[nv * nei], Nc);
        }
#endif
      }    // mult_zpb, mult_zmb, done

      int Nxyzv = Nxyv * Nz;
      int ixyz  = site - it * Nxyzv;
      {   // mult_tpb, mult_tmb
        int    idir = 3;
        real_t *u   = u0 + VLEN * Ndf * Nstv * idir;

#ifdef RUN_HOP_TP
        if ((it < Nt - 1) || (do_comm[3] == 0)) {
          int it2 = (it + 1) % Nt;
          int nei = ixyz + Nxyzv * it2;
          mult_coarse_tpb(out,
                          &u[nv2 * site], &v1[nv * nei], Nc);
        }
#endif
#ifdef RUN_HOP_TM
        if ((it > 0) || (do_comm[3] == 0)) {
          int it2 = (it - 1 + Nt) % Nt;
          int nei = ixyz + Nxyzv * it2;
          mult_coarse_tmb(out,
                          &u[nv2 * nei], &v1[nv * nei], Nc);
        }
#endif
      } // mult_tpb, mult_tmb, done
    }   // site
  }


//====================================================================
  void mult_coarse_2(real_t *v2, real_t *u0, real_t *v1,
                     real_t *buf2_xp, real_t *buf2_xm,
                     real_t *buf2_yp, real_t *buf2_ym,
                     real_t *buf2_zp, real_t *buf2_zm,
                     real_t *buf2_tp, real_t *buf2_tm,
                     const int *Nsize, int ncol, const int *do_comm,
                     real_t *work,
                     std::vector<int>& list)
  {
    int ith, nth, is, ns;
    int Nstv = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];
    int Nxv  = Nsize[0];
    int Nyv  = Nsize[1];
    int Nz   = Nsize[2];
    int Nt   = Nsize[3];
    int Nc   = ncol;
    int Nvc  = 2 * ncol; // 2 for complex
    int Nc2  = ncol * ncol;
    int Ndf  = 2 * Nc2;  // 2 for complex

#ifdef USE_QXS_ACLE
    svbool_t pg1_xp, pg2_xp, pg1_xm, pg2_xm;
    svbool_t pg1_yp, pg2_yp, pg1_ym, pg2_ym;
    set_predicate_xp(pg1_xp, pg2_xp);
    set_predicate_xm(pg1_xm, pg2_xm);
    set_predicate_yp(pg1_yp, pg2_yp);
    set_predicate_ym(pg1_ym, pg2_ym);
    svint_t svidx_xp, svidx_xm, svidx_yp, svidx_ym;
    set_index_xp(svidx_xp);
    set_index_xm(svidx_xm);
    set_index_yp(svidx_yp);
    set_index_ym(svidx_ym);
#endif

    int nv  = VLEN * Nvc;
    int nv2 = VLEN * Ndf;

    for (int i = 0; i < list.size(); i++) {
      int    site = list[i];
      real_t *out = &v2[nv * site];

      const int ix   = site % Nxv;
      const int iyzt = site / Nxv;

      if (do_comm[0] == 1) {
        int    idir = 0;
        int    ibf  = VLENY * Nvc * iyzt;
        real_t *u   = u0 + nv2 * Nstv * idir;
#ifdef RUN_HOP_XP
        if (ix == Nxv - 1) {
#ifdef USE_QXS_ACLE
          mult_coarse_xp2(pg1_xp, pg2_xp, svidx_xp,
                          out, &u[nv2 * site],
                          &v1[nv * site], &buf2_xp[ibf], Nc, work);
#else
          mult_coarse_xp2(out, &u[nv2 * site],
                          &v1[nv * site], &buf2_xp[ibf], Nc, work);
#endif
        }
#endif
#ifdef RUN_HOP_XM
        if (ix == 0) {
#ifdef USE_QXS_ACLE
          mult_coarse_xm2(pg1_xm, pg2_xm, svidx_xm,
                          out, &u[nv2 * site],
                          &v1[nv * site], &buf2_xm[ibf], Nc);
#else
          mult_coarse_xm2(out, &u[nv2 * site],
                          &v1[nv * site], &buf2_xm[ibf], Nc);
#endif
        }
#endif
      }   // do_comm[0] == 1


      const int iy  = iyzt % Nyv;
      const int izt = iyzt / Nyv;

      if (do_comm[1] == 1) {
        int    idir = 1;
        int    ixzt = ix + Nxv * izt;
        int    ibf  = VLENX * Nvc * ixzt;
        real_t *u   = u0 + nv2 * Nstv * idir;
#ifdef RUN_HOP_YP
        if (iy == Nyv - 1) {
#ifdef USE_QXS_ACLE
          mult_coarse_yp2(pg1_yp, pg2_yp, svidx_yp,
                          out,
                          &u[nv2 * site],
                          &v1[nv * site], &buf2_yp[ibf], Nc, work);
#else
          mult_coarse_yp2(out,
                          &u[nv2 * site],
                          &v1[nv * site], &buf2_yp[ibf], Nc, work);
#endif
        }
#endif
#ifdef RUN_HOP_YM
        if (iy == 0) {
#ifdef USE_QXS_ACLE
          mult_coarse_ym2(pg1_ym, pg2_ym, svidx_ym,
                          out,
                          &u[nv2 * site],
                          &v1[nv * site], &buf2_ym[ibf], Nc);
#else
          mult_coarse_ym2(out,
                          &u[nv2 * site],
                          &v1[nv * site], &buf2_ym[ibf], Nc);
#endif
        }
#endif
      }   // do_comm[1] == 1


      const int ixy  = ix + Nxv * iy;
      const int iz   = izt % Nz;
      const int it   = izt / Nz;
      const int Nxyv = Nxv * Nyv;

      if (do_comm[2] == 1) {
        int    idir = 2;
        int    ixyt = ixy + Nxyv * it;
        real_t *u   = u0 + nv2 * Nstv * idir;
#ifdef RUN_HOP_ZP
        if (iz == Nz - 1) {
          mult_coarse_zp2(out,
                          &u[nv2 * site], &buf2_zp[nv * ixyt], Nc);
        }
#endif
#ifdef RUN_HOP_ZM
        if (iz == 0) {
          mult_coarse_zm2(out,
                          &buf2_zm[nv * ixyt], Nc);
        }
#endif
      }   // do_comm[2] == 1

      if (do_comm[3] == 1) {
        int    idir = 3;
        int    ixyz = ixy + Nxyv * iz;
        real_t *u   = u0 + nv2 * Nstv * idir;
#ifdef RUN_HOP_TP
        if (it == Nt - 1) {
          mult_coarse_tp2(out,
                          &u[nv2 * site], &buf2_tp[nv * ixyz], Nc);
        }
#endif
#ifdef RUN_HOP_TM
        if (it == 0) {
          mult_coarse_tm2(out,
                          &buf2_tm[nv * ixyz], Nc);
        }
#endif
      } // do_comm[3] == 1
    }   // site
  }
}

#endif
//============================================================END=====
