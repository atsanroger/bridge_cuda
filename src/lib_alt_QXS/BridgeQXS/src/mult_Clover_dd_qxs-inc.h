/*!
      @file    mult_Clover_dd_qxs-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#ifndef MULT_CLOVER_DD_QXS_INCLUDED
#define MULT_CLOVER_DD_QXS_INCLUDED

#include "mult_common_th-inc.h"

//====================================================================
void BridgeQXS::mult_clover_dd_dirac(real_t *v2, real_t *up,
                                     real_t *ct, real_t *v1,
                                     real_t kappa, int *bc,
                                     int *Nsize, int *block_size,
                                     int ieo)
{
  // lattice size in units of SIMD vector
  int Nxv  = Nsize[0];
  int Nyv  = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Nyv * Nz * Nt;
  int Nst  = Nstv * VLEN;

  // size of block in units of SIMD vector
  int Bxv   = block_size[0];
  int Byv   = block_size[1];
  int Bz    = block_size[2];
  int Bt    = block_size[3];
  int Bsize = Bxv * Byv * Bz * Bt;

  // numbers of blocks
  int NBx    = Nsize[0] / block_size[0];
  int NBy    = Nsize[1] / block_size[1];
  int NBz    = Nsize[2] / block_size[2];
  int NBt    = Nsize[3] / block_size[3];
  int Nblock = NBx * NBy * NBz * NBt;

  svbool_t pg1_xp, pg2_xp, pg1_xm, pg2_xm;
  svbool_t pg1_yp, pg2_yp, pg1_ym, pg2_ym;
  set_predicate_xp(pg1_xp, pg2_xp);
  set_predicate_xm(pg1_xm, pg2_xm);
  set_predicate_yp(pg1_yp, pg2_yp);
  set_predicate_ym(pg1_ym, pg2_ym);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Bsize);

  for (int block = 0; block < Nblock; ++block) {
    int ibx = block % NBx;
    int iby = (block / NBx) % NBy;
    int ibz = (block / (NBx * NBy)) % NBz;
    int ibt = block / (NBx * NBy * NBz);
    int jeo = (ieo + ibx + iby + ibz + ibt) % 2;
    //    if(jeo == 1) continue;
    if ((ieo > -1) && (jeo == 1)) continue;

    for (int bsite = is; bsite < ns; ++bsite) {
      int kx   = bsite % Bxv;
      int ix   = kx + Bxv * ibx;
      int kyzt = bsite / Bxv;
      int ky   = kyzt % Byv;
      int iy   = ky + Byv * iby;
      int kzt  = kyzt / Byv;
      int kz   = kzt % Bz;
      int iz   = kz + Bz * ibz;
      int kt   = kzt / Bz;
      int it   = kt + Bt * ibt;
      int site = ix + Nxv * (iy + Nyv * (iz + Nz * it));

      Vsimd_t v2v[NVCD];
      clear_vec(v2v, NVCD);

      real_t zL[VLEN * NVCD];
      real_t uL[VLEN * NDF];

      if (kx < Bxv - 1) {
        real_t *u  = &up[NDF * Nst * 0];
        int    nei = site + 1;
        mult_wilson_xpb(pg1_xp, pg2_xp, v2v, &u[VLEN * NDF * site],
                        &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
      } else {  // kx = Bxv-1
        real_t *u  = &up[NDF * Nst * 0];
        int    nei = site;
        mult_wilson_xpb(pg1_xp, pg2_xp, v2v, &u[VLEN * NDF * site],
                        &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
      }

      if (kx > 0) {
        real_t *u  = &up[NDF * Nst * 0];
        int    nei = site - 1;
        mult_wilson_xmb(pg1_xm, pg2_xm, v2v,
                        &u[VLEN * NDF * site], &u[VLEN * NDF * nei],
                        &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
      } else {   // kx = 0
        real_t *u  = &up[NDF * Nst * 0];
        int    nei = site + Bxv - 1;
        mult_wilson_xmb(pg1_xm, pg2_xm, v2v,
                        &u[VLEN * NDF * site], &u[VLEN * NDF * nei],
                        &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
      }

      if (ky < Byv - 1) {
        int    nei = site + Nxv;
        real_t *u  = &up[NDF * Nst * 1];
        mult_wilson_ypb(pg1_yp, pg2_yp, v2v,
                        &u[VLEN * NDF * site],
                        &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
      } else {  // ky = Byv-1
        int    nei = site;
        real_t *u  = &up[NDF * Nst * 1];
        mult_wilson_ypb(pg1_yp, pg2_yp, v2v,
                        &u[VLEN * NDF * site],
                        &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
      }

      if (ky > 0) {
        int    nei = site - Nxv;
        real_t *u  = &up[NDF * Nst * 1];
        mult_wilson_ymb(pg1_ym, pg2_ym, v2v,
                        &u[VLEN * NDF * site], &u[VLEN * NDF * nei],
                        &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
      } else {  // ky = 0
        int    nei = site + Nxv * (Byv - 1);
        real_t *u  = &up[NDF * Nst * 1];
        mult_wilson_ymb(pg1_ym, pg2_ym, v2v,
                        &u[VLEN * NDF * site], &u[VLEN * NDF * nei],
                        &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
      }

      if (kz < Bz - 1) {
        int    nei = site + Nxv * Nyv;
        real_t *u  = &up[NDF * Nst * 2];
        mult_wilson_zpb(v2v, &u[VLEN * NDF * site], &v1[VLEN * NVCD * nei]);
      }

      if (kz > 0) {
        int    nei = site - Nxv * Nyv;
        real_t *u  = &up[NDF * Nst * 2];
        mult_wilson_zmb(v2v, &u[VLEN * NDF * nei], &v1[VLEN * NVCD * nei]);
      }

      if (kt < Bt - 1) {
        int    nei = site + Nxv * Nyv * Nz;
        real_t *u  = &up[NDF * Nst * 3];
        mult_wilson_tpb_dirac(v2v, &u[VLEN * NDF * site],
                              &v1[VLEN * NVCD * nei]);
      }

      if (kt > 0) {
        int    nei = site - Nxv * Nyv * Nz;
        real_t *u  = &up[NDF * Nst * 3];
        mult_wilson_tmb_dirac(v2v, &u[VLEN * NDF * nei],
                              &v1[VLEN * NVCD * nei]);
      }

      mult_clover_csw_aypx(&v2[VLEN * NVCD * site], -kappa, v2v,
                           &ct[VLEN * NDF * 2 * ND * site], &v1[VLEN * NVCD * site]);
    }
  }
}


//====================================================================
void BridgeQXS::mult_clover_dd_dirac_chrot(
  real_t *v2, real_t *up,
  real_t *ct, real_t *v1,
  real_t kappa, int *bc,
  int *Nsize, int *block_size,
  int ieo)
{
  // lattice size in units of SIMD vector
  int Nxv  = Nsize[0];
  int Nyv  = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Nyv * Nz * Nt;
  int Nst  = Nstv * VLEN;

  // size of block in units of SIMD vector
  int Bxv   = block_size[0];
  int Byv   = block_size[1];
  int Bz    = block_size[2];
  int Bt    = block_size[3];
  int Bsize = Bxv * Byv * Bz * Bt;

  // numbers of blocks
  int NBx    = Nsize[0] / block_size[0];
  int NBy    = Nsize[1] / block_size[1];
  int NBz    = Nsize[2] / block_size[2];
  int NBt    = Nsize[3] / block_size[3];
  int Nblock = NBx * NBy * NBz * NBt;

  svbool_t pg1_xp, pg2_xp, pg1_xm, pg2_xm;
  svbool_t pg1_yp, pg2_yp, pg1_ym, pg2_ym;
  set_predicate_xp(pg1_xp, pg2_xp);
  set_predicate_xm(pg1_xm, pg2_xm);
  set_predicate_yp(pg1_yp, pg2_yp);
  set_predicate_ym(pg1_ym, pg2_ym);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Bsize);

  for (int block = 0; block < Nblock; ++block) {
    int ibx = block % NBx;
    int iby = (block / NBx) % NBy;
    int ibz = (block / (NBx * NBy)) % NBz;
    int ibt = block / (NBx * NBy * NBz);
    int jeo = (ieo + ibx + iby + ibz + ibt) % 2;
    //    if(jeo == 1) continue;
    if ((ieo > -1) && (jeo == 1)) continue;

    for (int bsite = is; bsite < ns; ++bsite) {
      int kx   = bsite % Bxv;
      int ix   = kx + Bxv * ibx;
      int kyzt = bsite / Bxv;
      int ky   = kyzt % Byv;
      int iy   = ky + Byv * iby;
      int kzt  = kyzt / Byv;
      int kz   = kzt % Bz;
      int iz   = kz + Bz * ibz;
      int kt   = kzt / Bz;
      int it   = kt + Bt * ibt;
      int site = ix + Nxv * (iy + Nyv * (iz + Nz * it));

      Vsimd_t v2v[NVCD];
      clear_vec(v2v, NVCD);

      real_t zL[VLEN * NVCD];
      real_t uL[VLEN * NDF];

      if (kx < Bxv - 1) {
        real_t *u  = &up[NDF * Nst * 0];
        int    nei = site + 1;
        mult_wilson_xpb(pg1_xp, pg2_xp, v2v, &u[VLEN * NDF * site],
                        &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
      } else {  // kx = Bxv-1
        real_t *u  = &up[NDF * Nst * 0];
        int    nei = site;
        mult_wilson_xpb(pg1_xp, pg2_xp, v2v, &u[VLEN * NDF * site],
                        &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
      }

      if (kx > 0) {
        real_t *u  = &up[NDF * Nst * 0];
        int    nei = site - 1;
        mult_wilson_xmb(pg1_xm, pg2_xm, v2v,
                        &u[VLEN * NDF * site], &u[VLEN * NDF * nei],
                        &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
      } else {   // kx = 0
        real_t *u  = &up[NDF * Nst * 0];
        int    nei = site + Bxv - 1;
        mult_wilson_xmb(pg1_xm, pg2_xm, v2v,
                        &u[VLEN * NDF * site], &u[VLEN * NDF * nei],
                        &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
      }

      if (ky < Byv - 1) {
        int    nei = site + Nxv;
        real_t *u  = &up[NDF * Nst * 1];
        mult_wilson_ypb(pg1_yp, pg2_yp, v2v,
                        &u[VLEN * NDF * site],
                        &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
      } else {  // ky = Byv-1
        int    nei = site;
        real_t *u  = &up[NDF * Nst * 1];
        mult_wilson_ypb(pg1_yp, pg2_yp, v2v,
                        &u[VLEN * NDF * site],
                        &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
      }

      if (ky > 0) {
        int    nei = site - Nxv;
        real_t *u  = &up[NDF * Nst * 1];
        mult_wilson_ymb(pg1_ym, pg2_ym, v2v,
                        &u[VLEN * NDF * site], &u[VLEN * NDF * nei],
                        &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
      } else {  // ky = 0
        int    nei = site + Nxv * (Byv - 1);
        real_t *u  = &up[NDF * Nst * 1];
        mult_wilson_ymb(pg1_ym, pg2_ym, v2v,
                        &u[VLEN * NDF * site], &u[VLEN * NDF * nei],
                        &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
      }

      if (kz < Bz - 1) {
        int    nei = site + Nxv * Nyv;
        real_t *u  = &up[NDF * Nst * 2];
        mult_wilson_zpb(v2v, &u[VLEN * NDF * site], &v1[VLEN * NVCD * nei]);
      }

      if (kz > 0) {
        int    nei = site - Nxv * Nyv;
        real_t *u  = &up[NDF * Nst * 2];
        mult_wilson_zmb(v2v, &u[VLEN * NDF * nei], &v1[VLEN * NVCD * nei]);
      }

      if (kt < Bt - 1) {
        int    nei = site + Nxv * Nyv * Nz;
        real_t *u  = &up[NDF * Nst * 3];
        mult_wilson_tpb_dirac(v2v, &u[VLEN * NDF * site],
                              &v1[VLEN * NVCD * nei]);
      }

      if (kt > 0) {
        int    nei = site - Nxv * Nyv * Nz;
        real_t *u  = &up[NDF * Nst * 3];
        mult_wilson_tmb_dirac(v2v, &u[VLEN * NDF * nei],
                              &v1[VLEN * NVCD * nei]);
      }

      mult_clover_csw_aypx_chrot(
        &v2[VLEN * NVCD * site], -kappa, &v2v[0].v[0],
        &ct[VLEN * NDF * 2 * ND * site], &v1[VLEN * NVCD * site]);
    }
  }
}


//============================================================END=====
#endif
