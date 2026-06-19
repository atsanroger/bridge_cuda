/*!
      @file    mult_Doainwall_5din_qxs-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2581 $
*/

#ifndef MULT_DOMAINWALL_5DIN_DD_QXS_INCLUDED
#define MULT_DOMAINWALL_5DIN_DD_QXS_INCLUDED

#include "mult_common_th-inc.h"


//====================================================================
//void BridgeQXS::mult_domainwall_5din_dd_5dir_dirac(
void mult_domainwall_5din_dd_5dir_dirac_try(
  real_t *vp, real_t *yp, real_t *wp,
  real_t mq, real_t M0, int Ns, int *bc,
  real_t *b, real_t *c,
  int *Nsize, int *block_size,
  int ieo)
{
  int Nxv  = Nsize[0];
  int Nyv  = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Nyv * Nz * Nt;
  int Nst  = Nstv * VLEN;

  int Nxy  = Nxv * Nyv;
  int Nxyz = Nxv * Nyv * Nz;

  int Nin4 = VLEN * NVCD;
  int Nin5 = Nin4 * Ns;

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
  int Nblock_total = NBx * NBy * NBz * NBt;

  int Nblock = Nblock_total;
  int Nblock_eo = Nblock_total % 2;
  if( ieo > -1 ) {
    Nblock = (Nblock_total + 1-ieo)/2;
    // total=9, ieo=0: Nblock = 5
    // total=9, ieo=1: Nblock = 4
    // total=8, ieo=0: Nblock = 4
    // total=8, ieo=1: Nblock = 4
  }
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Bsize*Nblock);
  for (int idx = 0; idx < ns; ++idx) {
    int bsite = idx % Bsize;
    int block = idx / Bsize;
    if( ieo > -1 ){
      block = 2*block + ieo;
    }
    int ibx = block % NBx;
    int iby = (block / NBx) % NBy;
    int ibz = (block / (NBx * NBy)) % NBz;
    int ibt = block / (NBx * NBy * NBz);
    if( ieo > -1 ) {
      if( ieo !=  (ibx + iby + ibz + ibt) % 2 ) {
        printf("hoge! ieo=%d, (ibx,iby,ibz,ibt)=(%d,%d,%d,%d)\n",
               ieo, ibx, iby, ibz, ibt);
      }
    }
    // Nblock = 4, 2dim
    // 0 1 0 1
    // 1 0 1 0
    //   ieo=0:
    //     block = 0,2,4,6
    //     ibx = 0,2,0,2
    //     iby = 0,0,1,1
    //     yeo = 0 0 0 0
    //   ieo=1:
    //     block = 1,3,5,7
    //     ibx = 1,3,1,3
    //     iby = 0,0,1,1
    //     yeo = 0 0 0 0


    // 0 1 0
    // 1 0 1
    // 0 1 0
    //   ieo=0: Nblock = 5
    //     block = 0,2,4,6,8
    //     ibx = 0,2,1,0,2
    //     iby = 0,0,1,2,2
    //     yeo = 0 0 1 0 0
    //   ieo=1:
    //     block = 1,3,5,7
    //     ibx = 1,0,2,1
    //     iby = 0,1,1,2
    //     yeo = 0 1 1 0

    // 0 1 0
    // 1 0 1
    //   ieo=0: Nblock = 3
    //     block = 0,2,4
    //     ibx = 0,2,1
    //     iby = 0,0,1
    //     yeo = 0 0 1
    //   ieo=1:
    //     block = 1,3,5
    //     ibx = 1,0,2
    //     iby = 0,1,1
    //     yeo = 0 1 1

      int kx   = bsite % Bxv;
      int kyzt = bsite / Bxv;
      int ky   = kyzt % Byv;
      int kzt  = kyzt / Byv;
      int kz   = kzt % Bz;
      int kt   = kzt / Bz;

      int ix   = kx + Bxv * ibx;
      int iy   = ky + Byv * iby;
      int iz   = kz + Bz * ibz;
      int it   = kt + Bt * ibt;

      int site = ix + Nxv * (iy + Nyv * (iz + Nz * it));

      real_t *vp2 = &vp[Nin5 * site];
      real_t *wp2 = &wp[Nin5 * site];
      real_t *yp2 = &yp[Nin5 * site];

      for (int is = 0; is < Ns; ++is) {
        svbool_t pg = set_predicate();

        real_t Fup   = 0.5;
        if (is == Ns - 1) Fup = -0.5 * mq;
        real_t Fdn   = 0.5;
        if (is == 0) Fdn = -0.5 * mq;

        int    is_up = (is + 1) % Ns;
        int    is_dn = (is - 1 + Ns) % Ns;
        real_t   FF1    = b[is] * (4.0 - M0) + 1.0;
        real_t   FF2    = c[is] * (4.0 - M0) - 1.0;

        for (int ic = 0; ic < NC; ++ic) {
          svreal_t vt1r, vt1i, vt2r, vt2i, vt3r, vt3i, vt4r, vt4i;

          set_aPm5_dirac_vec(vt1r, vt1i, vt2r, vt2i,
                             vt3r, vt3i, vt4r, vt4i,
                             Fup, wp2, is_up, ic);

          add_aPp5_dirac_vec(vt1r, vt1i, vt2r, vt2i,
                             vt3r, vt3i, vt4r, vt4i,
                             Fdn, wp2, is_dn, ic);

          int      offset = 2 * ND * ic + NVCD * is;
          dw_5dir_axpy(pg, vp2, yp2, wp2, FF1, FF2, b[is], c[is], vt1r, 0 + offset);
          dw_5dir_axpy(pg, vp2, yp2, wp2, FF1, FF2, b[is], c[is], vt1i, 1 + offset);
          dw_5dir_axpy(pg, vp2, yp2, wp2, FF1, FF2, b[is], c[is], vt2r, 2 + offset);
          dw_5dir_axpy(pg, vp2, yp2, wp2, FF1, FF2, b[is], c[is], vt2i, 3 + offset);
          dw_5dir_axpy(pg, vp2, yp2, wp2, FF1, FF2, b[is], c[is], vt3r, 4 + offset);
          dw_5dir_axpy(pg, vp2, yp2, wp2, FF1, FF2, b[is], c[is], vt3i, 5 + offset);
          dw_5dir_axpy(pg, vp2, yp2, wp2, FF1, FF2, b[is], c[is], vt4r, 6 + offset);
          dw_5dir_axpy(pg, vp2, yp2, wp2, FF1, FF2, b[is], c[is], vt4i, 7 + offset);

        } // ic
      } //is
  } // idx
}

//====================================================================
void BridgeQXS::mult_domainwall_5din_dd_5dir_dirac(
  real_t *vp, real_t *yp, real_t *wp,
  real_t mq, real_t M0, int Ns, int *bc,
  real_t *b, real_t *c,
  int *Nsize, int *block_size,
  int ieo)
{
  int Nxv  = Nsize[0];
  int Nyv  = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Nyv * Nz * Nt;
  int Nst  = Nstv * VLEN;

  int Nxy  = Nxv * Nyv;
  int Nxyz = Nxv * Nyv * Nz;

  int Nin4 = VLEN * NVCD;
  int Nin5 = Nin4 * Ns;

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
      int kyzt = bsite / Bxv;
      int ky   = kyzt % Byv;
      int kzt  = kyzt / Byv;
      int kz   = kzt % Bz;
      int kt   = kzt / Bz;

      int ix   = kx + Bxv * ibx;
      int iy   = ky + Byv * iby;
      int iz   = kz + Bz * ibz;
      int it   = kt + Bt * ibt;

      int site = ix + Nxv * (iy + Nyv * (iz + Nz * it));

      real_t *vp2 = &vp[Nin5 * site];
      real_t *wp2 = &wp[Nin5 * site];
      real_t *yp2 = &yp[Nin5 * site];

      for (int is = 0; is < Ns; ++is) {
        svbool_t pg = set_predicate();

        real_t Fup   = 0.5;
        if (is == Ns - 1) Fup = -0.5 * mq;
        real_t Fdn   = 0.5;
        if (is == 0) Fdn = -0.5 * mq;

        int    is_up = (is + 1) % Ns;
        int    is_dn = (is - 1 + Ns) % Ns;
        real_t   FF1    = b[is] * (4.0 - M0) + 1.0;
        real_t   FF2    = c[is] * (4.0 - M0) - 1.0;

        for (int ic = 0; ic < NC; ++ic) {
          svreal_t vt1r, vt1i, vt2r, vt2i, vt3r, vt3i, vt4r, vt4i;

          set_aPm5_dirac_vec(vt1r, vt1i, vt2r, vt2i,
                             vt3r, vt3i, vt4r, vt4i,
                             Fup, wp2, is_up, ic);

          add_aPp5_dirac_vec(vt1r, vt1i, vt2r, vt2i,
                             vt3r, vt3i, vt4r, vt4i,
                             Fdn, wp2, is_dn, ic);

          int      offset = 2 * ND * ic + NVCD * is;
          dw_5dir_axpy(pg, vp2, yp2, wp2, FF1, FF2, b[is], c[is], vt1r, 0 + offset);
          dw_5dir_axpy(pg, vp2, yp2, wp2, FF1, FF2, b[is], c[is], vt1i, 1 + offset);
          dw_5dir_axpy(pg, vp2, yp2, wp2, FF1, FF2, b[is], c[is], vt2r, 2 + offset);
          dw_5dir_axpy(pg, vp2, yp2, wp2, FF1, FF2, b[is], c[is], vt2i, 3 + offset);
          dw_5dir_axpy(pg, vp2, yp2, wp2, FF1, FF2, b[is], c[is], vt3r, 4 + offset);
          dw_5dir_axpy(pg, vp2, yp2, wp2, FF1, FF2, b[is], c[is], vt3i, 5 + offset);
          dw_5dir_axpy(pg, vp2, yp2, wp2, FF1, FF2, b[is], c[is], vt4r, 6 + offset);
          dw_5dir_axpy(pg, vp2, yp2, wp2, FF1, FF2, b[is], c[is], vt4i, 7 + offset);

        } // ic
      } //is
    } // bsite
  } // block
}


//====================================================================
void BridgeQXS::mult_domainwall_5din_dd_5dirdag_dirac(
  real_t *vp, real_t *yp, real_t *wp,
  real_t mq, real_t M0, int Ns, int *bc,
  real_t *b, real_t *c,
  int *Nsize, int *block_size,
  int ieo)
{
  int Nxv  = Nsize[0];
  int Nyv  = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Nyv * Nz * Nt;
  int Nst  = Nstv * VLEN;

  int Nxy  = Nxv * Nyv;
  int Nxyz = Nxv * Nyv * Nz;

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

  int Nin4 = VLEN * NVCD;
  int Nin5 = Nin4 * Ns;

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
      int kyzt = bsite / Bxv;
      int ky   = kyzt % Byv;
      int kzt  = kyzt / Byv;
      int kz   = kzt % Bz;
      int kt   = kzt / Bz;

      int ix   = kx + Bxv * ibx;
      int iy   = ky + Byv * iby;
      int iz   = kz + Bz * ibz;
      int it   = kt + Bt * ibt;

      int site = ix + Nxv * (iy + Nyv * (iz + Nz * it));

      real_t *vp2 = &vp[Nin5 * site];
      real_t *wp2 = &wp[Nin5 * site];
      real_t *yp2 = &yp[Nin5 * site];


      for (int is = 0; is < Ns; ++is) {
        svbool_t pg = set_predicate();

        real_t Fup = 0.5;
        if (is == Ns - 1) Fup = -0.5 * mq;
        real_t Fdn = 0.5;
        if (is == 0) Fdn = -0.5 * mq;

        int is_up = (is + 1) % Ns;
        int is_dn = (is - 1 + Ns) % Ns;

        real_t FF1   = b[is] * (4.0 - M0) + 1.0;
        real_t a1    = -0.5 * b[is];
        real_t FF2up = c[is_up] * (4.0 - M0) - 1.0;
        real_t aup   = -0.5 * c[is_up];
        real_t FF2dn = c[is_dn] * (4.0 - M0) - 1.0;
        real_t adn   = -0.5 * c[is_dn];

        for (int ic = 0; ic < NC; ++ic) {
          svreal_t vt1r, vt1i, vt2r, vt2i, vt3r, vt3i, vt4r, vt4i;

          int    index = 2 * ND * ic + NVCD * is;
          dw_5dir_dag(pg, vt3r, vt3i, vt4r, vt4i, vt1r, vt1i, vt2r, vt2i,
                      wp2, yp2, -FF1, -a1, index);

          svreal_t xt1r, xt1i, xt2r, xt2i, xt3r, xt3i, xt4r, xt4i;
          int index_up = 2 * ND * ic + NVCD * is_up;
          dw_5dir_dag(pg, xt1r, xt1i, xt2r, xt2i, xt3r, xt3i, xt4r, xt4i,
                      wp2, yp2, FF2up, aup, index_up);

          add_aPp5_dirac_vec(pg, vt1r, vt1i, vt2r, vt2i, vt3r, vt3i, vt4r, vt4i,
                             Fup, xt1r, xt1i, xt2r, xt2i, xt3r, xt3i, xt4r, xt4i);

          int index_dn = 2 * ND * ic + NVCD * is_dn;
          dw_5dir_dag(pg, xt1r, xt1i, xt2r, xt2i, xt3r, xt3i, xt4r, xt4i,
                      wp2, yp2, FF2dn, adn, index_dn);

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
        } // ic
      } // is
    } // bsite
  } // block
}


//====================================================================
void BridgeQXS::mult_domainwall_5din_dd_mult_gm5_dirac(real_t *vp, real_t *wp,
                                                    int Ns, int *Nsize,
                                                    int *block_size, int ieo)
{
  int Nxv  = Nsize[0];
  int Nyv  = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Nyv * Nz * Nt;
  int Nst  = Nstv * VLEN;

  int Nxy  = Nxv * Nyv;
  int Nxyz = Nxv * Nyv * Nz;

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

  int Nin4 = VLEN * NVCD;
  int Nin5 = Nin4 * Ns;

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
      int kyzt = bsite / Bxv;
      int ky   = kyzt % Byv;
      int kzt  = kyzt / Byv;
      int kz   = kzt % Bz;
      int kt   = kzt / Bz;

      int ix   = kx + Bxv * ibx;
      int iy   = ky + Byv * iby;
      int iz   = kz + Bz * ibz;
      int it   = kt + Bt * ibt;

      int site = ix + Nxv * (iy + Nyv * (iz + Nz * it));

      real_t *vp2 = &vp[Nin5 * site];
      real_t *wp2 = &wp[Nin5 * site];

      svbool_t pg = set_predicate();
      for (int is = 0; is < Ns; ++is) {
        for (int ic = 0; ic < NC; ++ic) {
          svreal_t vt1r, vt1i, vt2r, vt2i, vt3r, vt3i, vt4r, vt4i;
          int      index = 2 * ND * ic + NVCD * is;
          load_vec(pg, vt3r, &wp2[VLEN * (0 + index)]);
          load_vec(pg, vt3i, &wp2[VLEN * (1 + index)]);
          load_vec(pg, vt4r, &wp2[VLEN * (2 + index)]);
          load_vec(pg, vt4i, &wp2[VLEN * (3 + index)]);
          load_vec(pg, vt1r, &wp2[VLEN * (4 + index)]);
          load_vec(pg, vt1i, &wp2[VLEN * (5 + index)]);
          load_vec(pg, vt2r, &wp2[VLEN * (6 + index)]);
          load_vec(pg, vt2i, &wp2[VLEN * (7 + index)]);

          flip_sign(pg, vt3r);
          flip_sign(pg, vt3i);
          flip_sign(pg, vt4r);
          flip_sign(pg, vt4i);
          flip_sign(pg, vt1r);  // vt1r := -vt1r
          flip_sign(pg, vt1i);
          flip_sign(pg, vt2r);
          flip_sign(pg, vt2i);
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
}


//====================================================================
void BridgeQXS::mult_domainwall_5din_dd_hopb_dirac(
  real_t *vp, real_t *up, real_t *wp,
  real_t mq, real_t M0, int Ns, int *bc,
  real_t *b, real_t *c,
  int *Nsize, int *block_size,
  int ieo)
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
      int kyzt = bsite / Bxv;
      int ky   = kyzt % Byv;
      int kzt  = kyzt / Byv;
      int kz   = kzt % Bz;
      int kt   = kzt / Bz;

      int ix   = kx + Bxv * ibx;
      int iy   = ky + Byv * iby;
      int iz   = kz + Bz * ibz;
      int it   = kt + Bt * ibt;

      int site = ix + Nxv * (iy + Nyv * (iz + Nz * it));

      real_t *wp2 = &wp[Nin5 * site];
      real_t *vp2 = &vp[Nin5 * site];

      Vsimd_t vL[NVCD * Ns];

      load_vec(vL, vp2, NVCD * Ns);

      {
        int idir = 0;
        int nei = site + 1;
        if ( kx == Bxv - 1 ) { nei = site; } // nei can be any, as the link is 0
        real_t *u   = &up[VLEN * NDF * site + NvU * idir];
        real_t *wpn = &wp[Nin5*nei];
        for (int is = 0; is < Ns; ++is) {
          mult_wilson_xpb(pg1_xp, pg2_xp, &vL[NVCD * is].v[0], u,
                          &wp2[Nin4 * is], &wpn[Nin4 * is]);
        }
      }

      {
        int idir = 0;
        int nei = site - 1;
        if ( kx == 0 ) { nei = site + Bxv - 1; } // nei can be any, as the link is 0
        // nei = site - 1 + (Bvx-kx)/Bvx
        real_t *u   = &up[VLEN * NDF * site + NvU * idir];
        real_t *un  = &up[VLEN * NDF * nei  + NvU * idir];
        real_t *wpn = &wp[Nin5*nei];
        if (kx < Bxv - 1) { wpn += Nin5; }
        for (int is = 0; is < Ns; ++is) {
          mult_wilson_xmb(pg1_xm, pg2_xm, &vL[NVCD * is].v[0], u, un,
                          &wp2[Nin4 * is], &wpn[Nin4 * is]);
        }
      }

      {
        int idir = 1;
        int nei = site + Nxv;
        if ( ky == Byv - 1 ) { nei = site; } // nei can be any, as the link is 0
        // nei = site + (1- (ky+1)/Bvy )*Nvx;
        real_t *u   = &up[VLEN * NDF * site + NvU * idir];
        real_t *wpn = &wp[Nin5*nei];
        for (int is = 0; is < Ns; ++is) {
          mult_wilson_ypb(pg1_yp, pg2_yp, &vL[NVCD * is].v[0], u,
                          &wp2[Nin4 * is], &wpn[Nin4 * is]);
        }
      }

      {
        int idir = 1;
        int nei = site - Nxv;
        if ( ky == 0 ) { nei = site + Nxv * (Byv - 1); } // nei can be any, as the link is 0
        // nei = site - (1- (Bvy-ky)/Bvy )*Nvx;
        real_t *u   = &up[VLEN * NDF * site + NvU * idir];
        real_t *un  = &up[VLEN * NDF * nei  + NvU * idir];
        real_t *wpn = &wp[Nin5*nei];
        for (int is = 0; is < Ns; ++is) {
          mult_wilson_ymb(pg1_ym, pg2_ym, &vL[NVCD * is].v[0], u, un,
                          &wp2[Nin4 * is], &wpn[Nin4 * is]);
        }
      }

      {
        int idir = 2;
        if ( kz < Bz - 1 ) {
          int    nei = site + Nxv * Nyv;
          real_t *u   = &up[VLEN * NDF * site + NvU * idir];
          real_t *wpn = &wp[Nin5*nei];
          for (int is = 0; is < Ns; ++is) {
            mult_wilson_zpb(&vL[NVCD * is].v[0], u, &wpn[Nin4 * is]);
          }
        }
        if ( kz > 0 ) {
          int    nei = site - Nxv * Nyv;
          real_t *un  = &up[VLEN * NDF * nei  + NvU * idir];
          real_t *wpn = &wp[Nin5*nei];
          for (int is = 0; is < Ns; ++is) {
            mult_wilson_zmb(&vL[NVCD * is].v[0], un, &wpn[Nin4 * is]);
          }
        }
      }

      {
        int idir = 3;
        if (kt < Bt - 1) {
          int    nei = site + Nxv * Nyv * Nz;
          real_t *u   = &up[VLEN * NDF * site + NvU * idir];
          real_t *wpn = &wp[Nin5*nei];
          for (int is = 0; is < Ns; ++is) {
            mult_wilson_tpb_dirac(&vL[NVCD * is].v[0], u, &wpn[Nin4 * is]);
          }
        }

        if (kt > 0) {
          int    nei = site - Nxv * Nyv * Nz;
          real_t *un  = &up[VLEN * NDF * nei + NvU * idir];
          real_t *wpn = &wp[Nin5 * nei];
          for (int is = 0; is < Ns; ++is) {
            mult_wilson_tmb_dirac(&vL[NVCD * is].v[0], un, &wpn[Nin4 * is]);
          }
        }

      }
      save_vec(vp2, vL, NVCD * Ns);
    }
  }
}

//============================================================END=====

#endif
