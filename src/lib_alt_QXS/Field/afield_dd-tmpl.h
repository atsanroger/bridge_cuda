/*!
        @file    afield_dd-inc.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2024-05-07 13:24:28 #$
        @version $LastChangedRevision: 2602 $
*/

#ifndef QXS_AFIELD_DD_INC_INCLUDED
#define QXS_AFIELD_DD_INC_INCLUDED

#include <cstdlib>

#include "complexTraits.h"
#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_QXS/inline/afield_th-inc.h"

//====================================================================
template<typename INDEX, typename AFIELD>
void block_dotc_eo_impl(typename AFIELD::complex_t *out,
                        const AFIELD& v, const AFIELD& w,
                        const int ieo, const INDEX& block_index)
{
  typedef typename AFIELD::real_t      real_t;
  typedef typename AFIELD::complex_t   complex_t;

  real_t *vp = const_cast<AFIELD *>(&v)->ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

  int Nin  = v.nin();
  int Nex  = v.nex();
  int Nstv = v.nvol() / VLEN;

  int Nxv = block_index.fine_lattice_size(0) / VLENX;
  int Nyv = block_index.fine_lattice_size(1) / VLENY;
  int Nz  = block_index.fine_lattice_size(2);

  int Nblock = block_index.coarse_nvol();
  int NBx    = block_index.coarse_lattice_size(0);
  int NBy    = block_index.coarse_lattice_size(1);
  int NBz    = block_index.coarse_lattice_size(2);
  int NBt    = block_index.coarse_lattice_size(3);

  int Bsize = block_index.block_nvol() / VLEN;
  int Bxv   = block_index.block_size(0) / VLENX;
  int Byv   = block_index.block_size(1) / VLENY;
  int Bz    = block_index.block_size(2);
  int Bt    = block_index.block_size(3);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nblock);

  int ieo_skip = 1 - ieo;
  for (int block = is; block < ns; ++block) {
    if (block_index.block_eo(block) == ieo_skip) {
      out[block] = cmplx(real_t(0.0), real_t(0.0));
      continue;
    }

    int ibx = block % NBx;
    int iby = (block / NBx) % NBy;
    int ibz = (block / (NBx * NBy)) % NBz;
    int ibt = block / (NBx * NBy * NBz);

    svbool_t pg = set_predicate();
    svreal_t ytr, yti;
    set_vec(pg, ytr, 0.0);
    set_vec(pg, yti, 0.0);

    for (int bsite = 0; bsite < Bsize; ++bsite) {
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

      svbool_t pg = set_predicate();
      svreal_t xtr, xti;
      set_vec(pg, xtr, 0.0);
      set_vec(pg, xti, 0.0);

      for (int ex = 0; ex < Nex; ++ex) {
        for (int in2 = 0; in2 < Nin / 2; ++in2) {
          int inr = 2 * in2;
          int ini = 2 * in2 + 1;

          svreal_t vtr, vti, wtr, wti;
          load_vec(pg, vtr, &vp[VLEN * (inr + Nin * (site + Nstv * ex))]);
          load_vec(pg, vti, &vp[VLEN * (ini + Nin * (site + Nstv * ex))]);
          load_vec(pg, wtr, &wp[VLEN * (inr + Nin * (site + Nstv * ex))]);
          load_vec(pg, wti, &wp[VLEN * (ini + Nin * (site + Nstv * ex))]);
          add_dot_vec(pg, xtr, vtr, wtr);
          add_dot_vec(pg, xtr, vti, wti);
          add_dot_vec(pg, xti, vtr, wti);
          sub_dot_vec(pg, xti, vti, wtr);
        }
      }

      add_vec(pg, ytr, xtr);
      add_vec(pg, yti, xti);
    }

    real_t atr, ati;
    reduce_vec(pg, atr, ytr);
    reduce_vec(pg, ati, yti);
    out[block] = cmplx(atr, ati);
  }
}
//====================================================================
template<typename INDEX, typename AFIELD>
void block_norm2_eo_impl(typename AFIELD::real_t *out,
                    const AFIELD& v, const int ieo,
                    const INDEX& block_index)
{
  typedef typename AFIELD::real_t real_t;

  real_t *vp = const_cast<AFIELD *>(&v)->ptr(0);

  int Nin  = v.nin();
  int Nex  = v.nex();
  int Nstv = v.nvol() / VLEN;

  int Nxv = block_index.fine_lattice_size(0) / VLENX;
  int Nyv = block_index.fine_lattice_size(1) / VLENY;
  int Nz  = block_index.fine_lattice_size(2);

  int Nblock = block_index.coarse_nvol();
  int NBx    = block_index.coarse_lattice_size(0);
  int NBy    = block_index.coarse_lattice_size(1);
  int NBz    = block_index.coarse_lattice_size(2);
  int NBt    = block_index.coarse_lattice_size(3);

  int Bsize = block_index.block_nvol() / VLEN;
  int Bxv   = block_index.block_size(0) / VLENX;
  int Byv   = block_index.block_size(1) / VLENY;
  int Bz    = block_index.block_size(2);
  int Bt    = block_index.block_size(3);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nblock);
  int ieo_skip = 1 - ieo;

  for (int block = is; block < ns; ++block) {
    if (block_index.block_eo(block) == ieo_skip) {
      // out[block] = real_t(0.0);
      continue;
    }

    int ibx = block % NBx;
    int iby = (block / NBx) % NBy;
    int ibz = (block / (NBx * NBy)) % NBz;
    int ibt = block / (NBx * NBy * NBz);

    svbool_t pg = set_predicate();
    svreal_t yt;
    set_vec(pg, yt, 0.0);

    for (int bsite = 0; bsite < Bsize; ++bsite) {
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

      svreal_t xt, vt;
      set_vec(pg, xt, 0.0);

      for (int ex = 0; ex < Nex; ++ex) {
        for (int in = 0; in < Nin; ++in) {
          load_vec(pg, vt, &vp[VLEN * (in + Nin * (site + Nstv * ex))]);
          add_dot_vec(pg, xt, vt, vt);
        }
      }
      add_vec(pg, yt, xt);
    }
    real_t a;
    reduce_vec(pg, a, yt);
    out[block] = a;
  }
}


//====================================================================
template<typename INDEX, typename AFIELD>
void block_scal_eo_impl(AFIELD& v, const typename AFIELD::real_t *a,
                        const int ieo, const INDEX& block_index)
{
  typedef typename AFIELD::real_t real_t;

  real_t *vp = v.ptr(0);

  int Nin  = v.nin();
  int Nex  = v.nex();
  int Nstv = v.nvol() / VLEN;

  int Nxv = block_index.fine_lattice_size(0) / VLENX;
  int Nyv = block_index.fine_lattice_size(1) / VLENY;
  int Nz  = block_index.fine_lattice_size(2);

  int Nblock = block_index.coarse_nvol();
  int NBx    = block_index.coarse_lattice_size(0);
  int NBy    = block_index.coarse_lattice_size(1);
  int NBz    = block_index.coarse_lattice_size(2);
  int NBt    = block_index.coarse_lattice_size(3);

  int Bsize = block_index.block_nvol() / VLEN;
  int Bxv   = block_index.block_size(0) / VLENX;
  int Byv   = block_index.block_size(1) / VLENY;
  int Bz    = block_index.block_size(2);
  int Bt    = block_index.block_size(3);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nblock);
  int ieo_skip = 1 - ieo;

  svbool_t pg = set_predicate();

  for (int block = is; block < ns; ++block) {
    if (block_index.block_eo(block) == ieo_skip) continue;
    int ibx = block % NBx;
    int iby = (block / NBx) % NBy;
    int ibz = (block / (NBx * NBy)) % NBz;
    int ibt = block / (NBx * NBy * NBz);

    real_t at = a[block];
    for (int bsite = 0; bsite < Bsize; ++bsite) {
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

      for (int ex = 0; ex < Nex; ++ex) {
        for (int in = 0; in < Nin; ++in) {
          svreal_t vt;
          load_vec(pg, vt, &vp[VLEN * (in + Nin * (site + Nstv * ex))]);
          scal_vec(pg, vt, at);
          save_vec(pg, &vp[VLEN * (in + Nin * (site + Nstv * ex))], vt);
        }
      }
    }
  }
}


//====================================================================
template<typename INDEX, typename AFIELD>
void block_scal_eo_impl(AFIELD& v, const typename AFIELD::complex_t *a,
                        const int ieo, const INDEX& block_index)
{
  typedef typename AFIELD::real_t    real_t;
  typedef typename AFIELD::complex_t complex_t;

  real_t *vp = v.ptr(0);

  int Nin  = v.nin();
  int Nex  = v.nex();
  int Nstv = v.nvol() / VLEN;

  int Nxv = block_index.fine_lattice_size(0) / VLENX;
  int Nyv = block_index.fine_lattice_size(1) / VLENY;
  int Nz  = block_index.fine_lattice_size(2);

  int Nblock = block_index.coarse_nvol();
  int NBx    = block_index.coarse_lattice_size(0);
  int NBy    = block_index.coarse_lattice_size(1);
  int NBz    = block_index.coarse_lattice_size(2);
  int NBt    = block_index.coarse_lattice_size(3);

  int Bsize = block_index.block_nvol() / VLEN;
  int Bxv   = block_index.block_size(0) / VLENX;
  int Byv   = block_index.block_size(1) / VLENY;
  int Bz    = block_index.block_size(2);
  int Bt    = block_index.block_size(3);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nblock);
  int ieo_skip = 1 - ieo;

  svbool_t pg = set_predicate();

  for (int block = is; block < ns; ++block) {
    if (block_index.block_eo(block) == ieo_skip) continue;
    int ibx = block % NBx;
    int iby = (block / NBx) % NBy;
    int ibz = (block / (NBx * NBy)) % NBz;
    int ibt = block / (NBx * NBy * NBz);

    for (int bsite = 0; bsite < Bsize; ++bsite) {
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

      complex_t at  = a[block];
      real_t    atr = real(at);
      real_t    ati = imag(at);

      for (int ex = 0; ex < Nex; ++ex) {
        for (int in2 = 0; in2 < Nin / 2; ++in2) {
          int      inr = 2 * in2;
          int      ini = 2 * in2 + 1;
          svreal_t vtr, vti, wtr, wti;
          load_vec(pg, wtr, &vp[VLEN * (inr + Nin * (site + Nstv * ex))]);
          load_vec(pg, wti, &vp[VLEN * (ini + Nin * (site + Nstv * ex))]);
          set_vec(pg, vtr, atr, wtr);
          axpy_vec(pg, vtr, -ati, wti);
          set_vec(pg, vti, atr, wti);
          axpy_vec(pg, vti, ati, wtr);
          save_vec(pg, &vp[VLEN * (inr + Nin * (site + Nstv * ex))], vtr);
          save_vec(pg, &vp[VLEN * (ini + Nin * (site + Nstv * ex))], vti);
        }
      }
    }
  }
}


//====================================================================
template<typename INDEX, typename AFIELD>
void block_axpy_eo_impl(AFIELD& v, const typename AFIELD::real_t *a,
                        const AFIELD& w, const int ieo,
                        const typename AFIELD::real_t fac,
                        const INDEX& block_index)
{
  typedef typename AFIELD::real_t real_t;

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

  int Nin  = v.nin();
  int Nex  = v.nex();
  int Nstv = v.nvol() / VLEN;

  int Nxv = block_index.fine_lattice_size(0) / VLENX;
  int Nyv = block_index.fine_lattice_size(1) / VLENY;
  int Nz  = block_index.fine_lattice_size(2);

  int Nblock = block_index.coarse_nvol();
  int NBx    = block_index.coarse_lattice_size(0);
  int NBy    = block_index.coarse_lattice_size(1);
  int NBz    = block_index.coarse_lattice_size(2);
  int NBt    = block_index.coarse_lattice_size(3);

  int Bsize = block_index.block_nvol() / VLEN;
  int Bxv   = block_index.block_size(0) / VLENX;
  int Byv   = block_index.block_size(1) / VLENY;
  int Bz    = block_index.block_size(2);
  int Bt    = block_index.block_size(3);

  svbool_t pg = set_predicate();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nblock);
  int ieo_skip = 1 - ieo;

  for (int block = is; block < ns; ++block) {
    if (block_index.block_eo(block) == ieo_skip) continue;
    int ibx = block % NBx;
    int iby = (block / NBx) % NBy;
    int ibz = (block / (NBx * NBy)) % NBz;
    int ibt = block / (NBx * NBy * NBz);

    for (int bsite = 0; bsite < Bsize; ++bsite) {
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

      real_t at = fac * a[block];

      for (int ex = 0; ex < Nex; ++ex) {
        for (int in = 0; in < Nin; ++in) {
          svreal_t wt, vt;
          load_vec(pg, wt, &wp[VLEN * (in + Nin * (site + Nstv * ex))]);
          load_vec(pg, vt, &vp[VLEN * (in + Nin * (site + Nstv * ex))]);
          axpy_vec(pg, vt, at, wt);
          save_vec(pg, &vp[VLEN * (in + Nin * (site + Nstv * ex))], vt);
        }
      }
    }
  }
}


//====================================================================
template<typename INDEX, typename AFIELD>
void block_axpy_eo_impl(AFIELD& v, const typename AFIELD::complex_t *a,
                        const AFIELD& w, const int ieo,
                        const typename AFIELD::real_t fac,
                        const INDEX& block_index)
{
  typedef typename AFIELD::real_t    real_t;
  typedef typename AFIELD::complex_t complex_t;

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

  int Nin  = v.nin();
  int Nex  = v.nex();
  int Nstv = v.nvol() / VLEN;

  int Nxv = block_index.fine_lattice_size(0) / VLENX;
  int Nyv = block_index.fine_lattice_size(1) / VLENY;
  int Nz  = block_index.fine_lattice_size(2);

  int Nblock = block_index.coarse_nvol();
  int NBx    = block_index.coarse_lattice_size(0);
  int NBy    = block_index.coarse_lattice_size(1);
  int NBz    = block_index.coarse_lattice_size(2);
  int NBt    = block_index.coarse_lattice_size(3);

  int Bsize = block_index.block_nvol() / VLEN;
  int Bxv   = block_index.block_size(0) / VLENX;
  int Byv   = block_index.block_size(1) / VLENY;
  int Bz    = block_index.block_size(2);
  int Bt    = block_index.block_size(3);

  svbool_t pg = set_predicate();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nblock);
  int ieo_skip = 1 - ieo;

  for (int block = is; block < ns; ++block) {
    if (block_index.block_eo(block) == ieo_skip) continue;
    int ibx = block % NBx;
    int iby = (block / NBx) % NBy;
    int ibz = (block / (NBx * NBy)) % NBz;
    int ibt = block / (NBx * NBy * NBz);

    for (int bsite = 0; bsite < Bsize; ++bsite) {
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

      complex_t at  = a[block];
      real_t    atr = fac * real(at);
      real_t    ati = fac * imag(at);

      for (int ex = 0; ex < Nex; ++ex) {
        for (int in2 = 0; in2 < Nin / 2; ++in2) {
          int      inr = 2 * in2;
          int      ini = 2 * in2 + 1;
          svreal_t wtr, wti, vtr, vti;
          load_vec(pg, wtr, &wp[VLEN * (inr + Nin * (site + Nstv * ex))]);
          load_vec(pg, wti, &wp[VLEN * (ini + Nin * (site + Nstv * ex))]);
          load_vec(pg, vtr, &vp[VLEN * (inr + Nin * (site + Nstv * ex))]);
          load_vec(pg, vti, &vp[VLEN * (ini + Nin * (site + Nstv * ex))]);
          axpy_vec(pg, vtr, atr, wtr);
          axpy_vec(pg, vtr, -ati, wti);
          axpy_vec(pg, vti, atr, wti);
          axpy_vec(pg, vti, ati, wtr);
          save_vec(pg, &vp[VLEN * (inr + Nin * (site + Nstv * ex))], vtr);
          save_vec(pg, &vp[VLEN * (ini + Nin * (site + Nstv * ex))], vti);
        }
      }
    }
  }
}

//============================================================END=====
#endif
