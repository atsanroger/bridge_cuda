/*!
        @file    afield_Gauge-inc.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-10-17 13:51:32 #$
        @version $LastChangedRevision: 2550 $
*/

#ifndef QXS_AFIELD_GAUGE_INC_INCLUDED
#define QXS_AFIELD_GAUGE_INC_INCLUDED

#include <cstdlib> 

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/aindex_lex.h"
#include "lib_alt_QXS/inline/afield_th-inc.h"

#include "lib_alt_QXS/Field/afield_Gauge_index-inc.h"

//namespace {
//  inline int idxGr(int ic1, int ic2){ return     2*(ic1 + NC * ic2); }
//  inline int idxGi(int ic1, int ic2){ return 1 + 2*(ic1 + NC * ic2); }
//}

namespace QXS_Gauge{

//====================================================================
template <typename REALTYPE>
void set_boundary(AField<REALTYPE, QXS>& ulex,
                  const std::vector<int>& boundary)
{
  typedef REALTYPE real_t;
  
#pragma omp barrier

  int Nx = CommonParameters::Nx();
  int Ny = CommonParameters::Ny();
  int Nz = CommonParameters::Nz();
  int Nt = CommonParameters::Nt();
  int Ndim = CommonParameters::Ndim();
  int Nvol = Nx * Ny * Nz * Nt;

  int Nin = ulex.nin();

  if(!ulex.check_size(Nin, Nvol, Ndim)){
    vout.crucial("set_boundary: wrong size of input field\n");
    exit(EXIT_FAILURE);
  }

  AIndex_lex<real_t, QXS> index;

  int mu = 0;
  int ipex = Communicator::ipe(mu);
  int npex = Communicator::npe(mu);

  if(boundary[mu] != 1 && ipex == npex-1){

    real_t bc = real_t(boundary[mu]);
    int Nyzt  = Ny * Nz * Nt;

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, Nyzt);

    for(int iyzt = is; iyzt < ns; ++iyzt){
      int site = Nx-1 + Nx * iyzt;
      for(int in = 0; in < Nin; ++in){
        int idx = index.idx_G(in, site, mu);
        real_t uv = ulex.cmp(idx);
        uv = uv * bc;
        ulex.set(idx, uv);
      }
    }
  }

  mu = 1;
  int ipey = Communicator::ipe(mu);
  int npey = Communicator::npe(mu);

  if(boundary[mu] != 1 && ipey == npey-1){

    real_t bc = real_t(boundary[mu]);
    int Nxzt  = Nx * Nz * Nt;

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, Nxzt);

    for(int ixzt = is; ixzt < ns; ++ixzt){
      int ix  = ixzt % Nx;
      int izt = ixzt / Nx;
      int site = ix + Nx * (Ny-1 + Ny * izt);
      for(int in = 0; in < Nin; ++in){
        int idx = index.idx_G(in, site, mu);
        real_t uv = ulex.cmp(idx);
        uv = uv * bc;
        ulex.set(idx, uv);
      }
    }
  }

  mu = 2;
  int ipez = Communicator::ipe(mu);
  int npez = Communicator::npe(mu);

  if(boundary[mu] != 1 && ipez == npez-1){

    real_t bc = real_t(boundary[mu]);
    int Nxy   = Nx * Ny;
    int Nxyt  =  Nxy * Nt;

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, Nxyt);

    for(int ixyt = is; ixyt < ns; ++ixyt){
      int ixy = ixyt % Nxy;
      int it  = ixyt / Nxy;
      int site = ixy + Nxy * (Nz-1 + Nz * it);
      for(int in = 0; in < Nin; ++in){
        int idx = index.idx_G(in, site, mu);
        real_t uv = ulex.cmp(idx);
        uv = uv * bc;
        ulex.set(idx, uv);
      }
    }
  }

  mu = 3;
  int ipet = Communicator::ipe(mu);
  int npet = Communicator::npe(mu);

  if(boundary[mu] != 1 && ipet == npet-1){

    real_t bc = real_t(boundary[mu]);
    int Nxyz  = Nx * Ny * Nz;

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, Nxyz);

    for(int ixyz = is; ixyz < ns; ++ixyz){
      int site = ixyz + Nxyz * (Nt-1);
      for(int in = 0; in < Nin; ++in){
        int idx = index.idx_G(in, site, mu);
        real_t uv = ulex.cmp(idx);
        uv = uv * bc;
        ulex.set(idx, uv);
      }
    }
  }

#pragma omp barrier

}

//====================================================================
template <typename REALTYPE>
void mult_Gnn(AField<REALTYPE,QXS>& u,       const int exu,
              const AField<REALTYPE,QXS>& v, const int exv,
              const AField<REALTYPE,QXS>& w, const int exw)
{
#pragma omp barrier

  int Nst  = w.nvol();
  int Nstv = Nst/VLEN;

  typedef AField<REALTYPE,QXS> AFIELD;

  REALTYPE* up = u.ptr(NDF * Nst * exu);
  REALTYPE* vp = const_cast<AFIELD*>(&v)->ptr(NDF * Nst * exv);
  REALTYPE* wp = const_cast<AFIELD*>(&w)->ptr(NDF * Nst * exw);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstv);

  svbool_t pg = set_predicate();

  for(int site = is; site < ns; ++site){
    REALTYPE* upt = &up[VLEN * NDF * site];
    REALTYPE* vpt = &vp[VLEN * NDF * site];
    REALTYPE* wpt = &wp[VLEN * NDF * site];

    svreal_t w00r, w00i, w01r, w01i, w02r, w02i;
    load_vec(pg, w00r, &wpt[VLEN * idxGr(0,0)]);
    load_vec(pg, w00i, &wpt[VLEN * idxGi(0,0)]);
    load_vec(pg, w01r, &wpt[VLEN * idxGr(0,1)]);
    load_vec(pg, w01i, &wpt[VLEN * idxGi(0,1)]);
    load_vec(pg, w02r, &wpt[VLEN * idxGr(0,2)]);
    load_vec(pg, w02i, &wpt[VLEN * idxGi(0,2)]);

    svreal_t w10r, w10i, w11r, w11i, w12r, w12i;
    load_vec(pg, w10r, &wpt[VLEN * idxGr(1,0)]);
    load_vec(pg, w10i, &wpt[VLEN * idxGi(1,0)]);
    load_vec(pg, w11r, &wpt[VLEN * idxGr(1,1)]);
    load_vec(pg, w11i, &wpt[VLEN * idxGi(1,1)]);
    load_vec(pg, w12r, &wpt[VLEN * idxGr(1,2)]);
    load_vec(pg, w12i, &wpt[VLEN * idxGi(1,2)]);

    svreal_t w20r, w20i, w21r, w21i, w22r, w22i;
    load_vec(pg, w20r, &wpt[VLEN * idxGr(2,0)]);
    load_vec(pg, w20i, &wpt[VLEN * idxGi(2,0)]);
    load_vec(pg, w21r, &wpt[VLEN * idxGr(2,1)]);
    load_vec(pg, w21i, &wpt[VLEN * idxGi(2,1)]);
    load_vec(pg, w22r, &wpt[VLEN * idxGr(2,2)]);
    load_vec(pg, w22i, &wpt[VLEN * idxGi(2,2)]);

    for(int ic1 = 0; ic1 < NC; ++ic1){

      svreal_t v0r, v0i, v1r, v1i, v2r, v2i;
      load_vec(pg, v0r, &vpt[VLEN * idxGr(ic1,0)]);
      load_vec(pg, v0i, &vpt[VLEN * idxGi(ic1,0)]);
      load_vec(pg, v1r, &vpt[VLEN * idxGr(ic1,1)]);
      load_vec(pg, v1i, &vpt[VLEN * idxGi(ic1,1)]);
      load_vec(pg, v2r, &vpt[VLEN * idxGr(ic1,2)]);
      load_vec(pg, v2i, &vpt[VLEN * idxGi(ic1,2)]);

      svreal_t ut0r, ut0i, ut1r, ut1i, ut2r, ut2i;
      mul_vec( pg, ut0r, v0r, w00r);
      mul_vec( pg, ut0i, v0r, w00i);
      mul_vec( pg, ut1r, v0r, w01r);
      mul_vec( pg, ut1i, v0r, w01i);
      mul_vec( pg, ut2r, v0r, w02r);
      mul_vec( pg, ut2i, v0r, w02i);

      ymax_vec(pg, ut0r, v0i, w00i);
      axpy_vec(pg, ut0i, v0i, w00r);
      ymax_vec(pg, ut1r, v0i, w01i);
      axpy_vec(pg, ut1i, v0i, w01r);
      ymax_vec(pg, ut2r, v0i, w02i);
      axpy_vec(pg, ut2i, v0i, w02r);

      axpy_vec(pg, ut0r, v1r, w10r);
      axpy_vec(pg, ut0i, v1r, w10i);
      axpy_vec(pg, ut1r, v1r, w11r);
      axpy_vec(pg, ut1i, v1r, w11i);
      axpy_vec(pg, ut2r, v1r, w12r);
      axpy_vec(pg, ut2i, v1r, w12i);

      ymax_vec(pg, ut0r, v1i, w10i);
      axpy_vec(pg, ut0i, v1i, w10r);
      ymax_vec(pg, ut1r, v1i, w11i);
      axpy_vec(pg, ut1i, v1i, w11r);
      ymax_vec(pg, ut2r, v1i, w12i);
      axpy_vec(pg, ut2i, v1i, w12r);

      axpy_vec(pg, ut0r, v2r, w20r);
      axpy_vec(pg, ut0i, v2r, w20i);
      axpy_vec(pg, ut1r, v2r, w21r);
      axpy_vec(pg, ut1i, v2r, w21i);
      axpy_vec(pg, ut2r, v2r, w22r);
      axpy_vec(pg, ut2i, v2r, w22i);

      ymax_vec(pg, ut0r, v2i, w20i);
      axpy_vec(pg, ut0i, v2i, w20r);
      ymax_vec(pg, ut1r, v2i, w21i);
      axpy_vec(pg, ut1i, v2i, w21r);
      ymax_vec(pg, ut2r, v2i, w22i);
      axpy_vec(pg, ut2i, v2i, w22r);

      save_vec(pg, &upt[VLEN * idxGr(ic1,0)], ut0r);
      save_vec(pg, &upt[VLEN * idxGi(ic1,0)], ut0i);
      save_vec(pg, &upt[VLEN * idxGr(ic1,1)], ut1r);
      save_vec(pg, &upt[VLEN * idxGi(ic1,1)], ut1i);
      save_vec(pg, &upt[VLEN * idxGr(ic1,2)], ut2r);
      save_vec(pg, &upt[VLEN * idxGi(ic1,2)], ut2i);
    }

  }

#pragma omp barrier
}

//====================================================================
template <typename REALTYPE>
void mult_Gnd(AField<REALTYPE,QXS>& u,       const int exu,
              const AField<REALTYPE,QXS>& v, const int exv,
              const AField<REALTYPE,QXS>& w, const int exw)
{
#pragma omp barrier

  int Nst  = w.nvol();
  int Nstv = Nst/VLEN;

  typedef AField<REALTYPE,QXS> AFIELD;

  REALTYPE* up = u.ptr(NDF * Nst * exu);
  REALTYPE* vp = const_cast<AFIELD*>(&v)->ptr(NDF * Nst * exv);
  REALTYPE* wp = const_cast<AFIELD*>(&w)->ptr(NDF * Nst * exw);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstv);

  svbool_t pg = set_predicate();

  for(int site = is; site < ns; ++site){
    REALTYPE* upt = &up[VLEN * NDF * site];
    REALTYPE* vpt = &vp[VLEN * NDF * site];
    REALTYPE* wpt = &wp[VLEN * NDF * site];

    svreal_t w00r, w00i, w01r, w01i, w02r, w02i;
    load_vec(pg, w00r, &wpt[VLEN * idxGr(0,0)]);
    load_vec(pg, w00i, &wpt[VLEN * idxGi(0,0)]);
    load_vec(pg, w01r, &wpt[VLEN * idxGr(0,1)]);
    load_vec(pg, w01i, &wpt[VLEN * idxGi(0,1)]);
    load_vec(pg, w02r, &wpt[VLEN * idxGr(0,2)]);
    load_vec(pg, w02i, &wpt[VLEN * idxGi(0,2)]);

    svreal_t w10r, w10i, w11r, w11i, w12r, w12i;
    load_vec(pg, w10r, &wpt[VLEN * idxGr(1,0)]);
    load_vec(pg, w10i, &wpt[VLEN * idxGi(1,0)]);
    load_vec(pg, w11r, &wpt[VLEN * idxGr(1,1)]);
    load_vec(pg, w11i, &wpt[VLEN * idxGi(1,1)]);
    load_vec(pg, w12r, &wpt[VLEN * idxGr(1,2)]);
    load_vec(pg, w12i, &wpt[VLEN * idxGi(1,2)]);

    svreal_t w20r, w20i, w21r, w21i, w22r, w22i;
    load_vec(pg, w20r, &wpt[VLEN * idxGr(2,0)]);
    load_vec(pg, w20i, &wpt[VLEN * idxGi(2,0)]);
    load_vec(pg, w21r, &wpt[VLEN * idxGr(2,1)]);
    load_vec(pg, w21i, &wpt[VLEN * idxGi(2,1)]);
    load_vec(pg, w22r, &wpt[VLEN * idxGr(2,2)]);
    load_vec(pg, w22i, &wpt[VLEN * idxGi(2,2)]);

    for(int ic1 = 0; ic1 < NC; ++ic1){

      svreal_t v0r, v0i, v1r, v1i, v2r, v2i;
      load_vec(pg, v0r, &vpt[VLEN * idxGr(ic1,0)]);
      load_vec(pg, v0i, &vpt[VLEN * idxGi(ic1,0)]);
      load_vec(pg, v1r, &vpt[VLEN * idxGr(ic1,1)]);
      load_vec(pg, v1i, &vpt[VLEN * idxGi(ic1,1)]);
      load_vec(pg, v2r, &vpt[VLEN * idxGr(ic1,2)]);
      load_vec(pg, v2i, &vpt[VLEN * idxGi(ic1,2)]);

      svreal_t ut0r, ut0i, ut1r, ut1i, ut2r, ut2i;
      mul_vec( pg, ut0r, v0r, w00r);
      mul_vec( pg, ut0i, v0i, w00r);
      mul_vec( pg, ut1r, v0r, w10r);
      mul_vec( pg, ut1i, v0i, w10r);
      mul_vec( pg, ut2r, v0r, w20r);
      mul_vec( pg, ut2i, v0i, w20r);

      axpy_vec(pg, ut0r, v0i, w00i);
      ymax_vec(pg, ut0i, v0r, w00i);
      axpy_vec(pg, ut1r, v0i, w10i);
      ymax_vec(pg, ut1i, v0r, w10i);
      axpy_vec(pg, ut2r, v0i, w20i);
      ymax_vec(pg, ut2i, v0r, w20i);

      axpy_vec(pg, ut0r, v1r, w01r);
      axpy_vec(pg, ut0i, v1i, w01r);
      axpy_vec(pg, ut1r, v1r, w11r);
      axpy_vec(pg, ut1i, v1i, w11r);
      axpy_vec(pg, ut2r, v1r, w21r);
      axpy_vec(pg, ut2i, v1i, w21r);

      axpy_vec(pg, ut0r, v1i, w01i);
      ymax_vec(pg, ut0i, v1r, w01i);
      axpy_vec(pg, ut1r, v1i, w11i);
      ymax_vec(pg, ut1i, v1r, w11i);
      axpy_vec(pg, ut2r, v1i, w21i);
      ymax_vec(pg, ut2i, v1r, w21i);

      axpy_vec(pg, ut0r, v2r, w02r);
      axpy_vec(pg, ut0i, v2i, w02r);
      axpy_vec(pg, ut1r, v2r, w12r);
      axpy_vec(pg, ut1i, v2i, w12r);
      axpy_vec(pg, ut2r, v2r, w22r);
      axpy_vec(pg, ut2i, v2i, w22r);

      axpy_vec(pg, ut0r, v2i, w02i);
      ymax_vec(pg, ut0i, v2r, w02i);
      axpy_vec(pg, ut1r, v2i, w12i);
      ymax_vec(pg, ut1i, v2r, w12i);
      axpy_vec(pg, ut2r, v2i, w22i);
      ymax_vec(pg, ut2i, v2r, w22i);

      save_vec(pg, &upt[VLEN * idxGr(ic1,0)], ut0r);
      save_vec(pg, &upt[VLEN * idxGi(ic1,0)], ut0i);
      save_vec(pg, &upt[VLEN * idxGr(ic1,1)], ut1r);
      save_vec(pg, &upt[VLEN * idxGi(ic1,1)], ut1i);
      save_vec(pg, &upt[VLEN * idxGr(ic1,2)], ut2r);
      save_vec(pg, &upt[VLEN * idxGi(ic1,2)], ut2i);
    }

  }

#pragma omp barrier
}

//====================================================================
template <typename REALTYPE>
void mult_Gdn(AField<REALTYPE,QXS>& u,       const int exu,
              const AField<REALTYPE,QXS>& v, const int exv,
              const AField<REALTYPE,QXS>& w, const int exw)
{
#pragma omp barrier

  int Nst  = w.nvol();
  int Nstv = Nst/VLEN;

  typedef AField<REALTYPE,QXS> AFIELD;

  REALTYPE* up = u.ptr(NDF * Nst * exu);
  REALTYPE* vp = const_cast<AFIELD*>(&v)->ptr(NDF * Nst * exv);
  REALTYPE* wp = const_cast<AFIELD*>(&w)->ptr(NDF * Nst * exw);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstv);

  svbool_t pg = set_predicate();

  for(int site = is; site < ns; ++site){
    REALTYPE* upt = &up[VLEN * NDF * site];
    REALTYPE* vpt = &vp[VLEN * NDF * site];
    REALTYPE* wpt = &wp[VLEN * NDF * site];

    svreal_t w00r, w00i, w01r, w01i, w02r, w02i;
    load_vec(pg, w00r, &wpt[VLEN * idxGr(0,0)]);
    load_vec(pg, w00i, &wpt[VLEN * idxGi(0,0)]);
    load_vec(pg, w01r, &wpt[VLEN * idxGr(0,1)]);
    load_vec(pg, w01i, &wpt[VLEN * idxGi(0,1)]);
    load_vec(pg, w02r, &wpt[VLEN * idxGr(0,2)]);
    load_vec(pg, w02i, &wpt[VLEN * idxGi(0,2)]);

    svreal_t w10r, w10i, w11r, w11i, w12r, w12i;
    load_vec(pg, w10r, &wpt[VLEN * idxGr(1,0)]);
    load_vec(pg, w10i, &wpt[VLEN * idxGi(1,0)]);
    load_vec(pg, w11r, &wpt[VLEN * idxGr(1,1)]);
    load_vec(pg, w11i, &wpt[VLEN * idxGi(1,1)]);
    load_vec(pg, w12r, &wpt[VLEN * idxGr(1,2)]);
    load_vec(pg, w12i, &wpt[VLEN * idxGi(1,2)]);

    svreal_t w20r, w20i, w21r, w21i, w22r, w22i;
    load_vec(pg, w20r, &wpt[VLEN * idxGr(2,0)]);
    load_vec(pg, w20i, &wpt[VLEN * idxGi(2,0)]);
    load_vec(pg, w21r, &wpt[VLEN * idxGr(2,1)]);
    load_vec(pg, w21i, &wpt[VLEN * idxGi(2,1)]);
    load_vec(pg, w22r, &wpt[VLEN * idxGr(2,2)]);
    load_vec(pg, w22i, &wpt[VLEN * idxGi(2,2)]);

    for(int ic1 = 0; ic1 < NC; ++ic1){

      svreal_t v0r, v0i, v1r, v1i, v2r, v2i;
      load_vec(pg, v0r, &vpt[VLEN * idxGr(0,ic1)]);
      load_vec(pg, v0i, &vpt[VLEN * idxGi(0,ic1)]);
      load_vec(pg, v1r, &vpt[VLEN * idxGr(1,ic1)]);
      load_vec(pg, v1i, &vpt[VLEN * idxGi(1,ic1)]);
      load_vec(pg, v2r, &vpt[VLEN * idxGr(2,ic1)]);
      load_vec(pg, v2i, &vpt[VLEN * idxGi(2,ic1)]);

      svreal_t ut0r, ut0i, ut1r, ut1i, ut2r, ut2i;
      mul_vec( pg, ut0r, v0r, w00r);
      mul_vec( pg, ut0i, v0r, w00i);
      mul_vec( pg, ut1r, v0r, w01r);
      mul_vec( pg, ut1i, v0r, w01i);
      mul_vec( pg, ut2r, v0r, w02r);
      mul_vec( pg, ut2i, v0r, w02i);

      axpy_vec(pg, ut0r, v0i, w00i);
      ymax_vec(pg, ut0i, v0i, w00r);
      axpy_vec(pg, ut1r, v0i, w01i);
      ymax_vec(pg, ut1i, v0i, w01r);
      axpy_vec(pg, ut2r, v0i, w02i);
      ymax_vec(pg, ut2i, v0i, w02r);

      axpy_vec(pg, ut0r, v1r, w10r);
      axpy_vec(pg, ut0i, v1r, w10i);
      axpy_vec(pg, ut1r, v1r, w11r);
      axpy_vec(pg, ut1i, v1r, w11i);
      axpy_vec(pg, ut2r, v1r, w12r);
      axpy_vec(pg, ut2i, v1r, w12i);

      axpy_vec(pg, ut0r, v1i, w10i);
      ymax_vec(pg, ut0i, v1i, w10r);
      axpy_vec(pg, ut1r, v1i, w11i);
      ymax_vec(pg, ut1i, v1i, w11r);
      axpy_vec(pg, ut2r, v1i, w12i);
      ymax_vec(pg, ut2i, v1i, w12r);

      axpy_vec(pg, ut0r, v2r, w20r);
      axpy_vec(pg, ut0i, v2r, w20i);
      axpy_vec(pg, ut1r, v2r, w21r);
      axpy_vec(pg, ut1i, v2r, w21i);
      axpy_vec(pg, ut2r, v2r, w22r);
      axpy_vec(pg, ut2i, v2r, w22i);

      axpy_vec(pg, ut0r, v2i, w20i);
      ymax_vec(pg, ut0i, v2i, w20r);
      axpy_vec(pg, ut1r, v2i, w21i);
      ymax_vec(pg, ut1i, v2i, w21r);
      axpy_vec(pg, ut2r, v2i, w22i);
      ymax_vec(pg, ut2i, v2i, w22r);

      save_vec(pg, &upt[VLEN * idxGr(ic1,0)], ut0r);
      save_vec(pg, &upt[VLEN * idxGi(ic1,0)], ut0i);
      save_vec(pg, &upt[VLEN * idxGr(ic1,1)], ut1r);
      save_vec(pg, &upt[VLEN * idxGi(ic1,1)], ut1i);
      save_vec(pg, &upt[VLEN * idxGr(ic1,2)], ut2r);
      save_vec(pg, &upt[VLEN * idxGi(ic1,2)], ut2i);
    }

  }

#pragma omp barrier
}

//====================================================================
template <typename REALTYPE>
void at_G(AField<REALTYPE,QXS>& u, const int ex)
{  // anti-Hermitian traceless

#pragma omp barrier

  int Nst  = u.nvol();
  int Nstv = Nst/VLEN;

  typedef AField<REALTYPE,QXS> AFIELD;

  REALTYPE* up = u.ptr(NDF * Nst * ex);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstv);

  svbool_t pg = set_predicate();
  svreal_t zero;
  set_vec(pg, zero, 0.0);

  for(int site = is; site < ns; ++site){

    REALTYPE* upt = &up[VLEN * NDF * site];

    svreal_t tri;
    set_vec(pg, tri, real_t(0.0));

    for(int ic1 = 0; ic1 < NC; ++ic1){
      for(int ic2 = ic1; ic2 < NC; ++ic2){
        if(ic1 != ic2){
          svreal_t ur, ui, unr, uni, utr, uti;
          load_vec(pg, unr, &upt[VLEN * idxGr(ic1,ic2)]);
          load_vec(pg, uni, &upt[VLEN * idxGi(ic1,ic2)]);
          load_vec(pg, utr, &upt[VLEN * idxGr(ic2,ic1)]);
          load_vec(pg, uti, &upt[VLEN * idxGi(ic2,ic1)]);
          sub_vec(pg, ur, unr, utr);
          add_vec(pg, ui, uni, uti);
          scal_vec(pg, ur, real_t(0.5));
          scal_vec(pg, ui, real_t(0.5));
          save_vec(pg, &upt[VLEN * idxGr(ic1,ic2)], ur);
          save_vec(pg, &upt[VLEN * idxGi(ic1,ic2)], ui);
          scal_vec(pg, ur, real_t(-1.0));
          save_vec(pg, &upt[VLEN * idxGr(ic2,ic1)], ur);
          save_vec(pg, &upt[VLEN * idxGi(ic2,ic1)], ui);
        }else{
          svreal_t uti;
          load_vec(pg, uti, &upt[VLEN * idxGi(ic1,ic1)]);
          add_vec(pg, tri, uti);
        }
      }
    }

    svreal_t uti;
    for(int ic1 = 0; ic1 < NC; ++ic1){
      load_vec(pg, uti, &upt[VLEN * idxGi(ic1,ic1)]);
      axpy_vec(pg, uti, -1.0/3.0, tri);
      save_vec(pg, &upt[VLEN * idxGr(ic1,ic1)], zero);
      save_vec(pg, &upt[VLEN * idxGi(ic1,ic1)], uti);
    }

  }

#pragma omp barrier

}

} // namespace QXS_Gauge

//============================================================END=====
#endif
