/*!
        @file    afield_dd-inc.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2024-04-01 12:13:09 #$
        @version $LastChangedRevision: 2595 $
*/

#ifndef ACCEL_AFIELD_DD_INC_INCLUDED
#define ACCEL_AFIELD_DD_INC_INCLUDED

#include <cstdlib> 

#include "lib_alt_Accel/Field/aindex_block_lex.h"
#include "lib_alt_Accel/BridgeACC/bridgeACC_AField_dd.h"

//====================================================================
template <typename INDEX, typename AFIELD>
  void block_dotc(typename AFIELD::complex_t *out,
                  const AFIELD& v, const AFIELD& w,
                  const INDEX& block_index)
{
  block_dotc_eo(out, v, w, -1, block_index);
}

//====================================================================
template <typename INDEX, typename AFIELD>
  void block_dotc_eo(typename AFIELD::complex_t *out,
                     const AFIELD& v, const AFIELD& w,
                     const int ieo, const INDEX& block_index)
{
  // Note that *out is assumed to be thread global array.

  typedef typename AFIELD::real_t real_t;
  typedef typename AFIELD::complex_t complex_t;

  int Nin = v.nin();
  int Nex = v.nex();
  int Nst = v.nvol();

  real_t *vp = const_cast<AFIELD*>(&v)->ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);

  real_t *buf_red = const_cast<INDEX*>(&block_index)->ptr_buf_red(0);
  real_t *buf_out = const_cast<INDEX*>(&block_index)->ptr_buf_out(0);

#pragma omp barrier

#pragma omp master
 {
   int Ndim = CommonParameters::Ndim();
   int Nsize[4], block_size[4];
   for(int mu = 0; mu < Ndim; ++mu){
     Nsize[mu]      = block_index.fine_lattice_size(mu);
     block_size[mu] = block_index.block_size(mu);
   }
   int Nblock = block_index.coarse_nvol();

   int Ieo_org = block_index.block_eo(0);
   int jeo = (ieo + Ieo_org) % 2;
   if(ieo < 0) jeo = ieo;

   BridgeACC::block_dotc_eo(buf_red, vp, wp, Nin, Nex,
                            //Nsize, block_size, ieo);
                            Nsize, block_size, jeo);

   BridgeACC::reduce_block(buf_out, buf_red, Nsize, block_size, jeo);

   const_cast<INDEX*>(&block_index)->buf_out_update_host();

   for(int block = 0; block < Nblock; ++block){
     out[block] = cmplx(buf_out[2*block], buf_out[2*block+1]);
   }

 }

#pragma omp barrier

}

//====================================================================
template <typename INDEX, typename AFIELD>
void block_norm2(typename AFIELD::real_t *out, const AFIELD& v,
                 const INDEX& block_index)
{
  block_norm2_eo(out, v, -1, block_index);
}

//====================================================================
template <typename INDEX, typename AFIELD>
void block_norm2_eo(typename AFIELD::real_t *out, const AFIELD& v,
                    const int ieo, const INDEX& block_index)
{
  typedef typename AFIELD::real_t real_t;
  typedef typename AFIELD::complex_t complex_t;

  int Nin = v.nin();
  int Nex = v.nex();
  int Nst = v.nvol();

  real_t *vp = const_cast<AFIELD*>(&v)->ptr(0);

  real_t *buf_red = const_cast<INDEX*>(&block_index)->ptr_buf_red(0);
  real_t *buf_out = const_cast<INDEX*>(&block_index)->ptr_buf_out(0);

#pragma omp barrier

#pragma omp master
 {
   int Ndim = CommonParameters::Ndim();
   int Nsize[4], block_size[4];
   for(int mu = 0; mu < Ndim; ++mu){
     Nsize[mu]      = block_index.fine_lattice_size(mu);
     block_size[mu] = block_index.block_size(mu);
   }
   int Nblock = block_index.coarse_nvol();

   int Ieo_org = block_index.block_eo(0);
   int jeo = (ieo + Ieo_org) % 2;
   if(ieo < 0) jeo = ieo;

   BridgeACC::block_norm2_eo(buf_red, vp, Nin, Nex,
                             Nsize, block_size, jeo);

   //  BridgeACC::reduce_block(buf_out, buf_red, Nblock, Nst, jeo);
   BridgeACC::reduce_block(buf_out, buf_red, Nsize, block_size, jeo);

   const_cast<INDEX*>(&block_index)->buf_out_update_host();

   for(int block = 0; block < Nblock; ++block){
     out[block] = buf_out[2*block];
   }

 }

#pragma omp barrier

}

//====================================================================
template <typename INDEX, typename AFIELD>
void block_axpy(AFIELD& v, const typename AFIELD::complex_t *a,
                const AFIELD& w, const typename AFIELD::real_t fac,
                const INDEX& block_index)
{
  block_axpy_eo(v, a, w, -1, fac, block_index);
}

//====================================================================
template <typename INDEX, typename AFIELD>
  void block_axpy_eo(AFIELD& v, const typename AFIELD::complex_t *a,
                     const AFIELD& w, const int ieo,
                     const typename AFIELD::real_t fac,
                     const INDEX& block_index)
{
  typedef typename AFIELD::real_t real_t;
  typedef typename AFIELD::complex_t complex_t;

  int Nin = v.nin();
  int Nex = v.nex();
  int Nst = v.nvol();

  real_t *vp = const_cast<AFIELD*>(&v)->ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);

  real_t *buf_out = const_cast<INDEX*>(&block_index)->ptr_buf_out(0);

#pragma omp barrier

#pragma omp master
 {
   int Ndim = CommonParameters::Ndim();
   int Nsize[4], block_size[4];
   for(int mu = 0; mu < Ndim; ++mu){
     Nsize[mu]      = block_index.fine_lattice_size(mu);
     block_size[mu] = block_index.block_size(mu);
   }
   int Nblock = block_index.coarse_nvol();

   int Ieo_org = block_index.block_eo(0);
   int jeo = (ieo + Ieo_org) % 2;
   if(ieo < 0) jeo = ieo;
   
   for(int block = 0; block < Nblock; ++block){
     buf_out[2*block]   = real(a[block]);
     buf_out[2*block+1] = imag(a[block]);
   }

   const_cast<INDEX*>(&block_index)->buf_out_update_device();

   BridgeACC::block_axpy_eo(vp, buf_out, wp, fac,
   //                       Nin, Nex, Nsize, block_size, ieo);
                            Nin, Nex, Nsize, block_size, jeo);

 }

#pragma omp barrier

}

//====================================================================
template <typename INDEX, typename AFIELD>
void block_axpy(AFIELD& v, const typename AFIELD::real_t *a,
                const AFIELD& w, const typename AFIELD::real_t fac,
                const INDEX& block_index)
{
  block_axpy_eo(v, a, w, -1, fac, block_index);
}

//====================================================================
template <typename INDEX, typename AFIELD>
  void block_axpy_eo(AFIELD& v, const typename AFIELD::real_t *a,
                     const AFIELD& w, const int ieo,
                     const typename AFIELD::real_t fac,
                     const INDEX& block_index)
{
  typedef typename AFIELD::real_t real_t;
  typedef typename AFIELD::complex_t complex_t;

  int Nin = v.nin();
  int Nex = v.nex();
  int Nst = v.nvol();

  real_t *vp = const_cast<AFIELD*>(&v)->ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);

  real_t *buf_out = const_cast<INDEX*>(&block_index)->ptr_buf_out(0);

#pragma omp barrier

#pragma omp master
 {
   int Ndim = CommonParameters::Ndim();
   int Nsize[4], block_size[4];
   for(int mu = 0; mu < Ndim; ++mu){
     Nsize[mu]      = block_index.fine_lattice_size(mu);
     block_size[mu] = block_index.block_size(mu);
   }
   int Nblock = block_index.coarse_nvol();

   int Ieo_org = block_index.block_eo(0);
   int jeo = (ieo + Ieo_org) % 2;
   if(ieo < 0) jeo = ieo;

   for(int block = 0; block < Nblock; ++block){
     buf_out[2*block]   = a[block];
     buf_out[2*block+1] = 0.0;
   }

   const_cast<INDEX*>(&block_index)->buf_out_update_device();

   BridgeACC::block_axpy_eo(vp, buf_out, wp, fac,
			    // Nin, Nex, Nsize, block_size, ieo);
                            Nin, Nex, Nsize, block_size, jeo);

 }

#pragma omp barrier

}

//====================================================================
template <typename INDEX, typename AFIELD>
void block_scal(const AFIELD& v, const typename AFIELD::real_t *a,
                const INDEX& block_index)
{
  block_scal_eo(v, a, -1, block_index);
}

//====================================================================
template <typename INDEX, typename AFIELD>
  void block_scal_eo(const AFIELD& v,
                     const typename AFIELD::real_t *a,
                     const int ieo, const INDEX& block_index)
{
  typedef typename AFIELD::real_t real_t;
  typedef typename AFIELD::complex_t complex_t;

  int Nin = v.nin();
  int Nex = v.nex();
  int Nst = v.nvol();

  real_t *vp = const_cast<AFIELD*>(&v)->ptr(0);

  real_t *buf_out = const_cast<INDEX*>(&block_index)->ptr_buf_out(0);

#pragma omp barrier

#pragma omp master
 {
   int Ndim = CommonParameters::Ndim();
   int Nsize[4], block_size[4];
   for(int mu = 0; mu < Ndim; ++mu){
     Nsize[mu]      = block_index.fine_lattice_size(mu);
     block_size[mu] = block_index.block_size(mu);
   }
   int Nblock = block_index.coarse_nvol();

   int Ieo_org = block_index.block_eo(0);
   int jeo = (ieo + Ieo_org) % 2;
   if(ieo < 0) jeo = ieo;

   for(int block = 0; block < Nblock; ++block){
     buf_out[2*block]   = a[block];
     buf_out[2*block+1] = 0.0;
   }

   const_cast<INDEX*>(&block_index)->buf_out_update_device();

   BridgeACC::block_scal_eo(vp, buf_out, Nin, Nex,
                            Nsize, block_size, jeo);

 }

#pragma omp barrier

}

//====================================================================
template <typename INDEX, typename AFIELD>
void block_scal(const AFIELD& v, const typename AFIELD::complex_t *a,
                const INDEX& block_index)
{
  block_scal_eo(v, a, -1, block_index);
}

//====================================================================
template <typename INDEX, typename AFIELD>
  void block_scal_eo(const AFIELD& v,
                     const typename AFIELD::complex_t *a,
                     const int ieo, const INDEX& block_index)
{
  typedef typename AFIELD::real_t real_t;
  typedef typename AFIELD::complex_t complex_t;

  int Nin = v.nin();
  int Nex = v.nex();
  int Nst = v.nvol();

  real_t *vp = const_cast<AFIELD*>(&v)->ptr(0);

  real_t *buf_out = const_cast<INDEX*>(&block_index)->ptr_buf_out(0);

#pragma omp barrier

#pragma omp master
 {
   int Ndim = CommonParameters::Ndim();
   int Nsize[4], block_size[4];
   for(int mu = 0; mu < Ndim; ++mu){
     Nsize[mu]      = block_index.fine_lattice_size(mu);
     block_size[mu] = block_index.block_size(mu);
   }
   int Nblock = block_index.coarse_nvol();

   int Ieo_org = block_index.block_eo(0);
   int jeo = (ieo + Ieo_org) % 2;
   if(ieo < 0) jeo = ieo;

   for(int block = 0; block < Nblock; ++block){
     buf_out[2*block]   = real(a[block]);
     buf_out[2*block+1] = imag(a[block]);
   }

   const_cast<INDEX*>(&block_index)->buf_out_update_device();

   BridgeACC::block_scal_eo(vp, buf_out, Nin, Nex,
                            Nsize, block_size, jeo);

 }

#pragma omp barrier

}

//============================================================END=====
#endif
