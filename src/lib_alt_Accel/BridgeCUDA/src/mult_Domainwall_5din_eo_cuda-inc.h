/*!
      @file    mult_Doainwall_5din_eo_cuda-inc.h
      @brief
      @author  Wei Lun Chen (wlchen)
               $LastChangedBy: wlchen $
      @date    $LastChangedDate: 2025-04-24 13:51:53 #$
      @version $LastChangedRevision: 2615 $
*/

#ifndef MULT_DOMAINWALL_5DIN_EO_ACC_INCLUDED
#define MULT_DOMAINWALL_5DIN_EO_ACC_INCLUDED

#include "../inline/constant_memory_inline.h"

//====================================================================
__global__ 
void mult_domainwall_5din_ee_5dir_dirac_kernel(
  real_t * __restrict__ vp, real_t * __restrict__ wp, 
  real_t mq, real_t M0, int Ns, real_t alpha,
  int Nst_pad)
{
  const int Nin5     = NVCD * Ns;
  const int ist      = blockIdx.x * blockDim.x + threadIdx.x;
  const int GridSize = blockDim.x * gridDim.x;
  const real_t* b_con  = const_b;
  const real_t* c_con  = const_c;

  for( int idx = ist; idx < Nst_pad * NVC; idx += GridSize){

    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int ivc     = idx2_wp % NVC;
    int idx_out = idx2_wp / NVC;
    int site    = idx_in  + NWP * idx_out;

    for (int is = 0; is < Ns; ++is) {

      real_t FF1 = b_con[is] * (4.0 - M0) + 1.0;
      real_t FF2 = c_con[is] * (4.0 - M0) - 1.0;  
      real_t vt1, vt2, vt3, vt4; 
      real_t wt1, wt2, wt3, wt4;

      int is_up = (is + 1) % Ns;
      real_t Fup   = 0.5 * FF2 * alpha;
      if (is == Ns-1) Fup = -0.5 * mq * FF2;

      wt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD * is_up), site)];
      wt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD * is_up), site)];
      wt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD * is_up), site)];
      wt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD * is_up), site)];

      vt1 = Fup * (wt1 - wt3);
      vt2 = Fup * (wt2 - wt4);
      vt3 = Fup * (wt3 - wt1);
      vt4 = Fup * (wt4 - wt2);

      int is_dn = (is-1 + Ns) % Ns;
      real_t Fdn   = 0.5 * FF2 * alpha;
      if (is == 0) Fdn = -0.5 * mq * FF2;

      wt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD * is_dn), site)];
      wt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD * is_dn), site)];
      wt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD * is_dn), site)];
      wt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD * is_dn), site)];

      vt1 += Fdn * (wt1 + wt3);
      vt2 += Fdn * (wt2 + wt4);
      vt3 += Fdn * (wt3 + wt1);
      vt4 += Fdn * (wt4 + wt2);

      if(is == 0){
        real_t fac1 = FF1 * 0.5 * ( 1.0 + alpha);
        real_t fac2 = FF1 * 0.5 * (-1.0 + alpha);
        real_t xt1, xt2, xt3, xt4;
        xt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)];
        xt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)];
        xt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)];
        xt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)];
        wt1 = fac1 * xt1 + fac2 * xt3;
        wt2 = fac1 * xt2 + fac2 * xt4;
        wt3 = fac1 * xt3 + fac2 * xt1;
        wt4 = fac1 * xt4 + fac2 * xt2;
      }else if(is == Ns-1){
        real_t fac1 = FF1 * 0.5 * (1.0 + alpha);
        real_t fac2 = FF1 * 0.5 * (1.0 - alpha);
        real_t xt1, xt2, xt3, xt4;
        xt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)];
        xt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)];
        xt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)];
        xt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)];
        wt1 = fac1 * xt1 + fac2 * xt3;
        wt2 = fac1 * xt2 + fac2 * xt4;
        wt3 = fac1 * xt3 + fac2 * xt1;
        wt4 = fac1 * xt4 + fac2 * xt2;
      }else{
        real_t f1 = FF1 * alpha;
        wt1 = f1 * wp[IDX2(Nin5, (ID1 + ivc + NVCD*is), site)];
        wt2 = f1 * wp[IDX2(Nin5, (ID2 + ivc + NVCD*is), site)];
        wt3 = f1 * wp[IDX2(Nin5, (ID3 + ivc + NVCD*is), site)];
        wt4 = f1 * wp[IDX2(Nin5, (ID4 + ivc + NVCD*is), site)];
      }

      vp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)] = wt1 + vt1;
      vp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)] = wt2 + vt2;
      vp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)] = wt3 + vt3;
      vp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)] = wt4 + vt4;
    }
  }
}

void mult_domainwall_5din_ee_5dir_dirac(
  real_t *vp, real_t *wp, real_t mq, real_t M0, int Ns, 
  real_t *b, real_t *c, real_t alpha, int *Nsize)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;
  int Nst_pad = ceil_nwp(Nst);

  real_t* vp_dev = (real_t*)dev_ptr(vp);
  real_t* wp_dev = (real_t*)dev_ptr(wp);

  int blockSize = VECTOR_LENGTH;
  int nSM = getNumSMs();
  dim3 threadsPerBlock(2 * blockSize);
  dim3 numBlocks((Nst_pad * NVC + threadsPerBlock.x - 1) / threadsPerBlock.x);
  int shared_mem_size = 0;
  mult_domainwall_5din_ee_5dir_dirac_kernel<<<numBlocks, threadsPerBlock, shared_mem_size>>>(
      vp_dev, wp_dev, mq, M0, Ns, alpha, Nst_pad);

  CHECK(cudaDeviceSynchronize());
}

//====================================================================
__global__ 
void mult_domainwall_5din_eo_5dir_dirac_kernel(
  real_t * __restrict__ yp, real_t * __restrict__ wp, 
  real_t mq, int Ns, real_t alpha, 
  int Nx, int Ny, int Nz, int Nt, int Nst_pad)
{
  const int Nin5     = NVCD * Ns;
  const int Nst      = Nx *Ny *Nz *Nt;
  const int ist      = blockIdx.x * blockDim.x + threadIdx.x;
  const int GridSize = blockDim.x * gridDim.x;

  for( int idx = ist; idx < Nst_pad * NVC; idx += GridSize){
    
    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int ivc     = idx2_wp % NVC;
    int idx_out = idx2_wp / NVC;
    int site    = idx_in + NWP * idx_out;

    for (int is = 0; is < Ns; ++is) {

      real_t vt1, vt2, vt3, vt4;
      real_t wt1, wt2, wt3, wt4;    

      int is_up           = (is+1) % Ns;
      real_t Fup          = 0.5 * const_c[is] * alpha;
      if (is == Ns-1) Fup = -0.5 * mq * const_c[is];

      wt1  = wp[IDX2(Nin5, (ID1 + ivc + NVCD* is_up), site)];
      wt2  = wp[IDX2(Nin5, (ID2 + ivc + NVCD* is_up), site)];
      wt3  = wp[IDX2(Nin5, (ID3 + ivc + NVCD* is_up), site)];
      wt4  = wp[IDX2(Nin5, (ID4 + ivc + NVCD* is_up), site)];

      vt1  = Fup * (wt1 - wt3);
      vt2  = Fup * (wt2 - wt4);
      vt3  = Fup * (wt3 - wt1);
      vt4  = Fup * (wt4 - wt2);

      int is_dn = (is-1 + Ns) % Ns;
      real_t Fdn   = 0.5 * const_c[is] * alpha;
      if (is == 0) Fdn = -0.5 * mq * const_c[is];

      wt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD* is_dn), site)];
      wt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD* is_dn), site)];
      wt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD* is_dn), site)];
      wt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD* is_dn), site)];

      vt1 += Fdn * (wt1 + wt3);
      vt2 += Fdn * (wt2 + wt4);
      vt3 += Fdn * (wt3 + wt1);
      vt4 += Fdn * (wt4 + wt2);

      if(is == 0){
        real_t b1 = const_b[is] * 0.5 * ( 1.0 + alpha);
        real_t b2 = const_b[is] * 0.5 * (-1.0 + alpha);
        real_t yt1, yt2, yt3, yt4;
        yt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD*is), site)];
        yt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD*is), site)];
        yt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD*is), site)];
        yt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD*is), site)];
        wt1 = b1 * yt1 + b2 * yt3;
        wt2 = b1 * yt2 + b2 * yt4;
        wt3 = b1 * yt3 + b2 * yt1;
        wt4 = b1 * yt4 + b2 * yt2;
      }else if(is == Ns-1){
        real_t b1 = const_b[is] * 0.5 * (1.0 + alpha);
        real_t b2 = const_b[is] * 0.5 * (1.0 - alpha);
        real_t yt1, yt2, yt3, yt4;
        yt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD*is), site)];
        yt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD*is), site)];
        yt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD*is), site)];
        yt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD*is), site)];
        wt1 = b1 * yt1 + b2 * yt3;
        wt2 = b1 * yt2 + b2 * yt4;
        wt3 = b1 * yt3 + b2 * yt1;
        wt4 = b1 * yt4 + b2 * yt2;
      }else{
        real_t bb = const_b[is] * alpha;
        wt1 = bb * wp[IDX2(Nin5, (ID1 + ivc + NVCD*is), site)];
        wt2 = bb * wp[IDX2(Nin5, (ID2 + ivc + NVCD*is), site)];
        wt3 = bb * wp[IDX2(Nin5, (ID3 + ivc + NVCD*is), site)];
        wt4 = bb * wp[IDX2(Nin5, (ID4 + ivc + NVCD*is), site)];
      }

      yp[IDX2(Nin5, (ID1 + ivc + NVCD*is), site)] = -0.5 * (wt1 + vt1);
      yp[IDX2(Nin5, (ID2 + ivc + NVCD*is), site)] = -0.5 * (wt2 + vt2);
      yp[IDX2(Nin5, (ID3 + ivc + NVCD*is), site)] = -0.5 * (wt3 + vt3);
      yp[IDX2(Nin5, (ID4 + ivc + NVCD*is), site)] = -0.5 * (wt4 + vt4);

    }
  }
}

void mult_domainwall_5din_eo_5dir_dirac(
  real_t *yp, real_t *wp, real_t mq, real_t M0, int Ns, 
  real_t *b, real_t *c, real_t alpha, int *Nsize)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;
  int Nst_pad = ceil_nwp(Nst);

  real_t* yp_dev = (real_t*)dev_ptr(yp);
  real_t* wp_dev = (real_t*)dev_ptr(wp);

  int GridSize = (Nst_pad * NVC + VECTOR_LENGTH -1) / VECTOR_LENGTH;
  mult_domainwall_5din_eo_5dir_dirac_kernel<<<GridSize, VECTOR_LENGTH>>>(
      yp_dev, wp_dev, mq, Ns, 
      alpha, 
      Nx, Ny, Nz, Nt, Nst_pad);

  CHECK(cudaDeviceSynchronize());
}

//====================================================================
__global__ 
void mult_domainwall_5din_ee_5dirdag_dirac_kernel(
  real_t * __restrict__ vp, real_t * __restrict__ wp, 
  real_t mq, real_t M0, 
  int Ns, real_t *b, real_t *c, real_t alpha, 
  int Nx, int Ny, int Nz, int Nt, int Nst_pad)
{
  const int Nin5     = NVCD * Ns;
  const int ist      = blockIdx.x * blockDim.x + threadIdx.x;
  const int GridSize = blockDim.x * gridDim.x;

  const real_t* b_con  = const_b;
  const real_t* c_con  = const_c;

  for( int idx = ist; idx < Nst_pad * NVC; idx += GridSize){

    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int ivc     = idx2_wp % NVC;
    int idx_out = idx2_wp / NVC;
    int site    = idx_in  + NWP * idx_out;

    real_t vt1,vt2,vt3,vt4;
    real_t xt1,xt2,xt3,xt4;

    for (int is = 0; is < Ns; ++is) {

      real_t B1 = b_con[is] * (4.0 - M0) + 1.0;
      real_t a1 = -0.5 * b_con[is];

      if(is == 0){
	      real_t fac1 = B1 * 0.5 * ( 1.0 + alpha);
        real_t fac2 = B1 * 0.5 * (-1.0 + alpha);
        real_t wt1, wt2, wt3, wt4;
        wt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)];
        wt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)];
        wt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)];
        wt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)];
        vt1 = fac1 * wt1 + fac2 * wt3;
        vt2 = fac1 * wt2 + fac2 * wt4;
        vt3 = fac1 * wt3 + fac2 * wt1;
        vt4 = fac1 * wt4 + fac2 * wt2;
      }else if(is == Ns-1){
	      real_t fac1 = B1 * 0.5 * (1.0 + alpha);
        real_t fac2 = B1 * 0.5 * (1.0 - alpha);
        real_t wt1, wt2, wt3, wt4;
        wt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)];
        wt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)];
        wt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)];
        wt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)];
        vt1 = fac1 * wt1 + fac2 * wt3;
        vt2 = fac1 * wt2 + fac2 * wt4;
        vt3 = fac1 * wt3 + fac2 * wt1;
        vt4 = fac1 * wt4 + fac2 * wt2;
      }else{
        real_t bb = B1 * alpha;
        vt1 = bb * wp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)];
        vt2 = bb * wp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)];
        vt3 = bb * wp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)];
        vt4 = bb * wp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)];
      }

      int is_up = (is+1) % Ns;
      xt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD * is_up), site)];
      xt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD * is_up), site)];
      xt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD * is_up), site)];
      xt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD * is_up), site)];

      real_t Fup = 0.5 * (c_con[is_up] * (4.0 - M0) - 1.0);
      if (is == Ns-1){
	      Fup *= -mq;
      }else{
        Fup *= alpha;
      }

      vt1 += Fup * (xt1 + xt3);
      vt2 += Fup * (xt2 + xt4);
      vt3 += Fup * (xt3 + xt1);
      vt4 += Fup * (xt4 + xt2);

      int is_dn = (is-1 + Ns) % Ns;
      xt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD * is_dn), site)];
      xt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD * is_dn), site)];
      xt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD * is_dn), site)];
      xt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD * is_dn), site)];

      real_t Fdn = 0.5 * (c_con[is_dn] * (4.0 - M0) - 1.0);
      if (is == 0){
        Fdn *= -mq;
      }else{
        Fdn *= alpha;
      }

      vt1 += Fdn * (xt1 - xt3);
      vt2 += Fdn * (xt2 - xt4);
      vt3 += Fdn * (xt3 - xt1);
      vt4 += Fdn * (xt4 - xt2);

      vp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)] = vt1;
      vp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)] = vt2;
      vp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)] = vt3;
      vp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)] = vt4;

    }
  }
}

void mult_domainwall_5din_ee_5dirdag_dirac(
    real_t *vp, real_t *wp, real_t mq, real_t M0, int Ns, real_t *b, real_t *c, 
    real_t alpha, int *Nsize) 
{
    int Nx = Nsize[0];
    int Ny = Nsize[1];
    int Nz = Nsize[2];
    int Nt = Nsize[3];
    int Nst = Nx * Ny * Nz * Nt;
    int Nst_pad = ceil_nwp(Nst);

    real_t* vp_dev = (real_t*)dev_ptr(vp);
    real_t* wp_dev = (real_t*)dev_ptr(wp);

    real_t* b_dev;
    real_t* c_dev; 

    int blockSize = VECTOR_LENGTH; 
    dim3 threadsPerBlock(2* blockSize);
    dim3 numBlocks((Nst_pad * NVC + threadsPerBlock.x - 1) / threadsPerBlock.x);

    mult_domainwall_5din_ee_5dirdag_dirac_kernel<<<numBlocks, threadsPerBlock>>>(
        vp_dev, wp_dev, mq, M0, Ns, b_dev, c_dev, alpha, Nx, Ny, Nz, Nt, Nst_pad);

    CHECK(cudaDeviceSynchronize());
}


//====================================================================
__global__ 
void mult_domainwall_5din_eo_5dirdag_dirac_kernel(
    real_t * __restrict__ vp, real_t * __restrict__ yp, 
    real_t mq, real_t M0, int Ns, 
    real_t *b, real_t *c, real_t alpha, 
    int Nx, int Ny, int Nz, int Nt, int Nst_pad)
{
  const int Nin5     = NVCD * Ns;
  const int ist      = blockIdx.x * blockDim.x + threadIdx.x;
  const int GridSize = blockDim.x * gridDim.x;

  const real_t* b_con  = const_b;
  const real_t* c_con  = const_c;

  for( int idx = ist; idx < Nst_pad * NVC; idx += GridSize){

    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int ivc     = idx2_wp % NVC;
    int idx_out = idx2_wp / NVC;
    int site    = idx_in + NWP * idx_out;

    real_t vt1, vt2, vt3, vt4;
    real_t yt1, yt2, yt3, yt4;

    for (int is = 0; is < Ns; ++is) {

      yt1  = yp[IDX2(Nin5, (ID1 + ivc + NVCD*is), site)];
      yt2  = yp[IDX2(Nin5, (ID2 + ivc + NVCD*is), site)];
      yt3  = yp[IDX2(Nin5, (ID3 + ivc + NVCD*is), site)];
      yt4  = yp[IDX2(Nin5, (ID4 + ivc + NVCD*is), site)];
      if(is == 0){
        real_t b1 = -0.5 * b_con[is] * 0.5 * ( 1.0 + alpha);
        real_t b2 = -0.5 * b_con[is] * 0.5 * (-1.0 + alpha);
        vt1  = b1 * yt3 + b2 * yt1;
        vt2  = b1 * yt4 + b2 * yt2;
        vt3  = b1 * yt1 + b2 * yt3;
        vt4  = b1 * yt2 + b2 * yt4;
      }else if(is == Ns-1){
        real_t b1 = -0.5 * b_con[is] * 0.5 * (1.0 + alpha);
        real_t b2 = -0.5 * b_con[is] * 0.5 * (1.0 - alpha);
        vt1  = b1 * yt3 + b2 * yt1;
        vt2  = b1 * yt4 + b2 * yt2;
        vt3  = b1 * yt1 + b2 * yt3;
        vt4  = b1 * yt2 + b2 * yt4;
      }else{
        real_t bb = -0.5 * b_con[is] * alpha;
        vt1  = bb * yt3;
        vt2  = bb * yt4;
        vt3  = bb * yt1;
        vt4  = bb * yt2;
      }

      int is_up = (is+1) % Ns;
      yt1  = yp[IDX2(Nin5, (ID1 + ivc + NVCD*is_up), site)];
      yt2  = yp[IDX2(Nin5, (ID2 + ivc + NVCD*is_up), site)];
      yt3  = yp[IDX2(Nin5, (ID3 + ivc + NVCD*is_up), site)];
      yt4  = yp[IDX2(Nin5, (ID4 + ivc + NVCD*is_up), site)];

      real_t Fup = -0.5 * c_con[is_up] * 0.5 * alpha;
      if (is == Ns-1) Fup = -0.5 * c_con[is_up] * (-0.5) * mq;

      vt1  += Fup * (yt3 + yt1);
      vt2  += Fup * (yt4 + yt2);
      vt3  += Fup * (yt1 + yt3);
      vt4  += Fup * (yt2 + yt4);

      int is_dn = (is-1 + Ns) % Ns;
      yt1  = yp[IDX2(Nin5, (ID1 + ivc + NVCD*is_dn), site)];
      yt2  = yp[IDX2(Nin5, (ID2 + ivc + NVCD*is_dn), site)];
      yt3  = yp[IDX2(Nin5, (ID3 + ivc + NVCD*is_dn), site)];
      yt4  = yp[IDX2(Nin5, (ID4 + ivc + NVCD*is_dn), site)];


      real_t Fdn = -0.5 * c_con[is_dn] * (0.5) * alpha;
      if (is == 0) Fdn = -0.5 * c_con[is_dn] * (-0.5) * mq;

      vt1  += Fdn * (yt3 - yt1);
      vt2  += Fdn * (yt4 - yt2);
      vt3  += Fdn * (yt1 - yt3);
      vt4  += Fdn * (yt2 - yt4);

      vp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)] = vt1;
      vp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)] = vt2;
      vp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)] = vt3;
      vp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)] = vt4;
    } 
  }
}

void mult_domainwall_5din_eo_5dirdag_dirac(
  real_t *vp, real_t *yp, real_t mq, real_t M0, int Ns, 
  real_t *b, real_t *c, real_t alpha, int *Nsize)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;
  int Nst_pad = ceil_nwp(Nst);

  real_t* vp_dev = (real_t*)dev_ptr(vp);
  real_t* yp_dev = (real_t*)dev_ptr(yp);

  real_t* b_dev;
  real_t* c_dev; 

  int blockSize = 2* VECTOR_LENGTH; 
  dim3 threadsPerBlock(blockSize);
  dim3 numBlocks((Nst_pad * NVC + threadsPerBlock.x - 1) / threadsPerBlock.x);
  mult_domainwall_5din_eo_5dirdag_dirac_kernel<<<numBlocks, threadsPerBlock>>>(
      vp_dev, yp_dev, mq, M0, Ns, b_dev, c_dev, alpha, Nx, Ny, Nz, Nt, Nst_pad);

  CHECK(cudaDeviceSynchronize());
}

//====================================================================
__global__
void mult_domainwall_5din_eo_hopb_dirac_4D_kernel(
  real_t * __restrict__ vp, real_t * __restrict__ up, real_t * __restrict__ wp, 
  int Ns, 
  int bc_x, int bc_y, int bc_z, int bc_t,
  int Nx, int Ny, int Nz, int Nt, 
  int ieo, int jeo, 
  int do_comm_x, int do_comm_y, int do_comm_z, int do_comm_t,
  int Nst_pad, int jgm5)
{
  const int Nin5 = NVCD * Ns;
  const int Nst  = Nx * Ny * Nz * Nt;
  const int site = blockIdx.x * blockDim.x + threadIdx.x;

  real_t ut[NDF];

  if (site < Nst) {
    int Nxy  = Nx  * Ny;
    int Nxyz = Nxy * Nz;
    int ix   = site % Nx;
    int iyzt = site / Nx;
    int ixy  = site % Nxy;
    int iy   = iyzt % Ny;
    int izt  = site / Nxy;
    int iz   = izt % Nz;
    int it   = izt / Nz;
    int ixyz = site % Nxyz;
    int keo  = (jeo + iy + iz + it) % 2;
    int idir;

    for(int is = 0; is < Ns; ++is) {
      real_t vL[NVCD];
      for(int ivcd = 0; ivcd < NVCD; ++ivcd) {
          vL[ivcd] = 0.0;
      }
      idir = 0;
      // mult_xp
      {
        int ix2 = (ix + keo) % Nx;
        int nei = ix2 + Nx * iyzt;
        real_t bc2 = 1.0;
        if(ix == Nx-1 && keo == 1) bc2 = bc_x;

          load_u(ut, up, site + Nst * (ieo + 2*idir));

          real_t wt[NVCD];
          if(jgm5 == 0){
            for(int ivcd = 0; ivcd < NVCD; ++ivcd){
              wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
            }
          }else{
            for(int ivc = 0; ivc < NVC; ++ivc){
              wt[ivc+ID1] = bc2 * wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), nei)];
              wt[ivc+ID2] = bc2 * wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), nei)];
              wt[ivc+ID3] = bc2 * wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), nei)];
              wt[ivc+ID4] = bc2 * wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), nei)];
            }
          }
        mult_wilson_xpb(vL, ut, wt);
      }
      // mult_xm
      {
        int ix2 = (ix - 1 + keo + Nx) % Nx;
        int nei  = ix2 + Nx * iyzt;
        real_t bc2 = 1.0;
        if(ix == 0 && keo == 0) bc2 = bc_x;

          load_u(ut, up, nei + Nst * (1-ieo + 2*idir));
          real_t wt[NVCD];
          if(jgm5 == 0){
            for(int ivcd = 0; ivcd < NVCD; ++ivcd){
              wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
            }
          }else{
            for(int ivc = 0; ivc < NVC; ++ivc){
              wt[ivc+ID1] = bc2 * wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), nei)];
              wt[ivc+ID2] = bc2 * wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), nei)];
              wt[ivc+ID3] = bc2 * wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), nei)];
              wt[ivc+ID4] = bc2 * wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), nei)];
            }
          }
        mult_wilson_xmb(vL, ut, wt);
      }

      idir = 1;

      if ((iy < Ny-1) || (do_comm_y == 0)) {
        int iy2 = (iy + 1) % Ny;
        int nei = ix + Nx * (iy2 + Ny * izt);
        real_t bc2 = 1.0;
        if(iy == Ny-1) bc2 = bc_y;

        load_u(ut, up, site + Nst * (ieo + 2*idir));
        real_t wt[NVCD];
        if(jgm5 == 0){
          for(int ivcd = 0; ivcd < NVCD; ++ivcd){
            wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
          }
        }else{
          for(int ivc = 0; ivc < NVC; ++ivc){
            wt[ivc+ID1] = bc2 * wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), nei)];
            wt[ivc+ID2] = bc2 * wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), nei)];
            wt[ivc+ID3] = bc2 * wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), nei)];
            wt[ivc+ID4] = bc2 * wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), nei)];
          }
        }
        mult_wilson_ypb(vL, ut, wt);
      }
      if ((iy > 0) || (do_comm_y == 0)) {
        int iy2 = (iy - 1 + Ny) % Ny;
        int nei = ix + Nx * (iy2 + Ny * izt);
        real_t bc2 = 1.0;
        if(iy == 0) bc2 = bc_y;

          load_u(ut, up, nei + Nst * (1-ieo + 2*idir));
          real_t wt[NVCD];
        if(jgm5 == 0){
          for(int ivcd = 0; ivcd < NVCD; ++ivcd){
            wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
          }
        }else{
          for(int ivc = 0; ivc < NVC; ++ivc){
            wt[ivc+ID1] = bc2 * wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), nei)];
            wt[ivc+ID2] = bc2 * wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), nei)];
            wt[ivc+ID3] = bc2 * wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), nei)];
            wt[ivc+ID4] = bc2 * wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), nei)];
          }
        }
        mult_wilson_ymb(vL, ut, wt);
      }

      idir = 2;

      if ((iz < Nz-1) || (do_comm_z == 0)) {
        int iz2 = (iz + 1) % Nz;
        int nei = ixy + Nxy * (iz2 + Nz * it);
        real_t bc2 = 1.0;
        if(iz == Nz-1) bc2 = bc_z;

          load_u(ut, up, site + Nst * (ieo + 2*idir));
          real_t wt[NVCD];
          if(jgm5 == 0){
            for(int ivcd = 0; ivcd < NVCD; ++ivcd){
              wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
            }
          }else{
            for(int ivc = 0; ivc < NVC; ++ivc){
              wt[ivc+ID1] = bc2 * wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), nei)];
              wt[ivc+ID2] = bc2 * wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), nei)];
              wt[ivc+ID3] = bc2 * wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), nei)];
              wt[ivc+ID4] = bc2 * wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), nei)];
            }
          }
        mult_wilson_zpb(vL, ut, wt);
      }

      if ((iz > 0) || (do_comm_z == 0)) {
        int iz2 = (iz - 1 + Nz) % Nz;
        int nei = ixy + Nxy * (iz2 + Nz * it);
        real_t bc2 = 1.0;
        if(iz == 0) bc2 = bc_z;

          load_u(ut, up, nei + Nst * (1-ieo + 2*idir));
          real_t wt[NVCD];
          if(jgm5 == 0){
            for(int ivcd = 0; ivcd < NVCD; ++ivcd){
              wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
            }
          }else{
            for(int ivc = 0; ivc < NVC; ++ivc){
              wt[ivc+ID1] = bc2 * wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), nei)];
              wt[ivc+ID2] = bc2 * wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), nei)];
              wt[ivc+ID3] = bc2 * wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), nei)];
              wt[ivc+ID4] = bc2 * wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), nei)];
            }
          }
        mult_wilson_zmb(vL, ut, wt);
      }

      idir = 3;

      if ((it < Nt-1) || (do_comm_t == 0)) {
        int it2 = (it + 1) % Nt;
        int nei = ixyz + Nxyz * it2;
        real_t bc2 = 1.0;
        if(it == Nt-1) bc2 = bc_t;

          load_u(ut, up, site + Nst * (ieo + 2*idir));
          real_t wt[NVCD];
          if(jgm5 == 0){
            for(int ivcd = 0; ivcd < NVCD; ++ivcd){
              wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
            }
          }else{
            for(int ivc = 0; ivc < NVC; ++ivc){
              wt[ivc+ID1] = bc2 * wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), nei)];
              wt[ivc+ID2] = bc2 * wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), nei)];
              wt[ivc+ID3] = bc2 * wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), nei)];
              wt[ivc+ID4] = bc2 * wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), nei)];
            }
          }
        mult_wilson_tpb_dirac(vL, ut, wt);
      }

      if ((it > 0) || (do_comm_t == 0)) {
        int it2 = (it - 1 + Nt) % Nt;
        int nei = ixyz + Nxyz * it2;
        real_t bc2 = 1.0;
        if(it == 0) bc2 = bc_t;

        load_u(ut, up, nei + Nst * (1-ieo + 2*idir));
        real_t wt[NVCD];
        if(jgm5 == 0){
          for(int ivcd = 0; ivcd < NVCD; ++ivcd){
            wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
          }
        }else{
          for(int ivc = 0; ivc < NVC; ++ivc){
            wt[ivc+ID1] = bc2 * wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), nei)];
            wt[ivc+ID2] = bc2 * wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), nei)];
            wt[ivc+ID3] = bc2 * wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), nei)];
            wt[ivc+ID4] = bc2 * wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), nei)];
          }
        }
        mult_wilson_tmb_dirac(vL, ut, wt);
      }

      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vL[ivcd];
      }
    }
  }
}

void mult_domainwall_5din_eo_hopb_dirac_4d(
    real_t *vp, real_t *up, real_t *wp, int Ns, int *bc, 
    int *Nsize, int *do_comm, int ieo, int jeo, int jgm5)
{
    int Nx = Nsize[0];
    int Ny = Nsize[1];
    int Nz = Nsize[2];
    int Nt = Nsize[3];
    int Nst = Nx * Ny * Nz * Nt;
    int Nst_pad = ceil_nwp(Nst);

    real_t* vp_dev = (real_t*)dev_ptr(vp);
    real_t* up_dev = (real_t*)dev_ptr(up);
    real_t* wp_dev = (real_t*)dev_ptr(wp);

    int bc_x = bc[0];
    int bc_y = bc[1];
    int bc_z = bc[2];
    int bc_t = bc[3];
    int do_comm_x = do_comm[0];
    int do_comm_y = do_comm[1];
    int do_comm_z = do_comm[2];
    int do_comm_t = do_comm[3];

    int blockSize  = VECTOR_LENGTH; 
    int gridSize   = (Nst_pad * Ns + blockSize - 1) / blockSize;
    int sharedsize = NDF * blockSize * sizeof(real_t);

    mult_domainwall_5din_eo_hopb_dirac_4D_kernel<<<gridSize, blockSize, sharedsize>>>(
        vp_dev, up_dev, wp_dev, Ns, 
        bc_x, bc_y, bc_z, bc_t,
        Nx, Ny, Nz, Nt, 
        ieo, jeo, 
        do_comm_x, do_comm_y, do_comm_z, do_comm_t, Nst_pad, jgm5);

    cudaDeviceSynchronize();
}

//====================================================================
// v2 of hopb_dirac_5D: block-level gauge link sharing via __shared__.
//   Thread block: (NWP, NS_PER_BLOCK).  NS_PER_BLOCK warps share the same
//   32-site block but process different is.  is_local==0 loads 12 gauge
//   doubles once per direction into u_sh, all warps read from u_sh.
//   Requires SU3_3RD_ROW_RECONST and Ns % NS_PER_BLOCK == 0.
//   Toggle off by removing -DUSE_5D_HOPB_V2 (falls back to v1 launcher).
//
// Tunables to benchmark:
//   NS_PER_BLOCK   2 / 4 / 8 / 16     (must divide Ns; smaller → more occupancy,
//                                      bigger → more u-sharing)
//   LB_BLOCKS_PER_SM   1 / 2 / 3      (2nd arg to __launch_bounds__; higher
//                                      forces lower REG/thread)
#ifndef NS_PER_BLOCK
#define NS_PER_BLOCK 8
#endif
#ifndef LB_BLOCKS_PER_SM
#define LB_BLOCKS_PER_SM 2
#endif

#define U_UP_LOAD_AND_SYNC(isg) do { \
  __syncthreads(); \
  if (is_local == 0) { \
    u_sh[ 0*NWP + lane] = u_up[IDX2_G_R(0,0,isg)]; \
    u_sh[ 1*NWP + lane] = u_up[IDX2_G_I(0,0,isg)]; \
    u_sh[ 2*NWP + lane] = u_up[IDX2_G_R(1,0,isg)]; \
    u_sh[ 3*NWP + lane] = u_up[IDX2_G_I(1,0,isg)]; \
    u_sh[ 4*NWP + lane] = u_up[IDX2_G_R(2,0,isg)]; \
    u_sh[ 5*NWP + lane] = u_up[IDX2_G_I(2,0,isg)]; \
    u_sh[ 6*NWP + lane] = u_up[IDX2_G_R(0,1,isg)]; \
    u_sh[ 7*NWP + lane] = u_up[IDX2_G_I(0,1,isg)]; \
    u_sh[ 8*NWP + lane] = u_up[IDX2_G_R(1,1,isg)]; \
    u_sh[ 9*NWP + lane] = u_up[IDX2_G_I(1,1,isg)]; \
    u_sh[10*NWP + lane] = u_up[IDX2_G_R(2,1,isg)]; \
    u_sh[11*NWP + lane] = u_up[IDX2_G_I(2,1,isg)]; \
  } \
  __syncthreads(); \
} while(0)

#define U_DN_LOAD_AND_SYNC(isg) do { \
  __syncthreads(); \
  if (is_local == 0) { \
    u_sh[ 0*NWP + lane] = u_dn[IDX2_G_R(0,0,isg)]; \
    u_sh[ 1*NWP + lane] = u_dn[IDX2_G_I(0,0,isg)]; \
    u_sh[ 2*NWP + lane] = u_dn[IDX2_G_R(0,1,isg)]; \
    u_sh[ 3*NWP + lane] = u_dn[IDX2_G_I(0,1,isg)]; \
    u_sh[ 4*NWP + lane] = u_dn[IDX2_G_R(0,2,isg)]; \
    u_sh[ 5*NWP + lane] = u_dn[IDX2_G_I(0,2,isg)]; \
    u_sh[ 6*NWP + lane] = u_dn[IDX2_G_R(1,0,isg)]; \
    u_sh[ 7*NWP + lane] = u_dn[IDX2_G_I(1,0,isg)]; \
    u_sh[ 8*NWP + lane] = u_dn[IDX2_G_R(1,1,isg)]; \
    u_sh[ 9*NWP + lane] = u_dn[IDX2_G_I(1,1,isg)]; \
    u_sh[10*NWP + lane] = u_dn[IDX2_G_R(1,2,isg)]; \
    u_sh[11*NWP + lane] = u_dn[IDX2_G_I(1,2,isg)]; \
  } \
  __syncthreads(); \
} while(0)

__global__ __launch_bounds__(NWP * NS_PER_BLOCK, LB_BLOCKS_PER_SM)
void mult_domainwall_5din_eo_hopb_dirac_5D_v2_kernel(
  real_t * __restrict__ vp, const real_t * __restrict__ up, const real_t * __restrict__ wp,
  int Ns,
  int bc_x, int bc_y, int bc_z, int bc_t,
  int Nx, int Ny, int Nz, int Nt,
  int ieo, int jeo,
  int do_comm_x, int do_comm_y, int do_comm_z, int do_comm_t,
  int Nst_pad, int jgm5){

  const int Nxy   = Nx * Ny;
  const int Nst   = Nx * Ny * Nz * Nt;
  const int Nxyz  = Nx * Ny * Nz;

  const int lane     = threadIdx.x;
  const int is_local = threadIdx.y;
  const int idx_out  = blockIdx.x;
  const int is_base  = blockIdx.y * NS_PER_BLOCK;
  const int is       = is_base + is_local;
  const int site     = lane + NWP * idx_out;

  const real_t * u_up = up;
  const real_t * u_dn = up;

  __shared__ real_t u_sh[12 * NWP];

  int ix   = site % Nx;
  int iyzt = site / Nx;
  int iy   = iyzt % Ny;
  int izt  = site / Nxy;
  int iz   = izt  % Nz;
  int it   = izt  / Nz;
  int ixyz = site % Nxyz;
  int keo  = (jeo + iy + iz + it) % 2;
  int idir;

  real_t u_0, u_1, u_2, u_3, u_4, u_5;
  real_t u_6, u_7, u_8, u_9, u10, u11;
  real_t u12, u13, u14, u15, u16, u17;
  real_t vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5;
  real_t vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5;
  real_t wt1r, wt1i, wt2r, wt2i;

  real_t v2_01, v2_11, v2_21, v2_31, v2_41, v2_51;
  real_t v2_02, v2_12, v2_22, v2_32, v2_42, v2_52;
  real_t v2_03, v2_13, v2_23, v2_33, v2_43, v2_53;
  real_t v2_04, v2_14, v2_24, v2_34, v2_44, v2_54;

  const real_t * v1 = wp;
  real_t * v2 = vp;

  #include "inc/mult_Domainwall_eo_cuda_xyz_shmem-inc.h"
  #include "inc/mult_Domainwall_eo_cuda_t_dirac_shmem-inc.h"

  v2[IDX2_SP_5D_R(0,0,is,Ns,site)] = v2_01;
  v2[IDX2_SP_5D_I(0,0,is,Ns,site)] = v2_11;
  v2[IDX2_SP_5D_R(1,0,is,Ns,site)] = v2_21;
  v2[IDX2_SP_5D_I(1,0,is,Ns,site)] = v2_31;
  v2[IDX2_SP_5D_R(2,0,is,Ns,site)] = v2_41;
  v2[IDX2_SP_5D_I(2,0,is,Ns,site)] = v2_51;

  v2[IDX2_SP_5D_R(0,1,is,Ns,site)] = v2_02;
  v2[IDX2_SP_5D_I(0,1,is,Ns,site)] = v2_12;
  v2[IDX2_SP_5D_R(1,1,is,Ns,site)] = v2_22;
  v2[IDX2_SP_5D_I(1,1,is,Ns,site)] = v2_32;
  v2[IDX2_SP_5D_R(2,1,is,Ns,site)] = v2_42;
  v2[IDX2_SP_5D_I(2,1,is,Ns,site)] = v2_52;

  v2[IDX2_SP_5D_R(0,2,is,Ns,site)] = v2_03;
  v2[IDX2_SP_5D_I(0,2,is,Ns,site)] = v2_13;
  v2[IDX2_SP_5D_R(1,2,is,Ns,site)] = v2_23;
  v2[IDX2_SP_5D_I(1,2,is,Ns,site)] = v2_33;
  v2[IDX2_SP_5D_R(2,2,is,Ns,site)] = v2_43;
  v2[IDX2_SP_5D_I(2,2,is,Ns,site)] = v2_53;

  v2[IDX2_SP_5D_R(0,3,is,Ns,site)] = v2_04;
  v2[IDX2_SP_5D_I(0,3,is,Ns,site)] = v2_14;
  v2[IDX2_SP_5D_R(1,3,is,Ns,site)] = v2_24;
  v2[IDX2_SP_5D_I(1,3,is,Ns,site)] = v2_34;
  v2[IDX2_SP_5D_R(2,3,is,Ns,site)] = v2_44;
  v2[IDX2_SP_5D_I(2,3,is,Ns,site)] = v2_54;
}

//====================================================================
__global__
void mult_domainwall_5din_eo_hopb_dirac_5D_kernel(
  real_t * __restrict__ vp, const real_t * __restrict__ up, const real_t * __restrict__ wp,
  int Ns,
  int bc_x, int bc_y, int bc_z, int bc_t,
  int Nx, int Ny, int Nz, int Nt, 
  int ieo, int jeo, 
  int do_comm_x, int do_comm_y, int do_comm_z, int do_comm_t,
  int Nst_pad, int jgm5){

  const int Nxy      = Nx * Ny;
  const int Nst      = Nx * Ny * Nz * Nt;
  const int Nxyz     = Nx * Ny * Nz;
  const int ist      = blockIdx.x * blockDim.x + threadIdx.x;
  const int GridSize = blockDim.x * gridDim.x;

  const real_t * u_up = up;
  const real_t * u_dn = up;

  for( int idx = ist; idx < Nst_pad * Ns; idx += GridSize){

    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int is      = idx2_wp % Ns;
    int idx_out = idx2_wp / Ns;
    int site    = idx_in + NWP*idx_out;

    int ix   = site % Nx;
    int iyzt = site / Nx;
    int ixy  = site % Nxy;
    int iy   = iyzt % Ny;
    int izt  = site / Nxy;
    int iz   = izt  % Nz;
    int it   = izt  / Nz;
    int ixyz = site % Nxyz;
    int keo  = (jeo + iy + iz + it) % 2;
    int idir;

    real_t u_0, u_1, u_2, u_3, u_4, u_5;
    real_t u_6, u_7, u_8, u_9, u10, u11;
    real_t u12, u13, u14, u15, u16, u17;
    real_t vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5;
    real_t vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5;
    real_t wt1r, wt1i, wt2r, wt2i;

    real_t v2_01, v2_11, v2_21, v2_31, v2_41, v2_51;
    real_t v2_02, v2_12, v2_22, v2_32, v2_42, v2_52;
    real_t v2_03, v2_13, v2_23, v2_33, v2_43, v2_53;
    real_t v2_04, v2_14, v2_24, v2_34, v2_44, v2_54;

    const real_t * v1 = wp;
    real_t * v2 = vp;

    #include "inc/mult_Domainwall_eo_openacc2_xyz-inc.h"
    #include "inc/mult_Domainwall_eo_openacc2_t_dirac-inc.h"

    v2[IDX2_SP_5D_R(0,0,is,Ns,site)] = v2_01;
    v2[IDX2_SP_5D_I(0,0,is,Ns,site)] = v2_11;
    v2[IDX2_SP_5D_R(1,0,is,Ns,site)] = v2_21;
    v2[IDX2_SP_5D_I(1,0,is,Ns,site)] = v2_31;
    v2[IDX2_SP_5D_R(2,0,is,Ns,site)] = v2_41;
    v2[IDX2_SP_5D_I(2,0,is,Ns,site)] = v2_51;

    v2[IDX2_SP_5D_R(0,1,is,Ns,site)] = v2_02;
    v2[IDX2_SP_5D_I(0,1,is,Ns,site)] = v2_12;
    v2[IDX2_SP_5D_R(1,1,is,Ns,site)] = v2_22;
    v2[IDX2_SP_5D_I(1,1,is,Ns,site)] = v2_32;
    v2[IDX2_SP_5D_R(2,1,is,Ns,site)] = v2_42;
    v2[IDX2_SP_5D_I(2,1,is,Ns,site)] = v2_52;

    v2[IDX2_SP_5D_R(0,2,is,Ns,site)] = v2_03;
    v2[IDX2_SP_5D_I(0,2,is,Ns,site)] = v2_13;
    v2[IDX2_SP_5D_R(1,2,is,Ns,site)] = v2_23;
    v2[IDX2_SP_5D_I(1,2,is,Ns,site)] = v2_33;
    v2[IDX2_SP_5D_R(2,2,is,Ns,site)] = v2_43;
    v2[IDX2_SP_5D_I(2,2,is,Ns,site)] = v2_53;

    v2[IDX2_SP_5D_R(0,3,is,Ns,site)] = v2_04;
    v2[IDX2_SP_5D_I(0,3,is,Ns,site)] = v2_14;
    v2[IDX2_SP_5D_R(1,3,is,Ns,site)] = v2_24;
    v2[IDX2_SP_5D_I(1,3,is,Ns,site)] = v2_34;
    v2[IDX2_SP_5D_R(2,3,is,Ns,site)] = v2_44;
    v2[IDX2_SP_5D_I(2,3,is,Ns,site)] = v2_54;
  }
}

void mult_domainwall_5din_eo_hopb_dirac_5d(
    real_t *vp, real_t *up, real_t *wp, int Ns, int *bc,
    int *Nsize, int *do_comm, int ieo, int jeo, int jgm5, bool recon)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;
  int Nst_pad = ceil_nwp(Nst);

  real_t* vp_dev = (real_t*)dev_ptr(vp);
  real_t* up_dev = (real_t*)dev_ptr(up);
  real_t* wp_dev = (real_t*)dev_ptr(wp);

  int bc_x = bc[0];
  int bc_y = bc[1];
  int bc_z = bc[2];
  int bc_t = bc[3];
  int do_comm_x = do_comm[0];
  int do_comm_y = do_comm[1];
  int do_comm_z = do_comm[2];
  int do_comm_t = do_comm[3];

  // su3_reconstruction (YAML) -> v2 kernel: caches 12 gauge reals in shared memory
  // and rebuilds column 2 (SU3_3RD_ROW_RECONST). recon=false uses v1 (full link).
  if (recon && Ns % NS_PER_BLOCK == 0 && (Nst_pad % NWP == 0)) {
    dim3 block(NWP, NS_PER_BLOCK);
    dim3 grid (Nst_pad / NWP, Ns / NS_PER_BLOCK);
    mult_domainwall_5din_eo_hopb_dirac_5D_v2_kernel<<<grid, block>>>(
        vp_dev, up_dev, wp_dev, Ns,
        bc_x, bc_y, bc_z, bc_t,
        Nx, Ny, Nz, Nt,
        ieo, jeo,
        do_comm_x, do_comm_y, do_comm_z, do_comm_t,
        Nst_pad, jgm5);
    cudaDeviceSynchronize();
    return;
  }

  int blockSize = VECTOR_LENGTH;
  int gridSize  = (Nst_pad * Ns + blockSize - 1) / blockSize;

  mult_domainwall_5din_eo_hopb_dirac_5D_kernel<<<gridSize, blockSize>>>(
      vp_dev, up_dev, wp_dev, Ns,
      bc_x, bc_y, bc_z, bc_t,
      Nx, Ny, Nz, Nt,
      ieo, jeo,
      do_comm_x, do_comm_y, do_comm_z, do_comm_t,
      Nst_pad, jgm5);

  cudaDeviceSynchronize();
}
//====================================================================

__global__  
void mult_domainwall_5din_eo_hop1x_dirac_dev(
  real_t * __restrict__ buf_xp, real_t * __restrict__ buf_xm,
  real_t * __restrict__ up, real_t * __restrict__ wp,
  int Ns, int bc, 
  int Nx, int Ny, int Nz, int Nt, 
  int ieo, int jeo, int jgm5){

  const int iyzt   = blockIdx.x * blockDim.x + threadIdx.x;
  const int Nyzt   = Ny   * Nz  * Nt     ;
  const int Nin5   = NVCD * Ns           ;
  const int Nin5bd = NVC  * ND2 * Ns     ;
  const int Nst    = Nx   * Ny  * Nz * Nt;
  int idir         = 0                   ;
  real_t ut[NDF];

  if ( iyzt < Nyzt ){

    int iy  = iyzt % Ny;
    int iz  = (iyzt/Ny) % Nz;
    int it  = iyzt / (Ny * Nz);
    int keo = (jeo + iy + iz + it) % 2;

    if(keo == 1){
      int ix = 0;
      int iyzt2 = iyzt/2;
      int site = ix + Nx * iyzt;
      real_t bc2 = bc;
      for (int is = 0; is < Ns; ++is) {
          real_t wt[NVCD], vt[NVC * ND2];
          if(jgm5 == 0){
            for(int ivcd = 0; ivcd < NVCD; ++ivcd){
              wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
            }
          }else{
            for(int ivc = 0; ivc < NVC; ++ivc){
              wt[ivc+ID1] = wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), site)];
              wt[ivc+ID2] = wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), site)];
              wt[ivc+ID3] = wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), site)];
              wt[ivc+ID4] = wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), site)];
            }
          }
          mult_wilson_xp1(vt, wt);
          for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
              buf_xp[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), iyzt2)] = bc2 * vt[ivcd];
          }
      }
    }

    if(keo == 0){
      int iyzt2 = iyzt/2;
      int ix = Nx-1;
      int site = ix + Nx * iyzt;
      load_u(ut, up, site + Nst * (1- ieo + 2*idir));
      for (int is = 0; is < Ns; ++is) {
        real_t wt[NVCD], vt[NVC * ND2];
        if(jgm5 == 0){
          for(int ivcd = 0; ivcd < NVCD; ++ivcd){
            wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
          }
        }else{
          for(int ivc = 0; ivc < NVC; ++ivc){
            wt[ivc+ID1] = wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), site)];
            wt[ivc+ID2] = wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), site)];
            wt[ivc+ID3] = wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), site)];
            wt[ivc+ID4] = wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), site)];
          }
        }
        mult_wilson_xm1(vt, ut, wt);
        for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
            buf_xm[IDX2(Nin5bd, (ivcd + NVC*ND2*is), iyzt2)] = vt[ivcd];
        }
      }
    }
  }
}

__global__  
void mult_domainwall_5din_eo_hop1y_dirac_dev(
  real_t * __restrict__ buf_yp, real_t * __restrict__ buf_ym,
  real_t * __restrict__ up, real_t * __restrict__ wp,
  int Ns, int bc, int Nx, int Ny, int Nz, int Nt, 
  int ieo, int jeo, int jgm5){

  const int site_thread  = blockIdx.x * blockDim.x + threadIdx.x;
  const int Nxzt         = Nx * Nz * Nt;
  const int ixzt         = site_thread; 
  const int Nst          = Nx * Ny * Nz * Nt;
  const int Nin5         = NVCD * Ns;
  const int Nin5bd       = NVC  * ND2 * Ns;
  int idir = 1;
  real_t ut[NDF];
  if ( ixzt < Nxzt ){

    int iy = 0;
    int ix  = ixzt % Nx;
    int izt = ixzt / Nx;
    int site = ix + Nx * (iy + Ny * izt);
    real_t bc2 = bc;

    for (int is = 0; is < Ns; ++is) {
      real_t wt[NVCD], vt[NVC * ND2];
      if(jgm5 == 0){
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
        }
      }else{
        for(int ivc = 0; ivc < NVC; ++ivc){
          wt[ivc+ID1] = wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), site)];
          wt[ivc+ID2] = wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), site)];
          wt[ivc+ID3] = wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), site)];
          wt[ivc+ID4] = wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), site)];
        }
      }
      mult_wilson_yp1(vt, wt);
      for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
          buf_yp[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixzt)] = bc2 * vt[ivcd];
      }
    }

    iy = Ny-1;
    site = ix + Nx * (iy + Ny * izt);

    load_u(ut, up, site + Nst * (1-ieo + 2*idir));
    for (int is = 0; is < Ns; ++is) {
      real_t wt[NVCD], vt[NVC * ND2];

      if(jgm5 == 0){
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
        }
      }else{
        for(int ivc = 0; ivc < NVC; ++ivc){
          wt[ivc+ID1] = wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), site)];
          wt[ivc+ID2] = wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), site)];
          wt[ivc+ID3] = wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), site)];
          wt[ivc+ID4] = wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), site)];
        }
      }
      mult_wilson_ym1(vt, ut, wt);
      for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
          buf_ym[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixzt)] = vt[ivcd];
      }
    }
  }
}

__global__  
void mult_domainwall_5din_eo_hop1z_dirac_dev(
  real_t * __restrict__ buf_zp, real_t * __restrict__ buf_zm,
  real_t * __restrict__ up, real_t * __restrict__ wp,
  int Ns, int bc, int Nx, int Ny, int Nz, int Nt, 
  int ieo, int jeo, int jgm5){

  const int ixyt   = blockIdx.x * blockDim.x + threadIdx.x;
  const int Nst    = Nx * Ny * Nz * Nt;
  const int Nxy    = Nx * Ny;
  const int Nxyt   = Nxy * Nt;
  const int Nin5   = NVCD * Ns;
  const int Nin5bd = NVC * ND2 * Ns;

  int idir = 2;
  real_t ut[NDF];

  if ( ixyt < Nxyt ){

    int iz = 0;
    int ixy = ixyt % Nxy;
    int it  = ixyt / Nxy;
    int site = ixy + Nxy * (iz + Nz * it);
    real_t bc2 = bc;
    for (int is = 0; is < Ns; ++is) {
      real_t wt[NVCD], vt[NVC * ND2];
      if(jgm5 == 0){
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
        }
      }else{
        for(int ivc = 0; ivc < NVC; ++ivc){
          wt[ivc+ID1] = wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), site)];
          wt[ivc+ID2] = wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), site)];
          wt[ivc+ID3] = wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), site)];
          wt[ivc+ID4] = wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), site)];
        }
      }
      mult_wilson_zp1(vt, wt);
      for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
          buf_zp[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixyt)] = bc2 * vt[ivcd];
      }
    }

    iz = Nz-1;
    site = ixy + Nxy * (iz + Nz * it);
    load_u(ut, up, site + Nst * (1-ieo + 2*idir));
    for (int is = 0; is < Ns; ++is) {
      real_t wt[NVCD], vt[NVC * ND2];
      if(jgm5 == 0){
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
        }
      }else{
        for(int ivc = 0; ivc < NVC; ++ivc){
          wt[ivc+ID1] = wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), site)];
          wt[ivc+ID2] = wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), site)];
          wt[ivc+ID3] = wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), site)];
          wt[ivc+ID4] = wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), site)];
        }
      }
      mult_wilson_zm1(vt, ut, wt);
      for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
          buf_zm[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixyt)] = vt[ivcd];
      }
    }
  }
}

__global__  void mult_domainwall_5din_eo_hop1t_dirac_dev(
  real_t * __restrict__ buf_tp, real_t * __restrict__ buf_tm,
  real_t * __restrict__ up, real_t * __restrict__ wp,
  int Ns, int bc, 
  int Nx, int Ny, int Nz, int Nt, 
  int ieo, int jeo, int jgm5){

  const int site_thread  = blockIdx.x * blockDim.x + threadIdx.x;
  const int ixyz         = site_thread;
  const int Nst          = Nx * Ny * Nz * Nt;
  const int Nxyz         = Nx * Ny * Nz;
  const int Nin5         = NVCD * Ns;
  const int Nin5bd       = NVC  * ND2 * Ns;
  int idir = 3;
  real_t ut[NDF];
  if ( ixyz < Nxyz ){

    int it = 0;
    int site = ixyz + Nxyz * it;
    real_t bc2 = bc;

    for (int is = 0; is < Ns; ++is) {
      real_t wt[NVCD], vt[NVC * ND2];
      if(jgm5 == 0){
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
        }
      }else{
        for(int ivc = 0; ivc < NVC; ++ivc){
          wt[ivc+ID1] = wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), site)];
          wt[ivc+ID2] = wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), site)];
          wt[ivc+ID3] = wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), site)];
          wt[ivc+ID4] = wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), site)];
        }
      }
      mult_wilson_tp1_dirac(vt, wt);
      for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
        buf_tp[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixyz)] = bc2 * vt[ivcd];
      }
    }

    it = Nt-1;
    site = ixyz + Nxyz * it;
    load_u(ut, up, site + Nst * (1-ieo + 2*idir));

    for (int is = 0; is < Ns; ++is) {
      real_t wt[NVCD], vt[NVC * ND2];

      if(jgm5 == 0){
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
        }
      }else{
        for(int ivc = 0; ivc < NVC; ++ivc){
          wt[ivc+ID1] = wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), site)];
          wt[ivc+ID2] = wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), site)];
          wt[ivc+ID3] = wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), site)];
          wt[ivc+ID4] = wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), site)];
        }
      }
      mult_wilson_tm1_dirac(vt, ut, wt);
      for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
          buf_tm[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixyz)] = vt[ivcd];
      }
    }
  }
}

__global__ void mult_domainwall_5din_eo_hop1x_5D_dirac_dev(
  real_t * __restrict__ buf_xp, real_t *__restrict__ buf_xm,
  real_t * __restrict__ up, real_t *__restrict__ wp,
  int Ns, int bc, 
  int Nx, int Ny, int Nz, int Nt, 
  int ieo, int jeo, int jgm5){

  const int Nyzt   = Ny * Nz * Nt;
  const int Nin5   = NVCD * Ns;
  const int Nin5bd = NVC * ND2 * Ns;
  const int Nst    = Nx * Ny * Nz * Nt;

  const int ist      = blockIdx.x * blockDim.x + threadIdx.x;
  const int gridSize = gridDim.x * blockDim.x;

  int idir   = 0;
  real_t ut[NDF];

  for (int idx = ist; idx < Nyzt * Ns; idx += gridSize){

    int is   = idx % Ns;
    int iyzt = idx / Ns;
    int iy   = iyzt % Ny;
    int iz   = (iyzt / Ny) % Nz;
    int it   = iyzt / (Ny * Nz);
    int keo  = (jeo + iy + iz + it) % 2;

    if (keo == 1) {
      int ix    = 0;
      int iyzt2 = iyzt / 2;
      int site  = ix + Nx * iyzt;

      real_t bc2 = bc;
      real_t wt[NVCD], vt[NVC * ND2];

      if(jgm5 == 0){
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
        }
      }else{
        for(int ivc = 0; ivc < NVC; ++ivc){
          wt[ivc+ID1] = wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), site)];
          wt[ivc+ID2] = wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), site)];
          wt[ivc+ID3] = wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), site)];
          wt[ivc+ID4] = wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), site)];
        }
      }
      mult_wilson_xp1(vt, wt);
      for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
        buf_xp[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), iyzt2)] = bc2 * vt[ivcd];
      }
    }

    if (keo == 0) {
      int iyzt2 = iyzt / 2;
      int ix = Nx - 1;
      int site = ix + Nx * iyzt;
      load_u(ut, up, site + Nst * (1 - ieo + 2 * idir));
      real_t wt[NVCD], vt[NVC * ND2];

      if(jgm5 == 0){
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
        }
      }else{
        for(int ivc = 0; ivc < NVC; ++ivc){
          wt[ivc+ID1] = wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), site)];
          wt[ivc+ID2] = wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), site)];
          wt[ivc+ID3] = wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), site)];
          wt[ivc+ID4] = wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), site)];
        }
      }
      mult_wilson_xm1(vt, ut, wt);
      for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
        buf_xm[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), iyzt2)] = vt[ivcd];
      }
    }
  }
}

__global__ void mult_domainwall_5din_eo_hop1y_5D_dirac_dev(
  real_t * __restrict__ buf_yp, real_t * __restrict__ buf_ym,
  real_t * __restrict__ up, real_t * __restrict__ wp,
  int Ns, int bc, 
  int Nx, int Ny, int Nz, int Nt, 
  int ieo, int jeo, int jgm5) 
{
  const int Nxzt   = Nx * Nz * Nt;
  const int Nin5   = NVCD * Ns;
  const int Nin5bd = NVC * ND2 * Ns;
  const int Nst    = Nx * Ny * Nz * Nt;

  const int ist      = blockIdx.x * blockDim.x + threadIdx.x;
  const int gridSize = gridDim.x *blockDim.x;

  int idir = 1;
  real_t ut[NDF];

  for (int idx = ist; idx < Nxzt * Ns; idx += gridSize){

    int is   = idx % Ns;
    int ixzt = idx / Ns;
    int ix   = ixzt % Nx;
    int izt  = ixzt / Nx;
    int site = ix + Nx * (0 + Ny * izt);

    real_t bc2 = bc;
    real_t wt[NVCD], vt[NVC * ND2];

    if(jgm5 == 0){
      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
      }
    }else{
      for(int ivc = 0; ivc < NVC; ++ivc){
        wt[ivc+ID1] = wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), site)];
        wt[ivc+ID2] = wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), site)];
        wt[ivc+ID3] = wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), site)];
        wt[ivc+ID4] = wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), site)];
      }
    }
    mult_wilson_yp1(vt, wt);
    for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
      buf_yp[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixzt)] = bc2 * vt[ivcd];
    }

    site = ix + Nx * (Ny - 1 + Ny * izt);
    load_u(ut, up, site + Nst * (1 - ieo + 2 * idir));
    if(jgm5 == 0){
      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
      }
    }else{
      for(int ivc = 0; ivc < NVC; ++ivc){
        wt[ivc+ID1] = wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), site)];
        wt[ivc+ID2] = wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), site)];
        wt[ivc+ID3] = wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), site)];
        wt[ivc+ID4] = wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), site)];
      }
    }
    mult_wilson_ym1(vt, ut, wt);
    for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
      buf_ym[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixzt)] = vt[ivcd];
    }
  }
}

__global__ void mult_domainwall_5din_eo_hop1z_5D_dirac_dev(
  real_t *__restrict__ buf_zp, real_t *__restrict__ buf_zm,
  real_t *__restrict__ up, real_t *__restrict__ wp,
  int Ns, int bc, 
  int Nx, int Ny, int Nz, int Nt, 
  int ieo, int jeo, int jgm5) 
{
  const int Nxy    = Nx * Ny;
  const int Nxyt   = Nxy * Nt;
  const int Nin5   = NVCD * Ns;
  const int Nin5bd = NVC * ND2 * Ns;
  const int Nst    = Nx * Ny * Nz * Nt;

  const int ist      = blockIdx.x * blockDim.x + threadIdx.x;
  const int GridSize = gridDim.x * blockDim.x;

  const int idir = 2;
  real_t ut[NDF];

  for (int idx = ist; idx < Nxyt * Ns; idx += GridSize){

    int is   = idx % Ns;
    int ixyt = idx / Ns;
    int ixy  = ixyt % Nxy;
    int it   = ixyt / Nxy;
    int site = ixy + Nxy * (0 + Nz * it);

    real_t bc2 = bc;
    real_t wt[NVCD], vt[NVC * ND2];

    if(jgm5 == 0){
      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
      }
    }else{
      for(int ivc = 0; ivc < NVC; ++ivc){
        wt[ivc+ID1] = wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), site)];
        wt[ivc+ID2] = wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), site)];
        wt[ivc+ID3] = wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), site)];
        wt[ivc+ID4] = wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), site)];
      }
    }
    mult_wilson_zp1(vt, wt);
    for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
        buf_zp[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixyt)] = bc2 * vt[ivcd];
    }

    site = ixy + Nxy * (Nz - 1 + Nz * it);
    load_u(ut, up, site + Nst * (1 - ieo + 2 * idir));
    if(jgm5 == 0){
      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
      }
    }else{
      for(int ivc = 0; ivc < NVC; ++ivc){
        wt[ivc+ID1] = wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), site)];
        wt[ivc+ID2] = wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), site)];
        wt[ivc+ID3] = wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), site)];
        wt[ivc+ID4] = wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), site)];
      }
    }
    mult_wilson_zm1(vt, ut, wt);
    for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
      buf_zm[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixyt)] = vt[ivcd];
    }
  }
}

__global__ void mult_domainwall_5din_eo_hop1t_5D_dirac_dev(
  real_t *__restrict__ buf_tp, real_t *__restrict__ buf_tm,
  real_t *__restrict__ up, real_t *__restrict__ wp,
  int Ns, int bc, 
  int Nx, int Ny, int Nz, int Nt, int ieo, int jeo, int jgm5){

  const int Nxyz     = Nx * Ny * Nz;
  const int Nin5     = NVCD * Ns;
  const int Nin5bd   = NVC * ND2 * Ns;
  const int Nst      = Nx * Ny * Nz * Nt;
  const int idir     = 3;

  const int ist      = blockIdx.x * blockDim.x + threadIdx.x;
  const int GridSize = blockDim.x * gridDim.x;

  real_t ut[NDF];

  for (int idx = ist; idx < Nxyz * Ns; idx += GridSize){

    int is   = idx  % Ns;
    int ixyz = idx  / Ns;
    int site = ixyz + Nxyz * 0;

    real_t bc2 = bc;
    real_t wt[NVCD], vt[NVC * ND2];

    if(jgm5 == 0){
      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
      }
    }else{
      for(int ivc = 0; ivc < NVC; ++ivc){
        wt[ivc+ID1] = wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), site)];
        wt[ivc+ID2] = wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), site)];
        wt[ivc+ID3] = wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), site)];
        wt[ivc+ID4] = wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), site)];
      }
    }
    mult_wilson_tp1_dirac(vt, wt);
    for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
        buf_tp[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixyz)] = bc2 * vt[ivcd];
    }

    site = ixyz + Nxyz * (Nt - 1);
    load_u(ut, up, site + Nst * (1 - ieo + 2 * idir));
    if(jgm5 == 0){
      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
      }
    }else{
      for(int ivc = 0; ivc < NVC; ++ivc){
        wt[ivc+ID1] = wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), site)];
        wt[ivc+ID2] = wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), site)];
        wt[ivc+ID3] = wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), site)];
        wt[ivc+ID4] = wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), site)];
      }
    }
    mult_wilson_tm1_dirac(vt, ut, wt);
    for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
      buf_tm[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixyz)] = vt[ivcd];
    }
  }
}

void mult_domainwall_5din_eo_hop1_dirac(
  real_t * buf1_xp, real_t * buf1_xm,
  real_t * buf1_yp, real_t * buf1_ym,
  real_t * buf1_zp, real_t * buf1_zm,
  real_t * buf1_tp, real_t * buf1_tm,
  real_t * up, real_t * wp,
  int Ns, int *bc, int *Nsize, int *do_comm,
  int ieo, int jeo, int jgm5){

  int Nx     = Nsize[0];
  int Ny     = Nsize[1];
  int Nz     = Nsize[2];
  int Nt     = Nsize[3];
  int Nin5bd = NVC * Ns * ND2;

  int bc_x = bc[0];
  int bc_y = bc[1];
  int bc_z = bc[2];
  int bc_t = bc[3];

  real_t * up_dev = (real_t*)dev_ptr(up);
  real_t * wp_dev = (real_t*)dev_ptr(wp);

  if( do_comm[0] > 0 ){

    int Nstx =  Ny * Nz * Nt;

    int blockSize  =  VECTOR_LENGTH;
    int sharedsize = NDF * sizeof(real_t);
    int gridSize5D =  (Nstx * Ns + blockSize - 1)/ blockSize;
    mult_domainwall_5din_eo_hop1x_5D_dirac_dev<<<gridSize5D, blockSize, sharedsize>>>(
    buf1_xp, buf1_xm,
    up_dev, wp_dev,
    Ns, bc_x, Nx, Ny, Nz, Nt, ieo, jeo, jgm5);

  }

  if( do_comm[1] > 0){

    int Nsty =  Nx * Nz * Nt;

    int blockSize  = VECTOR_LENGTH;
    int sharedsize = NDF * sizeof(real_t);
    int gridSize5D = (Nsty * Ns + blockSize - 1)/ blockSize;
    mult_domainwall_5din_eo_hop1y_5D_dirac_dev<<<gridSize5D, blockSize, sharedsize>>>(
    buf1_yp, buf1_ym,
    up_dev, wp_dev,
    Ns, bc_y, Nx, Ny, Nz, Nt, ieo, jeo, jgm5);

  }

  if( do_comm[2] > 0 ){

    int Nstz =  Nx * Ny * Nt;

    int blockSize  = VECTOR_LENGTH;
    int sharedsize = NDF* sizeof(real_t);
    int gridSize5D = (Nstz * Ns + blockSize - 1)/ blockSize;
    mult_domainwall_5din_eo_hop1z_5D_dirac_dev<<<gridSize5D, blockSize, sharedsize>>>(
    buf1_zp, buf1_zm,
    up_dev, wp_dev,
    Ns, bc_z, Nx, Ny, Nz, Nt, ieo, jeo, jgm5);

  }

  if( do_comm[3] > 0 ){

    int Nstt =  Nx * Ny * Nz;

    int blockSize = VECTOR_LENGTH;
    int sharedsize = NDF * sizeof(real_t);
    int gridSize5D  = (Nstt * Ns + blockSize - 1)/ blockSize;
    mult_domainwall_5din_eo_hop1t_5D_dirac_dev<<<gridSize5D, blockSize, sharedsize>>>(
    buf1_tp, buf1_tm,
    up_dev, wp_dev,
    Ns, bc_t, Nx, Ny, Nz, Nt, ieo, jeo, jgm5);

  }
  CHECK(cudaDeviceSynchronize());

};

//==============================================================
__global__ void mult_domainwall_5din_eo_hop2_dirac_dev(
  real_t *__restrict__ vp, real_t *__restrict__ up, real_t *__restrict__ wp,
  real_t *__restrict__ buf_xp, real_t *__restrict__ buf_xm,
  real_t *__restrict__ buf_yp, real_t *__restrict__ buf_ym,
  real_t *__restrict__ buf_zp, real_t *__restrict__ buf_zm,
  real_t *__restrict__ buf_tp, real_t *__restrict__ buf_tm,
  int Ns,
  int bc_x, int bc_y, int bc_z, int bc_t,
  int do_comm_x, int do_comm_y, int do_comm_z, int do_comm_t,
  int Nx, int Ny, int Nz, int Nt, int ieo, int jeo){

  const int ist      = blockIdx.x * blockDim.x + threadIdx.x;
  const int GridSize = gridDim.x * blockDim.x;

  const int Nxy    = Nx * Ny;
  const int Nxyz   = Nx   * Ny  * Nz     ;
  const int Nin5bd = (NVCD / 2) * Ns     ;
  const int Nst    = Nx   * Ny  * Nz * Nt;
  const int Nin5   = NVCD * Ns           ;
  real_t ut[NDF];
  for (int site = ist; site < Nst; site += GridSize){

    int ix   = site % Nx;
    int iyzt = site / Nx;
    int ixy  = site % Nxy;
    int iy   = iyzt % Ny;
    int izt  = site / Nxy;
    int iz   = izt % Nz;
    int it   = izt / Nz;
    int ixyz = site % Nxyz;
    int keo  = (jeo + iy + iz + it) % 2;
    int idir;

    for (int is = 0; is < Ns; ++is)
    {

      real_t vL[NVCD];

      for (int ivcd = 0; ivcd < NVCD; ++ivcd)
      {
        vL[ivcd] = 0.0;
      }

      int opr_any = 0;

      idir = 0;
      if (do_comm_x > 0)
      {

        if (ix == Nx - 1 && keo == 1)
        {
          load_u(ut, up, site + Nst * (ieo + 2 * idir));
          real_t wt[NVC * ND2];
          int iyzt2 = iyzt / 2;
          for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd)
          {
            wt[ivcd] = buf_xp[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), iyzt2)];
          }
          mult_wilson_xp2(vL, ut, wt);
          ++opr_any;
        }

        if (ix == 0 && keo == 0)
        {
          real_t bc2 = bc_x;
          int iyzt2 = iyzt / 2;
          real_t wt[NVC * ND2];
          for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd)
          {
            wt[ivcd] = bc2 * buf_xm[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), iyzt2)];
          }
          mult_wilson_xm2(vL, wt);
          ++opr_any;
        }
      }

      idir = 1;
      if (do_comm_y > 0)
      {
        int ixzt = ix + Nx * izt;
        if (iy == Ny - 1)
        {
          load_u(ut, up, site + Nst * (ieo + 2 * idir));
          real_t wt[NVC * ND2];
          for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd)
          {
            wt[ivcd] = buf_yp[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixzt)];
          }
          mult_wilson_yp2(vL, ut, wt);
          ++opr_any;
        }

        if (iy == 0)
        {
          real_t bc2 = bc_y;
          real_t wt[NVC * ND2];
          for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd)
          {
            wt[ivcd] = bc2 * buf_ym[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixzt)];
          }
          mult_wilson_ym2(vL, wt);
          ++opr_any;
        }
      }

      idir = 2;
      if (do_comm_z > 0)
      {
        int ixyt = ixy + Nxy * it;
        if (iz == Nz - 1)
        {
          load_u(ut, up, site + Nst * (ieo + 2 * idir));
          real_t wt[NVC * ND2];
          for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd)
          {
            wt[ivcd] = buf_zp[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixyt)];
          }
          mult_wilson_zp2(vL, ut, wt);
          ++opr_any;
        }

        if (iz == 0)
        {
          real_t bc2 = bc_z;
          real_t wt[NVC * ND2];
          for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd)
          {
            wt[ivcd] = bc2 * buf_zm[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixyt)];
          }
          mult_wilson_zm2(vL, wt);
          ++opr_any;
        }
      }

      idir = 3;
      if (do_comm_t > 0)
      {
        if (it == Nt - 1)
        {
          load_u(ut, up, site + Nst * (ieo + 2 * idir));
          real_t wt[NVC * ND2];
          for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd)
          {
            wt[ivcd] = buf_tp[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixyz)];
          }
          mult_wilson_tp2_dirac(vL, ut, wt);
          ++opr_any;
        }

        if (it == 0)
        {
          real_t bc2 = bc_t;
          real_t wt[NVC * ND2];
          for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd)
          {
            wt[ivcd] = bc2 * buf_tm[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixyz)];
          }
          mult_wilson_tm2_dirac(vL, wt);
          ++opr_any;
        }
      }
      if (opr_any > 0)
      {
        for (int ivcd = 0; ivcd < NVCD; ++ivcd)
        {
          vp[IDX2(Nin5, (ivcd + NVCD * is), site)] += vL[ivcd];
        }
      }
    } // is loop
  }
}

__global__ void mult_domainwall_5din_eo_hop2_5D_dirac_dev(
  real_t *__restrict__ vp, real_t *__restrict__ up, real_t *__restrict__ wp,
  real_t *__restrict__ buf_xp, real_t *__restrict__ buf_xm,
  real_t *__restrict__ buf_yp, real_t *__restrict__ buf_ym,
  real_t *__restrict__ buf_zp, real_t *__restrict__ buf_zm,
  real_t *__restrict__ buf_tp, real_t *__restrict__ buf_tm,
  int Ns,
  int bc_x, int bc_y, int bc_z, int bc_t,
  int do_comm_x, int do_comm_y, int do_comm_z, int do_comm_t,
  int Nx, int Ny, int Nz, int Nt, int ieo, int jeo, int Nst_pad)
{
  const int Nxy    = Nx * Ny;
  const int Nxyz   = Nx * Ny * Nz;
  const int Nin5bd = (NVCD / 2) * Ns;
  const int Nst    = Nx * Ny * Nz * Nt;
  const int Nin5   = NVCD * Ns;

  const int ist      = blockIdx.x * blockDim.x + threadIdx.x;
  const int GridSize = blockDim.x * gridDim.x;

  real_t ut[NDF];

  for (int idx = ist; idx < Nst_pad * Ns; idx += GridSize){

    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int is      = idx2_wp % Ns;
    int idx_out = idx2_wp / Ns;
    int site    = idx_in + NWP * idx_out;
    int ix      = site % Nx;
    int iyzt    = site / Nx;
    int ixy     = site % Nxy;
    int iy      = iyzt % Ny;
    int izt     = site / Nxy;
    int iz      = izt  % Nz;
    int it      = izt  / Nz;
    int ixyz    = site % Nxyz;
    int keo     = (jeo + iy + iz + it) % 2;
    int idir;

    real_t vL[NVCD];
    for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
      vL[ivcd] = 0.0;
    }

    int opr_any = 0;

    idir = 0;
    if (do_comm_x > 0) {
      if (ix == Nx - 1 && keo == 1) {
        load_u(ut, up, site + Nst * (ieo + 2 * idir));
        real_t wt[NVC * ND2];
        int iyzt2 = iyzt / 2;
        for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
          wt[ivcd] = buf_xp[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), iyzt2)];
        }
        mult_wilson_xp2(vL, ut, wt);
        ++opr_any;
      }

      if (ix == 0 && keo == 0) {
        real_t bc2 = bc_x;
        int iyzt2 = iyzt / 2;
        real_t wt[NVC * ND2];
        for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
          wt[ivcd] = bc2 * buf_xm[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), iyzt2)];
        }
        mult_wilson_xm2(vL, wt);
        ++opr_any;
      }
    }

    idir = 1;
    if (do_comm_y > 0) {
      int ixzt = ix + Nx * izt;
      if (iy == Ny - 1) {
        load_u(ut, up, site + Nst * (ieo + 2 * idir));
        real_t wt[NVC * ND2];
        for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
          wt[ivcd] = buf_yp[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixzt)];
        }
        mult_wilson_yp2(vL, ut, wt);
        ++opr_any;
      }

      if (iy == 0) {
        real_t bc2 = bc_y;
        real_t wt[NVC * ND2];
        for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
          wt[ivcd] = bc2 * buf_ym[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixzt)];
        }
        mult_wilson_ym2(vL, wt);
        ++opr_any;
      }
    }

    idir = 2;
    if (do_comm_z > 0) {
      int ixyt = ixy + Nxy * it;
      if (iz == Nz - 1) {
        load_u(ut, up, site + Nst * (ieo + 2 * idir));
        real_t wt[NVC * ND2];
        for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
          wt[ivcd] = buf_zp[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixyt)];
        }
        mult_wilson_zp2(vL, ut, wt);
        ++opr_any;
      }

      if (iz == 0) {
        real_t bc2 = bc_z;
        real_t wt[NVC * ND2];
        for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
          wt[ivcd] = bc2 * buf_zm[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixyt)];
        }
        mult_wilson_zm2(vL, wt);
        ++opr_any;
      }
    }

    idir = 3;
    if (do_comm_t > 0) {
      if (it == Nt - 1) {
        load_u(ut, up, site + Nst * (ieo + 2 * idir));
        real_t wt[NVC * ND2];
        for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
          wt[ivcd] = buf_tp[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixyz)];
        }
        mult_wilson_tp2_dirac(vL, ut, wt);
        ++opr_any;
      }

      if (it == 0) {
        real_t bc2 = bc_t;
        real_t wt[NVC * ND2];
        for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
          wt[ivcd] = bc2 * buf_tm[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixyz)];
        }
        mult_wilson_tm2_dirac(vL, wt);
        ++opr_any;
      }
    }

    if (opr_any > 0) {
      for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
        vp[IDX2(Nin5, (ivcd + NVCD * is), site)] += vL[ivcd];
      }
    }
  }
}


void mult_domainwall_5din_eo_hop2_dirac(
  real_t * vp     , real_t * up, real_t * wp,
  real_t * buf2_xp, real_t * buf2_xm,
  real_t * buf2_yp, real_t * buf2_ym,
  real_t * buf2_zp, real_t * buf2_zm,
  real_t * buf2_tp, real_t * buf2_tm,
  int Ns, int *bc, int *Nsize, int *do_comm, int ieo, int jeo)
{
  int Nx      = Nsize[0];
  int Ny      = Nsize[1];
  int Nz      = Nsize[2];
  int Nt      = Nsize[3];
  int Nst     = Nx * Ny * Nz * Nt;
  int Nst_pad = ceil_nwp(Nst);  
  int Nin5bd  = (NVCD / 2) * Ns;

  real_t * vp_dev = (real_t*)dev_ptr(vp);
  real_t * up_dev = (real_t*)dev_ptr(up);
  real_t * wp_dev = (real_t*)dev_ptr(wp);

  int blockSize  = VECTOR_LENGTH;
  int sharedsize = NDF * sizeof(real_t);
  int gridSize5D = (Nst_pad * Ns + blockSize - 1)/ blockSize;
  mult_domainwall_5din_eo_hop2_5D_dirac_dev<<<gridSize5D, blockSize, sharedsize>>>(
    vp_dev, up_dev, wp_dev,
    buf2_xp, buf2_xm,
    buf2_yp, buf2_ym,
    buf2_zp, buf2_zm,
    buf2_tp, buf2_tm,
    Ns,
    bc[0], bc[1], bc[2], bc[3],
    do_comm[0], do_comm[1], do_comm[2], do_comm[3],
    Nx, Ny, Nz, Nt, ieo, jeo, Nst_pad);

  CHECK(cudaDeviceSynchronize());
}


#endif
//=================================
