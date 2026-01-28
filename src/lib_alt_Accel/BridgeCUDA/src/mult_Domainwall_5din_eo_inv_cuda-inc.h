/*!
      @file    mult_Doainwall_5din_eo_inv_cuda-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: wlchen $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2621 $
*/

#ifndef MULT_DOMAINWALL_5DIN_EO_INV_ACC_INCLUDED
#define MULT_DOMAINWALL_5DIN_EO_INV_ACC_INCLUDED

#include "../inline/constant_memory_inline.h"

//====================================================================
__global__
void mult_domainwall_5din_ee_inv_dirac_4d_kernel(
                                   real_t * vp, real_t * wp,
                                   int jd, int Ns, real_t * mat_inv,
                                   int Nx, int Ny, int Nz, int Nt)
{
  int Nst  = Nx * Ny * Nz * Nt;
  int Nin5 = NVCD * Ns;

  //int site = blockIdx.x * blockDim.x + threadIdx.x;
  //if (site < Nst) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < Nst * NC) {

    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int ic      = idx2_wp % NC;
    int idx_out = idx2_wp / NC;
    int site    = idx_in + NWP * idx_out;

    for (int is1 = 0; is1 < Ns; ++is1) {
      // for (int ic = 0; ic < NC; ++ic) {
      {
        for (int id = 0; id < ND2; ++id) {

          real_t vt1r = 0.0;
          real_t vt1i = 0.0;
          real_t vt2r = 0.0;
          real_t vt2i = 0.0;
          int idx1 = id + ND2 * is1;
          for (int is2 = 0; is2 < Ns; ++is2) {
            int idx2, ivcd1, ivcd2;
            real_t mat1, mat2;

	    if(jd == 1){
              idx2 = 0 + ND2 * is2;
              mat1 = mat_inv[idx1 + (ND2 * Ns) * idx2];
              idx2 = 1 + ND2 * is2;
              mat2 = mat_inv[idx1 + (ND2 * Ns) * idx2];
            }else{
              idx2 = 0 + ND2 * is2;
              mat1 = mat_inv[idx2 + (ND2 * Ns) * idx1];
              idx2 = 1 + ND2 * is2;
              mat2 = mat_inv[idx2 + (ND2 * Ns) * idx1];
            }	    

            ivcd1 = 0 + 2 * (ic + NC * 0);
            ivcd2 = 0 + 2 * (ic + NC * 2);
            vt1r +=   mat1 * wp[IDX2(Nin5, (ivcd1 + NVCD * is2), site)]
                    + mat2 * wp[IDX2(Nin5, (ivcd2 + NVCD * is2), site)];

            ivcd1 = 1 + 2 * (ic + NC * 0);
            ivcd2 = 1 + 2 * (ic + NC * 2);
            vt1i +=   mat1 * wp[IDX2(Nin5, (ivcd1 + NVCD * is2), site)]
                    + mat2 * wp[IDX2(Nin5, (ivcd2 + NVCD * is2), site)];

            ivcd1 = 0 + 2 * (ic + NC * 1);
            ivcd2 = 0 + 2 * (ic + NC * 3);
            vt2r +=   mat1 * wp[IDX2(Nin5, (ivcd1 + NVCD * is2), site)]
                    + mat2 * wp[IDX2(Nin5, (ivcd2 + NVCD * is2), site)];

            ivcd1 = 1 + 2 * (ic + NC * 1);
            ivcd2 = 1 + 2 * (ic + NC * 3);
            vt2i +=   mat1 * wp[IDX2(Nin5, (ivcd1 + NVCD * is2), site)]
                    + mat2 * wp[IDX2(Nin5, (ivcd2 + NVCD * is2), site)];
	  }
          int id1 = 2 * id;
          int id2 = 2 * id + 1;
          vp[IDX2(Nin5, (0+2*(ic+NC*id1) + NVCD * is1), site)] = vt1r;
          vp[IDX2(Nin5, (1+2*(ic+NC*id1) + NVCD * is1), site)] = vt1i;
          vp[IDX2(Nin5, (0+2*(ic+NC*id2) + NVCD * is1), site)] = vt2r;
          vp[IDX2(Nin5, (1+2*(ic+NC*id2) + NVCD * is1), site)] = vt2i;
	}
      }
    }
  }

}

void mult_domainwall_5din_ee_inv_dirac_4d(
                        real_t * vp, real_t * wp, int jd,
                        int Ns, real_t * mat_inv, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  real_t* vp_dev = (real_t*)dev_ptr(vp);
  real_t* wp_dev = (real_t*)dev_ptr(wp);

  int mat_size = ND2 * Ns;
  size_t mat_size2 = mat_size * mat_size * sizeof(real_t);
  real_t* mat_inv_dev;
  CHECK(cudaMalloc((void **)&mat_inv_dev, mat_size2));
  CHECK(cudaMemcpy(mat_inv_dev, mat_inv, mat_size2, cudaMemcpyHostToDevice));

  int blockSize = VECTOR_LENGTH;
  int gridSize  = (Nst * NC + blockSize - 1) / blockSize;

  mult_domainwall_5din_ee_inv_dirac_4d_kernel<<<gridSize, blockSize>>>(
                 vp_dev, wp_dev, jd, Ns, mat_inv_dev, Nx, Ny, Nz, Nt);

  CHECK(cudaDeviceSynchronize());

  CHECK(cudaFree(mat_inv_dev));

}

//====================================================================
__global__
void mult_domainwall_5din_ee_inv_dirac_5d_kernel(
                                   real_t * vp, real_t * wp,
                                   int jd, int Ns, real_t * mat_inv,
                                   int Nx, int Ny, int Nz, int Nt)
// still not 5d implementation
{
  int Nst  = Nx * Ny * Nz * Nt;
  int Nin5 = NVCD * Ns;
  int Nst_pad = CEIL_NWP(Nst);
  
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (idx < Nst_pad * NC * Ns) {

    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int ic      = idx2_wp % NC;
    int is1     = (idx2_wp/NC) % Ns;
    int idx_out = idx2_wp/(NC * Ns);
    int site = idx_in + NWP * idx_out;

      {
        for (int id = 0; id < ND2; ++id) {

          real_t vt1r = 0.0;
          real_t vt1i = 0.0;
          real_t vt2r = 0.0;
          real_t vt2i = 0.0;
          int idx1 = id + ND2 * is1;
          for (int is2 = 0; is2 < Ns; ++is2) {
            int idx2, ivcd1, ivcd2;
            real_t mat1, mat2;

	    if(jd == 1){
              idx2 = 0 + ND2 * is2;
              mat1 = mat_inv[idx1 + (ND2 * Ns) * idx2];
              idx2 = 1 + ND2 * is2;
              mat2 = mat_inv[idx1 + (ND2 * Ns) * idx2];
            }else{
              idx2 = 0 + ND2 * is2;
              mat1 = mat_inv[idx2 + (ND2 * Ns) * idx1];
              idx2 = 1 + ND2 * is2;
              mat2 = mat_inv[idx2 + (ND2 * Ns) * idx1];
            }	    

            ivcd1 = 0 + 2 * (ic + NC * 0);
            ivcd2 = 0 + 2 * (ic + NC * 2);
            vt1r +=   mat1 * wp[IDX2(Nin5, (ivcd1 + NVCD * is2), site)]
                    + mat2 * wp[IDX2(Nin5, (ivcd2 + NVCD * is2), site)];

            ivcd1 = 1 + 2 * (ic + NC * 0);
            ivcd2 = 1 + 2 * (ic + NC * 2);
            vt1i +=   mat1 * wp[IDX2(Nin5, (ivcd1 + NVCD * is2), site)]
                    + mat2 * wp[IDX2(Nin5, (ivcd2 + NVCD * is2), site)];

            ivcd1 = 0 + 2 * (ic + NC * 1);
            ivcd2 = 0 + 2 * (ic + NC * 3);
            vt2r +=   mat1 * wp[IDX2(Nin5, (ivcd1 + NVCD * is2), site)]
                    + mat2 * wp[IDX2(Nin5, (ivcd2 + NVCD * is2), site)];

            ivcd1 = 1 + 2 * (ic + NC * 1);
            ivcd2 = 1 + 2 * (ic + NC * 3);
            vt2i +=   mat1 * wp[IDX2(Nin5, (ivcd1 + NVCD * is2), site)]
                    + mat2 * wp[IDX2(Nin5, (ivcd2 + NVCD * is2), site)];
	  }
          int id1 = 2 * id;
          int id2 = 2 * id + 1;
          vp[IDX2(Nin5, (0+2*(ic+NC*id1) + NVCD * is1), site)] = vt1r;
          vp[IDX2(Nin5, (1+2*(ic+NC*id1) + NVCD * is1), site)] = vt1i;
          vp[IDX2(Nin5, (0+2*(ic+NC*id2) + NVCD * is1), site)] = vt2r;
          vp[IDX2(Nin5, (1+2*(ic+NC*id2) + NVCD * is1), site)] = vt2i;
	}
      }
      // }
  }

}

__global__
void mult_domainwall_5din_ee_inv_dirac_5d_kernel_opt(
                                   real_t * vp, real_t * wp,
                                   int jd, int Ns, real_t * __restrict__ mat_inv,
                                   int Nx, int Ny, int Nz, int Nt)
{
  const int Nst     = Nx * Ny * Nz * Nt;
  const int Nin5    = NVCD * Ns;
  const int Nst_pad = CEIL_NWP(Nst);
  const int idx     = blockIdx.x * blockDim.x + threadIdx.x;
  //__shared__ real_t shared_wp[8];
  __shared__ real_t mat1, mat2;
  __shared__ real_t wp1r, wp1i, wp2r, wp2i, wp3r, wp3i, wp4r, wp4i;

  if (idx < Nst_pad * NC * Ns) {

    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int ic      = idx2_wp % NC;
    int is1     = (idx2_wp/NC) % Ns;
    int idx_out = idx2_wp/(NC * Ns);
    int site    = idx_in + NWP * idx_out;

    for (int id = 0; id < ND2; ++id) {

      real_t vt1r = 0.0;
      real_t vt1i = 0.0;
      real_t vt2r = 0.0;
      real_t vt2i = 0.0;
      int idx1 = id + ND2 * is1;
      for (int is2 = 0; is2 < Ns; ++is2) {
        int idx2, ivcd1;

	      if(jd == 1){
          idx2 = 0 + ND2 * is2;
          mat1 = mat_inv[idx1 + (ND2 * Ns) * idx2];
          idx2 = 1 + ND2 * is2;
          mat2 = mat_inv[idx1 + (ND2 * Ns) * idx2];
        }else{
          idx2 = 0 + ND2 * is2;
          mat1 = mat_inv[idx2 + (ND2 * Ns) * idx1];
          idx2 = 1 + ND2 * is2;
          mat2 = mat_inv[idx2 + (ND2 * Ns) * idx1];
        }	    

        int ivcd1r = 0 + 2 * (ic + NC * 0);
        int ivcd1i = 0 + 2 * (ic + NC * 2);
        int ivcd2r = 1 + 2 * (ic + NC * 0);
        int ivcd2i = 1 + 2 * (ic + NC * 2);
        int ivcd3r = 0 + 2 * (ic + NC * 1);
        int ivcd3i = 0 + 2 * (ic + NC * 3);
        int ivcd4r = 1 + 2 * (ic + NC * 1);
        int ivcd4i = 1 + 2 * (ic + NC * 3);
        wp1r = wp[IDX2(Nin5, (ivcd1r + NVCD * is2), site)];
        wp1i = wp[IDX2(Nin5, (ivcd1i + NVCD * is2), site)];
        wp2r = wp[IDX2(Nin5, (ivcd2r + NVCD * is2), site)];
        wp2i = wp[IDX2(Nin5, (ivcd2i + NVCD * is2), site)];
        wp3r = wp[IDX2(Nin5, (ivcd3r + NVCD * is2), site)];
        wp3i = wp[IDX2(Nin5, (ivcd3i + NVCD * is2), site)];
        wp4r = wp[IDX2(Nin5, (ivcd4r + NVCD * is2), site)];
        wp4i = wp[IDX2(Nin5, (ivcd4i + NVCD * is2), site)];
        vt1r +=  mat1 * wp1r + mat2 * wp1i;
        vt1i +=  mat1 * wp2r + mat2 * wp2i;
        vt2r +=  mat1 * wp3r + mat2 * wp3i;
        vt2i +=  mat1 * wp4r + mat2 * wp4i;
	    }
      int id1 = 2 * id;
      int id2 = 2 * id + 1;
      vp[IDX2(Nin5, (0+2*(ic+NC*id1) + NVCD * is1), site)] = vt1r;
      vp[IDX2(Nin5, (1+2*(ic+NC*id1) + NVCD * is1), site)] = vt1i;
      vp[IDX2(Nin5, (0+2*(ic+NC*id2) + NVCD * is1), site)] = vt2r;
      vp[IDX2(Nin5, (1+2*(ic+NC*id2) + NVCD * is1), site)] = vt2i;
	  }
  }
}

void mult_domainwall_5din_ee_inv_dirac_5d(
                        real_t * vp, real_t * wp, int jd,
                        int Ns, real_t * mat_inv, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;
  int Nst_pad = ceil_nwp(Nst);

  real_t* vp_dev = (real_t*)dev_ptr(vp);
  real_t* wp_dev = (real_t*)dev_ptr(wp);

  int mat_size = ND2 * Ns;
  size_t mat_size2 = mat_size * mat_size * sizeof(real_t);

  real_t* mat_inv_dev;
  CHECK(cudaMalloc((void **)&mat_inv_dev, mat_size2));
  CHECK(cudaMemcpy(mat_inv_dev, mat_inv, mat_size2, cudaMemcpyHostToDevice));

  int blockSize = VECTOR_LENGTH;
  int gridSize  = (Nst_pad * NC * Ns + blockSize - 1) / blockSize;

  mult_domainwall_5din_ee_inv_dirac_5d_kernel_opt<<<gridSize, blockSize>>>(
                 vp_dev, wp_dev, jd, Ns, mat_inv_dev, Nx, Ny, Nz, Nt);

  CHECK(cudaDeviceSynchronize());

  CHECK(cudaFree(mat_inv_dev));
}

__global__ void mult_domainwall_5din_ee_LUinv_dirac_kernel(
                        real_t* __restrict__ vp, real_t* __restrict__ wp,
                        int Ns, int Nst_pad, real_t alpha)
{
  const int Nin5     = NVCD * Ns;
  const int ist      = blockIdx.x * blockDim.x + threadIdx.x;
  const int GridSize = blockDim.x * gridDim.x;

  const double* e_con     = const_e;
  const double* f_con     = const_f;
  const double* dpinv_con = const_dpinv;
  const double* dm_con    = const_dm; 
  
  for( int idx = ist; idx < Nst_pad * NVC; idx += GridSize){

    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int ivc     = idx2_wp % NVC;
    int idx_out = idx2_wp / NVC;
    int site    = idx_in  + NWP * idx_out;

    real_t vt[ND], yt[ND], xt[ND];

    int is = 0;
    for(int id = 0; id < ND; ++id){
      int ivcd = ivc + NVC * id;
      vt[id] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
    }
    for(int id = 0; id < ND; ++id){
      int ivcd = ivc + NVC * id;
      vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vt[id];
    }
    for(int id = 0; id < ND; ++id){
      yt[id] = e_con[0] * vt[id];
    }

    for (int is = 1; is < Ns-1; ++is) {

      for(int id = 0; id < ND; ++id){
        xt[id] = vt[id];
      }

      for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
        vt[id] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
      }

      real_t a = real_t(0.5) * dm_con[is] * dpinv_con[is-1];

      vt[0] += a * (xt[0] + xt[2]);
      vt[1] += a * (xt[1] + xt[3]);
      vt[2] += a * (xt[2] + xt[0]);
      vt[3] += a * (xt[3] + xt[1]);

      for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
	      vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vt[id];
      }

      for(int id = 0; id < ND; ++id){
        yt[id] += e_con[is] * vt[id];
      }

    }

    is = Ns-1;

    for(int id = 0; id < ND; ++id){
      xt[id] = vt[id];
    }

    for(int id = 0; id < ND; ++id){
      int ivcd = ivc + NVC * id;
      vt[id] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
    }

    real_t a = real_t(0.5) * dm_con[is] * dpinv_con[is-1];

    vt[0] += a * (xt[0] + xt[2]);
    vt[1] += a * (xt[1] + xt[3]);
    vt[2] += a * (xt[2] + xt[0]);
    vt[3] += a * (xt[3] + xt[1]);

    vt[0] += -0.5 * (yt[0] - yt[2]);
    vt[1] += -0.5 * (yt[1] - yt[3]);
    vt[2] += -0.5 * (yt[2] - yt[0]);
    vt[3] += -0.5 * (yt[3] - yt[1]);

    for(int id = 0; id < ND; ++id){
      int ivcd = ivc + NVC * id;
      vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vt[id];
    }
    // L_inv completed

    is = Ns-1;

    a = dpinv_con[Ns-1];
    real_t f1 = 0.5 *  (1.0 + alpha);
    real_t f2 = 0.5 * (-1.0 + alpha);
    real_t vt1, vt2, vt3, vt4;

    vt1 = vp[IDX2(Nin5, (ID1 + ivc + NVCD*is), site)];
    vt2 = vp[IDX2(Nin5, (ID2 + ivc + NVCD*is), site)];
    vt3 = vp[IDX2(Nin5, (ID3 + ivc + NVCD*is), site)];
    vt4 = vp[IDX2(Nin5, (ID4 + ivc + NVCD*is), site)];
    vt[0] = a * (f1 * vt1 + f2 * vt3); 
    vt[1] = a * (f1 * vt2 + f2 * vt4); 
    vt[2] = a * (f1 * vt3 + f2 * vt1); 
    vt[3] = a * (f1 * vt4 + f2 * vt2); 

    for(int id = 0; id < ND; ++id){
      int ivcd = ivc + NVC * id;
      vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vt[id];
    }
      
    yt[0] = 0.5 * (vt[0] + vt[2]);
    yt[1] = 0.5 * (vt[1] + vt[3]);
    yt[2] = 0.5 * (vt[2] + vt[0]);
    yt[3] = 0.5 * (vt[3] + vt[1]);

    for (int is = Ns-2; is >= 0; --is) {

      for(int id = 0; id < ND; ++id){
        xt[id] = vt[id];
      }

      for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
        vt[id] = vp[IDX2(Nin5, (ivcd + NVCD*is), site)];
      }

      real_t a = real_t(0.5) * dm_con[is];

      vt[0] += a * (xt[0] - xt[2]);
      vt[1] += a * (xt[1] - xt[3]);
      vt[2] += a * (xt[2] - xt[0]);
      vt[3] += a * (xt[3] - xt[1]);

      for(int id = 0; id < ND; ++id){
        vt[id] += - f_con[is] * yt[id];
      }

      real_t aa = dpinv_con[is];

      for(int id = 0; id < ND; ++id){
        vt[id] *= aa;
      }

      if(is == 0){
        real_t f1 = 0.5 * (1.0 + alpha);
        real_t f2 = 0.5 * (1.0 - alpha);
        vt1 = f1 * vt[0] + f2 * vt[2];
        vt2 = f1 * vt[1] + f2 * vt[3];
        vt3 = f1 * vt[2] + f2 * vt[0];
        vt4 = f1 * vt[3] + f2 * vt[1];
        vp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)] = vt1;
        vp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)] = vt2;
        vp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)] = vt3;
        vp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)] = vt4;
      }else{
        for(int id = 0; id < ND; ++id){
          int ivcd = ivc + NVC * id;
          vp[IDX2(Nin5, (ivcd + NVCD * is), site)] = vt[id];
        }
      }
    }
  } // idx loop
}

void mult_domainwall_5din_ee_LUinv_dirac(
                        real_t* vp, real_t* wp,
                        int Ns, int *Nsize,
                        real_t *e, real_t *f,
                        real_t *dpinv, real_t *dm, real_t alpha)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;
  int Nst_pad = ceil_nwp(Nst);

  real_t* vp_dev = (real_t*)dev_ptr(vp);
  real_t* wp_dev = (real_t*)dev_ptr(wp);

  int blockSize = 2 * VECTOR_LENGTH;
  int gridSize  = (Nst_pad * NVC)/ blockSize;
  //int nSM       = getNumSMs();

  mult_domainwall_5din_ee_LUinv_dirac_kernel<<<gridSize, blockSize>>>(
    vp_dev, wp_dev, Ns, Nst_pad, alpha);

  CHECK(cudaDeviceSynchronize());

}

__global__ void mult_domainwall_5din_ee_LUdaginv_dirac_kernel(
                        real_t * __restrict__ vp, real_t *__restrict__ wp,
                        int Ns, int Nst_pad, real_t alpha)
{
  const int Nin5     = NVCD * Ns;
  const int ist      = blockIdx.x * blockDim.x + threadIdx.x;
  const int GridSize = blockDim.x * gridDim.x;

  const double* e_con     = const_e;
  const double* f_con     = const_f;
  const double* dpinv_con = const_dpinv;
  const double* dm_con    = const_dm; 
  
  for( int idx = ist; idx < Nst_pad * NVC; idx += GridSize){

    real_t vt[ND], yt[ND], xt[ND];

    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int ivc     = idx2_wp % NVC;
    int idx_out = idx2_wp / NVC;
    int site    = idx_in + NWP * idx_out;

    int is = 0;

    real_t a = dpinv_con[0];
    real_t f1 = 0.5 * (1.0 + alpha);
    real_t f2 = 0.5 * (1.0 - alpha);
    {
      real_t vt1, vt2, vt3, vt4;
      vt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD*is), site)];
      vt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD*is), site)];
      vt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD*is), site)];
      vt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD*is), site)];
      vt[0] = a * (f1 * vt1 + f2 * vt3);
      vt[1] = a * (f1 * vt2 + f2 * vt4);
      vt[2] = a * (f1 * vt3 + f2 * vt1);
      vt[3] = a * (f1 * vt4 + f2 * vt2);
      }
  
    for(int id = 0; id < ND; ++id){
      int ivcd = ivc + NVC * id;
      vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vt[id];
      yt[id] = f_con[0] * vt[id];
    }

    for (int is = 1; is < Ns-1; ++is) {

      for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
        xt[id] = vt[id];
        vt[id] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
      }

      real_t a = real_t(0.5) * dm_con[is - 1];

      vt[0] += a * (xt[0] - xt[2]);
      vt[1] += a * (xt[1] - xt[3]);
      vt[2] += a * (xt[2] - xt[0]);
      vt[3] += a * (xt[3] - xt[1]);

      real_t aa = dpinv_con[is];

      for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
        vt[id] *= aa;
        vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vt[id];
        yt[id] += f_con[is] * vt[id];
      }

    }

    is = Ns-1;

    for(int id = 0; id < ND; ++id){
      xt[id] = vt[id];
    }

    for(int id = 0; id < ND; ++id){
      int ivcd = ivc + NVC * id;
      vt[id] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
    }

    a = real_t(0.5) * dm_con[is - 1];

    vt[0] += a * (xt[0] - xt[2]);
    vt[1] += a * (xt[1] - xt[3]);
    vt[2] += a * (xt[2] - xt[0]);
    vt[3] += a * (xt[3] - xt[1]);

    vt[0] += -0.5 * (yt[0] + yt[2]);
    vt[1] += -0.5 * (yt[1] + yt[3]);
    vt[2] += -0.5 * (yt[2] + yt[0]);
    vt[3] += -0.5 * (yt[3] + yt[1]);

    real_t aa = dpinv_con[is];

    for(int id = 0; id < ND; ++id){
      vt[id] *= aa;
    }

    real_t ff1 = 0.5 * ( 1.0 + alpha);
    real_t ff2 = 0.5 * (-1.0 + alpha);
    {
    real_t vt1, vt2, vt3, vt4;
    vt1 = ff1 * vt[0] + ff2 * vt[2];
    vt2 = ff1 * vt[1] + ff2 * vt[3];
    vt3 = ff1 * vt[2] + ff2 * vt[0];
    vt4 = ff1 * vt[3] + ff2 * vt[1];
    vp[IDX2(Nin5, (ID1 + ivc + NVCD*is), site)] = vt1;
    vp[IDX2(Nin5, (ID2 + ivc + NVCD*is), site)] = vt2;
    vp[IDX2(Nin5, (ID3 + ivc + NVCD*is), site)] = vt3;
    vp[IDX2(Nin5, (ID4 + ivc + NVCD*is), site)] = vt4;
    }
    // Udag_inv completed

    is = Ns-1;

    for(int id = 0; id < ND; ++id){
      int ivcd = ivc + NVC * id;
      vt[id] = vp[IDX2(Nin5, (ivcd + NVCD*is), site)];
      vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vt[id];
    }

    yt[0] = 0.5 * (vt[0] - vt[2]);
    yt[1] = 0.5 * (vt[1] - vt[3]);
    yt[2] = 0.5 * (vt[2] - vt[0]);
    yt[3] = 0.5 * (vt[3] - vt[1]);

    for (int is = Ns-2; is >= 0; --is) {

      for(int id = 0; id < ND; ++id){
        xt[id] = vt[id];
      }

      for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
        vt[id] = vp[IDX2(Nin5, (ivcd + NVCD*is), site)];
      }

      real_t a = real_t(0.5) * dm_con[is + 1] * dpinv_con[is];

      vt[0] += a * (xt[0] + xt[2]);
      vt[1] += a * (xt[1] + xt[3]);
      vt[2] += a * (xt[2] + xt[0]);
      vt[3] += a * (xt[3] + xt[1]);

      for(int id = 0; id < ND; ++id){
        vt[id] += -e_con[is] * yt[id];
      }
      
      for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
        vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vt[id];
      }
    }

  } // idx loop end

}

void mult_domainwall_5din_ee_LUdaginv_dirac(
                        real_t* vp, real_t* wp,
                        int Ns, int *Nsize,
                        real_t *e, real_t *f,
                        real_t *dpinv, real_t *dm, real_t alpha)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;
  int Nst_pad = ceil_nwp(Nst);

  real_t* vp_dev = (real_t*)dev_ptr(vp);
  real_t* wp_dev = (real_t*)dev_ptr(wp);

  int blockSize  = 2 * VECTOR_LENGTH;
  int gridSize   = (Nst_pad * NVC)/blockSize;
  //int nSM        = getNumSMs();

  mult_domainwall_5din_ee_LUdaginv_dirac_kernel<<<gridSize, blockSize>>>(
    vp_dev, wp_dev,
    Ns, Nst_pad, alpha);

  CHECK(cudaDeviceSynchronize());
}

#endif
//============================================================END=====

