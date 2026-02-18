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

  const real_t* e_con     = ConstantMemoryTraits<real_t>::e();
  const real_t* f_con     = ConstantMemoryTraits<real_t>::f();
  const real_t* dpinv_con = ConstantMemoryTraits<real_t>::dpinv();
  const real_t* dm_con    = ConstantMemoryTraits<real_t>::dm(); 
  
  for( int idx = ist; idx < Nst_pad * NVC; idx += GridSize){

    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int ivc     = idx2_wp % NVC;
    int idx_out = idx2_wp / NVC;
    // int site    = idx_in  + NWP * idx_out; // Removed as we use base_idx

    const int base_site_idx = idx_in + NWP * (ivc + Nin5 * idx_out);
    const int stride_id = NWP * NVC;
    const int stride_is = NWP * NVCD;

    real_t vt[ND], yt[ND], xt[ND];

    if (sizeof(real_t) == 4) {
      int cur_idx_base = base_site_idx; // is=0
      vt[0] = wp[cur_idx_base + stride_id * 0];
      vt[1] = wp[cur_idx_base + stride_id * 1];
      vt[2] = wp[cur_idx_base + stride_id * 2];
      vt[3] = wp[cur_idx_base + stride_id * 3];

      vp[cur_idx_base + stride_id * 0] = vt[0];
      vp[cur_idx_base + stride_id * 1] = vt[1];
      vp[cur_idx_base + stride_id * 2] = vt[2];
      vp[cur_idx_base + stride_id * 3] = vt[3];

      real_t e0 = e_con[0];
      yt[0] = e0 * vt[0]; yt[1] = e0 * vt[1]; yt[2] = e0 * vt[2]; yt[3] = e0 * vt[3];
    } else {
      for(int id = 0; id < ND; ++id){
        int cur_idx = base_site_idx + stride_id * id;
        vt[id] = wp[cur_idx];
        vp[cur_idx] = vt[id];
        yt[id] = e_con[0] * vt[id];
      }
    }

    if (sizeof(real_t) == 4) {
      for (int is = 1; is < Ns-1; ++is) {
        real_t a = real_t(0.5) * dm_con[is] * dpinv_con[is-1];
        real_t eis = e_con[is];
        xt[0] = vt[0]; xt[1] = vt[1]; xt[2] = vt[2]; xt[3] = vt[3];
        int is_offset = stride_is * is;
        vt[0] = wp[base_site_idx + stride_id * 0 + is_offset];
        vt[1] = wp[base_site_idx + stride_id * 1 + is_offset];
        vt[2] = wp[base_site_idx + stride_id * 2 + is_offset];
        vt[3] = wp[base_site_idx + stride_id * 3 + is_offset];
        vt[0] += a * (xt[0] + xt[2]); vt[1] += a * (xt[1] + xt[3]); vt[2] += a * (xt[2] + xt[0]); vt[3] += a * (xt[3] + xt[1]);
        vp[base_site_idx + stride_id * 0 + is_offset] = vt[0];
        vp[base_site_idx + stride_id * 1 + is_offset] = vt[1];
        vp[base_site_idx + stride_id * 2 + is_offset] = vt[2];
        vp[base_site_idx + stride_id * 3 + is_offset] = vt[3];
        yt[0] += eis * vt[0]; yt[1] += eis * vt[1]; yt[2] += eis * vt[2]; yt[3] += eis * vt[3];
      }
    } else {
      for (int is = 1; is < Ns-1; ++is) {
        real_t a = real_t(0.5) * dm_con[is] * dpinv_con[is-1];
        real_t eis = e_con[is];
        for(int id=0; id<ND; ++id) xt[id] = vt[id];
        for(int id = 0; id < ND; ++id) vt[id] = wp[base_site_idx + stride_id * id + stride_is * is];
        vt[0] += a * (xt[0] + xt[2]); vt[1] += a * (xt[1] + xt[3]); vt[2] += a * (xt[2] + xt[0]); vt[3] += a * (xt[3] + xt[1]);
        for(int id = 0; id < ND; ++id){
          int cur_idx = base_site_idx + stride_id * id + stride_is * is;
          vp[cur_idx] = vt[id];
          yt[id] += eis * vt[id];
        }
      }
    }

    if (sizeof(real_t) == 4) {
      int is = Ns-1;
      xt[0] = vt[0]; xt[1] = vt[1]; xt[2] = vt[2]; xt[3] = vt[3];
      int last_is_offset = stride_is * is;
      vt[0] = wp[base_site_idx + stride_id * 0 + last_is_offset];
      vt[1] = wp[base_site_idx + stride_id * 1 + last_is_offset];
      vt[2] = wp[base_site_idx + stride_id * 2 + last_is_offset];
      vt[3] = wp[base_site_idx + stride_id * 3 + last_is_offset];
      real_t a_last = real_t(0.5) * dm_con[is] * dpinv_con[is-1];
      vt[0] += a_last * (xt[0] + xt[2]) - real_t(0.5) * (yt[0] - yt[2]);
      vt[1] += a_last * (xt[1] + xt[3]) - real_t(0.5) * (yt[1] - yt[3]);
      vt[2] += a_last * (xt[2] + xt[0]) - real_t(0.5) * (yt[2] - yt[0]);
      vt[3] += a_last * (xt[3] + xt[1]) - real_t(0.5) * (yt[3] - yt[1]);
      vp[base_site_idx + stride_id * 0 + last_is_offset] = vt[0];
      vp[base_site_idx + stride_id * 1 + last_is_offset] = vt[1];
      vp[base_site_idx + stride_id * 2 + last_is_offset] = vt[2];
      vp[base_site_idx + stride_id * 3 + last_is_offset] = vt[3];
    } else {
      int is = Ns-1;
      for(int id = 0; id < ND; ++id) xt[id] = vt[id];
      for(int id = 0; id < ND; ++id) vt[id] = wp[base_site_idx + stride_id * id + stride_is * is];
      real_t a_last = real_t(0.5) * dm_con[is] * dpinv_con[is-1];
      vt[0] += a_last * (xt[0] + xt[2]) - real_t(0.5) * (yt[0] - yt[2]);
      vt[1] += a_last * (xt[1] + xt[3]) - real_t(0.5) * (yt[1] - yt[3]);
      vt[2] += a_last * (xt[2] + xt[0]) - real_t(0.5) * (yt[2] - yt[0]);
      vt[3] += a_last * (xt[3] + xt[1]) - real_t(0.5) * (yt[3] - yt[1]);
      for(int id = 0; id < ND; ++id) vp[base_site_idx + stride_id * id + stride_is * is] = vt[id];
    }
    // L_inv completed

    if (sizeof(real_t) == 4) {
      int is = Ns-1;
      real_t aa = dpinv_con[Ns-1];
      real_t f1 = real_t(0.5) * (real_t(1.0) + alpha);
      real_t f2 = real_t(0.5) * (real_t(-1.0) + alpha);
      int last_is_offset = stride_is * is;
      vt[0] = aa * (f1 * vp[base_site_idx + stride_is * is + NWP * ID1] + f2 * vp[base_site_idx + stride_is * is + NWP * ID3]); 
      vt[1] = aa * (f1 * vp[base_site_idx + stride_is * is + NWP * ID2] + f2 * vp[base_site_idx + stride_is * is + NWP * ID4]); 
      vt[2] = aa * (f1 * vp[base_site_idx + stride_is * is + NWP * ID3] + f2 * vp[base_site_idx + stride_is * is + NWP * ID1]); 
      vt[3] = aa * (f1 * vp[base_site_idx + stride_is * is + NWP * ID4] + f2 * vp[base_site_idx + stride_is * is + NWP * ID2]); 
      vp[base_site_idx + last_is_offset + stride_id * 0] = vt[0];
      vp[base_site_idx + last_is_offset + stride_id * 1] = vt[1];
      vp[base_site_idx + last_is_offset + stride_id * 2] = vt[2];
      vp[base_site_idx + last_is_offset + stride_id * 3] = vt[3];
      yt[0] = real_t(0.5) * (vt[0] + vt[2]); yt[1] = real_t(0.5) * (vt[1] + vt[3]); yt[2] = real_t(0.5) * (vt[2] + vt[0]); yt[3] = real_t(0.5) * (vt[3] + vt[1]);
    } else {
      int is = Ns-1;
      real_t aa = dpinv_con[Ns-1];
      real_t f1 = real_t(0.5) * (real_t(1.0) + alpha);
      real_t f2 = real_t(0.5) * (real_t(-1.0) + alpha);
      real_t vt_tmp[ND];
      vt_tmp[0] = aa * (f1 * vp[base_site_idx + stride_is * is + NWP * ID1] + f2 * vp[base_site_idx + stride_is * is + NWP * ID3]); 
      vt_tmp[1] = aa * (f1 * vp[base_site_idx + stride_is * is + NWP * ID2] + f2 * vp[base_site_idx + stride_is * is + NWP * ID4]); 
      vt_tmp[2] = aa * (f1 * vp[base_site_idx + stride_is * is + NWP * ID3] + f2 * vp[base_site_idx + stride_is * is + NWP * ID1]); 
      vt_tmp[3] = aa * (f1 * vp[base_site_idx + stride_is * is + NWP * ID4] + f2 * vp[base_site_idx + stride_is * is + NWP * ID2]); 
      for(int id = 0; id < ND; ++id){
        vt[id] = vt_tmp[id];
        vp[base_site_idx + stride_id * id + stride_is * is] = vt[id];
      }
      yt[0] = real_t(0.5) * (vt[0] + vt[2]); yt[1] = real_t(0.5) * (vt[1] + vt[3]); yt[2] = real_t(0.5) * (vt[2] + vt[0]); yt[3] = real_t(0.5) * (vt[3] + vt[1]);
    }

    if (sizeof(real_t) == 4) {
      for (int is = Ns-2; is >= 0; --is) {
        real_t a = real_t(0.5) * dm_con[is];
        real_t f_is = f_con[is];
        real_t aa_is = dpinv_con[is];
        xt[0] = vt[0]; xt[1] = vt[1]; xt[2] = vt[2]; xt[3] = vt[3];
        int is_offset = stride_is * is;
        vt[0] = vp[base_site_idx + stride_id * 0 + is_offset];
        vt[1] = vp[base_site_idx + stride_id * 1 + is_offset];
        vt[2] = vp[base_site_idx + stride_id * 2 + is_offset];
        vt[3] = vp[base_site_idx + stride_id * 3 + is_offset];
        vt[0] += a * (xt[0] - xt[2]) - f_is * yt[0];
        vt[1] += a * (xt[1] - xt[3]) - f_is * yt[1];
        vt[2] += a * (xt[2] - xt[0]) - f_is * yt[2];
        vt[3] += a * (xt[3] - xt[1]) - f_is * yt[3];
        vt[0] *= aa_is; vt[1] *= aa_is; vt[2] *= aa_is; vt[3] *= aa_is;
        if(is == 0){
          real_t ff1 = real_t(0.5) * (real_t(1.0) + alpha);
          real_t ff2 = real_t(0.5) * (real_t(1.0) - alpha);
          vp[base_site_idx + is_offset + NWP * ID1] = ff1 * vt[0] + ff2 * vt[2];
          vp[base_site_idx + is_offset + NWP * ID2] = ff1 * vt[1] + ff2 * vt[3];
          vp[base_site_idx + is_offset + NWP * ID3] = ff1 * vt[2] + ff2 * vt[0];
          vp[base_site_idx + is_offset + NWP * ID4] = ff1 * vt[3] + ff2 * vt[1];
        }else{
          vp[base_site_idx + is_offset + stride_id * 0] = vt[0];
          vp[base_site_idx + is_offset + stride_id * 1] = vt[1];
          vp[base_site_idx + is_offset + stride_id * 2] = vt[2];
          vp[base_site_idx + is_offset + stride_id * 3] = vt[3];
        }
      }
    } else {
      for (int is = Ns-2; is >= 0; --is) {
        real_t a = real_t(0.5) * dm_con[is];
        real_t f_is = f_con[is];
        real_t aa_is = dpinv_con[is];
        for(int id=0; id<ND; ++id) xt[id] = vt[id];
        for(int id = 0; id < ND; ++id) vt[id] = vp[base_site_idx + stride_id * id + stride_is * is];
        vt[0] += a * (xt[0] - xt[2]) - f_is * yt[0];
        vt[1] += a * (xt[1] - xt[3]) - f_is * yt[1];
        vt[2] += a * (xt[2] - xt[0]) - f_is * yt[2];
        vt[3] += a * (xt[3] - xt[1]) - f_is * yt[3];
        for(int id=0; id<ND; ++id) vt[id] *= aa_is;
        if(is == 0){
          real_t ff1 = real_t(0.5) * (real_t(1.0) + alpha);
          real_t ff2 = real_t(0.5) * (real_t(1.0) - alpha);
          vp[base_site_idx + stride_is * is + NWP * ID1] = ff1 * vt[0] + ff2 * vt[2];
          vp[base_site_idx + stride_is * is + NWP * ID2] = ff1 * vt[1] + ff2 * vt[3];
          vp[base_site_idx + stride_is * is + NWP * ID3] = ff1 * vt[2] + ff2 * vt[0];
          vp[base_site_idx + stride_is * is + NWP * ID4] = ff1 * vt[3] + ff2 * vt[1];
        }else{
          for(int id = 0; id < ND; ++id) vp[base_site_idx + stride_id * id + stride_is * is] = vt[id];
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

  const real_t* e_con     = ConstantMemoryTraits<real_t>::e();
  const real_t* f_con     = ConstantMemoryTraits<real_t>::f();
  const real_t* dpinv_con = ConstantMemoryTraits<real_t>::dpinv();
  const real_t* dm_con    = ConstantMemoryTraits<real_t>::dm(); 
  
  for( int idx = ist; idx < Nst_pad * NVC; idx += GridSize){

    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int ivc     = idx2_wp % NVC;
    int idx_out = idx2_wp / NVC;
    // int site    = idx_in + NWP * idx_out;

    const int base_site_idx = idx_in + NWP * (ivc + Nin5 * idx_out);
    const int stride_id = NWP * NVC;
    const int stride_is = NWP * NVCD;

    real_t vt[ND], yt[ND], xt[ND];

    int is = 0;

    if (sizeof(real_t) == 4) {
      real_t a0 = dpinv_con[0];
      real_t f1_0 = real_t(0.5) * (real_t(1.0) + alpha);
      real_t f2_0 = real_t(0.5) * (real_t(1.0) - alpha);
      vt[0] = a0 * (f1_0 * wp[base_site_idx + stride_is * is + NWP * ID1] + f2_0 * wp[base_site_idx + stride_is * is + NWP * ID3]);
      vt[1] = a0 * (f1_0 * wp[base_site_idx + stride_is * is + NWP * ID2] + f2_0 * wp[base_site_idx + stride_is * is + NWP * ID4]);
      vt[2] = a0 * (f1_0 * wp[base_site_idx + stride_is * is + NWP * ID3] + f2_0 * wp[base_site_idx + stride_is * is + NWP * ID1]);
      vt[3] = a0 * (f1_0 * wp[base_site_idx + stride_is * is + NWP * ID4] + f2_0 * wp[base_site_idx + stride_is * is + NWP * ID2]);
    
      real_t f0 = f_con[0];
      int is_offset_0 = stride_is * is;
      vp[base_site_idx + stride_id * 0 + is_offset_0] = vt[0];
      vp[base_site_idx + stride_id * 1 + is_offset_0] = vt[1];
      vp[base_site_idx + stride_id * 2 + is_offset_0] = vt[2];
      vp[base_site_idx + stride_id * 3 + is_offset_0] = vt[3];
      yt[0] = f0 * vt[0]; yt[1] = f0 * vt[1]; yt[2] = f0 * vt[2]; yt[3] = f0 * vt[3];
    } else {
      real_t a0 = dpinv_con[0];
      real_t f1_0 = real_t(0.5) * (real_t(1.0) + alpha);
      real_t f2_0 = real_t(0.5) * (real_t(1.0) - alpha);
      real_t vt_tmp[ND];
      vt_tmp[0] = a0 * (f1_0 * wp[base_site_idx + stride_is * is + NWP * ID1] + f2_0 * wp[base_site_idx + stride_is * is + NWP * ID3]);
      vt_tmp[1] = a0 * (f1_0 * wp[base_site_idx + stride_is * is + NWP * ID2] + f2_0 * wp[base_site_idx + stride_is * is + NWP * ID4]);
      vt_tmp[2] = a0 * (f1_0 * wp[base_site_idx + stride_is * is + NWP * ID3] + f2_0 * wp[base_site_idx + stride_is * is + NWP * ID1]);
      vt_tmp[3] = a0 * (f1_0 * wp[base_site_idx + stride_is * is + NWP * ID4] + f2_0 * wp[base_site_idx + stride_is * is + NWP * ID2]);
      for(int id=0; id<ND; ++id) vt[id] = vt_tmp[id];
      real_t f0 = f_con[0];
      for(int id=0; id<ND; ++id){
        int cur_idx = base_site_idx + stride_id * id + stride_is * is;
        vp[cur_idx] = vt[id];
        yt[id] = f0 * vt[id];
      }
    }

    if (sizeof(real_t) == 4) {
      for (is = 1; is < Ns-1; ++is) {
        real_t a = real_t(0.5) * dm_con[is - 1];
        real_t aa = dpinv_con[is];
        real_t f_is = f_con[is];
        xt[0] = vt[0]; xt[1] = vt[1]; xt[2] = vt[2]; xt[3] = vt[3];
        int is_offset = stride_is * is;
        vt[0] = wp[base_site_idx + stride_id * 0 + is_offset];
        vt[1] = wp[base_site_idx + stride_id * 1 + is_offset];
        vt[2] = wp[base_site_idx + stride_id * 2 + is_offset];
        vt[3] = wp[base_site_idx + stride_id * 3 + is_offset];
        vt[0] += a * (xt[0] - xt[2]); vt[1] += a * (xt[1] - xt[3]); vt[2] += a * (xt[2] - xt[0]); vt[3] += a * (xt[3] - xt[1]);
        vt[0] *= aa; vt[1] *= aa; vt[2] *= aa; vt[3] *= aa;
        vp[base_site_idx + stride_id * 0 + is_offset] = vt[0];
        vp[base_site_idx + stride_id * 1 + is_offset] = vt[1];
        vp[base_site_idx + stride_id * 2 + is_offset] = vt[2];
        vp[base_site_idx + stride_id * 3 + is_offset] = vt[3];
        yt[0] += f_is * vt[0]; yt[1] += f_is * vt[1]; yt[2] += f_is * vt[2]; yt[3] += f_is * vt[3];
      }
    } else {
      for (is = 1; is < Ns-1; ++is) {
        real_t a = real_t(0.5) * dm_con[is - 1];
        real_t aa = dpinv_con[is];
        real_t f_is = f_con[is];
        for(int id=0; id<ND; ++id) xt[id] = vt[id];
        for(int id = 0; id < ND; ++id) vt[id] = wp[base_site_idx + stride_id * id + stride_is * is];
        vt[0] += a * (xt[0] - xt[2]); vt[1] += a * (xt[1] - xt[3]); vt[2] += a * (xt[2] - xt[0]); vt[3] += a * (xt[3] - xt[1]);
        for(int id = 0; id < ND; ++id){
          int cur_idx = base_site_idx + stride_id * id + stride_is * is;
          vt[id] *= aa;
          vp[cur_idx] = vt[id];
          yt[id] += f_is * vt[id];
        }
      }
    }

    if (sizeof(real_t) == 4) {
      is = Ns-1;
      xt[0] = vt[0]; xt[1] = vt[1]; xt[2] = vt[2]; xt[3] = vt[3];
      int last_is_offset = stride_is * is;
      vt[0] = wp[base_site_idx + stride_id * 0 + last_is_offset];
      vt[1] = wp[base_site_idx + stride_id * 1 + last_is_offset];
      vt[2] = wp[base_site_idx + stride_id * 2 + last_is_offset];
      vt[3] = wp[base_site_idx + stride_id * 3 + last_is_offset];
      real_t a_last = real_t(0.5) * dm_con[is - 1];
      vt[0] += a_last * (xt[0] - xt[2]) - real_t(0.5) * (yt[0] + yt[2]);
      vt[1] += a_last * (xt[1] - xt[3]) - real_t(0.5) * (yt[1] + yt[3]);
      vt[2] += a_last * (xt[2] - xt[0]) - real_t(0.5) * (yt[2] + yt[0]);
      vt[3] += a_last * (xt[3] - xt[1]) - real_t(0.5) * (yt[3] + yt[1]);
      real_t aa_last = dpinv_con[is];
      vt[0] *= aa_last; vt[1] *= aa_last; vt[2] *= aa_last; vt[3] *= aa_last;
      real_t ff1 = real_t(0.5) * (real_t(1.0) + alpha);
      real_t ff2 = real_t(0.5) * (real_t(-1.0) + alpha);
      vp[base_site_idx + stride_is * is + NWP * ID1] = ff1 * vt[0] + ff2 * vt[2];
      vp[base_site_idx + stride_is * is + NWP * ID2] = ff1 * vt[1] + ff2 * vt[3];
      vp[base_site_idx + stride_is * is + NWP * ID3] = ff1 * vt[2] + ff2 * vt[0];
      vp[base_site_idx + stride_is * is + NWP * ID4] = ff1 * vt[3] + ff2 * vt[1];
    } else {
      is = Ns-1;
      for(int id=0; id<ND; ++id) xt[id] = vt[id];
      for(int id = 0; id < ND; ++id) vt[id] = wp[base_site_idx + stride_id * id + stride_is * is];
      real_t a_last = real_t(0.5) * dm_con[is - 1];
      vt[0] += a_last * (xt[0] - xt[2]) - real_t(0.5) * (yt[0] + yt[2]);
      vt[1] += a_last * (xt[1] - xt[3]) - real_t(0.5) * (yt[1] + yt[3]);
      vt[2] += a_last * (xt[2] - xt[0]) - real_t(0.5) * (yt[2] + yt[0]);
      vt[3] += a_last * (xt[3] - xt[1]) - real_t(0.5) * (yt[3] + yt[1]);
      real_t aa_last = dpinv_con[is];
      for(int id=0; id<ND; ++id) vt[id] *= aa_last;
      real_t ff1 = real_t(0.5) * (real_t(1.0) + alpha);
      real_t ff2 = real_t(0.5) * (real_t(-1.0) + alpha);
      vp[base_site_idx + stride_is * is + NWP * ID1] = ff1 * vt[0] + ff2 * vt[2];
      vp[base_site_idx + stride_is * is + NWP * ID2] = ff1 * vt[1] + ff2 * vt[3];
      vp[base_site_idx + stride_is * is + NWP * ID3] = ff1 * vt[2] + ff2 * vt[0];
      vp[base_site_idx + stride_is * is + NWP * ID4] = ff1 * vt[3] + ff2 * vt[1];
    }

    // Udag_inv completed

    if (sizeof(real_t) == 4) {
      is = Ns-1;
      yt[0] = real_t(0.5) * (vt[0] - vt[2]); yt[1] = real_t(0.5) * (vt[1] - vt[3]); yt[2] = real_t(0.5) * (vt[2] - vt[0]); yt[3] = real_t(0.5) * (vt[3] - vt[1]);
      for (is = Ns-2; is >= 0; --is) {
        real_t a = real_t(0.5) * dm_con[is + 1] * dpinv_con[is];
        real_t e_is = e_con[is];
        xt[0] = vt[0]; xt[1] = vt[1]; xt[2] = vt[2]; xt[3] = vt[3];
        int is_offset = stride_is * is;
        vt[0] = vp[base_site_idx + stride_id * 0 + is_offset];
        vt[1] = vp[base_site_idx + stride_id * 1 + is_offset];
        vt[2] = vp[base_site_idx + stride_id * 2 + is_offset];
        vt[3] = vp[base_site_idx + stride_id * 3 + is_offset];
        vt[0] += a * (xt[0] + xt[2]) - e_is * yt[0];
        vt[1] += a * (xt[1] + xt[3]) - e_is * yt[1];
        vt[2] += a * (xt[2] + xt[0]) - e_is * yt[2];
        vt[3] += a * (xt[3] + xt[1]) - e_is * yt[3];
        vp[base_site_idx + stride_id * 0 + is_offset] = vt[0];
        vp[base_site_idx + stride_id * 1 + is_offset] = vt[1];
        vp[base_site_idx + stride_id * 2 + is_offset] = vt[2];
        vp[base_site_idx + stride_id * 3 + is_offset] = vt[3];
      }
    } else {
      is = Ns-1;
      yt[0] = real_t(0.5) * (vt[0] - vt[2]); yt[1] = real_t(0.5) * (vt[1] - vt[3]); yt[2] = real_t(0.5) * (vt[2] - vt[0]); yt[3] = real_t(0.5) * (vt[3] - vt[1]);
      for (is = Ns-2; is >= 0; --is) {
        real_t a = real_t(0.5) * dm_con[is + 1] * dpinv_con[is];
        real_t e_is = e_con[is];
        for(int id=0; id<ND; ++id) xt[id] = vt[id];
        for(int id = 0; id < ND; ++id) vt[id] = vp[base_site_idx + stride_id * id + stride_is * is];
        vt[0] += a * (xt[0] + xt[2]) - e_is * yt[0];
        vt[1] += a * (xt[1] + xt[3]) - e_is * yt[1];
        vt[2] += a * (xt[2] + xt[0]) - e_is * yt[2];
        vt[3] += a * (xt[3] + xt[1]) - e_is * yt[3];
        for(int id = 0; id < ND; ++id) vp[base_site_idx + stride_id * id + stride_is * is] = vt[id];
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

