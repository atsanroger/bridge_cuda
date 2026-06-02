/*!
      @file    aindex_eo_alt_cuda-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-08-20 14:25:12 #$
      @version $LastChangedRevision: 2535 $
*/

//====================================================================
__global__ void split_kernel(real_t* ve, real_t* vo, real_t* w,
                             int ieo_origin, int nin,
                             int Nx, int Ny, int Nz, int Nt)
{
  int ist = blockIdx.x * blockDim.x + threadIdx.x;

  int ix = ist % Nx;
  int iy = (ist/Nx) % Ny;
  int iz = (ist/(Nx*Ny)) % Nz;
  int it = ist/(Nx*Ny*Nz);
  int ieo = (ix + iy + iz + it + ieo_origin) % 2;
  int ist2 = ist/2;
  if(ieo == 0){
    for(int in = 0; in < nin; ++in){
      ve[IDX2(nin, in, ist2)] = w[IDX2(nin, in, ist)];
    }
  }else{
    for(int in = 0; in < nin; ++in){
      vo[IDX2(nin, in, ist2)] = w[IDX2(nin, in, ist)];
    }
  }

}

//====================================================================
void split(real_t* ve, real_t* vo, real_t* w, int ieo_origin,
           int nin, int* Nsize, int nvol2_pad, int nvol_pad)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  real_t* ve_dev  = (real_t*)dev_ptr(ve);
  real_t* vo_dev  = (real_t*)dev_ptr(vo);
  real_t* w_dev   = (real_t*)dev_ptr(w);

  int nth = VECTOR_LENGTH;
  int nbl = Nst/nth;

  split_kernel<<<nbl,nth>>>(ve_dev, vo_dev, w_dev,
                            ieo_origin, nin, Nx, Ny, Nz, Nt);

  CHECK(cudaDeviceSynchronize());

}

//====================================================================
__global__ void merge_kernel(real_t* v, real_t* we, real_t* wo,
                             int ieo_origin, int nin,
                             int Nx, int Ny, int Nz, int Nt)
{
  int ist = blockIdx.x * blockDim.x + threadIdx.x;

  int ix = ist % Nx;
  int iy = (ist/Nx) % Ny;
  int iz = (ist/(Nx*Ny)) % Nz;
  int it = ist/(Nx*Ny*Nz);
  int ieo = (ix + iy + iz + it + ieo_origin) % 2;
  int ist2 = ist/2;
  if(ieo == 0){
    for(int in = 0; in < nin; ++in){
      v[IDX2(nin, in, ist)] = we[IDX2(nin, in, ist2)];
    }
  }else{
    for(int in = 0; in < nin; ++in){
      v[IDX2(nin, in, ist)] = wo[IDX2(nin, in, ist2)];
    }
  }

}

//====================================================================
void merge(real_t* v, real_t* we, real_t* wo, int ieo_origin,
           int nin, int* Nsize, int nvol2_pad, int nvol_pad)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  real_t* v_dev   = (real_t*)dev_ptr(v);
  real_t* we_dev  = (real_t*)dev_ptr(we);
  real_t* wo_dev  = (real_t*)dev_ptr(wo);

  int nth = VECTOR_LENGTH;
  int nbl = Nst/nth;

  merge_kernel<<<nbl,nth>>>(v_dev, we_dev, wo_dev,
                            ieo_origin, nin, Nx, Ny, Nz, Nt);

  CHECK(cudaDeviceSynchronize());

}

//====================================================================
// QDW-aware split/merge: the QDW field is double4-contiguous
// (4*IDX2(nin/4, c, site)+k), NOT the real NWP-strided IDX2(nin,...) layout.
// Move whole double4 units to preserve the layout.
//====================================================================
__global__ void split_kernel_qdw(real_t* ve, real_t* vo, real_t* w,
                                 int ieo_origin, int nin,
                                 int Nx, int Ny, int Nz, int Nt)
{
  int ist = blockIdx.x * blockDim.x + threadIdx.x;
  int ix = ist % Nx;
  int iy = (ist/Nx) % Ny;
  int iz = (ist/(Nx*Ny)) % Nz;
  int it = ist/(Nx*Ny*Nz);
  int ieo = (ix + iy + iz + it + ieo_origin) % 2;
  int ist2 = ist/2;
  int nin4 = nin/4;
  real_t* dst = (ieo == 0) ? ve : vo;
  for(int c = 0; c < nin4; ++c){
    int s = 4*IDX2(nin4, c, ist);
    int d = 4*IDX2(nin4, c, ist2);
    dst[d+0]=w[s+0]; dst[d+1]=w[s+1]; dst[d+2]=w[s+2]; dst[d+3]=w[s+3];
  }
}

void split_qdw(real_t* ve, real_t* vo, real_t* w, int ieo_origin,
               int nin, int* Nsize, int nvol2_pad, int nvol_pad)
{
  int Nx = Nsize[0], Ny = Nsize[1], Nz = Nsize[2], Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;
  real_t* ve_dev = (real_t*)dev_ptr(ve);
  real_t* vo_dev = (real_t*)dev_ptr(vo);
  real_t* w_dev  = (real_t*)dev_ptr(w);
  int nth = VECTOR_LENGTH;
  int nbl = Nst/nth;
  split_kernel_qdw<<<nbl,nth>>>(ve_dev, vo_dev, w_dev,
                                ieo_origin, nin, Nx, Ny, Nz, Nt);
  CHECK(cudaDeviceSynchronize());
}

__global__ void merge_kernel_qdw(real_t* v, real_t* we, real_t* wo,
                                 int ieo_origin, int nin,
                                 int Nx, int Ny, int Nz, int Nt)
{
  int ist = blockIdx.x * blockDim.x + threadIdx.x;
  int ix = ist % Nx;
  int iy = (ist/Nx) % Ny;
  int iz = (ist/(Nx*Ny)) % Nz;
  int it = ist/(Nx*Ny*Nz);
  int ieo = (ix + iy + iz + it + ieo_origin) % 2;
  int ist2 = ist/2;
  int nin4 = nin/4;
  real_t* src = (ieo == 0) ? we : wo;
  for(int c = 0; c < nin4; ++c){
    int dv = 4*IDX2(nin4, c, ist);
    int sh = 4*IDX2(nin4, c, ist2);
    v[dv+0]=src[sh+0]; v[dv+1]=src[sh+1]; v[dv+2]=src[sh+2]; v[dv+3]=src[sh+3];
  }
}

void merge_qdw(real_t* v, real_t* we, real_t* wo, int ieo_origin,
               int nin, int* Nsize, int nvol2_pad, int nvol_pad)
{
  int Nx = Nsize[0], Ny = Nsize[1], Nz = Nsize[2], Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;
  real_t* v_dev  = (real_t*)dev_ptr(v);
  real_t* we_dev = (real_t*)dev_ptr(we);
  real_t* wo_dev = (real_t*)dev_ptr(wo);
  int nth = VECTOR_LENGTH;
  int nbl = Nst/nth;
  merge_kernel_qdw<<<nbl,nth>>>(v_dev, we_dev, wo_dev,
                                ieo_origin, nin, Nx, Ny, Nz, Nt);
  CHECK(cudaDeviceSynchronize());
}

//====================================================================
// QTW-aware split/merge: the QTW field is 6-reals-per-cplx-contiguous
// (6*IDX2(nin/6, c, site)+k), parallel to QDW's 4-real packing.
//====================================================================
__global__ void split_kernel_qtw(real_t* ve, real_t* vo, real_t* w,
                                 int ieo_origin, int nin,
                                 int Nx, int Ny, int Nz, int Nt)
{
  int ist = blockIdx.x * blockDim.x + threadIdx.x;
  int ix = ist % Nx;
  int iy = (ist/Nx) % Ny;
  int iz = (ist/(Nx*Ny)) % Nz;
  int it = ist/(Nx*Ny*Nz);
  int ieo = (ix + iy + iz + it + ieo_origin) % 2;
  int ist2 = ist/2;
  int nin6 = nin/6;
  real_t* dst = (ieo == 0) ? ve : vo;
  for(int c = 0; c < nin6; ++c){
    int s = 6*IDX2(nin6, c, ist);
    int d = 6*IDX2(nin6, c, ist2);
    dst[d+0]=w[s+0]; dst[d+1]=w[s+1]; dst[d+2]=w[s+2];
    dst[d+3]=w[s+3]; dst[d+4]=w[s+4]; dst[d+5]=w[s+5];
  }
}

void split_qtw(real_t* ve, real_t* vo, real_t* w, int ieo_origin,
               int nin, int* Nsize, int nvol2_pad, int nvol_pad)
{
  int Nx = Nsize[0], Ny = Nsize[1], Nz = Nsize[2], Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;
  real_t* ve_dev = (real_t*)dev_ptr(ve);
  real_t* vo_dev = (real_t*)dev_ptr(vo);
  real_t* w_dev  = (real_t*)dev_ptr(w);
  int nth = VECTOR_LENGTH;
  int nbl = Nst/nth;
  split_kernel_qtw<<<nbl,nth>>>(ve_dev, vo_dev, w_dev,
                                ieo_origin, nin, Nx, Ny, Nz, Nt);
  CHECK(cudaDeviceSynchronize());
}

__global__ void merge_kernel_qtw(real_t* v, real_t* we, real_t* wo,
                                 int ieo_origin, int nin,
                                 int Nx, int Ny, int Nz, int Nt)
{
  int ist = blockIdx.x * blockDim.x + threadIdx.x;
  int ix = ist % Nx;
  int iy = (ist/Nx) % Ny;
  int iz = (ist/(Nx*Ny)) % Nz;
  int it = ist/(Nx*Ny*Nz);
  int ieo = (ix + iy + iz + it + ieo_origin) % 2;
  int ist2 = ist/2;
  int nin6 = nin/6;
  real_t* src = (ieo == 0) ? we : wo;
  for(int c = 0; c < nin6; ++c){
    int dv = 6*IDX2(nin6, c, ist);
    int sh = 6*IDX2(nin6, c, ist2);
    v[dv+0]=src[sh+0]; v[dv+1]=src[sh+1]; v[dv+2]=src[sh+2];
    v[dv+3]=src[sh+3]; v[dv+4]=src[sh+4]; v[dv+5]=src[sh+5];
  }
}

void merge_qtw(real_t* v, real_t* we, real_t* wo, int ieo_origin,
               int nin, int* Nsize, int nvol2_pad, int nvol_pad)
{
  int Nx = Nsize[0], Ny = Nsize[1], Nz = Nsize[2], Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;
  real_t* v_dev  = (real_t*)dev_ptr(v);
  real_t* we_dev = (real_t*)dev_ptr(we);
  real_t* wo_dev = (real_t*)dev_ptr(wo);
  int nth = VECTOR_LENGTH;
  int nbl = Nst/nth;
  merge_kernel_qtw<<<nbl,nth>>>(v_dev, we_dev, wo_dev,
                                ieo_origin, nin, Nx, Ny, Nz, Nt);
  CHECK(cudaDeviceSynchronize());
}

//============================================================END=====
