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

//============================================================END=====
