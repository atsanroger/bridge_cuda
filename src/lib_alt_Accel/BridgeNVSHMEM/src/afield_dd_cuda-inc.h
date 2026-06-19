/*!
      @file    afield_dd_cuda-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2024-04-01 12:13:09 #$
      @version $LastChangedRevision: 2595 $
*/

//====================================================================
__global__
void block_axpy_eo_kernel(real_t* vp, real_t* out, real_t* wp, 
                          real_t fac, int Nin, int Nex,
                          int Nx, int Ny, int Nz, int Nt,
                          int Bx, int By, int Bz, int Bt, int ieo)
{
  // In this function, out is an array of size (2 * Nblock),
  // assuming Nblock complex values.

  int Nst  = Nx * Ny * Nz * Nt;

  // size of block
  int Bsize = Bx * By * Bz * Bt;

  // numbers of blocks
  int NBx  = Nx/Bx;
  int NBy  = Ny/By;
  int NBz  = Nz/Bz;
  //int NBt  = Nt/Bt;  //(not used)

  int site2 = blockIdx.x * blockDim.x + threadIdx.x;

  int bsite = site2 % Bsize;
  int block = site2 / Bsize;

  int ibx = block % NBx;
  int iby = (block/NBx) % NBy;
  int ibz = (block/(NBx * NBy)) % NBz;
  int ibt = block/(NBx * NBy * NBz);

  int ieo_skip = 1 - ieo;
  int jeo = (ibx + iby + ibz + ibt) % 2;
  if(jeo == ieo_skip) return;

  int kx   = bsite % Bx;
  int ix   = kx + Bx * ibx;
  int kyzt = bsite/Bx;
  int ky   = kyzt % By;
  int iy   = ky + By * iby;
  int kzt  = kyzt/By;
  int kz   = kzt % Bz;
  int iz   = kz + Bz * ibz;
  int kt   = kzt/Bz;
  int it   = kt + Bt * ibt;
  int site = ix + Nx * (iy + Ny * (iz + Nz * it));

  real_t ar = out[2*block]   * fac;
  real_t ai = out[2*block+1] * fac;

  for(int ex = 0; ex < Nex; ++ex){
    int nv1 = Nin * Nst * ex;
    int nv2 = Nin * Nst * ex;
    for(int in2 = 0; in2 < Nin/2; ++in2){
      int inr = 2 * in2;
      int ini = 2 * in2 + 1;
      real_t vtr = vp[nv1 + IDX2(Nin, inr, site)];
      real_t vti = vp[nv1 + IDX2(Nin, ini, site)];
      real_t wtr = wp[nv2 + IDX2(Nin, inr, site)];
      real_t wti = wp[nv2 + IDX2(Nin, ini, site)];
      vtr += ar * wtr - ai * wti;
      vti += ar * wti + ai * wtr;
      vp[nv1 + IDX2(Nin, inr, site)] = vtr;
      vp[nv1 + IDX2(Nin, ini, site)] = vti;
    }
  }

}

//====================================================================
void block_axpy_eo(real_t* vp, real_t* out, real_t* wp,
                   real_t fac, int Nin, int Nex,
                   int *Nsize, int *block_size, int ieo)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst  = Nx * Ny * Nz * Nt;

  // size of block
  int Bx = block_size[0];
  int By = block_size[1];
  int Bz = block_size[2];
  int Bt = block_size[3];

  real_t* vp_dev  = (real_t*)dev_ptr(vp);
  real_t* out_dev = (real_t*)dev_ptr(out);
  real_t* wp_dev  = (real_t*)dev_ptr(wp);

  int nth = VECTOR_LENGTH;
  int nbl = Nst/nth;

  block_axpy_eo_kernel<<<nbl,nth>>>(vp_dev, out_dev, wp_dev,
                                    fac, Nin, Nex,
                                    Nx, Ny, Nz, Nt,
                                    Bx, By, Bz, Bt, ieo);

  CHECK(cudaDeviceSynchronize());

}

//====================================================================
__global__
void block_scal_eo_kernel(real_t* vp, real_t* out,
                          int Nin, int Nex,
                          int Nx, int Ny, int Nz, int Nt,
                          int Bx, int By, int Bz, int Bt, int ieo)
{
  // In this function, out is an array of size (2 * Nblock).

  int Nst  = Nx * Ny * Nz * Nt;

  // size of block
  int Bsize = Bx * By * Bz * Bt;

  // numbers of blocks
  int NBx  = Nx/Bx;
  int NBy  = Ny/By;
  int NBz  = Nz/Bz;
  // int NBt  = Nt/Bt;  //(not used)

  int site2 = blockIdx.x * blockDim.x + threadIdx.x;

  int bsite = site2 % Bsize;
  int block = site2 / Bsize;

  int ibx = block % NBx;
  int iby = (block/NBx) % NBy;
  int ibz = (block/(NBx * NBy)) % NBz;
  int ibt = block/(NBx * NBy * NBz);

  int ieo_skip = 1 - ieo;
  int jeo = (ibx + iby + ibz + ibt) % 2;
  if(jeo == ieo_skip) return;

  int kx   = bsite % Bx;
  int ix   = kx + Bx * ibx;
  int kyzt = bsite/Bx;
  int ky   = kyzt % By;
  int iy   = ky + By * iby;
  int kzt  = kyzt/By;
  int kz   = kzt % Bz;
  int iz   = kz + Bz * ibz;
  int kt   = kzt/Bz;
  int it   = kt + Bt * ibt;
  int site = ix + Nx * (iy + Ny * (iz + Nz * it));

  real_t ar = out[2*block];
  real_t ai = out[2*block+1];

  for(int ex = 0; ex < Nex; ++ex){
    int nv1 = Nin * Nst * ex;
    for(int in2 = 0; in2 < Nin/2; ++in2){
      int inr = 2 * in2;
      int ini = 2 * in2 + 1;
      real_t vtr = vp[nv1 + IDX2(Nin, inr, site)];
      real_t vti = vp[nv1 + IDX2(Nin, ini, site)];
      vp[nv1 + IDX2(Nin, inr, site)] = ar * vtr - ai * vti;
      vp[nv1 + IDX2(Nin, ini, site)] = ar * vti + ai * vtr;
    }
  }

}

//====================================================================
void block_scal_eo(real_t* vp, real_t* out, int Nin, int Nex,
                   int *Nsize, int *block_size, int ieo)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst  = Nx * Ny * Nz * Nt;

  // size of block
  int Bx = block_size[0];
  int By = block_size[1];
  int Bz = block_size[2];
  int Bt = block_size[3];

  real_t* vp_dev  = (real_t*)dev_ptr(vp);
  real_t* out_dev = (real_t*)dev_ptr(out);

  int nth = VECTOR_LENGTH;
  int nbl = Nst/nth;

  block_scal_eo_kernel<<<nbl,nth>>>(vp_dev, out_dev, Nin, Nex,
                                    Nx, Ny, Nz, Nt,
                                    Bx, By, Bz, Bt, ieo);

  CHECK(cudaDeviceSynchronize());

}

//====================================================================
__global__
void block_dotc_eo_kernel(real_t* out, real_t* vp, real_t* wp,
                          int Nin, int Nex,
                          int Nx, int Ny, int Nz, int Nt,
                          int Bx, int By, int Bz, int Bt, int ieo)
{
  // In this function, out is an array of size (2 * Nst).

  int Nst  = Nx * Ny * Nz * Nt;

  // size of block
  int Bsize = Bx * By * Bz * Bt;

  int NBx  = Nx/Bx;
  int NBy  = Ny/By;
  int NBz  = Nz/Bz;
  // int NBt  = Nt/Bt;  //(not used)

  int site2 = blockIdx.x * blockDim.x + threadIdx.x;

  int bsite = site2 % Bsize;
  int block = site2/Bsize;

  int ibx = block % NBx;
  int iby = (block/NBx) % NBy;
  int ibz = (block/(NBx * NBy)) % NBz;
  int ibt = block/(NBx * NBy * NBz);

  int ieo_skip = 1 - ieo;
  int jeo = (ibx + iby + ibz + ibt) % 2;
  if(jeo == ieo_skip) return;

  int kx   = bsite % Bx;
  int ix   = kx + Bx * ibx;
  int kyzt = bsite/Bx;
  int ky   = kyzt % By;
  int iy   = ky + By * iby;
  int kzt  = kyzt/By;
  int kz   = kzt % Bz;
  int iz   = kz + Bz * ibz;
  int kt   = kzt/Bz;
  int it   = kt + Bt * ibt;
  int site = ix + Nx * (iy + Ny * (iz + Nz * it));

  real_t xtr = 0.0;
  real_t xti = 0.0;

  for(int ex = 0; ex < Nex; ++ex){
    int nv1 = Nin * Nst * ex;
    int nv2 = Nin * Nst * ex;
    for(int in2 = 0; in2 < Nin/2; ++in2){
      int inr = 2 * in2;
      int ini = 2 * in2 + 1;
      real_t vtr = vp[nv1 + IDX2(Nin, inr, site)];
      real_t vti = vp[nv1 + IDX2(Nin, ini, site)];
      real_t wtr = wp[nv2 + IDX2(Nin, inr, site)];
      real_t wti = wp[nv2 + IDX2(Nin, ini, site)];
      xtr += vtr * wtr + vti * wti;
      xti += vtr * wti - vti * wtr;
    }
  }

  out[2 * site2]     = xtr;
  out[2 * site2 + 1] = xti;

}

//====================================================================
void block_dotc_eo(real_t* out, real_t* vp, real_t* wp,
                   int Nin, int Nex,
                   int *Nsize, int *block_size, int ieo)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst  = Nx * Ny * Nz * Nt;

  // size of block
  int Bx = block_size[0];
  int By = block_size[1];
  int Bz = block_size[2];
  int Bt = block_size[3];

  real_t* out_dev = (real_t*)dev_ptr(out);
  real_t* vp_dev  = (real_t*)dev_ptr(vp);
  real_t* wp_dev  = (real_t*)dev_ptr(wp);

  int nth = VECTOR_LENGTH;
  int nbl = Nst/nth;

  block_dotc_eo_kernel<<<nbl,nth>>>(out_dev, vp_dev, wp_dev,
                                    Nin, Nex,
                                    Nx, Ny, Nz, Nt,
                                    Bx, By, Bz, Bt, ieo);

  CHECK(cudaDeviceSynchronize());

}

//====================================================================
__global__ void block_norm2_eo_kernel(
                             real_t* out, real_t* vp,
                             int Nin, int Nex,
                             int Nx, int Ny, int Nz, int Nt,
                             int Bx, int By, int Bz, int Bt, int ieo)
{
  // In this function, out is an array of size (2 * Nst).

  int Nst  = Nx * Ny * Nz * Nt;

  // size of block
  int Bsize = Bx * By * Bz * Bt;

  int NBx  = Nx/Bx;
  int NBy  = Ny/By;
  int NBz  = Nz/Bz;
  // int NBt  = Nt/Bt;  //(not used)

  int site2 = blockIdx.x * blockDim.x + threadIdx.x;

  int bsite = site2 % Bsize;
  int block = site2/Bsize;

  int ibx = block % NBx;
  int iby = (block/NBx) % NBy;
  int ibz = (block/(NBx * NBy)) % NBz;
  int ibt = block/(NBx * NBy * NBz);

  int ieo_skip = 1 - ieo;
  int jeo = (ibx + iby + ibz + ibt) % 2;
  if(jeo == ieo_skip) return;

  int kx   = bsite % Bx;
  int ix   = kx + Bx * ibx;
  int kyzt = bsite/Bx;
  int ky   = kyzt % By;
  int iy   = ky + By * iby;
  int kzt  = kyzt/By;
  int kz   = kzt % Bz;
  int iz   = kz + Bz * ibz;
  int kt   = kzt/Bz;
  int it   = kt + Bt * ibt;
  int site = ix + Nx * (iy + Ny * (iz + Nz * it));

  real_t xtr = 0.0;

  for(int ex = 0; ex < Nex; ++ex){
    int nv1 = Nin * Nst * ex;
    for(int in2 = 0; in2 < Nin/2; ++in2){
      int inr = 2 * in2;
      int ini = 2 * in2 + 1;
      real_t vtr = vp[nv1 + IDX2(Nin, inr, site)];
      real_t vti = vp[nv1 + IDX2(Nin, ini, site)];
      xtr += vtr * vtr + vti * vti;
    }
  }

  out[2 * site2]     = xtr;
  out[2 * site2 + 1] = 0.0;

}

//====================================================================
void block_norm2_eo(real_t* out, real_t* vp,
                    int Nin, int Nex,
                    int *Nsize, int *block_size, int ieo)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst  = Nx * Ny * Nz * Nt;

  // size of block
  int Bx = block_size[0];
  int By = block_size[1];
  int Bz = block_size[2];
  int Bt = block_size[3];

  real_t* out_dev = (real_t*)dev_ptr(out);
  real_t* vp_dev  = (real_t*)dev_ptr(vp);

  int nth = VECTOR_LENGTH;
  int nbl = Nst/nth;

  block_norm2_eo_kernel<<<nbl,nth>>>(out_dev, vp_dev,
                                     Nin, Nex,
                                     Nx, Ny, Nz, Nt,
                                     Bx, By, Bz, Bt, ieo);
 
  CHECK(cudaDeviceSynchronize());

}

//====================================================================
__global__
void reduce_block_kernel(real_t* out, real_t* out_red,
                         int Nx, int Ny, int Nz, int Nt,
                         int Bx, int By, int Bz, int Bt, int ieo)
{
  // size of fine lattice
  int Nst  = Nx * Ny * Nz * Nt;

  // size of block
  int Bsize = Bx * By * Bz * Bt;

  // number of blocks
  int NBx  = Nx/Bx;
  int NBy  = Ny/By;
  int NBz  = Nz/Bz;
  int NBt  = Nt/Bt;
  int Nblock = NBx * NBy * NBz * NBt;

  int block = blockIdx.x * blockDim.x + threadIdx.x;

  int ibx = block % NBx;
  int iby = (block/NBx) % NBy;
  int ibz = (block/(NBx * NBy)) % NBz;
  int ibt = block/(NBx * NBy * NBz);

  int ieo_skip = 1 - ieo;
  int jeo = (ibx + iby + ibz + ibt) % 2;
  if(jeo == ieo_skip){
    out[2*block]   = 0.0;
    out[2*block+1] = 0.0;
  }else{
    real_t vtr = 0.0;
    real_t vti = 0.0;
    for(int bsite = 0; bsite < Bsize; ++ bsite){
      real_t *vp = &out_red[2 * Bsize * block];
      vtr += vp[2 * bsite];
      vti += vp[2 * bsite +1];
    }
    out[2*block]   = vtr;
    out[2*block+1] = vti;
  }

}

//====================================================================
void reduce_block(real_t* out, real_t* out_red,
                  int* Nsize, int* block_size, int ieo)
{
  real_t* out_dev = (real_t*)dev_ptr(out);
  real_t* red_dev = (real_t*)dev_ptr(out_red);

  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst  = Nx * Ny * Nz * Nt;

  // size of block
  int Bx = block_size[0];
  int By = block_size[1];
  int Bz = block_size[2];
  int Bt = block_size[3];

  // numbers of blocks
  int NBx  = Nx/Bx;
  int NBy  = Ny/By;
  int NBz  = Nz/Bz;
  int NBt  = Nt/Bt;
  int Nblock = NBx * NBy * NBz * NBt;

  int nth = VECTOR_LENGTH;
  int nbl = Nblock/nth;

  if(nbl * nth != Nblock){
    printf("unsupported value of Nblock = %d\n", Nblock);
    exit(1);
  }

  reduce_block_kernel<<<nbl,nth>>>(out_dev, red_dev,
                                   Nx, Ny, Nz, Nt,
                                   Bx, By, Bz, Bt, ieo);
 
  CHECK(cudaDeviceSynchronize());

}

//============================================================END=====
