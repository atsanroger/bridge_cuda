/*!
      @file    afield_dd_cuda-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2024-04-01 12:13:09 #$
      @version $LastChangedRevision: 2595 $
*/

//====================================================================
// Launch-overhead reduction for the block (domain-decomposed) field ops.
//
// These kernels run on the legacy default stream, so a block_dotc_eo ->
// reduce_block -> block_axpy_eo chain (and the fused block_orthogonalize Gram-
// Schmidt step) is already correctly ordered without a host-side
// cudaDeviceSynchronize() after each launch. The real barrier is provided where
// the host reads device data: the buf_out_update_host()/update_device() copies in
// the orchestration (afield_dd-inc.h) are synchronizing cudaMemcpy on the same
// stream. The per-kernel sync only forced the host to block, which on this WSL
// box (high cudaDeviceSynchronize latency) dominated the O(nvector) Gram-Schmidt
// setup cost. Set AFIELD_DD_SYNC_EACH_KERNEL=1 to restore it for debugging.
#ifndef AFIELD_DD_SYNC_EACH_KERNEL
#define AFIELD_DD_SYNC_EACH_KERNEL 0
#endif
static inline void afield_dd_kernel_sync()
{
#if AFIELD_DD_SYNC_EACH_KERNEL
  CHECK(cudaDeviceSynchronize());
#endif
}

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

   afield_dd_kernel_sync();

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

   afield_dd_kernel_sync();

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

   afield_dd_kernel_sync();

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
 
   afield_dd_kernel_sync();

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
 
   afield_dd_kernel_sync();

}

//====================================================================
// Batched block Gram-Schmidt (Classical GS, optionally repeated = CGS2).
//
// Orthogonalize a target vector against nbasis basis vectors *at once*, computing
// all per-block inner products <basis[j]|target> in one kernel and subtracting all
// components in one kernel. This collapses the launch count from O(nbasis) per
// target (the per-j block_dotc/block_axpy MGS loop) to O(1) per CGS pass, which is
// the win on launch-bound hardware. CGS is numerically weaker than MGS, so it is
// repeated ncgs times (ncgs=2 -> CGS2) to recover orthogonality. Coefficients stay
// on the device throughout. Only the all-block case (as used by the testvector
// Gram-Schmidt, ieo<0) is implemented. Per-(site,j) arithmetic mirrors the single-
// vector block_dotc_eo/block_axpy_eo kernels exactly, so each pair's result is
// identical to the unbatched path.

// per-site partials: part[ 2*(j*Nst + site2) + {0,1} ] = <basis[j]|target> at site
__global__
void block_dotc_multi_kernel(real_t* part, const real_t* tp,
                             real_t* const* basis, int nbasis,
                             int Nin, int Nex,
                             int Nx, int Ny, int Nz, int Nt,
                             int Bx, int By, int Bz, int Bt)
{
  int Nst   = Nx * Ny * Nz * Nt;
  int Bsize = Bx * By * Bz * Bt;
  int NBx = Nx/Bx, NBy = Ny/By, NBz = Nz/Bz;

  int site2 = blockIdx.x * blockDim.x + threadIdx.x;
  if(site2 >= Nst) return;

  int bsite = site2 % Bsize;
  int block = site2 / Bsize;
  int ibx = block % NBx;
  int iby = (block/NBx) % NBy;
  int ibz = (block/(NBx*NBy)) % NBz;
  int ibt = block/(NBx*NBy*NBz);
  int kx = bsite % Bx;      int ix = kx + Bx*ibx;
  int kyzt = bsite/Bx;      int ky = kyzt % By;  int iy = ky + By*iby;
  int kzt = kyzt/By;        int kz = kzt % Bz;   int iz = kz + Bz*ibz;
  int kt = kzt/Bz;          int it = kt + Bt*ibt;
  int site = ix + Nx*(iy + Ny*(iz + Nz*it));

  for(int j = 0; j < nbasis; ++j){
    const real_t* wp = basis[j];
    real_t xtr = 0.0, xti = 0.0;
    for(int ex = 0; ex < Nex; ++ex){
      int nv = Nin * Nst * ex;
      for(int in2 = 0; in2 < Nin/2; ++in2){
        int inr = 2*in2, ini = 2*in2 + 1;
        real_t vtr = wp[nv + IDX2(Nin, inr, site)];
        real_t vti = wp[nv + IDX2(Nin, ini, site)];
        real_t wtr = tp[nv + IDX2(Nin, inr, site)];
        real_t wti = tp[nv + IDX2(Nin, ini, site)];
        xtr += vtr*wtr + vti*wti;   // Re <basis[j]|target>
        xti += vtr*wti - vti*wtr;   // Im <basis[j]|target>
      }
    }
    part[2*(j*Nst + site2)]     = xtr;
    part[2*(j*Nst + site2) + 1] = xti;
  }
}

// reduce per (block,j): c[ 2*(j*Nblock + block) + {0,1} ]
__global__
void reduce_block_multi_kernel(real_t* c, const real_t* part, int nbasis,
                               int Nx, int Ny, int Nz, int Nt,
                               int Bx, int By, int Bz, int Bt)
{
  int Nst   = Nx * Ny * Nz * Nt;
  int Bsize = Bx * By * Bz * Bt;
  int Nblock = (Nx/Bx)*(Ny/By)*(Nz/Bz)*(Nt/Bt);

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if(idx >= Nblock*nbasis) return;
  int block = idx % Nblock;
  int j     = idx / Nblock;
  const real_t* vp = &part[2*(j*Nst + block*Bsize)];
  real_t vtr = 0.0, vti = 0.0;
  for(int bsite = 0; bsite < Bsize; ++bsite){
    vtr += vp[2*bsite];
    vti += vp[2*bsite + 1];
  }
  c[2*(j*Nblock + block)]     = vtr;
  c[2*(j*Nblock + block) + 1] = vti;
}

// target += fac * sum_j c[j][block] * basis[j]   (per site, sequential j accum)
__global__
void block_axpy_multi_kernel(real_t* tp, const real_t* c,
                             real_t* const* basis, int nbasis, real_t fac,
                             int Nin, int Nex,
                             int Nx, int Ny, int Nz, int Nt,
                             int Bx, int By, int Bz, int Bt)
{
  int Nst   = Nx * Ny * Nz * Nt;
  int Bsize = Bx * By * Bz * Bt;
  int NBx = Nx/Bx, NBy = Ny/By, NBz = Nz/Bz;
  int Nblock = NBx*NBy*NBz*(Nt/Bt);

  int site2 = blockIdx.x * blockDim.x + threadIdx.x;
  if(site2 >= Nst) return;

  int bsite = site2 % Bsize;
  int block = site2 / Bsize;
  int ibx = block % NBx;
  int iby = (block/NBx) % NBy;
  int ibz = (block/(NBx*NBy)) % NBz;
  int ibt = block/(NBx*NBy*NBz);
  int kx = bsite % Bx;      int ix = kx + Bx*ibx;
  int kyzt = bsite/Bx;      int ky = kyzt % By;  int iy = ky + By*iby;
  int kzt = kyzt/By;        int kz = kzt % Bz;   int iz = kz + Bz*ibz;
  int kt = kzt/Bz;          int it = kt + Bt*ibt;
  int site = ix + Nx*(iy + Ny*(iz + Nz*it));

  for(int j = 0; j < nbasis; ++j){
    real_t ar = c[2*(j*Nblock + block)]     * fac;
    real_t ai = c[2*(j*Nblock + block) + 1] * fac;
    const real_t* wp = basis[j];
    for(int ex = 0; ex < Nex; ++ex){
      int nv = Nin * Nst * ex;
      for(int in2 = 0; in2 < Nin/2; ++in2){
        int inr = 2*in2, ini = 2*in2 + 1;
        real_t wtr = wp[nv + IDX2(Nin, inr, site)];
        real_t wti = wp[nv + IDX2(Nin, ini, site)];
        real_t vtr = tp[nv + IDX2(Nin, inr, site)];
        real_t vti = tp[nv + IDX2(Nin, ini, site)];
        vtr += ar*wtr - ai*wti;
        vti += ar*wti + ai*wtr;
        tp[nv + IDX2(Nin, inr, site)] = vtr;
        tp[nv + IDX2(Nin, ini, site)] = vti;
      }
    }
  }
}

// CGS2 driver: orthogonalize tp against basis_host[0..nbasis-1], ncgs passes.
void block_orthogonalize_all_eo(real_t* tp, real_t* const* basis_host, int nbasis,
                                real_t fac, int Nin, int Nex,
                                int* Nsize, int* block_size, int ncgs)
{
  int Nx = Nsize[0], Ny = Nsize[1], Nz = Nsize[2], Nt = Nsize[3];
  int Nst = Nx*Ny*Nz*Nt;
  int Bx = block_size[0], By = block_size[1], Bz = block_size[2], Bt = block_size[3];
  int Nblock = (Nx/Bx)*(Ny/By)*(Nz/Bz)*(Nt/Bt);

  real_t* tp_dev = (real_t*)dev_ptr(tp);

  // device scratch, grown on demand and kept for the program lifetime
  static real_t*  part_dev  = nullptr;  // [2*Nst*nbasis]
  static real_t*  c_dev     = nullptr;  // [2*Nblock*nbasis]
  static real_t** basis_dev = nullptr;  // [nbasis] device array of device pointers
  static int cap = 0;
  if(nbasis > cap){
    if(part_dev)  cudaFree(part_dev);
    if(c_dev)     cudaFree(c_dev);
    if(basis_dev) cudaFree(basis_dev);
    CHECK(cudaMalloc((void**)&part_dev,  sizeof(real_t) * 2 * Nst * nbasis));
    CHECK(cudaMalloc((void**)&c_dev,     sizeof(real_t) * 2 * Nblock * nbasis));
    CHECK(cudaMalloc((void**)&basis_dev, sizeof(real_t*) * nbasis));
    cap = nbasis;
  }

  // map basis host base pointers to device pointers, upload the pointer array
  real_t* bptr_host[256];
  if(nbasis > 256){ printf("block_orthogonalize_all_eo: nbasis=%d > 256\n", nbasis); exit(1); }
  for(int j = 0; j < nbasis; ++j) bptr_host[j] = (real_t*)dev_ptr(basis_host[j]);
  CHECK(cudaMemcpy(basis_dev, bptr_host, sizeof(real_t*)*nbasis, cudaMemcpyHostToDevice));

  int nth = VECTOR_LENGTH;
  int nbl_site = (Nst + nth - 1)/nth;
  int nbl_bj   = (Nblock*nbasis + nth - 1)/nth;

  // Classical block Gram-Schmidt, ncgs passes (ncgs=1 here; see the note in
  // MultiGrid_Domainwall-tmpl.h -- a 2nd reorthogonalization pass is invalid for
  // the two-chirality projector basis as it exactly undoes the first).
  for(int pass = 0; pass < ncgs; ++pass){
    block_dotc_multi_kernel<<<nbl_site, nth>>>(part_dev, tp_dev, basis_dev, nbasis,
                                               Nin, Nex, Nx,Ny,Nz,Nt, Bx,By,Bz,Bt);
    reduce_block_multi_kernel<<<nbl_bj, nth>>>(c_dev, part_dev, nbasis,
                                               Nx,Ny,Nz,Nt, Bx,By,Bz,Bt);
    block_axpy_multi_kernel<<<nbl_site, nth>>>(tp_dev, c_dev, basis_dev, nbasis, fac,
                                               Nin, Nex, Nx,Ny,Nz,Nt, Bx,By,Bz,Bt);
  }
  afield_dd_kernel_sync();
}

//============================================================END=====
