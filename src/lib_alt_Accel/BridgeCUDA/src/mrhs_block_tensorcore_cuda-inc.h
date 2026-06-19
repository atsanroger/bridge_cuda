/*!
      @file    mrhs_block_tensorcore_cuda-inc.h
      @brief   Multiple-RHS (MRHS) batched block restriction / prolongation,
               cast as small dense GEMMs and run on the tensor cores.
      @author  Wei-Lun Chen (wlchen)

      This kernel is the supercomputer-target form of the two-grid transfer
      operators used by the PV-preconditioned domain-wall AMG solver.  It is the
      intersection of three ingredients:

        (1) MRHS batching  -- Boyle, "Multiple right hand side multigrid for
            domain wall fermions" (arXiv:2409.03904).  A block of N_rhs sources
            is restricted/prolonged at once, so the per-block transfer is no
            longer a mat-vec but a dense matrix-matrix product
                 C[nbasis x N_rhs] = B^H[nbasis x K] * R[K x N_rhs]   (restrict)
                 F[K x N_rhs]      = B  [K x nbasis] * C[nbasis x N_rhs] (prolong)
            with the contraction K = (block volume) x Nin running over the fine
            degrees of freedom inside one aggregate.

        (2) Tensor cores  -- Tu et al., "Accelerating High-Order Finite Element
            Simulations at Extreme Scale with FP64 Tensor Cores".  The transfer
            decomposes into a batch of small dense GEMMs of order O(10); the
            CUDA-core form is shared-memory-bandwidth bound (flop/byte ~0.1), so
            we gather each aggregate to a packed staging tile (the paper's
            L-vector -> E-vector step) and contract it with warp-level MMA tiles.

        (3) Coarse precision is FP32, by design.  The V-cycle is only a
            preconditioner: its accuracy sets the convergence rate, never the
            final solution (the outer FGMRES on A carries the real precision), so
            the coarse transfer needs nothing beyond FP32 -- and emphatically not
            the sub-FP32 multiword/eps^3 arithmetic, which is a fine-solve tool
            (the QDW/QTW CG) and is deliberately out of scope here.  The only
            subtlety is the hardware: this card's tensor cores accumulate in
            TF32 (~10-bit mantissa), i.e. *below* FP32.  For a preconditioner
            that is normally fine (cf. half-precision coarse grids in QUDA), so
            the default tensor-core path is a single TF32 pass.  If TF32 ever
            proves too coarse, MRHS_TC_FP32_SPLIT enables a 3-pass float-float
            split (hi*hi + hi*lo + lo*hi) that reconstructs full FP32 -- still
            FP32, never FP64, just recovered on TF32 units.

      The complex fine field is stored as interleaved (re,im) in the existing
      NWP-tiled IDX2 layout, identical to block_dotc_multi_kernel /
      block_axpy_multi_kernel in afield_dd_cuda-inc.h, so the reference path
      below is bit-compatible with the single-RHS solver and serves as the
      correctness oracle for the tensor-core path.
*/

#ifndef MRHS_BLOCK_TENSORCORE_CUDA_INCLUDED
#define MRHS_BLOCK_TENSORCORE_CUDA_INCLUDED

namespace {

//====================================================================
// Reference (CUDA-core) MRHS restriction.  Bit-compatible with the proven
// single-RHS block_dotc_multi_kernel: C[j,r,block] = sum_{site in block,in}
// conj(basis[j])_in * rhs[r]_in.  rhs_base[] is N_rhs device pointers into the
// same fine layout; coarse[] is laid out 2*(r + N_rhs*(j + nbasis*block)).
__global__
void block_restrict_mrhs_ref_dev(real_t* coarse,
                                 const real_t* tp_unused,      // kept for symmetry
                                 real_t* const* basis, int nbasis,
                                 real_t* const* rhs,   int nrhs,
                                 int Nin, int Nex,
                                 int Nx, int Ny, int Nz, int Nt,
                                 int Bx, int By, int Bz, int Bt)
{
  (void)tp_unused;
  int Nst    = Nx * Ny * Nz * Nt;
  int Bsize  = Bx * By * Bz * Bt;
  int NBx = Nx/Bx, NBy = Ny/By, NBz = Nz/Bz, NBt = Nt/Bt;
  int Nblock = NBx*NBy*NBz*NBt;

  // one warp-ish thread per (block, j, r)
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if(idx >= Nblock*nbasis*nrhs) return;
  int r     = idx % nrhs;
  int j     = (idx / nrhs) % nbasis;
  int block = idx / (nrhs*nbasis);

  int ibx = block % NBx;
  int iby = (block/NBx) % NBy;
  int ibz = (block/(NBx*NBy)) % NBz;
  int ibt = block/(NBx*NBy*NBz);

  const real_t* bp = basis[j];
  const real_t* rp = rhs[r];
  real_t cr = 0.0, ci = 0.0;
  for(int bsite = 0; bsite < Bsize; ++bsite){
    int kx = bsite % Bx;            
    int ix = kx + Bx*ibx;
    int kyzt = bsite/Bx;            
    int ky = kyzt % By;  
    int iy = ky + By*iby;
    int kzt = kyzt/By;             
    int kz = kzt % Bz;   
    int iz = kz + Bz*ibz;
    int kt = kzt/Bz;              
    int it = kt + Bt*ibt;
    int site = ix + Nx*(iy + Ny*(iz + Nz*it));
    for(int ex = 0; ex < Nex; ++ex){
      int nv = Nin * Nst * ex;
      for(int in2 = 0; in2 < Nin/2; ++in2){
        int inr = 2*in2, ini = 2*in2 + 1;
        real_t br = bp[nv + IDX2(Nin, inr, site)];
        real_t bi = bp[nv + IDX2(Nin, ini, site)];
        real_t xr = rp[nv + IDX2(Nin, inr, site)];
        real_t xi = rp[nv + IDX2(Nin, ini, site)];
        cr += br*xr + bi*xi;   // Re conj(b).x
        ci += br*xi - bi*xr;   // Im conj(b).x
      }
    }
  }
  coarse[2*(r + nrhs*(j + nbasis*block))]     = cr;
  coarse[2*(r + nrhs*(j + nbasis*block)) + 1] = ci;
}

// Reference MRHS prolongation: fine[r] += fac * sum_j basis[j] * C[j,r,block].
__global__
void block_prolong_mrhs_ref_dev(real_t* const* fine, int nrhs,
                                const real_t* coarse,
                                real_t* const* basis, int nbasis, real_t fac,
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

  for(int r = 0; r < nrhs; ++r){
    real_t* fp = fine[r];
    for(int ex = 0; ex < Nex; ++ex){
      int nv = Nin * Nst * ex;
      for(int in2 = 0; in2 < Nin/2; ++in2){
        int inr = 2*in2, ini = 2*in2 + 1;
        real_t vr = fp[nv + IDX2(Nin, inr, site)];
        real_t vi = fp[nv + IDX2(Nin, ini, site)];
        for(int j = 0; j < nbasis; ++j){
          real_t cr = coarse[2*(r + nrhs*(j + nbasis*block))]     * fac;
          real_t ci = coarse[2*(r + nrhs*(j + nbasis*block)) + 1] * fac;
          const real_t* bp = basis[j];
          real_t br = bp[nv + IDX2(Nin, inr, site)];
          real_t bi = bp[nv + IDX2(Nin, ini, site)];
          vr += br*cr - bi*ci;
          vi += br*ci + bi*cr;
        }
        fp[nv + IDX2(Nin, inr, site)] = vr;
        fp[nv + IDX2(Nin, ini, site)] = vi;
      }
    }
  }
}

// Reference (CUDA-core) MRHS coarse Dirac operator: the FP32 oracle / fallback
// for mult_coarse_mrhs_tc_dev.  Bit-compatible with accum_mult_u: one thread per
// (site, a, r); V2[a,r] = V1[a,r] + sum_8dir bc * sum_b U_dir[a,b] V1_nei[b,r].
__global__
void mult_coarse_mrhs_ref_dev(real_t* const* v2, real_t* const* v1,
                              const real_t* u, int nrhs, int Ncol,
                              int Nx, int Ny, int Nz, int Nt,
                              int bcx, int bcy, int bcz, int bct)
{
  int Nst = Nx*Ny*Nz*Nt;
  int Ncol2 = 2*Ncol, Ndf = 2*Ncol*Ncol, Ngst1 = Ndf*Nst;
  int Nxy = Nx*Ny, Nxyz = Nx*Ny*Nz;

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if(idx >= Nst*Ncol*nrhs) return;
  int r = idx % nrhs; 
  int a = (idx / nrhs) % Ncol; 
  int site = idx / (nrhs*Ncol);
  int ix = site%Nx, iy = (site/Nx)%Ny, iz = (site/Nxy)%Nz, it = site/Nxyz;

  float accr = v1[r][IDX2(Ncol2, 2*a,   site)];
  float acci = v1[r][IDX2(Ncol2, 2*a+1, site)];
  for(int d = 0; d < 8; ++d){
    int uoff; int nei; float bc2 = 1.0f;
    switch(d){
      case 0: uoff=Ngst1*0; nei=(ix==Nx-1)?(site-(Nx-1)):(site+1);       if(ix==Nx-1)bc2=bcx; break;
      case 1: uoff=Ngst1*4; nei=(ix==0)?(site+Nx-1):(site-1);            if(ix==0)   bc2=bcx; break;
      case 2: uoff=Ngst1*1; nei=(iy==Ny-1)?(site-Nx*(Ny-1)):(site+Nx);   if(iy==Ny-1)bc2=bcy; break;
      case 3: uoff=Ngst1*5; nei=(iy==0)?(site+Nx*(Ny-1)):(site-Nx);      if(iy==0)   bc2=bcy; break;
      case 4: uoff=Ngst1*2; nei=(iz==Nz-1)?(site-Nxy*(Nz-1)):(site+Nxy); if(iz==Nz-1)bc2=bcz; break;
      case 5: uoff=Ngst1*6; nei=(iz==0)?(site+Nxy*(Nz-1)):(site-Nxy);    if(iz==0)   bc2=bcz; break;
      case 6: uoff=Ngst1*3; nei=(it==Nt-1)?(site-Nxyz*(Nt-1)):(site+Nxyz);if(it==Nt-1)bc2=bct; break;
      default:uoff=Ngst1*7; nei=(it==0)?(site+Nxyz*(Nt-1)):(site-Nxyz);  if(it==0)   bc2=bct; break;
    }
    for(int b = 0; b < Ncol; ++b){
      int in_r = Ncol2*b + 2*a;
      float mr = u[uoff + IDX2(Ndf, in_r,   site)];
      float mi = u[uoff + IDX2(Ndf, in_r+1, site)];
      float vr = v1[r][IDX2(Ncol2, 2*b,   nei)];
      float vi = v1[r][IDX2(Ncol2, 2*b+1, nei)];
      accr += bc2 * (mr*vr - mi*vi);
      acci += bc2 * (mr*vi + mi*vr);
    }
  }
  v2[r][IDX2(Ncol2, 2*a,   site)] = accr;
  v2[r][IDX2(Ncol2, 2*a+1, site)] = acci;
}

// Reference (CUDA-core) block inner product G[i,j] = <V_i|W_j>, one thread per
// (i,j) reducing over the whole volume.  FP32 oracle / fallback for block_inner.
__global__
void block_inner_ref_dev(real_t* G, real_t* const* V, real_t* const* W,
                         int nv, int Nin, int Nex, int Nst)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if(idx >= nv*nv) return;
  int j = idx % nv, i = idx / nv;
  double gr = 0.0, gi = 0.0;
  for(int ex = 0; ex < Nex; ++ex){
    int nvo = Nin*Nst*ex;
    for(int site = 0; site < Nst; ++site){
      for(int in2 = 0; in2 < Nin/2; ++in2){
        float vr = V[i][nvo + IDX2(Nin, 2*in2, site)], vi = V[i][nvo + IDX2(Nin, 2*in2+1, site)];
        float wr = W[j][nvo + IDX2(Nin, 2*in2, site)], wi = W[j][nvo + IDX2(Nin, 2*in2+1, site)];
        gr += (double)vr*wr + (double)vi*wi;
        gi += (double)vr*wi - (double)vi*wr;
      }
    }
  }
  G[2*(j + nv*i)]     = gr;
  G[2*(j + nv*i) + 1] = gi;
}

// Reference block update W_j -= sum_i V_i G[i,j], one thread per (dof,j).
__global__
void block_update_ref_dev(real_t* const* W, real_t* const* V, const real_t* G,
                          int nv, int Nin, int Nex, int Nst)
{
  int Ktot = Nex*Nst*(Nin/2);
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if(idx >= Ktot*nv) return;
  int j = idx % nv, kg = idx / nv;
  int in2 = kg % (Nin/2); int t1 = kg / (Nin/2);
  int site = t1 % Nst; int ex = t1 / Nst; int nvo = Nin*Nst*ex;
  float wr = W[j][nvo + IDX2(Nin, 2*in2, site)], wi = W[j][nvo + IDX2(Nin, 2*in2+1, site)];
  for(int i = 0; i < nv; ++i){
    float gr = G[2*(j + nv*i)], gi = G[2*(j + nv*i) + 1];
    float vr = V[i][nvo + IDX2(Nin, 2*in2, site)], vi = V[i][nvo + IDX2(Nin, 2*in2+1, site)];
    wr -= vr*gr - vi*gi;
    wi -= vr*gi + vi*gr;
  }
  W[j][nvo + IDX2(Nin, 2*in2,   site)] = wr;
  W[j][nvo + IDX2(Nin, 2*in2+1, site)] = wi;
}

//====================================================================
// Tensor-core MRHS restriction.
//
// One warp (one CUDA block) owns one aggregate.  Step 1 gathers the aggregate's
// basis and rhs columns out of the NWP-tiled IDX2 layout into row-major shared
// staging tiles (the paper's L->E-vector step).  Step 2 contracts them with the
// m16n8k8 TF32 tensor-core instruction, issued directly as inline PTX.  Step 3
// writes the (nbasis x N_rhs) complex coarse block.  TF32-accumulate-in-FP32 is
// ample for a preconditioner; this is the supercomputer-target transfer kernel.

#ifdef USE_TENSORCORE_MRHS

// One MMA tile shape: M16 x N8 x K8, the FP32-accumulate TF32 tensor-core op.
// nbasis padded to 16 (M), N_rhs handled 8 at a time (N), contraction tiled by 8.
static const int MMA_M = 16;
static const int MMA_N = 8;
static const int MMA_K = 8;

// FP32 -> TF32, rounded, packed in the low bits of a .b32 register.  This is the
// only "precision" step: the tensor core consumes TF32, accumulates in FP32.
__device__ __forceinline__ unsigned f2tf32(float x)
{
  unsigned r;
  asm("cvt.rna.tf32.f32 %0, %1;" : "=r"(r) : "f"(x));
  return r;
}

// The instruction itself -- programmed directly, the way Tu et al. drive the
// FP64 DMMA in the tsunami code, here the sm_80 TF32 variant:
//   D[16x8] = A[16x8] * B[8x8] + C[16x8]
// A: 4 tf32 regs/thread, B: 2 tf32 regs/thread, C/D: 4 f32 regs/thread.
__device__ __forceinline__
void mma_m16n8k8(float d[4], const unsigned a[4], const unsigned b[2], const float c[4])
{
  asm volatile(
    "mma.sync.aligned.m16n8k8.row.col.f32.tf32.tf32.f32 "
    "{%0,%1,%2,%3}, {%4,%5,%6,%7}, {%8,%9}, {%10,%11,%12,%13};\n"
    : "=f"(d[0]), "=f"(d[1]), "=f"(d[2]), "=f"(d[3])
    : "r"(a[0]), "r"(a[1]), "r"(a[2]), "r"(a[3]),
      "r"(b[0]), "r"(b[1]),
      "f"(c[0]), "f"(c[1]), "f"(c[2]), "f"(c[3]));
}

// Accumulate one K-chunk c[4] += sgn * A[16x8] * B[8x8], A and B in row-major
// shared memory.  The lane->element maps are the canonical m16n8k8 fragment
// layouts (groupID = lane/4, threadID-in-group = lane%4).
__device__ __forceinline__
void tc_chunk(float c[4], const float* Ash, const float* Bsh, float sgn)
{
  int lane = threadIdx.x & 31;
  int gid  = lane >> 2;     // groupID         0..7
  int tig  = lane & 3;      // threadID-in-grp 0..3
  unsigned a[4], b[2];
  // A (16x8 row-major): a0=(gid,tig) a1=(gid+8,tig) a2=(gid,tig+4) a3=(gid+8,tig+4)
  a[0] = f2tf32(sgn * Ash[(gid    )*MMA_K + (tig    )]);
  a[1] = f2tf32(sgn * Ash[(gid + 8)*MMA_K + (tig    )]);
  a[2] = f2tf32(sgn * Ash[(gid    )*MMA_K + (tig + 4)]);
  a[3] = f2tf32(sgn * Ash[(gid + 8)*MMA_K + (tig + 4)]);
  // B (8x8): b0=(row tig, col gid)  b1=(row tig+4, col gid)
  b[0] = f2tf32(Bsh[(tig    )*MMA_N + gid]);
  b[1] = f2tf32(Bsh[(tig + 4)*MMA_N + gid]);
  mma_m16n8k8(c, a, b, c);
}

//====================================================================
// Tensor-core MRHS restriction: C[j,r,block] = sum_{k in aggregate} conj(B[j,k]) R[k,r]
// run as a batch of m16n8k8 TF32 GEMMs.  One warp (one CUDA block) per aggregate;
// N_rhs handled in tiles of 8.  Complex via the real expansion
//   C_re = Bre.Rre + Bim.Rim ,  C_im = Bre.Rim - Bim.Rre .
__global__
void block_restrict_mrhs_tc_dev(real_t* coarse,
                                real_t* const* basis, int nbasis,
                                real_t* const* rhs,   int nrhs,
                                int Nin, int Nex,
                                int Nx, int Ny, int Nz, int Nt,
                                int Bx, int By, int Bz, int Bt)
{
  int Nst   = Nx * Ny * Nz * Nt;
  int Bsize = Bx * By * Bz * Bt;
  int NBx = Nx/Bx, NBy = Ny/By, NBz = Nz/Bz;
  int Ktot = Nex * Bsize * (Nin/2);            // complex contraction length

  int block = blockIdx.x;
  int ibx = block % NBx;
  int iby = (block/NBx) % NBy;
  int ibz = (block/(NBx*NBy)) % NBz;
  int ibt = block/(NBx*NBy*NBz);

  __shared__ float Bre[MMA_M*MMA_K], Bim[MMA_M*MMA_K];   // basis tiles (M x K)
  __shared__ float Rre[MMA_K*MMA_N], Rim[MMA_K*MMA_N];   // rhs tiles   (K x N)

  for(int rtile = 0; rtile < nrhs; rtile += MMA_N){       // N tiles of 8 RHS
    float c_re[4] = {0,0,0,0};
    float c_im[4] = {0,0,0,0};

    for(int kc = 0; kc < Ktot; kc += MMA_K){              // contraction tiles
      // --- gather the aggregate into row-major staging tiles (L->E vector) ---
      for(int idx = (threadIdx.x & 31); idx < MMA_M*MMA_K; idx += 32){
        int m = idx / MMA_K, kl = idx % MMA_K;
        int kg = kc + kl;
        float vr = 0.f, vi = 0.f;
        if(m < nbasis && kg < Ktot){
          int in2 = kg % (Nin/2); int t1 = kg / (Nin/2);
          int bsite = t1 % Bsize; int ex = t1 / Bsize;
          int kx = bsite % Bx;      int ix = kx + Bx*ibx;
          int kyzt = bsite/Bx;      int ky = kyzt % By;  int iy = ky + By*iby;
          int kzt = kyzt/By;        int kz = kzt % Bz;   int iz = kz + Bz*ibz;
          int kt = kzt/Bz;          int it = kt + Bt*ibt;
          int site = ix + Nx*(iy + Ny*(iz + Nz*it));
          int nv = Nin * Nst * ex;
          vr = basis[m][nv + IDX2(Nin, 2*in2,   site)];
          vi = basis[m][nv + IDX2(Nin, 2*in2+1, site)];
        }
        Bre[idx] = vr; Bim[idx] = vi;
      }
      for(int idx = (threadIdx.x & 31); idx < MMA_K*MMA_N; idx += 32){
        int kl = idx / MMA_N, r = idx % MMA_N;
        int kg = kc + kl; int rr = rtile + r;
        float vr = 0.f, vi = 0.f;
        if(rr < nrhs && kg < Ktot){
          int in2 = kg % (Nin/2); int t1 = kg / (Nin/2);
          int bsite = t1 % Bsize; int ex = t1 / Bsize;
          int kx = bsite % Bx;      int ix = kx + Bx*ibx;
          int kyzt = bsite/Bx;      int ky = kyzt % By;  int iy = ky + By*iby;
          int kzt = kyzt/By;        int kz = kzt % Bz;   int iz = kz + Bz*ibz;
          int kt = kzt/Bz;          int it = kt + Bt*ibt;
          int site = ix + Nx*(iy + Ny*(iz + Nz*it));
          int nv = Nin * Nst * ex;
          vr = rhs[rr][nv + IDX2(Nin, 2*in2,   site)];
          vi = rhs[rr][nv + IDX2(Nin, 2*in2+1, site)];
        }
        Rre[idx] = vr; Rim[idx] = vi;
      }
      __syncwarp();
      // --- the GEMM tiles: conj(B).R, real expansion ---
      tc_chunk(c_re, Bre, Rre, +1.0f);   // C_re += Bre.Rre
      tc_chunk(c_re, Bim, Rim, +1.0f);   // C_re += Bim.Rim
      tc_chunk(c_im, Bre, Rim, +1.0f);   // C_im += Bre.Rim
      tc_chunk(c_im, Bim, Rre, -1.0f);   // C_im -= Bim.Rre
      __syncwarp();
    }

    // --- scatter the 16x8 complex result (canonical C-fragment map) ---
    int lane = threadIdx.x & 31, gid = lane >> 2, tig = lane & 3;
    int jj[2] = { gid, gid + 8 };
    for(int half = 0; half < 2; ++half){
      int j = jj[half];
      for(int cc = 0; cc < 2; ++cc){
        int r = rtile + (tig*2 + cc);
        if(j < nbasis && r < nrhs){
          int o = 2*(r + nrhs*(j + nbasis*block));
          coarse[o]     = c_re[half*2 + cc];
          coarse[o + 1] = c_im[half*2 + cc];
        }
      }
    }
  }
}

//====================================================================
// Tensor-core MRHS prolongation: fine[k,r] += fac * sum_j B[j,k] C[j,r,block],
// the same m16n8k8 GEMM with the roles swapped -- contraction over the nbasis
// near-null directions (<=16, tiled by 8), output over the aggregate dof (tiled
// by 16) and scattered back to the fine field.  No conjugation here:
//   F_re = Bre.Cre - Bim.Cim ,  F_im = Bre.Cim + Bim.Cre .
__global__
void block_prolong_mrhs_tc_dev(real_t* const* fine, int nrhs,
                               const real_t* coarse,
                               real_t* const* basis, int nbasis, real_t fac,
                               int Nin, int Nex,
                               int Nx, int Ny, int Nz, int Nt,
                               int Bx, int By, int Bz, int Bt)
{
  int Nst   = Nx * Ny * Nz * Nt;
  int Bsize = Bx * By * Bz * Bt;
  int NBx = Nx/Bx, NBy = Ny/By, NBz = Nz/Bz;
  int Ktot = Nex * Bsize * (Nin/2);            // aggregate dof = output rows

  int block = blockIdx.x;
  int ibx = block % NBx;
  int iby = (block/NBx) % NBy;
  int ibz = (block/(NBx*NBy)) % NBz;
  int ibt = block/(NBx*NBy*NBz);

  __shared__ float Are[MMA_M*MMA_K], Aim[MMA_M*MMA_K];   // basis  (Mtile x nbasis-chunk)
  __shared__ float Cre[MMA_K*MMA_N], Cim[MMA_K*MMA_N];   // coarse (nbasis-chunk x N)

  // decode a fine dof index (0..Ktot) into (site,in2)
  auto decode = [&](int kg, int& site, int& in2){
    in2 = kg % (Nin/2); int t1 = kg / (Nin/2);
    int bsite = t1 % Bsize; int ex = t1 / Bsize;
    int kx = bsite % Bx;      int ix = kx + Bx*ibx;
    int kyzt = bsite/Bx;      int ky = kyzt % By;  int iy = ky + By*iby;
    int kzt = kyzt/By;        int kz = kzt % Bz;   int iz = kz + Bz*ibz;
    int kt = kzt/Bz;          int it = kt + Bt*ibt;
    site = ix + Nx*(iy + Ny*(iz + Nz*it));
    return ex;
  };

  for(int rtile = 0; rtile < nrhs; rtile += MMA_N){          // N tiles of 8 RHS
    for(int mt = 0; mt < Ktot; mt += MMA_M){                 // M tiles of 16 dof
      float c_re[4] = {0,0,0,0};
      float c_im[4] = {0,0,0,0};

      for(int kt = 0; kt < nbasis; kt += MMA_K){             // contract over nbasis
        // --- A = basis tile [16 dof x 8 basis] ---
        for(int idx = (threadIdx.x & 31); idx < MMA_M*MMA_K; idx += 32){
          int mrow = idx / MMA_K, jl = idx % MMA_K;
          int kg = mt + mrow; int j = kt + jl;
          float vr = 0.f, vi = 0.f;
          if(kg < Ktot && j < nbasis){
            int site, in2; int ex = decode(kg, site, in2);
            int nv = Nin * Nst * ex;
            vr = basis[j][nv + IDX2(Nin, 2*in2,   site)];
            vi = basis[j][nv + IDX2(Nin, 2*in2+1, site)];
          }
          Are[idx] = vr; Aim[idx] = vi;
        }
        // --- B = coarse tile [8 basis x 8 rhs] ---
        for(int idx = (threadIdx.x & 31); idx < MMA_K*MMA_N; idx += 32){
          int jl = idx / MMA_N, r = idx % MMA_N;
          int j = kt + jl; int rr = rtile + r;
          float vr = 0.f, vi = 0.f;
          if(j < nbasis && rr < nrhs){
            int o = 2*(rr + nrhs*(j + nbasis*block));
            vr = coarse[o]; vi = coarse[o + 1];
          }
          Cre[idx] = vr; Cim[idx] = vi;
        }
        __syncwarp();
        tc_chunk(c_re, Are, Cre, +1.0f);   // F_re += Bre.Cre
        tc_chunk(c_re, Aim, Cim, -1.0f);   // F_re -= Bim.Cim
        tc_chunk(c_im, Are, Cim, +1.0f);   // F_im += Bre.Cim
        tc_chunk(c_im, Aim, Cre, +1.0f);   // F_im += Bim.Cre
        __syncwarp();
      }

      // --- scatter: fine[k,r] += fac * F (canonical C-fragment map) ---
      int lane = threadIdx.x & 31, gid = lane >> 2, tig = lane & 3;
      int mm[2] = { mt + gid, mt + gid + 8 };
      for(int half = 0; half < 2; ++half){
        int kg = mm[half];
        if(kg >= Ktot) continue;
        int site, in2; int ex = decode(kg, site, in2);
        int nv = Nin * Nst * ex;
        for(int cc = 0; cc < 2; ++cc){
          int rr = rtile + (tig*2 + cc);
          if(rr < nrhs){
            real_t* fp = fine[rr];
            fp[nv + IDX2(Nin, 2*in2,   site)] += fac * c_re[half*2 + cc];
            fp[nv + IDX2(Nin, 2*in2+1, site)] += fac * c_im[half*2 + cc];
          }
        }
      }
    }
  }
}

//====================================================================
// Tensor-core MRHS coarse Dirac operator (single node, no halo):
//   V2[a,r] = V1[a,r] + sum_{8 dirs} bc * sum_b U_dir[a,b] V1_nei[b,r]
// The per-(site,direction) link is a dense Ncol x Ncol complex matrix; batching
// N_rhs makes each hop an (Ncol x Ncol)(Ncol x N_rhs) GEMM on the m16n8k8 tiles.
// Same complex expansion as prolongation (no conjugate).  One warp per site.
// Link element U_dir[a][b] real = u[uoff + IDX2(Ndf, Ncol2*b + 2*a, site)],
// neighbour vector V1_nei[b] = v1[r][IDX2(Ncol2, 2*b, nei)].  Mirrors the proven
// accum_mult_u arithmetic in mult_Domainwall_coarse_cuda-inc.h exactly.
__global__
void mult_coarse_mrhs_tc_dev(real_t* const* v2, real_t* const* v1,
                             const real_t* u, int nrhs, int Ncol,
                             int Nx, int Ny, int Nz, int Nt,
                             int bcx, int bcy, int bcz, int bct)
{
  int Nst   = Nx * Ny * Nz * Nt;
  int Ncol2 = 2 * Ncol;
  int Ndf   = 2 * Ncol * Ncol;
  int Ngst1 = Ndf * Nst;            // per-direction link stride
  int Nxy = Nx*Ny, Nxyz = Nx*Ny*Nz;

  int site = blockIdx.x;
  if(site >= Nst) return;
  int ix = site%Nx, iy = (site/Nx)%Ny, iz = (site/Nxy)%Nz, it = site/Nxyz;

  __shared__ float Mre[MMA_M*MMA_K], Mim[MMA_M*MMA_K];   // link tile  (a x b-chunk)
  __shared__ float Vre[MMA_K*MMA_N], Vim[MMA_K*MMA_N];   // nei block  (b-chunk x N)

  for(int rtile = 0; rtile < nrhs; rtile += MMA_N){
    float c_re[4] = {0,0,0,0};
    float c_im[4] = {0,0,0,0};

    for(int d = 0; d < 8; ++d){
      // direction d -> (link offset, neighbour site, boundary factor)
      int uoff; int nei; float bc2 = 1.0f;
      switch(d){
        case 0: uoff=Ngst1*0; nei=(ix==Nx-1)?(site-(Nx-1)):(site+1);      if(ix==Nx-1)bc2=bcx; break;
        case 1: uoff=Ngst1*4; nei=(ix==0)?(site+Nx-1):(site-1);           if(ix==0)   bc2=bcx; break;
        case 2: uoff=Ngst1*1; nei=(iy==Ny-1)?(site-Nx*(Ny-1)):(site+Nx);  if(iy==Ny-1)bc2=bcy; break;
        case 3: uoff=Ngst1*5; nei=(iy==0)?(site+Nx*(Ny-1)):(site-Nx);     if(iy==0)   bc2=bcy; break;
        case 4: uoff=Ngst1*2; nei=(iz==Nz-1)?(site-Nxy*(Nz-1)):(site+Nxy);if(iz==Nz-1)bc2=bcz; break;
        case 5: uoff=Ngst1*6; nei=(iz==0)?(site+Nxy*(Nz-1)):(site-Nxy);   if(iz==0)   bc2=bcz; break;
        case 6: uoff=Ngst1*3; nei=(it==Nt-1)?(site-Nxyz*(Nt-1)):(site+Nxyz);if(it==Nt-1)bc2=bct; break;
        default:uoff=Ngst1*7; nei=(it==0)?(site+Nxyz*(Nt-1)):(site-Nxyz); if(it==0)   bc2=bct; break;
      }

      for(int kt = 0; kt < Ncol; kt += MMA_K){             // contract over b
        for(int idx = (threadIdx.x & 31); idx < MMA_M*MMA_K; idx += 32){
          int a = idx / MMA_K, bl = idx % MMA_K; int b = kt + bl;
          float mr = 0.f, mi = 0.f;
          if(a < Ncol && b < Ncol){
            int in_r = Ncol2*b + 2*a;
            mr = bc2 * u[uoff + IDX2(Ndf, in_r,   site)];
            mi = bc2 * u[uoff + IDX2(Ndf, in_r+1, site)];
          }
          Mre[idx] = mr; Mim[idx] = mi;
        }
        for(int idx = (threadIdx.x & 31); idx < MMA_K*MMA_N; idx += 32){
          int bl = idx / MMA_N, r = idx % MMA_N; int b = kt + bl; int rr = rtile + r;
          float vr = 0.f, vi = 0.f;
          if(b < Ncol && rr < nrhs){
            vr = v1[rr][IDX2(Ncol2, 2*b,   nei)];
            vi = v1[rr][IDX2(Ncol2, 2*b+1, nei)];
          }
          Vre[idx] = vr; Vim[idx] = vi;
        }
        __syncwarp();
        tc_chunk(c_re, Mre, Vre, +1.0f);   // V2_re += Mre.Vre
        tc_chunk(c_re, Mim, Vim, -1.0f);   // V2_re -= Mim.Vim
        tc_chunk(c_im, Mre, Vim, +1.0f);   // V2_im += Mre.Vim
        tc_chunk(c_im, Mim, Vre, +1.0f);   // V2_im += Mim.Vre
        __syncwarp();
      }
    }

    // scatter: v2[a,r] = v1[site][a,r] + accumulated hops  (canonical C-map)
    int lane = threadIdx.x & 31, gid = lane >> 2, tig = lane & 3;
    int aa[2] = { gid, gid + 8 };
    for(int half = 0; half < 2; ++half){
      int a = aa[half];
      if(a >= Ncol) continue;
      for(int cc = 0; cc < 2; ++cc){
        int rr = rtile + (tig*2 + cc);
        if(rr < nrhs){
          float dr = v1[rr][IDX2(Ncol2, 2*a,   site)];
          float di = v1[rr][IDX2(Ncol2, 2*a+1, site)];
          v2[rr][IDX2(Ncol2, 2*a,   site)] = dr + c_re[half*2 + cc];
          v2[rr][IDX2(Ncol2, 2*a+1, site)] = di + c_im[half*2 + cc];
        }
      }
    }
  }
}

//====================================================================
// (C) Block-Krylov primitives -- the outer block-FGMRES orthogonalisation.
//
// For N_rhs right-hand sides the outer Arnoldi works on *blocks* of N_rhs full
// vectors.  Its two hot operations are dense GEMMs over the lattice volume:
//   block inner product   G[i,j] = <V_i | W_j>     (N_rhs x N_rhs Gram matrix)
//   block update          W_j   -= sum_i V_i G[i,j]
// i.e. the global analogue of restriction / prolongation with a single block =
// the whole lattice.  Same m16n8k8 tiles; the contraction is the full volume,
// parallelised over warps (one K-chunk each) with an atomic reduction into G.

// G[N_rhs x N_rhs] (complex) = V^H W, accumulated over the volume.  G must be
// pre-zeroed.  Layout G[2*(j_col + nv*i_row) + {0,1}] = <V_i | W_j>.
__global__
void block_inner_tc_dev(real_t* G, real_t* const* V, real_t* const* W,
                        int nv, int Nin, int Nex, int Nst)
{
  int Ktot = Nex * Nst * (Nin/2);
  int kc = blockIdx.x * MMA_K;
  if(kc >= Ktot) return;

  __shared__ float Vre[MMA_M*MMA_K], Vim[MMA_M*MMA_K];
  __shared__ float Wre[MMA_K*MMA_N], Wim[MMA_K*MMA_N];

  for(int rtile = 0; rtile < nv; rtile += MMA_N){
    float c_re[4] = {0,0,0,0};
    float c_im[4] = {0,0,0,0};
    for(int idx = (threadIdx.x & 31); idx < MMA_M*MMA_K; idx += 32){
      int m = idx / MMA_K, kl = idx % MMA_K; int kg = kc + kl;
      float vr = 0.f, vi = 0.f;
      if(m < nv && kg < Ktot){
        int in2 = kg % (Nin/2); int t1 = kg / (Nin/2);
        int site = t1 % Nst; int ex = t1 / Nst; int nvo = Nin*Nst*ex;
        vr = V[m][nvo + IDX2(Nin, 2*in2, site)];
        vi = V[m][nvo + IDX2(Nin, 2*in2+1, site)];
      }
      Vre[idx] = vr; Vim[idx] = vi;
    }
    for(int idx = (threadIdx.x & 31); idx < MMA_K*MMA_N; idx += 32){
      int kl = idx / MMA_N, r = idx % MMA_N; int kg = kc + kl; int rr = rtile + r;
      float vr = 0.f, vi = 0.f;
      if(rr < nv && kg < Ktot){
        int in2 = kg % (Nin/2); int t1 = kg / (Nin/2);
        int site = t1 % Nst; int ex = t1 / Nst; int nvo = Nin*Nst*ex;
        vr = W[rr][nvo + IDX2(Nin, 2*in2, site)];
        vi = W[rr][nvo + IDX2(Nin, 2*in2+1, site)];
      }
      Wre[idx] = vr; Wim[idx] = vi;
    }
    __syncwarp();
    tc_chunk(c_re, Vre, Wre, +1.0f);   // Re<V|W> = Vre.Wre + Vim.Wim
    tc_chunk(c_re, Vim, Wim, +1.0f);
    tc_chunk(c_im, Vre, Wim, +1.0f);   // Im<V|W> = Vre.Wim - Vim.Wre
    tc_chunk(c_im, Vim, Wre, -1.0f);
    __syncwarp();

    int lane = threadIdx.x & 31, gid = lane >> 2, tig = lane & 3;
    int ii[2] = { gid, gid + 8 };
    for(int half = 0; half < 2; ++half){
      int i = ii[half]; if(i >= nv) continue;
      for(int cc = 0; cc < 2; ++cc){
        int jc = rtile + (tig*2 + cc);
        if(jc < nv){
          int o = 2*(jc + nv*i);
          atomicAdd(&G[o],     c_re[half*2 + cc]);
          atomicAdd(&G[o + 1], c_im[half*2 + cc]);
        }
      }
    }
  }
}

// W_j -= sum_i V_i G[i,j], over the volume (one warp per output dof tile).
__global__
void block_update_tc_dev(real_t* const* W, real_t* const* V, const real_t* G,
                         int nv, int Nin, int Nex, int Nst)
{
  int Ktot = Nex * Nst * (Nin/2);
  int mt = blockIdx.x * MMA_M;
  if(mt >= Ktot) return;

  __shared__ float Vre[MMA_M*MMA_K], Vim[MMA_M*MMA_K];   // V tile (dof x nv-chunk)
  __shared__ float Gre[MMA_K*MMA_N], Gim[MMA_K*MMA_N];   // G tile (nv-chunk x N)

  for(int rtile = 0; rtile < nv; rtile += MMA_N){
    float c_re[4] = {0,0,0,0};
    float c_im[4] = {0,0,0,0};
    for(int kt = 0; kt < nv; kt += MMA_K){
      for(int idx = (threadIdx.x & 31); idx < MMA_M*MMA_K; idx += 32){
        int mrow = idx / MMA_K, il = idx % MMA_K; int kg = mt + mrow; int i = kt + il;
        float vr = 0.f, vi = 0.f;
        if(kg < Ktot && i < nv){
          int in2 = kg % (Nin/2); int t1 = kg / (Nin/2);
          int site = t1 % Nst; int ex = t1 / Nst; int nvo = Nin*Nst*ex;
          vr = V[i][nvo + IDX2(Nin, 2*in2, site)];
          vi = V[i][nvo + IDX2(Nin, 2*in2+1, site)];
        }
        Vre[idx] = vr; Vim[idx] = vi;
      }
      for(int idx = (threadIdx.x & 31); idx < MMA_K*MMA_N; idx += 32){
        int il = idx / MMA_N, r = idx % MMA_N; int i = kt + il; int jc = rtile + r;
        float gr = 0.f, gi = 0.f;
        if(i < nv && jc < nv){ int o = 2*(jc + nv*i); gr = G[o]; gi = G[o+1]; }
        Gre[idx] = gr; Gim[idx] = gi;
      }
      __syncwarp();
      tc_chunk(c_re, Vre, Gre, +1.0f);   // (V G)_re = Vre.Gre - Vim.Gim
      tc_chunk(c_re, Vim, Gim, -1.0f);
      tc_chunk(c_im, Vre, Gim, +1.0f);   // (V G)_im = Vre.Gim + Vim.Gre
      tc_chunk(c_im, Vim, Gre, +1.0f);
      __syncwarp();
    }
    int lane = threadIdx.x & 31, gid = lane >> 2, tig = lane & 3;
    int mm[2] = { mt + gid, mt + gid + 8 };
    for(int half = 0; half < 2; ++half){
      int kg = mm[half]; if(kg >= Ktot) continue;
      int in2 = kg % (Nin/2); int t1 = kg / (Nin/2);
      int site = t1 % Nst; int ex = t1 / Nst; int nvo = Nin*Nst*ex;
      for(int cc = 0; cc < 2; ++cc){
        int jc = rtile + (tig*2 + cc);
        if(jc < nv){
          real_t* wp = W[jc];
          wp[nvo + IDX2(Nin, 2*in2,   site)] -= c_re[half*2 + cc];
          wp[nvo + IDX2(Nin, 2*in2+1, site)] -= c_im[half*2 + cc];
        }
      }
    }
  }
}

#endif // USE_TENSORCORE_MRHS

//====================================================================
// Host drivers.  Both dispatch to the reference path by default; defining
// USE_TENSORCORE_MRHS at build time swaps in the tensor-core GEMM kernels.
// coarse_dev layout: 2*(r + nrhs*(j + nbasis*block)) interleaved (re,im).

// Runtime TF32(tensor-core) <-> FP32(CUDA-core reference) selector for the GEMM
// path.  Only has an effect when compiled with USE_TENSORCORE_MRHS (both kernel
// sets are then present in the TU); it lets the live solve A/B TF32 vs FP32
// transfer precision without a rebuild (verify_block_fgmres toggles it).  When
// the macro is OFF only the reference kernels exist, so the flag is forced false.
#ifdef USE_TENSORCORE_MRHS
static bool g_mrhs_use_tc = true;
#else
static bool g_mrhs_use_tc = false;
#endif
void mrhs_set_use_tc_dev(bool on) {
#ifdef USE_TENSORCORE_MRHS
  g_mrhs_use_tc = on;
#else
  (void)on;   // no tensor-core kernels compiled -- stays FP32 reference
#endif
}
bool mrhs_get_use_tc_dev() { return g_mrhs_use_tc; }

// fwd decl: cached host->device pointer-array upload (defined below).
static real_t** map_upload(real_t* const* host, int n, real_t**& dev, int& cap);

void block_restrict_mrhs(real_t* coarse_dev,
                         real_t* const* basis_host, int nbasis,
                         real_t* const* rhs_host,   int nrhs,
                         int Nin, int Nex, int* Nsize, int* block_size)
{
  int Nx = Nsize[0], Ny = Nsize[1], Nz = Nsize[2], Nt = Nsize[3];
  int Bx = block_size[0], By = block_size[1], Bz = block_size[2], Bt = block_size[3];
  int Nblock = (Nx/Bx)*(Ny/By)*(Nz/Bz)*(Nt/Bt);

  // map host base pointers to device pointers; map_upload caches & skips the
  // H2D when the (persistent) basis/rhs pointer arrays are unchanged.
  static real_t** basis_dev = nullptr; static real_t** rhs_dev = nullptr;
  static int capb = 0, capr = 0;
  map_upload(basis_host, nbasis, basis_dev, capb);
  map_upload(rhs_host,   nrhs,   rhs_dev,   capr);

#ifdef USE_TENSORCORE_MRHS
  if (g_mrhs_use_tc) {
    // one warp (32 threads) per aggregate; the m16n8k8 MMA tiles do the GEMM
    block_restrict_mrhs_tc_dev<<<Nblock, 32>>>(coarse_dev,
                                               basis_dev, nbasis, rhs_dev, nrhs,
                                               Nin, Nex, Nx,Ny,Nz,Nt, Bx,By,Bz,Bt);
  } else {
    int nth = VECTOR_LENGTH;
    int work = Nblock*nbasis*nrhs;
    int nbl  = (work + nth - 1)/nth;
    block_restrict_mrhs_ref_dev<<<nbl, nth>>>(coarse_dev, nullptr,
                                              basis_dev, nbasis, rhs_dev, nrhs,
                                              Nin, Nex, Nx,Ny,Nz,Nt, Bx,By,Bz,Bt);
  }
#else
  int nth = VECTOR_LENGTH;
  int work = Nblock*nbasis*nrhs;
  int nbl  = (work + nth - 1)/nth;
  block_restrict_mrhs_ref_dev<<<nbl, nth>>>(coarse_dev, nullptr,
                                            basis_dev, nbasis, rhs_dev, nrhs,
                                            Nin, Nex, Nx,Ny,Nz,Nt, Bx,By,Bz,Bt);
#endif
  afield_dd_kernel_sync();
}

void block_prolong_mrhs(real_t* const* fine_host, int nrhs,
                        const real_t* coarse_dev,
                        real_t* const* basis_host, int nbasis, real_t fac,
                        int Nin, int Nex, int* Nsize, int* block_size)
{
  int Nx = Nsize[0], Ny = Nsize[1], Nz = Nsize[2], Nt = Nsize[3];
  int Bx = block_size[0], By = block_size[1], Bz = block_size[2], Bt = block_size[3];

  // cached pointer-array upload (skips H2D when basis/fine pointers unchanged).
  static real_t** basis_dev = nullptr; static real_t** fine_dev = nullptr;
  static int capb = 0, capf = 0;
  map_upload(basis_host, nbasis, basis_dev, capb);
  map_upload(fine_host,  nrhs,   fine_dev,  capf);

#ifdef USE_TENSORCORE_MRHS
  if (g_mrhs_use_tc) {
    int Nblock = (Nx/Bx)*(Ny/By)*(Nz/Bz)*(Nt/Bt);
    block_prolong_mrhs_tc_dev<<<Nblock, 32>>>(fine_dev, nrhs, coarse_dev,
                                              basis_dev, nbasis, fac,
                                              Nin, Nex, Nx,Ny,Nz,Nt, Bx,By,Bz,Bt);
  } else {
    int nth = VECTOR_LENGTH;
    int nbl = (Nx*Ny*Nz*Nt + nth - 1)/nth;
    block_prolong_mrhs_ref_dev<<<nbl, nth>>>(fine_dev, nrhs, coarse_dev,
                                             basis_dev, nbasis, fac,
                                             Nin, Nex, Nx,Ny,Nz,Nt, Bx,By,Bz,Bt);
  }
#else
  int nth = VECTOR_LENGTH;
  int nbl = (Nx*Ny*Nz*Nt + nth - 1)/nth;
  block_prolong_mrhs_ref_dev<<<nbl, nth>>>(fine_dev, nrhs, coarse_dev,
                                           basis_dev, nbasis, fac,
                                           Nin, Nex, Nx,Ny,Nz,Nt, Bx,By,Bz,Bt);
#endif
  afield_dd_kernel_sync();
}

// MRHS coarse Dirac operator: V2[r] = D_coarse V1[r] for r = 0..nrhs-1.
// u_dev = coarse gauge (8 directions packed).  Ncol = coarse dof per site.
void mult_coarse_mrhs(real_t* const* v2_host, real_t* const* v1_host,
                      const real_t* u_dev, int nrhs, int Ncol,
                      int* Nsize, int* bc)
{
  int Nx = Nsize[0], Ny = Nsize[1], Nz = Nsize[2], Nt = Nsize[3];
  int Nst = Nx*Ny*Nz*Nt;

  // cached pointer-array upload (skips H2D when v2/v1 pointers unchanged).
  static real_t** v2_dev = nullptr; static real_t** v1_dev = nullptr;
  static int capv2 = 0, capv1 = 0;
  map_upload(v2_host, nrhs, v2_dev, capv2);
  map_upload(v1_host, nrhs, v1_dev, capv1);

#ifdef USE_TENSORCORE_MRHS
  if (g_mrhs_use_tc) {
    mult_coarse_mrhs_tc_dev<<<Nst, 32>>>(v2_dev, v1_dev, u_dev, nrhs, Ncol,
                                         Nx,Ny,Nz,Nt, bc[0],bc[1],bc[2],bc[3]);
  } else {
    int nth = VECTOR_LENGTH;
    int nbl = (Nst*Ncol*nrhs + nth - 1)/nth;
    mult_coarse_mrhs_ref_dev<<<nbl, nth>>>(v2_dev, v1_dev, u_dev, nrhs, Ncol,
                                           Nx,Ny,Nz,Nt, bc[0],bc[1],bc[2],bc[3]);
  }
#else
  int nth = VECTOR_LENGTH;
  int nbl = (Nst*Ncol*nrhs + nth - 1)/nth;
  mult_coarse_mrhs_ref_dev<<<nbl, nth>>>(v2_dev, v1_dev, u_dev, nrhs, Ncol,
                                         Nx,Ny,Nz,Nt, bc[0],bc[1],bc[2],bc[3]);
#endif
  afield_dd_kernel_sync();
}

// MRHS coarse clover-inverse preconditioner: V2[r] = Clov^-1 V1[r] (site-local
// Ncol x Ncol complex GEMM, NO hop, NO identity), r = 0..nrhs-1.  ct = Clov_inv
// in the per-site Ndf layout (same ct that mult_domainwall_coarse_prec reads):
// out[a] = sum_b ct[Ncol2*b + 2*a, site] * v1[b].  One thread per (site,a,r);
// FP32 reference only (a tiny site-local source-preconditioner, no TF32 needed).
__global__
void mult_coarse_prec_mrhs_ref_dev(real_t* const* v2, real_t* const* v1,
                                   const real_t* ct, int nrhs, int Ncol,
                                   int Nx, int Ny, int Nz, int Nt)
{
  int Nst = Nx*Ny*Nz*Nt;
  int Ncol2 = 2*Ncol, Ndf = 2*Ncol*Ncol;
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if(idx >= Nst*Ncol*nrhs) return;
  int r    = idx % nrhs;
  int a    = (idx / nrhs) % Ncol;
  int site = idx / (nrhs*Ncol);

  float accr = 0.0f, acci = 0.0f;
  for(int b = 0; b < Ncol; ++b){
    int in_r = Ncol2*b + 2*a;
    float mr = ct[IDX2(Ndf, in_r,   site)];
    float mi = ct[IDX2(Ndf, in_r+1, site)];
    float vr = v1[r][IDX2(Ncol2, 2*b,   site)];
    float vi = v1[r][IDX2(Ncol2, 2*b+1, site)];
    accr += mr*vr - mi*vi;
    acci += mr*vi + mi*vr;
  }
  v2[r][IDX2(Ncol2, 2*a,   site)] = accr;
  v2[r][IDX2(Ncol2, 2*a+1, site)] = acci;
}

// driver: v2[r] = Clov^-1 v1[r] for all nrhs, one launch.  ct_dev = Clov_inv raw
// device pointer (caller dev_ptr's it); v2_host/v1_host = coarse AField host bases.
void mult_coarse_prec_mrhs(real_t* const* v2_host, real_t* const* v1_host,
                           const real_t* ct_dev, int nrhs, int Ncol, int* Nsize)
{
  int Nx = Nsize[0], Ny = Nsize[1], Nz = Nsize[2], Nt = Nsize[3];
  int Nst = Nx*Ny*Nz*Nt;
  real_t* v2p[256]; real_t* v1p[256];
  if(nrhs > 256){ printf("mult_coarse_prec_mrhs: nrhs > 256\n"); exit(1); }
  for(int r = 0; r < nrhs; ++r){ v2p[r] = (real_t*)dev_ptr(v2_host[r]);
                                 v1p[r] = (real_t*)dev_ptr(v1_host[r]); }
  static real_t** v2d = nullptr; static real_t** v1d = nullptr; static int capp = 0;
  if(nrhs > capp){ if(v2d) cudaFree(v2d); if(v1d) cudaFree(v1d);
    CHECK(cudaMalloc((void**)&v2d, sizeof(real_t*)*nrhs));
    CHECK(cudaMalloc((void**)&v1d, sizeof(real_t*)*nrhs)); capp = nrhs; }
  CHECK(cudaMemcpy(v2d, v2p, sizeof(real_t*)*nrhs, cudaMemcpyHostToDevice));
  CHECK(cudaMemcpy(v1d, v1p, sizeof(real_t*)*nrhs, cudaMemcpyHostToDevice));
  int nth = VECTOR_LENGTH;
  int nbl = (Nst*Ncol*nrhs + nth - 1)/nth;
  mult_coarse_prec_mrhs_ref_dev<<<nbl, nth>>>(v2d, v1d, ct_dev, nrhs, Ncol,
                                              Nx,Ny,Nz,Nt);
  afield_dd_kernel_sync();
}

//====================================================================
// (C) Block-Krylov host drivers.  V_host/W_host = nv full vectors; Nsize is the
// fine 4d lattice, Nin/Nex the field's inner/outer extents (Nex carries Ns for
// 5d).  G_dev = 2*nv*nv complex reals (zeroed inside block_inner_mrhs).

// Upload a field-pointer array to the device, but SKIP the H2D copy when the
// pointer values are unchanged since the last call for this `dev` buffer.  The
// AMG Krylov / Z / scratch fields are persistent, so their device pointers are
// stable across iterations -- re-uploading them every block_inner/block_update/
// A-apply call is exactly the per-iteration host<->device traffic we must avoid.
// Cache is fixed static slots (no heap, no per-call allocation).
static real_t** map_upload(real_t* const* host, int n, real_t**& dev, int& cap)
{
  real_t* p[256];
  if(n > 256){ printf("block-Krylov: nv > 256\n"); exit(1); }
  for(int i = 0; i < n; ++i) p[i] = (real_t*)dev_ptr(host[i]);
  bool grew = false;
  if(n > cap){ if(dev) cudaFree(dev);
    CHECK(cudaMalloc((void**)&dev, sizeof(real_t*)*n)); cap = n; grew = true; }

  static const int NSLOT = 128;
  static real_t** keys[NSLOT] = {0};
  static real_t*  vals[NSLOT][256];
  static int      lens[NSLOT] = {0};
  int slot = -1, freeslot = -1;
  for(int c = 0; c < NSLOT; ++c){ if(keys[c] == dev){ slot = c; break; }
                                  if(freeslot < 0 && keys[c] == 0) freeslot = c; }
  bool same = (slot >= 0 && !grew && lens[slot] == n);
  for(int i = 0; same && i < n; ++i) if(vals[slot][i] != p[i]) same = false;
  if(!same){
    CHECK(cudaMemcpy(dev, p, sizeof(real_t*)*n, cudaMemcpyHostToDevice));
    if(slot < 0) slot = (freeslot >= 0) ? freeslot : 0;
    keys[slot] = dev; lens[slot] = n;
    for(int i = 0; i < n; ++i) vals[slot][i] = p[i];
  }
  return dev;
}

// G = V^H W  (nv x nv complex Gram matrix), reduced over the volume.
void block_inner_mrhs(real_t* G_dev, real_t* const* V_host, real_t* const* W_host,
                      int nv, int Nin, int Nex, int* Nsize)
{
  int Nst = Nsize[0]*Nsize[1]*Nsize[2]*Nsize[3];
  static real_t **Vd = nullptr, **Wd = nullptr; static int cV = 0, cW = 0;
  real_t** Vdev = map_upload(V_host, nv, Vd, cV);
  real_t** Wdev = map_upload(W_host, nv, Wd, cW);
  CHECK(cudaMemset(G_dev, 0, sizeof(real_t)*2*nv*nv));
#ifdef USE_TENSORCORE_MRHS
  if (g_mrhs_use_tc) {
    int Ktot = Nex*Nst*(Nin/2);
    int nbl = (Ktot + MMA_K - 1)/MMA_K;
    block_inner_tc_dev<<<nbl, 32>>>(G_dev, Vdev, Wdev, nv, Nin, Nex, Nst);
  } else {
    int nth = VECTOR_LENGTH; int nbl = (nv*nv + nth - 1)/nth;
    block_inner_ref_dev<<<nbl, nth>>>(G_dev, Vdev, Wdev, nv, Nin, Nex, Nst);
  }
#else
  int nth = VECTOR_LENGTH; int nbl = (nv*nv + nth - 1)/nth;
  block_inner_ref_dev<<<nbl, nth>>>(G_dev, Vdev, Wdev, nv, Nin, Nex, Nst);
#endif
  afield_dd_kernel_sync();
}

// W_j -= sum_i V_i G[i,j]  (block-Arnoldi orthogonalisation update).
void block_update_mrhs(real_t* const* W_host, real_t* const* V_host, const real_t* G_dev,
                       int nv, int Nin, int Nex, int* Nsize)
{
  int Nst = Nsize[0]*Nsize[1]*Nsize[2]*Nsize[3];
  int Ktot = Nex*Nst*(Nin/2);
  static real_t **Vd = nullptr, **Wd = nullptr; static int cV = 0, cW = 0;
  real_t** Vdev = map_upload(V_host, nv, Vd, cV);
  real_t** Wdev = map_upload(W_host, nv, Wd, cW);
#ifdef USE_TENSORCORE_MRHS
  if (g_mrhs_use_tc) {
    int nbl = (Ktot + MMA_M - 1)/MMA_M;
    block_update_tc_dev<<<nbl, 32>>>(Wdev, Vdev, G_dev, nv, Nin, Nex, Nst);
  } else {
    int nth = VECTOR_LENGTH; int nbl = (Ktot*nv + nth - 1)/nth;
    block_update_ref_dev<<<nbl, nth>>>(Wdev, Vdev, G_dev, nv, Nin, Nex, Nst);
  }
#else
  int nth = VECTOR_LENGTH; int nbl = (Ktot*nv + nth - 1)/nth;
  block_update_ref_dev<<<nbl, nth>>>(Wdev, Vdev, G_dev, nv, Nin, Nex, Nst);
#endif
  afield_dd_kernel_sync();
}

// ---------------------------------------------------------------------------
// On-device s x s complex dense LA for the block-Krylov solver, so the Gram /
// Cholesky / Hessenberg NEVER round-trip to the host inside the iteration.
// G (Hermitian PD, 2*s*s floats, G[i][j] at 2*(i*s+j) re / +1 im) -> upper L
// (G = L^H L) and negLi = -L^{-1}.  Run as ONE block of s threads cooperating in
// shared memory; mirrors the host choleskyU/invU exactly so results match.
//
//  - Cholesky (upper) is column-sequential (column j needs all k<j), but inside
//    each column the s-j off-diagonal entries are computed in parallel, one per
//    thread, with a __syncthreads() barrier between columns.
//  - The triangular inverse parallelises the other way: the columns of M are
//    independent, so thread j back-substitutes the whole of column j on its own
//    (no inter-thread dependence), no barriers needed inside the phase.
#ifndef BLK_DENSE_MAXS
#define BLK_DENSE_MAXS 32
#endif
// One thread block of s*s threads: thread (i,j) = tid = j*s + i owns matrix
// entry (i,j).  The s entries below a diagonal share the work of that column's
// dot-product reduction (each contributes one term via shared partials), so the
// inner sums are parallel too, not a per-thread serial loop.
__global__
void block_chol_inv_dev(real_t* Lout, real_t* negLi, const real_t* G, int s)
{
  __shared__ float Lr[BLK_DENSE_MAXS*BLK_DENSE_MAXS];
  __shared__ float Li[BLK_DENSE_MAXS*BLK_DENSE_MAXS];
  __shared__ float Mr[BLK_DENSE_MAXS*BLK_DENSE_MAXS];
  __shared__ float Mi[BLK_DENSE_MAXS*BLK_DENSE_MAXS];
  __shared__ float pr[BLK_DENSE_MAXS*BLK_DENSE_MAXS];   // per-(k,i) partial re
  __shared__ float pi[BLK_DENSE_MAXS*BLK_DENSE_MAXS];   // per-(k,i) partial im
  const int tid = threadIdx.x;
  const int ii  = tid % s;          // row    of this thread's entry
  const int jj  = tid / s;          // column of this thread's entry
  for (int t = tid; t < s*s; t += blockDim.x) { Lr[t]=0.f; Li[t]=0.f; Mr[t]=0.f; Mi[t]=0.f; }
  __syncthreads();

  // Cholesky (upper): G = L^H L, column by column.  For column j, thread (k,i)
  // with k<j and i>=j computes one product term of the (j,i) reduction.
  for (int j = 0; j < s; ++j) {
    // each thread (k=jj, i=ii) with k<j, i>=j stages its product term.
    if (jj < j && ii >= j) {
      const int kk = jj;
      if (ii == j) {                             // diagonal term: |L[k][j]|^2
        float a=Lr[kk*s+j], b=Li[kk*s+j];
        pr[kk*s+ii] = a*a + b*b; pi[kk*s+ii] = 0.f;
      } else {                                   // off-diag term: conj(L[k][j])*L[k][i]
        float ar=Lr[kk*s+j], ai=Li[kk*s+j], br=Lr[kk*s+ii], bi=Li[kk*s+ii];
        pr[kk*s+ii] = ar*br + ai*bi;
        pi[kk*s+ii] = ar*bi - ai*br;
      }
    }
    __syncthreads();
    // diagonal FIRST (single thread): L[j][j] = sqrt(G[j][j] - sum_k |L[k][j]|^2).
    if (jj == 0 && ii == j) {
      float sr = G[2*(j*s+j)];
      for (int k = 0; k < j; ++k) sr -= pr[k*s+j];
      if (sr < 1.0e-30f) sr = 1.0e-30f;
      Lr[j*s+j] = sqrtf(sr); Li[j*s+j] = 0.f;
    }
    __syncthreads();                              // L[j][j] visible before /Ljj below
    // off-diagonals i>j (one thread each) -- now safe to read L[j][j].
    if (jj == 0 && ii > j) {
      float sr = G[2*(j*s+ii)], si = G[2*(j*s+ii)+1];
      for (int k = 0; k < j; ++k) { sr -= pr[k*s+ii]; si -= pi[k*s+ii]; }
      float Ljj = Lr[j*s+j];
      Lr[j*s+ii] = sr/Ljj; Li[j*s+ii] = si/Ljj;
    }
    __syncthreads();
  }

  // inverse: M = L^{-1} (upper triangular); thread (0,j) owns column j of M.
  if (jj == 0 && ii < s) {
    const int j = ii;
    Mr[j*s+j] = 1.f/Lr[j*s+j]; Mi[j*s+j] = 0.f;
    for (int i = j-1; i >= 0; --i) {
      float sr = 0.f, si = 0.f;
      for (int k = i+1; k <= j; ++k) {           // += L[i][k] * M[k][j]
        float ar=Lr[i*s+k], ai=Li[i*s+k], br=Mr[k*s+j], bi=Mi[k*s+j];
        sr += ar*br - ai*bi;
        si += ar*bi + ai*br;
      }
      float d = Lr[i*s+i];                        // L[i][i] real
      Mr[i*s+j] = -sr/d; Mi[i*s+j] = -si/d;
    }
  }
  __syncthreads();

  for (int t = tid; t < s*s; t += blockDim.x) {
    Lout[2*t]   = Lr[t];  Lout[2*t+1]  = Li[t];
    negLi[2*t]  = -Mr[t]; negLi[2*t+1] = -Mi[t];
  }
}

// driver: all pointers are RAW device buffers (caller dev_alloc'd) -- no dev_ptr,
// no host transfer.  Lout/negLi/G are each 2*s*s floats.  ONE block (a single
// Gram per Arnoldi step, sequentially dependent), s*s threads inside it.
void block_chol_inv(real_t* Lout_dev, real_t* negLi_dev, const real_t* G_dev, int s)
{
  if (s > BLK_DENSE_MAXS) { printf("block_chol_inv: s=%d > BLK_DENSE_MAXS\n", s); exit(1); }
  block_chol_inv_dev<<<1, s*s>>>(Lout_dev, negLi_dev, G_dev, s);
  afield_dd_kernel_sync();
}

// ---------------------------------------------------------------------------
// Batched BLAS-1 over s columns.  Each field is nflt contiguous (re,im) floats;
// copy/axpy/scal are element-wise and norm2 is a sum of squares, all invariant
// under the NWP-tiled layout (a bijection of the nflt floats), so one flat
// kernel over all s columns is bit-identical to s per-column AField ops -- and
// replaces s kernel launches with ONE.  This is the block BLAS-1 layer the
// nested block-Krylov was missing (per-column copy/axpy/norm = the launch-bound
// floor on a small lattice).  All pointers are host base ptrs (dev_ptr-mapped).
__global__
void mrhs_copy_dev(real_t* const* y, real_t* const* x, int s, long nflt)
{
  const long stride = (long)gridDim.x * blockDim.x, tot = nflt * (long)s;
  for (long t = (long)blockIdx.x*blockDim.x + threadIdx.x; t < tot; t += stride) {
    int c = (int)(t / nflt); long i = t - (long)c*nflt; y[c][i] = x[c][i];
  }
}
__global__
void mrhs_axpy_dev(real_t* const* y, real_t* const* x, real_t a, int s, long nflt)
{
  const long stride = (long)gridDim.x * blockDim.x, tot = nflt * (long)s;
  for (long t = (long)blockIdx.x*blockDim.x + threadIdx.x; t < tot; t += stride) {
    int c = (int)(t / nflt); long i = t - (long)c*nflt; y[c][i] += a * x[c][i];
  }
}
__global__
void mrhs_scal_dev(real_t* const* y, real_t a, int s, long nflt)
{
  const long stride = (long)gridDim.x * blockDim.x, tot = nflt * (long)s;
  for (long t = (long)blockIdx.x*blockDim.x + threadIdx.x; t < tot; t += stride) {
    int c = (int)(t / nflt); long i = t - (long)c*nflt; y[c][i] *= a;
  }
}
// one block per column; block-stride sum of squares -> out[c] (double accum).
__global__
void mrhs_norm2_dev(double* out, real_t* const* x, int s, long nflt)
{
  const int c = blockIdx.x;
  if (c >= s) return;
  const real_t* xc = x[c];
  double loc = 0.0;
  for (long i = threadIdx.x; i < nflt; i += blockDim.x) { double v = xc[i]; loc += v*v; }
  __shared__ double sh[256];
  sh[threadIdx.x] = loc; __syncthreads();
  for (int k = blockDim.x/2; k > 0; k >>= 1) {
    if (threadIdx.x < k) sh[threadIdx.x] += sh[threadIdx.x+k];
    __syncthreads();
  }
  if (threadIdx.x == 0) out[c] = sh[0];
}

static int mrhs_blas1_blocks(long tot, int nth)
{ long nbl = (tot + nth - 1)/nth; return (int)(nbl > 65535 ? 65535 : (nbl < 1 ? 1 : nbl)); }

// y[c] = x[c]                       (batched copy, all s columns, one launch)
void mrhs_copy(real_t* const* y_host, real_t* const* x_host, int s, long nflt)
{
  static real_t **yd=nullptr,**xd=nullptr; static int cy=0,cx=0;
  real_t** y = map_upload(y_host, s, yd, cy);
  real_t** x = map_upload(x_host, s, xd, cx);
  int nth=256; mrhs_copy_dev<<<mrhs_blas1_blocks(nflt*(long)s,nth),nth>>>(y,x,s,nflt);
  afield_dd_kernel_sync();
}
// y[c] += a * x[c]                  (batched axpy, shared scalar a)
void mrhs_axpy(real_t* const* y_host, real_t* const* x_host, real_t a, int s, long nflt)
{
  static real_t **yd=nullptr,**xd=nullptr; static int cy=0,cx=0;
  real_t** y = map_upload(y_host, s, yd, cy);
  real_t** x = map_upload(x_host, s, xd, cx);
  int nth=256; mrhs_axpy_dev<<<mrhs_blas1_blocks(nflt*(long)s,nth),nth>>>(y,x,a,s,nflt);
  afield_dd_kernel_sync();
}
// y[c] *= a                         (batched scal, shared scalar a)
void mrhs_scal(real_t* const* y_host, real_t a, int s, long nflt)
{
  static real_t **yd=nullptr; static int cy=0;
  real_t** y = map_upload(y_host, s, yd, cy);
  int nth=256; mrhs_scal_dev<<<mrhs_blas1_blocks(nflt*(long)s,nth),nth>>>(y,a,s,nflt);
  afield_dd_kernel_sync();
}
// out_host[c] = ||x[c]||^2          (batched norm2, s norms -> host once)
void mrhs_norm2(double* out_host, real_t* const* x_host, int s, long nflt)
{
  static real_t **xd=nullptr; static int cx=0;
  static double* od=nullptr;  static int co=0;
  real_t** x = map_upload(x_host, s, xd, cx);
  if (s > co) { if(od) cudaFree(od); CHECK(cudaMalloc((void**)&od, sizeof(double)*s)); co=s; }
  mrhs_norm2_dev<<<s, 256>>>(od, x, s, nflt);
  CHECK(cudaMemcpy(out_host, od, sizeof(double)*s, cudaMemcpyDeviceToHost));
  afield_dd_kernel_sync();
}

} // namespace

#endif // MRHS_BLOCK_TENSORCORE_CUDA_INCLUDED
