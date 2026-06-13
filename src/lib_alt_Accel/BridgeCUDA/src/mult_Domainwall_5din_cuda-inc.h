/*!
      @file    mult_Domainwall_5din_cuda-inc.h
      @brief
      @author  Wei-Lun Chen (wlchen)
      @date    $LastChangedDate: 2026-06-09 16:42:53 #$
      @version $LastChangedRevision: 2606 $
*/

#ifndef MULT_DOMAINWALL_5DIN_ACC_INCLUDED
#define MULT_DOMAINWALL_5DIN_ACC_INCLUDED

#include "../inline/constant_memory_inline.h"

//====================================================================
// Launch-overhead reduction.
//
// The compute kernels in this file run on the legacy default stream (stream 0),
// so consecutive launches are implicitly ordered and their data dependencies are
// respected WITHOUT a host-side cudaDeviceSynchronize() between them. Every point
// where the host actually reads device data -- reductions (norm2/dotc) and the
// field->host convert/gramschmidt -- goes through a *synchronizing* cudaMemcpy,
// which drains the default-stream queue and provides the real barrier. The only
// mid-operator host read is the MPI halo path (chset_send/recv), guarded by
// do_comm_any>0; the syncs in the hop1/hop2 wrappers (which pack/consume the MPI
// buffers) are therefore KEPT for multi-rank correctness.
//
// A per-kernel cudaDeviceSynchronize() here only forces the host to BLOCK instead
// of queuing the next launch, which serializes launch+execute and starves the GPU
// (measured: GPU ~35% util, one host core pegged at 100% spinning) -- i.e. the
// domain-wall operator apply was launch-latency bound on small/medium local
// volumes. Dropping it lets the host run ahead and pipeline the launches.
//
// Set DW5DIN_SYNC_EACH_KERNEL=1 to restore the old per-kernel sync (e.g. to
// localize an async kernel error during debugging).
#ifndef DW5DIN_SYNC_EACH_KERNEL
#define DW5DIN_SYNC_EACH_KERNEL 0
#endif
static inline void dw5din_kernel_sync()
{
#if DW5DIN_SYNC_EACH_KERNEL
  CHECK(cudaDeviceSynchronize());
#endif
}

//====================================================================
__global__
void mult_domainwall_5din_5dir_dirac_kernel(
    real_t * __restrict__ vp, real_t * __restrict__ yp,
    real_t * __restrict__ wp,
    real_t mq, real_t M0, int Ns, real_t alpha, int Nst_pad) {

  // Forward 5d block of the (non-eo) Moebius domain-wall operator.
  // Produces vp = D_5 w (diagonal + 5d-Wilson hop) and the precond auxiliary
  // yp = -0.5 (b*self + c*hop) used by the 4d hopping.  Color/Re-Im (ivc) is
  // parallelized over threads (one (site,ivc) per thread, NWP-coalesced) and
  // b/c are read from __constant__ memory; matches the eo/dd kernel layout.
  //
  // Neff alpha (bulk-vs-boundary 5D scaling): the 5d hop gets 0.5*alpha in the
  // bulk (mass-boundary terms keep -0.5*mq), and the self block becomes
  // alpha*w in the bulk but a chiral (gamma5) mix 0.5(1+alpha) +/- 0.5(-+1+alpha)
  // at the s=0 / s=Ns-1 corners.  alpha = 1 reduces to the standard operator.

  const int Nin5     = NVCD * Ns;
  const int ist      = blockIdx.x * blockDim.x + threadIdx.x;
  const int GridSize = blockDim.x * gridDim.x;

  const real_t *b_con = const_b;
  const real_t *c_con = const_c;

  for (int idx = ist; idx < Nst_pad * NVC; idx += GridSize) {

    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int ivc     = idx2_wp % NVC;
    int idx_out = idx2_wp / NVC;
    int site    = idx_in + NWP * idx_out;

    for (int is = 0; is < Ns; ++is) {

      real_t B_is = b_con[is] * (4.0 - M0) + 1.0;
      real_t C_is = c_con[is] * (4.0 - M0) - 1.0;

      real_t vt1, vt2, vt3, vt4;
      real_t wt1, wt2, wt3, wt4;

      int is_up = (is + 1) % Ns;
      real_t Fup = 0.5 * alpha;
      if (is == Ns - 1) Fup = -0.5 * mq;
      wt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD * is_up), site)];
      wt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD * is_up), site)];
      wt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD * is_up), site)];
      wt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD * is_up), site)];
      vt1 = Fup * (wt1 - wt3);
      vt2 = Fup * (wt2 - wt4);
      vt3 = Fup * (wt3 - wt1);
      vt4 = Fup * (wt4 - wt2);

      int is_dn = (is - 1 + Ns) % Ns;
      real_t Fdn = 0.5 * alpha;
      if (is == 0) Fdn = -0.5 * mq;
      wt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD * is_dn), site)];
      wt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD * is_dn), site)];
      wt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD * is_dn), site)];
      wt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD * is_dn), site)];
      vt1 += Fdn * (wt1 + wt3);
      vt2 += Fdn * (wt2 + wt4);
      vt3 += Fdn * (wt3 + wt1);
      vt4 += Fdn * (wt4 + wt2);

      wt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)];
      wt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)];
      wt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)];
      wt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)];

      real_t st1, st2, st3, st4;
      if (is == 0) {
        real_t f1 = 0.5 * (1.0 + alpha);
        real_t f2 = 0.5 * (-1.0 + alpha);
        st1 = f1 * wt1 + f2 * wt3;
        st2 = f1 * wt2 + f2 * wt4;
        st3 = f1 * wt3 + f2 * wt1;
        st4 = f1 * wt4 + f2 * wt2;
      } else if (is == Ns - 1) {
        real_t f1 = 0.5 * (1.0 + alpha);
        real_t f2 = 0.5 * (1.0 - alpha);
        st1 = f1 * wt1 + f2 * wt3;
        st2 = f1 * wt2 + f2 * wt4;
        st3 = f1 * wt3 + f2 * wt1;
        st4 = f1 * wt4 + f2 * wt2;
      } else {
        st1 = alpha * wt1;
        st2 = alpha * wt2;
        st3 = alpha * wt3;
        st4 = alpha * wt4;
      }

      real_t b_is = b_con[is];
      real_t c_is = c_con[is];
      vp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)] = B_is * st1 + C_is * vt1;
      vp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)] = B_is * st2 + C_is * vt2;
      vp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)] = B_is * st3 + C_is * vt3;
      vp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)] = B_is * st4 + C_is * vt4;
      yp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)] = -0.5 * (b_is * st1 + c_is * vt1);
      yp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)] = -0.5 * (b_is * st2 + c_is * vt2);
      yp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)] = -0.5 * (b_is * st3 + c_is * vt3);
      yp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)] = -0.5 * (b_is * st4 + c_is * vt4);
    }
  }
}


void mult_domainwall_5din_5dir_dirac(
    real_t *vp, real_t *yp, real_t *wp,
    real_t mq, real_t M0, int Ns, real_t *b, real_t *c, int *Nsize,
    real_t alpha) {
  // b/c are read from __constant__ memory (const_b/const_c) inside the kernel;
  // the b/c pointer args are kept only for backward signature compatibility.
  (void)b; (void)c;

  int Nst = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];
  int Nst_pad = ceil_nwp(Nst);

  real_t* vp_dev = (real_t*)dev_ptr(vp);
  real_t* yp_dev = (real_t*)dev_ptr(yp);
  real_t* wp_dev = (real_t*)dev_ptr(wp);

  int blockSize = VECTOR_LENGTH;
  int gridSize  = (Nst_pad * NVC + blockSize - 1) / blockSize;
  mult_domainwall_5din_5dir_dirac_kernel<<<gridSize, blockSize>>>(
      vp_dev, yp_dev, wp_dev, mq, M0, Ns, alpha, Nst_pad);

  dw5din_kernel_sync();
}

//====================================================================
__global__
void mult_domainwall_5din_5dirdag_dirac_kernel(
    real_t * __restrict__ vp, real_t * __restrict__ yp,
    real_t * __restrict__ wp,
    real_t mq, real_t M0, int Ns, real_t alpha, int Nst_pad) {

  // Adjoint (dagger) of the forward 5d block.  Takes both wp and the precond
  // auxiliary yp as inputs.  Color/Re-Im parallel, const-memory b/c, same Neff
  // alpha congruence as the forward kernel (self block chirally mixed at the
  // s corners, 5d hop scaled by alpha in the bulk).  alpha = 1 -> standard.

  const int Nin5     = NVCD * Ns;
  const int ist      = blockIdx.x * blockDim.x + threadIdx.x;
  const int GridSize = blockDim.x * gridDim.x;

  const real_t *b_con = const_b;
  const real_t *c_con = const_c;

  for (int idx = ist; idx < Nst_pad * NVC; idx += GridSize) {

    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int ivc     = idx2_wp % NVC;
    int idx_out = idx2_wp / NVC;
    int site    = idx_in + NWP * idx_out;

    for (int is = 0; is < Ns; ++is) {

      real_t B1 = b_con[is] * (4.0 - M0) + 1.0;
      real_t a1 = -0.5 * b_con[is];

      real_t wt1, wt2, wt3, wt4;
      real_t yt1, yt2, yt3, yt4;

      // self block: S = B1*w + a1*gamma5(y), then alpha congruence.
      wt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)];
      wt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)];
      wt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)];
      wt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)];
      yt1 = yp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)];
      yt2 = yp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)];
      yt3 = yp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)];
      yt4 = yp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)];
      real_t s1 = B1 * wt1 + a1 * yt3;
      real_t s2 = B1 * wt2 + a1 * yt4;
      real_t s3 = B1 * wt3 + a1 * yt1;
      real_t s4 = B1 * wt4 + a1 * yt2;

      real_t vt1, vt2, vt3, vt4;
      if (is == 0) {
        real_t f1 = 0.5 * (1.0 + alpha);
        real_t f2 = 0.5 * (-1.0 + alpha);
        vt1 = f1 * s1 + f2 * s3;
        vt2 = f1 * s2 + f2 * s4;
        vt3 = f1 * s3 + f2 * s1;
        vt4 = f1 * s4 + f2 * s2;
      } else if (is == Ns - 1) {
        real_t f1 = 0.5 * (1.0 + alpha);
        real_t f2 = 0.5 * (1.0 - alpha);
        vt1 = f1 * s1 + f2 * s3;
        vt2 = f1 * s2 + f2 * s4;
        vt3 = f1 * s3 + f2 * s1;
        vt4 = f1 * s4 + f2 * s2;
      } else {
        vt1 = alpha * s1;
        vt2 = alpha * s2;
        vt3 = alpha * s3;
        vt4 = alpha * s4;
      }

      int is_up = (is + 1) % Ns;
      real_t C1  = c_con[is_up] * (4.0 - M0) - 1.0;
      real_t aup = -0.5 * c_con[is_up];
      wt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD * is_up), site)];
      wt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD * is_up), site)];
      wt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD * is_up), site)];
      wt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD * is_up), site)];
      yt1 = yp[IDX2(Nin5, (ID1 + ivc + NVCD * is_up), site)];
      yt2 = yp[IDX2(Nin5, (ID2 + ivc + NVCD * is_up), site)];
      yt3 = yp[IDX2(Nin5, (ID3 + ivc + NVCD * is_up), site)];
      yt4 = yp[IDX2(Nin5, (ID4 + ivc + NVCD * is_up), site)];
      real_t xu1 = C1 * wt1 + aup * yt3;
      real_t xu2 = C1 * wt2 + aup * yt4;
      real_t xu3 = C1 * wt3 + aup * yt1;
      real_t xu4 = C1 * wt4 + aup * yt2;
      real_t Fup = 0.5 * alpha;
      if (is == Ns - 1) Fup = -0.5 * mq;
      vt1 += Fup * (xu1 + xu3);
      vt2 += Fup * (xu2 + xu4);
      vt3 += Fup * (xu3 + xu1);
      vt4 += Fup * (xu4 + xu2);

      int is_dn = (is - 1 + Ns) % Ns;
      real_t C2  = c_con[is_dn] * (4.0 - M0) - 1.0;
      real_t adn = -0.5 * c_con[is_dn];
      wt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD * is_dn), site)];
      wt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD * is_dn), site)];
      wt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD * is_dn), site)];
      wt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD * is_dn), site)];
      yt1 = yp[IDX2(Nin5, (ID1 + ivc + NVCD * is_dn), site)];
      yt2 = yp[IDX2(Nin5, (ID2 + ivc + NVCD * is_dn), site)];
      yt3 = yp[IDX2(Nin5, (ID3 + ivc + NVCD * is_dn), site)];
      yt4 = yp[IDX2(Nin5, (ID4 + ivc + NVCD * is_dn), site)];
      real_t xd1 = C2 * wt1 + adn * yt3;
      real_t xd2 = C2 * wt2 + adn * yt4;
      real_t xd3 = C2 * wt3 + adn * yt1;
      real_t xd4 = C2 * wt4 + adn * yt2;
      real_t Fdn = 0.5 * alpha;
      if (is == 0) Fdn = -0.5 * mq;
      vt1 += Fdn * (xd1 - xd3);
      vt2 += Fdn * (xd2 - xd4);
      vt3 += Fdn * (xd3 - xd1);
      vt4 += Fdn * (xd4 - xd2);

      vp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)] = vt1;
      vp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)] = vt2;
      vp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)] = vt3;
      vp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)] = vt4;
    }
  }
}

void mult_domainwall_5din_5dirdag_dirac(
    real_t *vp, real_t *yp, real_t *wp,
    real_t mq, real_t M0, int Ns, real_t *b, real_t *c, int *Nsize,
    real_t alpha) {
  // b/c from __constant__ memory; pointer args kept for ABI compatibility.
  (void)b; (void)c;

  int Nst = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];
  int Nst_pad = ceil_nwp(Nst);

  real_t* vp_dev = (real_t*)dev_ptr(vp);
  real_t* yp_dev = (real_t*)dev_ptr(yp);
  real_t* wp_dev = (real_t*)dev_ptr(wp);

  int blockSize = VECTOR_LENGTH;
  int gridSize = (Nst_pad * NVC + blockSize - 1) / blockSize;

  mult_domainwall_5din_5dirdag_dirac_kernel<<<gridSize, blockSize>>>(
      vp_dev, yp_dev, wp_dev, mq, M0, Ns, alpha, Nst_pad);

  dw5din_kernel_sync();
}

//====================================================================
__global__ void mult_domainwall_5din_mult_gm5_dirac_kernel(
    real_t *vp, real_t *wp, int Ns, int Nst) {

    const int Nin5 = NVCD * Ns;
    const int site = blockIdx.x * blockDim.x + threadIdx.x;

    if (site < Nst) {
        for (int is = 0; is < Ns; ++is) {

            real_t vt[NVCD], wt[NVCD]; 

            for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
                wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
            }

            for (int ivc = 0; ivc < NVC; ++ivc) {
                vt[ID1 + ivc] = wt[ID3 + ivc];
                vt[ID2 + ivc] = wt[ID4 + ivc];
                vt[ID3 + ivc] = wt[ID1 + ivc];
                vt[ID4 + ivc] = wt[ID2 + ivc];
            }

            for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
                vp[IDX2(Nin5, (ivcd + NVCD * is), site)] = vt[ivcd];
            }
        }
    }
}

void mult_domainwall_5din_mult_gm5_dirac(
         real_t *vp, real_t *wp, int Ns, int *Nsize) {

    int Nx = Nsize[0];
    int Ny = Nsize[1];
    int Nz = Nsize[2];
    int Nt = Nsize[3];
    int Nst = Nx * Ny * Nz * Nt; 

    real_t* vp_dev = (real_t*)dev_ptr(vp);
    real_t* wp_dev = (real_t*)dev_ptr(wp);

    int blockSize = VECTOR_LENGTH;
    int gridSize = (Nst + blockSize - 1)/ blockSize;

    mult_domainwall_5din_mult_gm5_dirac_kernel<<<gridSize, blockSize>>>(
        vp_dev, wp_dev, Ns, Nst);

    dw5din_kernel_sync();
}


//====================================================================
__global__ void mult_domainwall_5din_mult_R_kernel(
         real_t *vp, real_t *wp, int Ns, int Nst) {

    int Nin5 = NVCD * Ns; 
    int site = blockIdx.x * blockDim.x + threadIdx.x; 

    if (site < Nst) {
        for (int is = 0; is < Ns; ++is) {

            real_t vt[NVCD]; 
            int isR = Ns - 1 - is; 

            for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
                vt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * isR), site)];
            }

            for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
                vp[IDX2(Nin5, (ivcd + NVCD * is), site)] = vt[ivcd];
            }
        }
    }
}

void mult_domainwall_5din_mult_R(
         real_t *vp, real_t *wp, int Ns, int *Nsize) {

    int Nx = Nsize[0];
    int Ny = Nsize[1];
    int Nz = Nsize[2];
    int Nt = Nsize[3];
    int Nst = Nx * Ny * Nz * Nt; 

    real_t* vp_dev = (real_t*)dev_ptr(vp);
    real_t* wp_dev = (real_t*)dev_ptr(wp);

    int blockSize = VECTOR_LENGTH;
    int gridSize = (Nst + blockSize - 1)/ blockSize;

    mult_domainwall_5din_mult_R_kernel<<<gridSize, blockSize>>>(
        vp_dev, wp_dev, Ns, Nst);

    dw5din_kernel_sync();
}
//----------------------------------------------------------------
__global__ void mult_domainwall_5din_mult_gm5R_dirac_kernel(
         real_t *vp, real_t *wp, int Ns, int Nx, int Ny, int Nz, int Nt) {
    
    int site = blockIdx.x * blockDim.x + threadIdx.x;
    int Nst = Nx * Ny * Nz * Nt;
    int Nin5 = NVCD * Ns;

    if (site < Nst) {
        for (int is = 0; is < Ns; ++is) {

            int isR = Ns - 1 - is;
            real_t vt[NVCD], wt[NVCD];

            for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
                wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * isR), site)];
            }

            for (int ivc = 0; ivc < NVC; ++ivc) {
                vt[ID1 + ivc] = wt[ID3 + ivc];
                vt[ID2 + ivc] = wt[ID4 + ivc];
                vt[ID3 + ivc] = wt[ID1 + ivc];
                vt[ID4 + ivc] = wt[ID2 + ivc];
            }

            for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
                vp[IDX2(Nin5, (ivcd + NVCD * is), site)] = vt[ivcd];
            }
        }
    }
}

void mult_domainwall_5din_mult_gm5R_dirac(
         real_t *vp, real_t *wp, int Ns, int *Nsize) {

    int Nx = Nsize[0];
    int Ny = Nsize[1];
    int Nz = Nsize[2];
    int Nt = Nsize[3];
    int Nst = Nx * Ny * Nz * Nt;

    real_t* vp_dev = (real_t*)dev_ptr(vp);
    real_t* wp_dev = (real_t*)dev_ptr(wp);

    int blockSize = VECTOR_LENGTH;
    int gridSize = (Nst + blockSize - 1)/ blockSize;

    mult_domainwall_5din_mult_gm5R_dirac_kernel<<<gridSize, blockSize>>>(
        vp_dev, wp_dev, Ns, Nx, Ny, Nz, Nt);

    dw5din_kernel_sync();
}

//====================================================================
__global__ void mult_domainwall_5din_clear_kernel(real_t *vp, int Ns, int Nst) {

    // Set each element of vp to 0.0
    
    int site = blockIdx.x * blockDim.x + threadIdx.x; 
    int Nin5 = NVCD * Ns;

    if (site < Nst) {
        for (int is = 0; is < Ns; ++is) {
            for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
                vp[IDX2(Nin5, (ivcd + NVCD * is), site)] = 0.0;
            }
        }
    }
}

void mult_domainwall_5din_clear(real_t *vp, int Ns, int *Nsize) {

    int Nx = Nsize[0];
    int Ny = Nsize[1];
    int Nz = Nsize[2];
    int Nt = Nsize[3];
    int Nst = Nx * Ny * Nz * Nt; 

    real_t* vp_dev = (real_t*)dev_ptr(vp);

    int blockSize = VECTOR_LENGTH;
    int gridSize = (Nst + blockSize - 1)/ blockSize;

    mult_domainwall_5din_clear_kernel<<<gridSize, blockSize>>>(vp_dev, Ns, Nst);

    dw5din_kernel_sync();
}

//====================================================================
__global__ 
void mult_domainwall_5din_hopb_4D_dirac_kernel_dir(
     real_t *vp, real_t *up, real_t *wp, int Ns, 
     int bc_x, int bc_y, int bc_z, int bc_t,
     int Nx, int Ny, int Nz, int Nt,
     int do_comm_x, int do_comm_y, int do_comm_z, int do_comm_t,
     int flag) 
{
    int Nxy   = Nx   * Ny;
    int Nxyz  = Nxy  * Nz;
    int Nst   = Nxyz * Nt;
    int Nin5  = NVCD * Ns; 
    int site = blockIdx.x * blockDim.x + threadIdx.x;
    extern __shared__ real_t sharedMemory[];
    real_t* ut = &sharedMemory[0];

    if (site < Nst) {

        int ix    = site % Nx;
        int iyzt  = site / Nx;
        int ixy   = site % Nxy;
        int iy    = iyzt % Ny;
        int izt   = site / Nxy;
        int iz    = izt  % Nz;
        int it    = izt  / Nz;
        int ixyz  = site % Nxyz;
        int idir;

        for(int is = 0; is < Ns; ++is) {
            real_t vL[NVCD];
            for(int ivcd = 0; ivcd < NVCD; ++ivcd) {
                vL[ivcd] = (flag == 0) ? 0.0 : vp[IDX2(Nin5, (ivcd + NVCD * is), site)];
            }

            // X direction computations
            idir = 0;

            // Positive X direction computation
            if ((ix < Nx - 1) || (do_comm_x == 0)) {
                int ix2 = (ix + 1) % Nx;
                int nei = ix2 + Nx * iyzt;

                real_t bc2 = (ix == Nx - 1) ? bc_x : 1.0;
                load_u(ut, up, site + Nst * idir);
                real_t wt[NVCD];

                for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                    wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD * is), nei)];
                }
                mult_wilson_xpb(vL, ut, wt);
            }

            // Negative X direction computation
            if ((ix > 0) || (do_comm_x == 0)) {
                int ix2 = (ix - 1 + Nx) % Nx;
                int nei = ix2 + Nx * iyzt;

                real_t bc2 = (ix == 0) ? bc_x : 1.0;
                load_u(ut, up, nei + Nst * idir);
                real_t wt[NVCD];

                for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                    wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD * is), nei)];
                }
                mult_wilson_xmb(vL, ut, wt);
            }

            // Y direction computations
            idir = 1;

            if ((iy < Ny-1) || (do_comm_y == 0)) {
              int iy2 = (iy + 1) % Ny;
              int nei = ix + Nx * (iy2 + Ny * izt);

              real_t bc2 = (iy == Ny-1) ? bc_y : 1.0;
              load_u(ut, up, site + Nst * idir);
              real_t wt[NVCD];

              for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                  wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
	            }
              mult_wilson_ypb(vL, ut, wt);
            }

            if ((iy > 0) || (do_comm_y == 0)) {
              int iy2 = (iy - 1 + Ny) % Ny;
              int nei = ix + Nx * (iy2 + Ny * izt);

              real_t bc2 = (iy == 0) ? bc_y : 1.0;
              load_u(ut, up, nei + Nst * idir);
              real_t wt[NVCD];

              for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                  wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
	            }
              mult_wilson_ymb(vL, ut, wt);
            }

            idir = 2;

            if ((iz < Nz-1) || (do_comm_z == 0)) {
              int iz2 = (iz + 1) % Nz;
              int nei = ixy + Nxy * (iz2 + Nz * it);

              real_t bc2 = (iz == Nz-1) ? bc_z : 1.0;
              load_u(ut, up, site + Nst * idir);
              real_t wt[NVCD];

              for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                  wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
	            }
              mult_wilson_zpb(vL, ut, wt);
            }

            if ((iz > 0) || (do_comm_z == 0)) {
              int iz2 = (iz - 1 + Nz) % Nz;
              int nei = ixy + Nxy * (iz2 + Nz * it);

              real_t bc2 = (iz == 0) ? bc_z : 1.0;
              load_u(ut, up, nei + Nst * idir);
              real_t wt[NVCD];

              for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                  wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
	          }
              mult_wilson_zmb(vL, ut, wt);
            }

            idir = 3;

            if ((it < Nt-1) || (do_comm_t == 0)) {
              int it2 = (it + 1) % Nt;
              int nei = ixyz + Nxyz * it2;

              real_t bc2 = (it == Nt-1) ? bc_t : 1.0;
              load_u(ut, up, site + Nst * idir);
              real_t wt[NVCD];

              for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                  wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
	          }
              mult_wilson_tpb_dirac(vL, ut, wt);
            }

            if ((it > 0) || (do_comm_t == 0)) {
              int it2 = (it - 1 + Nt) % Nt;
              int nei = ixyz + Nxyz * it2;

              real_t bc2 = (it == 0) ? bc_t : 1.0;
              load_u(ut, up, nei + Nst * idir);
              real_t wt[NVCD];

              for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
	          }
                mult_wilson_tmb_dirac(vL, ut, wt);
              }

            for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vL[ivcd];
            }
        } // is loop
    }
}

__global__ 
void mult_domainwall_5din_hopb_5D_dirac_kernel(
     real_t *vp, real_t *up, real_t *wp, int Ns, 
     int bc_x, int bc_y, int bc_z, int bc_t,
     int Nx, int Ny, int Nz, int Nt,
     int do_comm_x, int do_comm_y, int do_comm_z, int do_comm_t,
     int flag, int Nst_pad) 
{
    const int Nin5 = NVCD * Ns;
    const int Nxy  = Nx   * Ny;
    const int Nxyz = Nxy * Nz;
    const int Nst  = Nx * Ny * Nz * Nt;
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    extern __shared__ real_t sharedMemory[];
    real_t * u_up = up;
    real_t * u_dn = up;

    if (idx < Nst_pad * Ns) {

        int idx2_wp = idx / NWP;
        int idx_in  = idx % NWP;
        int is      = idx2_wp % Ns;
        int idx_out = idx2_wp / Ns;
        int site    = idx_in + NWP * idx_out;

        int ix   = site % Nx;
        int iyzt = site / Nx;
        int ixy  = site % Nxy;
        int iy   = iyzt % Ny;
        int izt  = site / Nxy;
        int iz   = izt  % Nz;
        int it   = izt  / Nz;
        int ixyz = site % Nxyz;
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

        #include "inc/mult_Domainwall_cuda_flag-inc.h"

        real_t * v1 = wp;
        real_t * v2 = vp;

        #include "inc/mult_Domainwall_cuda_hopb_dirac-inc.h"

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

void mult_domainwall_5din_hopb_dirac(
     real_t *vp, real_t *up, real_t *wp, int Ns, int *bc, int *Nsize, int *do_comm, int flag) {

    int Nx = Nsize[0];
    int Ny = Nsize[1];
    int Nz = Nsize[2];
    int Nt = Nsize[3];
    int Nst = Nx * Ny * Nz * Nt; 
    int Nst_pad = ceil_nwp(Nst);  

    int bc_x = bc[0];
    int bc_y = bc[1];
    int bc_z = bc[2];
    int bc_t = bc[3];
    int do_comm_x = do_comm[0];
    int do_comm_y = do_comm[1];
    int do_comm_z = do_comm[2];
    int do_comm_t = do_comm[3];

    real_t* vp_dev = (real_t*)dev_ptr(vp);
    real_t* up_dev = (real_t*)dev_ptr(up);
    real_t* wp_dev = (real_t*)dev_ptr(wp);

    int blockSize = VECTOR_LENGTH;
    int gridSize  = (Nst_pad * Ns + blockSize - 1) / blockSize;
    int sharedsize = NDF * blockSize * sizeof(real_t); 

    mult_domainwall_5din_hopb_5D_dirac_kernel<<<gridSize, blockSize, sharedsize>>>(
    vp_dev, up_dev, wp_dev, Ns,
    bc_x, bc_y, bc_z, bc_t,
    Nx, Ny, Nz, Nt,
    do_comm_x, do_comm_y, do_comm_z, do_comm_t,
    flag, Nst_pad); 
    
  dw5din_kernel_sync();
}
//====================================================================
__global__ void mult_domainwall_5din_L_inv_dirac_kernel(
    real_t * __restrict__ vp, real_t * __restrict__ wp, int Ns, int Nin5,
    const real_t * __restrict__ e, const real_t * __restrict__ dpinv,
    const real_t * __restrict__ dm, int Nst_pad) {

  // Forward substitution L^{-1} (serial over s, independent per color/Re-Im).
  // One thread per (site,ivc); carries the 4 spin components across the
  // s-recurrence.  e/dpinv/dm are per-operator device arrays (no cudaMalloc).
  const int ist      = blockIdx.x * blockDim.x + threadIdx.x;
  const int GridSize = blockDim.x * gridDim.x;

  for (int idx = ist; idx < Nst_pad * NVC; idx += GridSize) {
    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int ivc     = idx2_wp % NVC;
    int idx_out = idx2_wp / NVC;
    int site    = idx_in + NWP * idx_out;

    real_t vt1, vt2, vt3, vt4, yt1, yt2, yt3, yt4, xt1, xt2, xt3, xt4;

    vt1 = wp[IDX2(Nin5, (ID1 + ivc), site)];
    vt2 = wp[IDX2(Nin5, (ID2 + ivc), site)];
    vt3 = wp[IDX2(Nin5, (ID3 + ivc), site)];
    vt4 = wp[IDX2(Nin5, (ID4 + ivc), site)];
    vp[IDX2(Nin5, (ID1 + ivc), site)] = vt1;
    vp[IDX2(Nin5, (ID2 + ivc), site)] = vt2;
    vp[IDX2(Nin5, (ID3 + ivc), site)] = vt3;
    vp[IDX2(Nin5, (ID4 + ivc), site)] = vt4;
    real_t e0 = e[0];
    yt1 = e0 * vt1; yt2 = e0 * vt2; yt3 = e0 * vt3; yt4 = e0 * vt4;

    for (int is = 1; is < Ns - 1; ++is) {
      xt1 = vt1; xt2 = vt2; xt3 = vt3; xt4 = vt4;
      vt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)];
      vt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)];
      vt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)];
      vt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)];
      real_t a = real_t(0.5) * dm[is] * dpinv[is - 1];
      vt1 += a * (xt1 + xt3);
      vt2 += a * (xt2 + xt4);
      vt3 += a * (xt3 + xt1);
      vt4 += a * (xt4 + xt2);
      vp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)] = vt1;
      vp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)] = vt2;
      vp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)] = vt3;
      vp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)] = vt4;
      real_t eis = e[is];
      yt1 += eis * vt1; yt2 += eis * vt2; yt3 += eis * vt3; yt4 += eis * vt4;
    }

    {
      int is = Ns - 1;
      xt1 = vt1; xt2 = vt2; xt3 = vt3; xt4 = vt4;
      vt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)];
      vt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)];
      vt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)];
      vt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)];
      real_t a = real_t(0.5) * dm[is] * dpinv[is - 1];
      vt1 += a * (xt1 + xt3);
      vt2 += a * (xt2 + xt4);
      vt3 += a * (xt3 + xt1);
      vt4 += a * (xt4 + xt2);
      vt1 += -0.5 * (yt1 - yt3);
      vt2 += -0.5 * (yt2 - yt4);
      vt3 += -0.5 * (yt3 - yt1);
      vt4 += -0.5 * (yt4 - yt2);
      vp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)] = vt1;
      vp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)] = vt2;
      vp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)] = vt3;
      vp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)] = vt4;
    }
  }
}

void mult_domainwall_5din_L_inv_dirac(
    real_t *vp, real_t *wp, int Ns, int *Nsize, real_t *e, real_t *dpinv, real_t *dm) {
  // e/dpinv/dm are device-resident (mapped once by the operator); fetch the
  // device pointers -- no per-call cudaMalloc/cudaMemcpy/cudaFree.
  int Nin5 = NVCD * Ns;
  int Nst  = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];
  int Nst_pad = ceil_nwp(Nst);

  real_t* vp_dev    = (real_t*)dev_ptr(vp);
  real_t* wp_dev    = (real_t*)dev_ptr(wp);
  real_t* e_dev     = (real_t*)dev_ptr(e);
  real_t* dpinv_dev = (real_t*)dev_ptr(dpinv);
  real_t* dm_dev    = (real_t*)dev_ptr(dm);

  int blockSize = VECTOR_LENGTH;
  int gridSize  = (Nst_pad * NVC + blockSize - 1) / blockSize;
  mult_domainwall_5din_L_inv_dirac_kernel<<<gridSize, blockSize>>>(
      vp_dev, wp_dev, Ns, Nin5, e_dev, dpinv_dev, dm_dev, Nst_pad);

  dw5din_kernel_sync();
}

//====================================================================
__global__ void mult_domainwall_5din_U_inv_dirac_kernel(
    real_t * __restrict__ vp, real_t * __restrict__ wp, int Ns, int Nin5,
    const real_t * __restrict__ f, const real_t * __restrict__ dpinv,
    const real_t * __restrict__ dm, int Nst_pad) {

  // Backward substitution U^{-1}.  yt is fixed at the s=Ns-1 boundary value.
  const int ist      = blockIdx.x * blockDim.x + threadIdx.x;
  const int GridSize = blockDim.x * gridDim.x;

  for (int idx = ist; idx < Nst_pad * NVC; idx += GridSize) {
    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int ivc     = idx2_wp % NVC;
    int idx_out = idx2_wp / NVC;
    int site    = idx_in + NWP * idx_out;

    real_t vt1, vt2, vt3, vt4, yt1, yt2, yt3, yt4, xt1, xt2, xt3, xt4;

    int is0 = Ns - 1;
    real_t a0 = dpinv[Ns - 1];
    vt1 = a0 * wp[IDX2(Nin5, (ID1 + ivc + NVCD * is0), site)];
    vt2 = a0 * wp[IDX2(Nin5, (ID2 + ivc + NVCD * is0), site)];
    vt3 = a0 * wp[IDX2(Nin5, (ID3 + ivc + NVCD * is0), site)];
    vt4 = a0 * wp[IDX2(Nin5, (ID4 + ivc + NVCD * is0), site)];
    vp[IDX2(Nin5, (ID1 + ivc + NVCD * is0), site)] = vt1;
    vp[IDX2(Nin5, (ID2 + ivc + NVCD * is0), site)] = vt2;
    vp[IDX2(Nin5, (ID3 + ivc + NVCD * is0), site)] = vt3;
    vp[IDX2(Nin5, (ID4 + ivc + NVCD * is0), site)] = vt4;
    yt1 = 0.5 * (vt1 + vt3);
    yt2 = 0.5 * (vt2 + vt4);
    yt3 = 0.5 * (vt3 + vt1);
    yt4 = 0.5 * (vt4 + vt2);

    for (int is = Ns - 2; is >= 0; --is) {
      xt1 = vt1; xt2 = vt2; xt3 = vt3; xt4 = vt4;
      vt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)];
      vt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)];
      vt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)];
      vt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)];
      real_t a = real_t(0.5) * dm[is];
      vt1 += a * (xt1 - xt3);
      vt2 += a * (xt2 - xt4);
      vt3 += a * (xt3 - xt1);
      vt4 += a * (xt4 - xt2);
      real_t fis = f[is];
      vt1 += -fis * yt1; vt2 += -fis * yt2; vt3 += -fis * yt3; vt4 += -fis * yt4;
      real_t aa = dpinv[is];
      vt1 *= aa; vt2 *= aa; vt3 *= aa; vt4 *= aa;
      vp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)] = vt1;
      vp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)] = vt2;
      vp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)] = vt3;
      vp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)] = vt4;
    }
  }
}

void mult_domainwall_5din_U_inv_dirac(
    real_t *vp, real_t *wp, int Ns, int *Nsize, real_t *f, real_t *dpinv, real_t *dm) {
  int Nin5 = NVCD * Ns;
  int Nst  = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];
  int Nst_pad = ceil_nwp(Nst);

  real_t* vp_dev    = (real_t*)dev_ptr(vp);
  real_t* wp_dev    = (real_t*)dev_ptr(wp);
  real_t* f_dev     = (real_t*)dev_ptr(f);
  real_t* dpinv_dev = (real_t*)dev_ptr(dpinv);
  real_t* dm_dev    = (real_t*)dev_ptr(dm);

  int blockSize = VECTOR_LENGTH;
  int gridSize  = (Nst_pad * NVC + blockSize - 1) / blockSize;
  mult_domainwall_5din_U_inv_dirac_kernel<<<gridSize, blockSize>>>(
      vp_dev, wp_dev, Ns, Nin5, f_dev, dpinv_dev, dm_dev, Nst_pad);

  dw5din_kernel_sync();
}

//====================================================================
__global__ void mult_domainwall_5din_Ldag_inv_dirac_kernel(
    real_t * __restrict__ vp, real_t * __restrict__ wp, int Ns, int Nin5,
    const real_t * __restrict__ e, const real_t * __restrict__ dpinv,
    const real_t * __restrict__ dm, int Nst_pad) {

  const int ist      = blockIdx.x * blockDim.x + threadIdx.x;
  const int GridSize = blockDim.x * gridDim.x;

  for (int idx = ist; idx < Nst_pad * NVC; idx += GridSize) {
    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int ivc     = idx2_wp % NVC;
    int idx_out = idx2_wp / NVC;
    int site    = idx_in + NWP * idx_out;

    real_t vt1, vt2, vt3, vt4, yt1, yt2, yt3, yt4, xt1, xt2, xt3, xt4;

    int is0 = Ns - 1;
    vt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD * is0), site)];
    vt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD * is0), site)];
    vt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD * is0), site)];
    vt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD * is0), site)];
    vp[IDX2(Nin5, (ID1 + ivc + NVCD * is0), site)] = vt1;
    vp[IDX2(Nin5, (ID2 + ivc + NVCD * is0), site)] = vt2;
    vp[IDX2(Nin5, (ID3 + ivc + NVCD * is0), site)] = vt3;
    vp[IDX2(Nin5, (ID4 + ivc + NVCD * is0), site)] = vt4;
    yt1 = 0.5 * (vt1 - vt3);
    yt2 = 0.5 * (vt2 - vt4);
    yt3 = 0.5 * (vt3 - vt1);
    yt4 = 0.5 * (vt4 - vt2);

    for (int is = Ns - 2; is >= 0; --is) {
      xt1 = vt1; xt2 = vt2; xt3 = vt3; xt4 = vt4;
      vt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)];
      vt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)];
      vt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)];
      vt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)];
      real_t a = real_t(0.5) * dm[is + 1] * dpinv[is];
      vt1 += a * (xt1 + xt3);
      vt2 += a * (xt2 + xt4);
      vt3 += a * (xt3 + xt1);
      vt4 += a * (xt4 + xt2);
      real_t eis = e[is];
      vt1 += -eis * yt1; vt2 += -eis * yt2; vt3 += -eis * yt3; vt4 += -eis * yt4;
      vp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)] = vt1;
      vp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)] = vt2;
      vp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)] = vt3;
      vp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)] = vt4;
    }
  }
}

void mult_domainwall_5din_Ldag_inv_dirac(
    real_t *vp, real_t *wp, int Ns, int *Nsize, real_t *e, real_t *dpinv, real_t *dm) {
  int Nin5 = NVCD * Ns;
  int Nst  = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];
  int Nst_pad = ceil_nwp(Nst);

  real_t* vp_dev    = (real_t*)dev_ptr(vp);
  real_t* wp_dev    = (real_t*)dev_ptr(wp);
  real_t* e_dev     = (real_t*)dev_ptr(e);
  real_t* dpinv_dev = (real_t*)dev_ptr(dpinv);
  real_t* dm_dev    = (real_t*)dev_ptr(dm);

  int blockSize = VECTOR_LENGTH;
  int gridSize  = (Nst_pad * NVC + blockSize - 1) / blockSize;
  mult_domainwall_5din_Ldag_inv_dirac_kernel<<<gridSize, blockSize>>>(
      vp_dev, wp_dev, Ns, Nin5, e_dev, dpinv_dev, dm_dev, Nst_pad);

  dw5din_kernel_sync();
}

//====================================================================
__global__ void mult_domainwall_5din_Udag_inv_dirac_kernel(
    real_t * __restrict__ vp, real_t * __restrict__ wp, int Ns, int Nin5,
    const real_t * __restrict__ f, const real_t * __restrict__ dpinv,
    const real_t * __restrict__ dm, int Nst_pad) {

  // Forward substitution Udag^{-1}; yt accumulates over s.
  const int ist      = blockIdx.x * blockDim.x + threadIdx.x;
  const int GridSize = blockDim.x * gridDim.x;

  for (int idx = ist; idx < Nst_pad * NVC; idx += GridSize) {
    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int ivc     = idx2_wp % NVC;
    int idx_out = idx2_wp / NVC;
    int site    = idx_in + NWP * idx_out;

    real_t vt1, vt2, vt3, vt4, yt1, yt2, yt3, yt4, xt1, xt2, xt3, xt4;

    real_t a0 = dpinv[0];
    vt1 = a0 * wp[IDX2(Nin5, (ID1 + ivc), site)];
    vt2 = a0 * wp[IDX2(Nin5, (ID2 + ivc), site)];
    vt3 = a0 * wp[IDX2(Nin5, (ID3 + ivc), site)];
    vt4 = a0 * wp[IDX2(Nin5, (ID4 + ivc), site)];
    vp[IDX2(Nin5, (ID1 + ivc), site)] = vt1;
    vp[IDX2(Nin5, (ID2 + ivc), site)] = vt2;
    vp[IDX2(Nin5, (ID3 + ivc), site)] = vt3;
    vp[IDX2(Nin5, (ID4 + ivc), site)] = vt4;
    real_t f0 = f[0];
    yt1 = f0 * vt1; yt2 = f0 * vt2; yt3 = f0 * vt3; yt4 = f0 * vt4;

    for (int is = 1; is < Ns - 1; ++is) {
      xt1 = vt1; xt2 = vt2; xt3 = vt3; xt4 = vt4;
      vt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)];
      vt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)];
      vt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)];
      vt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)];
      real_t a = real_t(0.5) * dm[is - 1];
      vt1 += a * (xt1 - xt3);
      vt2 += a * (xt2 - xt4);
      vt3 += a * (xt3 - xt1);
      vt4 += a * (xt4 - xt2);
      real_t aa = dpinv[is];
      vt1 *= aa; vt2 *= aa; vt3 *= aa; vt4 *= aa;
      vp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)] = vt1;
      vp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)] = vt2;
      vp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)] = vt3;
      vp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)] = vt4;
      real_t fis = f[is];
      yt1 += fis * vt1; yt2 += fis * vt2; yt3 += fis * vt3; yt4 += fis * vt4;
    }

    {
      int is = Ns - 1;
      xt1 = vt1; xt2 = vt2; xt3 = vt3; xt4 = vt4;
      vt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)];
      vt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)];
      vt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)];
      vt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)];
      real_t a = real_t(0.5) * dm[Ns - 2];
      vt1 += a * (xt1 - xt3);
      vt2 += a * (xt2 - xt4);
      vt3 += a * (xt3 - xt1);
      vt4 += a * (xt4 - xt2);
      vt1 += -0.5 * (yt1 + yt3);
      vt2 += -0.5 * (yt2 + yt4);
      vt3 += -0.5 * (yt3 + yt1);
      vt4 += -0.5 * (yt4 + yt2);
      real_t aa = dpinv[Ns - 1];
      vt1 *= aa; vt2 *= aa; vt3 *= aa; vt4 *= aa;
      vp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)] = vt1;
      vp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)] = vt2;
      vp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)] = vt3;
      vp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)] = vt4;
    }
  }
}

void mult_domainwall_5din_Udag_inv_dirac(
    real_t *vp, real_t *wp, int Ns, int *Nsize, real_t *f, real_t *dpinv, real_t *dm) {
  int Nin5 = NVCD * Ns;
  int Nst  = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];
  int Nst_pad = ceil_nwp(Nst);

  real_t* vp_dev    = (real_t*)dev_ptr(vp);
  real_t* wp_dev    = (real_t*)dev_ptr(wp);
  real_t* f_dev     = (real_t*)dev_ptr(f);
  real_t* dpinv_dev = (real_t*)dev_ptr(dpinv);
  real_t* dm_dev    = (real_t*)dev_ptr(dm);

  int blockSize = VECTOR_LENGTH;
  int gridSize  = (Nst_pad * NVC + blockSize - 1) / blockSize;
  mult_domainwall_5din_Udag_inv_dirac_kernel<<<gridSize, blockSize>>>(
      vp_dev, wp_dev, Ns, Nin5, f_dev, dpinv_dev, dm_dev, Nst_pad);

  dw5din_kernel_sync();
}


//anchor...
__global__  void mult_domainwall_5din_hop1x_4D_dirac_dev(
        real_t * buf_xp, real_t * buf_xm,
        real_t * up, real_t * wp,
        int Ns, int bc, int Nx, int Ny, int Nz, int Nt){

        int iyzt         = blockIdx.x * blockDim.x + threadIdx.x;
        int Nyzt         = Ny   * Nz  * Nt     ;
        int Nin5         = NVCD * Ns           ;
        int Nin5bd       = NVC  * ND2 * Ns     ;
        int Nst          = Nx   * Ny  * Nz * Nt;
        int idir         = 0                   ;
        extern __shared__ real_t sharedMemory[];
        real_t* ut = &sharedMemory[0];

        if ( iyzt < Nyzt ){

            int ix   = 0;
            int site = ix + Nx * iyzt;
            real_t bc2 = bc;

            for (int is = 0; is < Ns; ++is) {
                real_t wt[NVCD], vt[NVC * ND2];
                for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                    wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
	            }
                mult_wilson_xp1(vt, wt);
                for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
                    buf_xp[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), iyzt)] = bc2 * vt[ivcd];
             	}
            }

            ix   = Nx - 1;
            site = ix + Nx * iyzt;
            load_u(ut, up, site + Nst * idir);

            for(int is = 0; is < Ns; ++is) {
                real_t wt[NVCD], vt[NVC * ND2];

                for(int ivcd = 0; ivcd < NVCD; ++ivcd) {
                    wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
	            }
                mult_wilson_xm1(vt, ut, wt);
                for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
                    buf_xm[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), iyzt)] = vt[ivcd];
	            }
            }
        }
    }

__global__  void mult_domainwall_5din_hop1y_4D_dirac_dev(
        real_t * buf_yp, real_t * buf_ym,
        real_t * up, real_t * wp,
        int Ns, int bc, int Nx, int Ny, int Nz, int Nt){

        int site_thread  = blockIdx.x * blockDim.x + threadIdx.x;
        int Nxzt         = Nx * Nz * Nt;
        int ixzt         = site_thread; 
        int Nst          = Nx * Ny * Nz * Nt;
        int Nin5         = NVCD * Ns;
        int Nin5bd       = NVC  * ND2 * Ns;
        extern __shared__ real_t sharedMemory[];
        real_t* ut = &sharedMemory[0];

        if ( ixzt < Nxzt ){

            int idir = 1;
            int iy   = 0;
            int ix   = ixzt % Nx;
            int izt  = ixzt / Nx;
            int site = ix + Nx * (iy + Ny * izt);
            real_t bc2 = bc;

            for(int is = 0; is < Ns; ++is){
                real_t wt[NVCD], vt[NVC * ND2];

                for(int ivcd = 0; ivcd < NVCD ; ++ivcd){
                    wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
                }
                mult_wilson_yp1(vt, wt);
                for(int ivcd = 0; ivcd < NVC * ND2 ; ++ivcd){
                    buf_yp[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixzt)] = bc2 * vt[ivcd];
                }
            }
            
            iy   = Ny - 1;            
            site = ix + Nx * (iy + Ny * izt);
            load_u(ut, up , site + Nst * idir);

            for(int is = 0; is < Ns ; ++is){

                real_t wt[NVCD], vt[NVC * ND2];

                for(int ivcd = 0 ; ivcd < NVCD; ++ivcd) {
                    wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
                }
                mult_wilson_ym1(vt , ut , wt);
                for(int ivcd = 0; ivcd < NVC * ND2 ; ++ivcd) {
                    buf_ym[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixzt)] =  vt[ivcd];
                }
            }
        }
    }

__global__  void mult_domainwall_5din_hop1z_4D_dirac_dev(
        real_t * buf_zp, real_t * buf_zm,
        real_t * up, real_t * wp,
        int Ns, int bc, int Nx, int Ny, int Nz, int Nt){

        int site_thread  = blockIdx.x * blockDim.x + threadIdx.x;

        int ixyt         = site_thread;
        int Nst          = Nx * Ny * Nz * Nt;
        int Nxy          = Nx * Ny;
        int Nxyt         = Nx * Ny * Nt;
        int Nin5         = NVCD * Ns;
        int Nin5bd       = NVC  * ND2 * Ns;
        extern __shared__ real_t sharedMemory[];
        real_t* ut = &sharedMemory[0];

        if ( ixyt < Nxyt ){

            int idir  = 2;
            int iz    = 0;
            int ixy   = ixyt % Nxy;
            int it    = ixyt / Nxy;
            int site  = ixy  + Nxy * (iz + Nz * it);
            real_t bc2 = bc;

            for(int is = 0; is < Ns; ++is){
                
                real_t wt[NVCD], vt[NVC * ND2];

                for(int ivcd = 0; ivcd < NVCD ; ++ivcd){
                    wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
                }
                mult_wilson_zp1(vt, wt);
                for(int ivcd = 0; ivcd < NVC * ND2 ; ++ivcd){
                    buf_zp[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixyt)] = bc2 * vt[ivcd];
                }
            }

            iz   = Nz - 1;
            site = ixy + Nxy * (iz + Nz * it);
            load_u(ut, up , site + Nst * idir);

            for(int is = 0; is < Ns ; ++is){

                real_t wt[NVCD], vt[NVC * ND2];

                for(int ivcd = 0 ; ivcd < NVCD; ++ivcd) {
                    wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
                }
                mult_wilson_zm1(vt , ut , wt);
                for(int ivcd = 0; ivcd <NVC * ND2 ; ++ivcd) {
                    buf_zm[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixyt)] =  vt[ivcd];
                }
            }
        }
    }

__global__  void mult_domainwall_5din_hop1t_4D_dirac_dev(
        real_t * buf_tp, real_t * buf_tm,
        real_t * up, real_t * wp,
        int Ns, int bc, int Nx, int Ny, int Nz, int Nt){

        int site_thread  = blockIdx.x * blockDim.x + threadIdx.x;
        int ixyz         = site_thread;
        int Nst          = Nx * Ny * Nz * Nt;
        int Nxyz         = Nx * Ny * Nz;
        int Nin5         = NVCD * Ns;
        int Nin5bd       = NVC  * ND2 * Ns;
        extern __shared__ real_t sharedMemory[];
        real_t* ut = &sharedMemory[0];

        if ( ixyz < Nxyz ){

            int idir   = 3;
            int it     = 0;
            int site   = ixyz + Nxyz * it;
            real_t bc2 = bc;

            for(int is = 0; is < Ns; ++is){

                real_t wt[NVCD], vt[NVC * ND2];

                for(int ivcd = 0; ivcd < NVCD ; ++ivcd){
                    wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
                    }   

                mult_wilson_tp1_dirac(vt, wt);

                for(int ivcd = 0; ivcd < NVC * ND2 ; ++ivcd){
                    buf_tp[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixyz)] = bc2 * vt[ivcd];
                    }
            }

            it    = Nt - 1;
            site  = ixyz + Nxyz * it;
            load_u(ut, up , site + Nst *idir);

            for(int is = 0; is < Ns ; ++is){

                real_t wt[NVCD], vt[NVC * ND2];

                for(int ivcd = 0 ; ivcd < NVCD; ++ivcd) {
                    wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
                }
                mult_wilson_tm1_dirac(vt , ut , wt);
                for(int ivcd = 0; ivcd < NVC * ND2 ; ++ivcd){
                    buf_tm[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixyz)] = vt[ivcd];
                }
            }
        }
    }

__global__ void mult_domainwall_5din_hop1x_5D_dirac_dev(
    real_t *buf_xp, real_t *buf_xm,
    real_t *up, real_t *wp,
    int Ns, int bc, int Nx, int Ny, int Nz, int Nt, int Nstx_pad) 
{
    int Nyzt   = Ny * Nz * Nt;
    int Nin5   = NVCD * Ns;
    int Nin5bd = NVC * ND2 * Ns;
    int Nst    = Nx * Ny * Nz * Nt;
    int idx    = blockIdx.x * blockDim.x + threadIdx.x;

    extern __shared__ real_t sharedMemory[];
    real_t* ut = &sharedMemory[0];

    if (idx < Nstx_pad * Ns) {

        int idx_wp  = idx / NWP;
        int idx_in  = idx % NWP;
        int is      = idx_wp % Ns;
        int idx_out = idx_wp / Ns;
        int iyzt    = idx_in + NWP * idx_out;

        {
            int ix = 0;
            int site = ix + Nx * iyzt;
            real_t bc2 = bc;

            real_t wt[NVCD], vt[NVC * ND2];
            for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
                wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
            }
            mult_wilson_xp1(vt, wt);
            for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
                buf_xp[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), iyzt)] = bc2 * vt[ivcd];
            }
        }
        {
            int ix = Nx - 1;
            int site = ix + Nx * iyzt;
            load_u(ut, up, site + Nst * 0);

            real_t wt[NVCD], vt[NVC * ND2];
            for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
                wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
            }
            mult_wilson_xm1(vt, ut, wt);
            for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
                buf_xm[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), iyzt)] = vt[ivcd];
            }
        }
    }
}

__global__ void mult_domainwall_5din_hop1y_5D_dirac_dev(
    real_t *buf_yp, real_t *buf_ym,
    real_t *up, real_t *wp,
    int Ns, int bc, int Nx, int Ny, int Nz, int Nt, int Nxzt_pad) 
{
    int Nxzt = Nx * Nz * Nt;
    int Nin5 = NVCD * Ns;
    int Nin5bd = NVC * ND2 * Ns;
    int Nst = Nx * Ny * Nz * Nt;
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    extern __shared__ real_t sharedMemory[];

    real_t* ut = &sharedMemory[0];
    real_t* wt = &sharedMemory[NVCD];

    if (idx < Nxzt_pad * Ns) {
 
        int idx_wp   = idx / NWP;
        int idx_in   = idx % NWP;
        int is       = idx_wp % Ns;
        int idx_out  = idx_wp / Ns;
        int ixzt     = idx_in + NWP * idx_out;

        {
            int iy   = 0;
            int ix   = ixzt % Nx;
            int izt  = ixzt / Nx;
            int site = ix + Nx * (iy + Ny * izt);

            real_t bc2 = bc;
            real_t vt[NVC * ND2];

            #pragma unroll
            for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
                wt[ivcd] = wp[IDX2_5D(NVCD, ivcd, is, Ns,site)];
            }
            mult_wilson_yp1(vt, wt);

            #pragma unroll
            for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
                buf_yp[IDX2_5D(NVC * ND2, ivcd, is, Ns, ixzt)]= bc2 * vt[ivcd];
            }
        }

        {
            int iy = Ny - 1;
            int ix = ixzt % Nx;
            int izt = ixzt / Nx;
            int site = ix + Nx * (iy + Ny * izt);
            load_u(ut, up, site + Nst * 1);

            real_t vt[NVC * ND2];
            #pragma unroll
            for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
                wt[ivcd] = wp[IDX2_5D(NVCD, ivcd, is, Ns,site)];
            }
            mult_wilson_ym1(vt, ut, wt);
            
            #pragma unroll
            for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
                buf_ym[IDX2_5D(NVC * ND2, ivcd, is, Ns, ixzt)] = vt[ivcd];
            }
        }
    }
}

__global__ void mult_domainwall_5din_hop1z_5D_dirac_dev(
    real_t *buf_zp, real_t *buf_zm,
    real_t *up, real_t *wp,
    int Ns, int bc, int Nx, int Ny, int Nz, int Nt, int Nstz_pad) 
{
    int Nxy = Nx * Ny;
    int Nxyt = Nx * Ny * Nt;
    int Nin5 = NVCD * Ns;
    int Nin5bd = NVC * ND2 * Ns;
    int Nst = Nx * Ny * Nz * Nt;
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    extern __shared__ real_t sharedMemory[];
    real_t* ut = &sharedMemory[0];

    if (idx < Nstz_pad * Ns) {

        int idx_wp   = idx / NWP;
        int idx_in   = idx % NWP;
        int is       = idx_wp % Ns;
        int idx_out  = idx_wp / Ns;
        int ixyt     = idx_in + NWP * idx_out;

        {
            int iz = 0;
            int ixy = ixyt % Nxy;
            int it = ixyt / Nxy;
            int site = ixy + Nxy * (iz + Nz * it);
            real_t bc2 = bc;

            real_t wt[NVCD], vt[NVC * ND2];
            for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
                wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
            }
            mult_wilson_zp1(vt, wt);
            for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
                buf_zp[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixyt)] = bc2 * vt[ivcd];
            }
        }

        {
            int iz = Nz - 1;
            int ixy = ixyt % Nxy;
            int it = ixyt / Nxy;
            int site = ixy + Nxy * (iz + Nz * it);
            load_u(ut, up, site + Nst * 2);

            real_t wt[NVCD], vt[NVC * ND2];
            for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
                wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
            }
            mult_wilson_zm1(vt, ut, wt);
            for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
                buf_zm[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixyt)] = vt[ivcd];
            }
        }
    }
}
__global__ void mult_domainwall_5din_hop1t_5D_dirac_dev(
    real_t *buf_tp, real_t *buf_tm,
    real_t *up, real_t *wp,
    int Ns, int bc, int Nx, int Ny, int Nz, int Nt, int Nstt_pad) 
{
    int Nxyz = Nx * Ny * Nz;
    int Nin5 = NVCD * Ns;
    int Nin5bd = NVC * ND2 * Ns;
    int Nst = Nx * Ny * Nz * Nt;
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    extern __shared__ real_t sharedMemory[];
    real_t* ut = &sharedMemory[0];

    if (idx < Nstt_pad * Ns) {

        int idx_wp   = idx / NWP;
        int idx_in   = idx % NWP;
        int is       = idx_wp % Ns;
        int idx_out  = idx_wp / Ns;
        int ixyz     = idx_in + NWP * idx_out;

        {
            int it = 0;
            int site = ixyz + Nxyz * it;
            real_t bc2 = bc;

            real_t wt[NVCD], vt[NVC * ND2];
            for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
                wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
            }
            mult_wilson_tp1_dirac(vt, wt);
            for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
                buf_tp[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixyz)] = bc2 * vt[ivcd];
            }
        }

        {
            int it = Nt - 1;
            int site = ixyz + Nxyz * it;
            load_u(ut, up, site + Nst * 3);

            real_t wt[NVCD], vt[NVC * ND2];
            for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
                wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
            }
            mult_wilson_tm1_dirac(vt, ut, wt);
            for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
                buf_tm[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixyz)] = vt[ivcd];
            }
        }
    }
}


//-----------------------------
void mult_domainwall_5din_hop1_dirac(
    real_t * buf1_xp, real_t * buf1_xm,
    real_t * buf1_yp, real_t * buf1_ym,
    real_t * buf1_zp, real_t * buf1_zm,
    real_t * buf1_tp, real_t * buf1_tm,
    real_t * up, real_t * wp,
    int Ns, int *bc, int *Nsize, int *do_comm){

    int Nx     = Nsize[0];
    int Ny     = Nsize[1];
    int Nz     = Nsize[2];
    int Nt     = Nsize[3];
    int Nst    = Nx * Ny * Nz * Nt;
    int Nin5bd = NVC * Ns * ND2;

    int bc_x = bc[0];
    int bc_y = bc[1];
    int bc_z = bc[2];
    int bc_t = bc[3];

    int do_comm_x = do_comm[0];
    int do_comm_y = do_comm[1];
    int do_comm_z = do_comm[2];
    int do_comm_t = do_comm[3];

    size_t size_bx = Nin5bd * Ny * Nz * Nt;
    size_t size_by = Nin5bd * Nx * Nz * Nt;
    size_t size_bz = Nin5bd * Ny * Nx * Nt;
    size_t size_bt = Nin5bd * Ny * Nz * Nx;

    real_t * up_dev = (real_t*)dev_ptr(up);
    real_t * wp_dev = (real_t*)dev_ptr(wp);

    //int Nst_pad = ceil_nwp(Nst);  

    if( do_comm_x > 0 ){

        int Nstx =  Ny * Nz * Nt;
        int Nstx_pad = ceil_nwp(Nstx);  

        // real_t * buf1_xp_dev = (real_t*)dev_ptr(buf1_xp);
        // real_t * buf1_xm_dev = (real_t*)dev_ptr(buf1_xm);

        int blockSize   =  VECTOR_LENGTH;
        int gridSize5D  =  (Nstx * Ns + blockSize - 1)/ blockSize;
        int sharedsize  = NDF * sizeof(real_t);

        mult_domainwall_5din_hop1x_5D_dirac_dev<<<gridSize5D, blockSize, sharedsize>>>(
        buf1_xp, buf1_xm,
        up_dev, wp_dev,
        Ns, bc_x, Nx, Ny, Nz, Nt, Nstx_pad);
    }

    if( do_comm_y > 0){

        int Nsty =  Nx * Nz * Nt;
        int Nsty_pad = ceil_nwp(Nsty);  

        // real_t * buf1_yp_dev = (real_t*)dev_ptr(buf1_yp);
        // real_t * buf1_ym_dev = (real_t*)dev_ptr(buf1_ym);

        int blockSize  = VECTOR_LENGTH;
        int gridSize5D =  (Nsty * Ns + blockSize - 1)/ blockSize;
        int sharedsize = NDF *  sizeof(real_t);
        mult_domainwall_5din_hop1y_5D_dirac_dev<<<gridSize5D, blockSize, sharedsize>>>(
        buf1_yp, buf1_ym,
        up_dev, wp_dev,
        Ns, bc_y, Nx, Ny, Nz, Nt, Nsty_pad);

    }

    if( do_comm_z > 0 ){

        int Nstz =  Nx * Ny * Nt;
        int Nstz_pad = ceil_nwp(Nstz);  

        // real_t * buf1_zp_dev = (real_t*)dev_ptr(buf1_zp);
        // real_t * buf1_zm_dev = (real_t*)dev_ptr(buf1_zm);

        int blockSize = VECTOR_LENGTH;
        int gridSize5D  =  (Nstz_pad * Ns + blockSize - 1)/ blockSize;
        int sharedsize = NDF * sizeof(real_t);
        mult_domainwall_5din_hop1z_5D_dirac_dev<<<gridSize5D, blockSize, sharedsize>>>(
        buf1_zp, buf1_zm,
        up_dev, wp_dev,
        Ns, bc_z, Nx, Ny, Nz, Nt, Nstz_pad);

    }

    if( do_comm_t > 0 ){

        int Nstt =  Nx * Ny * Nz;
        int Nstt_pad = ceil_nwp(Nstt);  

        // real_t * buf1_tp_dev = (real_t*)dev_ptr(buf1_tp);
        // real_t * buf1_tm_dev = (real_t*)dev_ptr(buf1_tm);

        int blockSize  = VECTOR_LENGTH;
        int sharedsize = NDF  * sizeof(real_t);
        int gridSize5D =  (Nstt_pad * Ns + blockSize - 1)/ blockSize;
        mult_domainwall_5din_hop1t_5D_dirac_dev<<<gridSize5D, blockSize, sharedsize>>>(
        buf1_tp, buf1_tm,
        up_dev, wp_dev,
        Ns, bc_t, Nx, Ny, Nz, Nt, Nstt_pad);

/*
        int gridSize4D  = (Nstt + blockSize - 1)/ blockSize;
        mult_domainwall_5din_hop1t_dirac_dev<<<gridSize4D, blockSize, sharedsize>>>(
        buf1_tp_dev, buf1_tm_dev,
        up_dev, wp_dev,
        Ns, bc_t, Nx, Ny, Nz, Nt);
*/
    }

    CHECK(cudaDeviceSynchronize());

    // cudaStream_t* cuda_stream = nullptr;
    // int nstream = 8;
    // CREATE_CUDA_STREAMS(&cuda_stream, nstream);

    // if (cuda_stream == nullptr) {
    //     fprintf(stderr, "Failed to create CUDA streams in hop1.\n");
    // }

    // if ( do_comm_x > 0) {
    //     update_host_asy(buf1_xp, 0 , size_bx, cuda_stream[0]);
    //     update_host_asy(buf1_xm, 0 , size_bx, cuda_stream[1]);
    // } 
        
    // if ( do_comm_y > 0) {
    //     update_host_asy(buf1_yp, 0 , size_by, cuda_stream[2]);
    //     update_host_asy(buf1_ym, 0 , size_by, cuda_stream[3]);
    // } 

    // if ( do_comm_z > 0) {
    //     update_host_asy(buf1_zp, 0 , size_bz, cuda_stream[4]);
    //     update_host_asy(buf1_zm, 0 , size_bz, cuda_stream[5]);
    // } 

    // if ( do_comm_t > 0) {
    //     update_host_asy(buf1_tp, 0 , size_bt, cuda_stream[6]);
    //     update_host_asy(buf1_tm, 0 , size_bt, cuda_stream[7]);
    // } 

    // DESTROY_CUDA_STREAMS(&cuda_stream, nstream);

}

__global__ void mult_domainwall_5din_hop2_4D_dirac_dev(
    real_t *vp, real_t *up, real_t *wp,
    real_t *buf_xp, real_t *buf_xm,
    real_t *buf_yp, real_t *buf_ym,
    real_t *buf_zp, real_t *buf_zm,
    real_t *buf_tp, real_t *buf_tm,
    int Ns,
    int bc_x, int bc_y, int bc_z, int bc_t,
    int do_comm_x, int do_comm_y, int do_comm_z, int do_comm_t, 
    int Nx, int Ny, int Nz, int Nt) {

    int site         = blockIdx.x * blockDim.x + threadIdx.x;
    int Nxy          = Nx   * Ny           ;
    int Nxyz         = Nx   * Ny  * Nz     ;
    int Nin5bd       = (NVCD / 2) * Ns     ;
    int Nst          = Nx   * Ny  * Nz * Nt;
    int Nin5         = NVCD * Ns           ;
    extern __shared__ real_t sharedMemory[];
    real_t* ut = &sharedMemory[0];

    if (site < Nst) {

        int ix   = site % Nx;
        int iyzt = site / Nx;
        int ixy  = site % Nxy;
        int iy   = iyzt % Ny;
        int izt  = site / Nxy;
        int iz   = izt  % Nz;
        int it   = izt  / Nz;
        int ixyz = site % Nxyz;
        int idir;

        for(int is = 0; is < Ns; ++is){

            real_t vL[NVCD];

            for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                vL[ivcd] = 0.0;
            }
            
            int opr_any = 0;

            idir = 0;
            if (do_comm_x > 0) {

                if (ix == Nx-1) {
                    load_u(ut, up, site + Nst * idir);
                    real_t wt[NVC * ND2];
                    for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
                        wt[ivcd] = buf_xp[IDX2(Nin5bd, (ivcd + NVC*ND2*is), iyzt)];
                    }
                    mult_wilson_xp2(vL, ut, wt);
                    ++opr_any;
                }

                if (ix == 0) {
                    real_t bc2 = bc_x;
                    real_t wt[NVC * ND2];
                    for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
                        wt[ivcd] = bc2 * buf_xm[IDX2(Nin5bd, (ivcd + NVC*ND2*is), iyzt)];
                    }
                    mult_wilson_xm2(vL, wt);
                    ++opr_any;
                }
            }

            idir = 1;
            if (do_comm_y > 0) {
                int ixzt = ix + Nx * izt;
                if (iy == Ny-1) {
                    load_u(ut, up, site + Nst * idir);
                    real_t wt[NVC * ND2];
                    for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
                        wt[ivcd] = buf_yp[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixzt)];
                    }
                    mult_wilson_yp2(vL, ut, wt);
                    ++opr_any;
                }

                if (iy == 0) {
                    real_t bc2 = bc_y;
                    real_t wt[NVC * ND2];
                    for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
                        wt[ivcd] = bc2 * buf_ym[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixzt)];
                    }
                    mult_wilson_ym2(vL, wt);
                    ++opr_any;
                }
            }

            idir = 2;
            if (do_comm_z > 0) {
                int ixyt = ixy + Nxy * it;
                if (iz == Nz-1) {
                    load_u(ut, up, site + Nst * idir);
                    real_t wt[NVC * ND2];
                    for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
                        wt[ivcd] = buf_zp[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixyt)];
                    }
                    mult_wilson_zp2(vL, ut, wt);
                    ++opr_any;
                }

                if (iz == 0) {
                    real_t bc2 = bc_z;
                    real_t wt[NVC * ND2];
                    for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
                        wt[ivcd] = bc2 * buf_zm[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixyt)];
                    }
                    mult_wilson_zm2(vL, wt);
                    ++opr_any;
                }
            }

            idir = 3;
            if (do_comm_t > 0) {
                if (it == Nt-1) {
                    load_u(ut, up, site + Nst * idir);
                    real_t wt[NVC * ND2];
                    for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
                        wt[ivcd] = buf_tp[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixyz)];
                    }
                    mult_wilson_tp2_dirac(vL, ut, wt);
                    ++opr_any;
                }

                if (it == 0) {
                    real_t bc2 = bc_t;
                    real_t wt[NVC * ND2];
                    for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
                        wt[ivcd] = bc2 * buf_tm[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixyz)];
                    }
                    mult_wilson_tm2_dirac(vL, wt);
                    ++opr_any;
                }
            }

            if (opr_any > 0) {
                for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                    vp[IDX2(Nin5, (ivcd + NVCD*is), site)] += vL[ivcd];
	            }
            }
        }   // is loop
    }
}

__global__ void mult_domainwall_5din_hop2_5D_dirac_dev(
    real_t *vp, real_t *up, real_t *wp,
    real_t *buf_xp, real_t *buf_xm,
    real_t *buf_yp, real_t *buf_ym,
    real_t *buf_zp, real_t *buf_zm,
    real_t *buf_tp, real_t *buf_tm,
    int Ns,
    int bc_x, int bc_y, int bc_z, int bc_t,
    int do_comm_x, int do_comm_y, int do_comm_z, int do_comm_t, 
    int Nx, int Ny, int Nz, int Nt, int Nst_pad) 
{
    int Nxy    = Nx * Ny;
    int Nxyz   = Nx * Ny * Nz;
    int Nin5bd = (NVCD / 2) * Ns;
    int Nst    = Nx * Ny * Nz * Nt;
    int Nin5   = NVCD * Ns;
    int idx    = blockIdx.x * blockDim.x + threadIdx.x;
    extern __shared__ real_t sharedMemory[];
    real_t* ut = &sharedMemory[0];

    if (idx < Nst * Ns) {

        int idx2_wp = idx / NWP;
        int idx_in  = idx % NWP;
        int is      = idx2_wp % Ns;
        int idx_out = idx2_wp / Ns;
        int site    = idx_in + NWP * idx_out;

        int ix   = site % Nx;
        int iyzt = site / Nx;
        int ixy  = site % Nxy;
        int iy   = iyzt % Ny;
        int izt  = site / Nxy;
        int iz   = izt % Nz;
        int it   = izt / Nz;
        int ixyz = site % Nxyz;

        real_t vL[NVCD];
        for(int ivcd = 0; ivcd < NVCD; ++ivcd) {
            vL[ivcd] = 0.0;
        }

        int opr_any = 0;
        int idir;

        idir = 0;
        if (do_comm_x > 0) {

            if (ix == Nx-1) {
                load_u(ut, up, site + Nst * idir);
                real_t wt[NVC * ND2];
                for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
                    wt[ivcd] = buf_xp[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), iyzt)];
                }
                mult_wilson_xp2(vL, ut, wt);
                ++opr_any;
            }

            if (ix == 0) {
                real_t bc2 = bc_x;
                real_t wt[NVC * ND2];
                for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
                    wt[ivcd] = bc2 * buf_xm[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), iyzt)];
                }
                mult_wilson_xm2(vL, wt);
                ++opr_any;
            }
        }

        idir = 1;
        if (do_comm_y > 0) {
            int ixzt = ix + Nx * izt;
            if (iy == Ny-1) {
                load_u(ut, up, site + Nst * idir);
                real_t wt[NVC * ND2];
                for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
                    wt[ivcd] = buf_yp[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixzt)];
                }
                mult_wilson_yp2(vL, ut, wt);
                ++opr_any;
            }

            if (iy == 0) {
                real_t bc2 = bc_y;
                real_t wt[NVC * ND2];
                for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
                    wt[ivcd] = bc2 * buf_ym[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixzt)];
                }
                mult_wilson_ym2(vL, wt);
                ++opr_any;
            }
        }

        idir = 2;
        if (do_comm_z > 0) {
            int ixyt = ixy + Nxy * it;
            if (iz == Nz-1) {
                load_u(ut, up, site + Nst * idir);
                real_t wt[NVC * ND2];
                for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
                    wt[ivcd] = buf_zp[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixyt)];
                }
                mult_wilson_zp2(vL, ut, wt);
                ++opr_any;
            }

            if (iz == 0) {
                real_t bc2 = bc_z;
                real_t wt[NVC * ND2];
                for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
                    wt[ivcd] = bc2 * buf_zm[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixyt)];
                }
                mult_wilson_zm2(vL, wt);
                ++opr_any;
            }
        }

        idir = 3;
        if (do_comm_t > 0) {
            if (it == Nt-1) {
                load_u(ut, up, site + Nst * idir);
                real_t wt[NVC * ND2];
                for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
                    wt[ivcd] = buf_tp[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixyz)];
                }
                mult_wilson_tp2_dirac(vL, ut, wt);
                ++opr_any;
            }

            if (it == 0) {
                real_t bc2 = bc_t;
                real_t wt[NVC * ND2];
                for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
                    wt[ivcd] = bc2 * buf_tm[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixyz)];
                }
                mult_wilson_tm2_dirac(vL, wt);
                ++opr_any;
            }
        }

        if (opr_any > 0) {
            for(int ivcd = 0; ivcd < NVCD; ++ivcd) {
                vp[IDX2(Nin5, (ivcd + NVCD * is), site)] += vL[ivcd];
            }
        }
    }
}

void mult_domainwall_5din_hop2_dirac(
        real_t * vp     , real_t * up, real_t * wp,
        real_t * buf2_xp, real_t * buf2_xm,
        real_t * buf2_yp, real_t * buf2_ym,
        real_t * buf2_zp, real_t * buf2_zm,
        real_t * buf2_tp, real_t * buf2_tm,
        int Ns, int *bc, int *Nsize, int *do_comm){

        int Nx     = Nsize[0];
        int Ny     = Nsize[1];
        int Nz     = Nsize[2];
        int Nt     = Nsize[3];
        int Nst    = Nx * Ny * Nz * Nt;
        int Nin5bd = (NVCD / 2) * Ns;
        int Nst_pad = ceil_nwp(Nst);  

        int bc_x = bc[0];
        int bc_y = bc[1];
        int bc_z = bc[2];
        int bc_t = bc[3];
        int do_comm_x = do_comm[0];
        int do_comm_y = do_comm[1];
        int do_comm_z = do_comm[2];
        int do_comm_t = do_comm[3];

        real_t * vp_dev = (real_t*)dev_ptr(vp);
        real_t * up_dev = (real_t*)dev_ptr(up);
        real_t * wp_dev = (real_t*)dev_ptr(wp);

        int blockSize  = VECTOR_LENGTH;

        int sharedsize = NDF * sizeof(real_t);
        int gridSize5D   = (Nst_pad * Ns + blockSize - 1)/ blockSize;
        mult_domainwall_5din_hop2_5D_dirac_dev<<<gridSize5D, VECTOR_LENGTH, sharedsize>>>(
            vp_dev,  up_dev,  wp_dev,
            buf2_xp, buf2_xm,
            buf2_yp, buf2_ym,
            buf2_zp, buf2_zm,
            buf2_tp, buf2_tm,
            Ns,  
            bc_x, bc_y, bc_z, bc_t,
            do_comm_x, do_comm_y, do_comm_z, do_comm_t,
            Nx, Ny, Nz, Nt, Nst_pad);

        CHECK(cudaDeviceSynchronize());
}

#endif
//============================================================END=====