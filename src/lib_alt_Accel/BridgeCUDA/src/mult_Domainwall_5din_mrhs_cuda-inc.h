/*!
      @file    mult_Domainwall_5din_mrhs_cuda-inc.h
      @brief   MRHS (multiple-RHS) batched fine Domainwall_5din D operator.

      Single-GPU (no-comm) path of the non-eo Moebius domain-wall D:
            D = 5dir_dirac  +  hopb_dirac
      batched over nrhs columns in ONE launch each, fully on device (no D2H/H2D).
      The per-(site,is) arithmetic is reused VERBATIM from the production
      single-RHS kernels (mult_Domainwall_5din_cuda-inc.h + the shared
      mult_Domainwall_cuda_{flag,hopb_dirac}-inc.h includes); the only change is
      an inner rhs loop and b/c read from this module's own __constant__ memory
      (uploaded once via mrhs_set_moebius_bc; mass-independent so one upload
      serves both the mq and PV operators).

      NOTE (consolidation TODO): the proper home is inside the BridgeACC library
      (bridgeACC_Domainwall_5din_float.cu) where the production const_b/const_c
      and dev_ptr resolve natively and the module-local copy can be dropped.  Kept
      dev8-local for now -- relocate at the main-project merge.

      This is the "make it work" milestone for batching the whole AMG fine
      operator: correct + single-launch + on-device.  Link-hoisting (loading each
      gauge link once and reusing across is AND rhs to raise arithmetic
      intensity) is a follow-on optimisation verified against this baseline.

      @author  Wei-Lun Chen (wlchen)
*/

#include "inline/mult_Wilson_cuda_inline-inc.h"   // MULT_GXr/i, MULT_GDXr/i

namespace {  // anonymous: lives inside the BridgeACC TU (dev_ptr/macros visible)

// Moebius b/c in this module's OWN __constant__ memory.  The production
// const_b_float lives in libbridgecuda's device module and is NOT reachable
// across the static-library device-link boundary (nvlink undefined reference),
// so we keep a dev8-local copy and upload it ONCE via mrhs_set_moebius_bc.
__constant__ real_t mrhs_cb[NS_MAX];
__constant__ real_t mrhs_cc[NS_MAX];

void mrhs_set_moebius_bc_dev(const real_t* b, const real_t* c, int Ns)
{
  CHECK(cudaMemcpyToSymbol(mrhs_cb, b, sizeof(real_t) * Ns));
  CHECK(cudaMemcpyToSymbol(mrhs_cc, c, sizeof(real_t) * Ns));
}

//====================================================================
// MRHS 5d-direction Dirac block.  For each rhs r:
//   vp_r = D5 w_r   (diagonal + 5d-Wilson hop, Neff-alpha scaled)
//   yp_r = -0.5 (b self + c hop)   (auxiliary consumed by hopb)
// One thread per (site,ivc); inner loop over rhs and is.  No gauge => batching
// only amortises launch overhead, but keeps the operator on-device + single-launch.
//====================================================================
__global__
void mult_dw5din_5dir_mrhs_dev(
    real_t* const* __restrict__ vp_arr,
    real_t* const* __restrict__ yp_arr,
    real_t* const* __restrict__ wp_arr,
    int nrhs, real_t mq, real_t M0, int Ns, real_t alpha,
    int Nst_pad)
{
  // b/c from this module's __constant__ memory (uploaded once, mass-independent).
  const real_t* b_arr = mrhs_cb;
  const real_t* c_arr = mrhs_cc;
  const int Nin5     = NVCD * Ns;
  const int ist      = blockIdx.x * blockDim.x + threadIdx.x;
  const int GridSize = blockDim.x * gridDim.x;

  for (int idx = ist; idx < Nst_pad * NVC; idx += GridSize) {

    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int ivc     = idx2_wp % NVC;
    int idx_out = idx2_wp / NVC;
    int site    = idx_in + NWP * idx_out;

    for (int r = 0; r < nrhs; ++r) {
      real_t* __restrict__ vp = vp_arr[r];
      real_t* __restrict__ yp = yp_arr[r];
      real_t* __restrict__ wp = wp_arr[r];

      for (int is = 0; is < Ns; ++is) {

        real_t B_is = b_arr[is] * (4.0 - M0) + 1.0;
        real_t C_is = c_arr[is] * (4.0 - M0) - 1.0;

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

        real_t b_is = b_arr[is];
        real_t c_is = c_arr[is];
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
}

//====================================================================
// MRHS 4d bulk hopping (hopb).  For each rhs r: vp_r += U-hop( yp_r ).
// Reuses the shared per-(site,is) flag-init + hopb_dirac arithmetic includes
// verbatim inside an rhs loop.  Single-GPU: do_comm_* are all 0 so the
// boundary branches use the periodic neighbour (no halo).
//====================================================================
__global__
void mult_dw5din_hopb_mrhs_dev(
    real_t* const* __restrict__ vp_arr, real_t* __restrict__ up,
    real_t* const* __restrict__ wp_arr, int nrhs, int Ns,
    int bc_x, int bc_y, int bc_z, int bc_t,
    int Nx, int Ny, int Nz, int Nt,
    int do_comm_x, int do_comm_y, int do_comm_z, int do_comm_t,
    int flag, int Nst_pad)
{
  const int Nin5 = NVCD * Ns;
  const int Nxy  = Nx   * Ny;
  const int Nxyz = Nxy * Nz;
  const int Nst  = Nx * Ny * Nz * Nt;
  const int idx  = blockIdx.x * blockDim.x + threadIdx.x;
  real_t* u_up = up;
  real_t* u_dn = up;

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

    for (int r = 0; r < nrhs; ++r) {

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

      real_t* __restrict__ vp = vp_arr[r];
      real_t* __restrict__ wp = wp_arr[r];

      #include "inc/mult_Domainwall_cuda_flag-inc.h"

      real_t* v1 = wp;
      real_t* v2 = vp;

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
}

//====================================================================
// Host driver: v_r = D w_r for r=0..nrhs-1, single-GPU no-comm path.
// u_field_host = host base pointer of the fine gauge AField (dev_ptr'd here).
// b_host/c_host = Moebius b/c arrays (length Ns), uploaded to device.
//====================================================================
static real_t*  g_dw_yp_buf = nullptr;   // nrhs * Nin5 * Nst_pad scratch (yp)
static long     g_dw_yp_cap = 0;

void fineD_mrhs(real_t* const* v_host, real_t* const* w_host, real_t* u_field_host,
                int nrhs, real_t mq, real_t M0, int Ns, real_t alpha,
                int* Nsize, int* bc, int* do_comm)
{
  int Nst     = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];
  int Nst_pad = ceil_nwp(Nst);
  int Nin5    = NVCD * Ns;

  real_t* u_dev = (real_t*)dev_ptr(u_field_host);

  // device pointer arrays for v (out) and w (in)
  static real_t **vd = nullptr, **wd = nullptr; static int cv = 0, cw = 0;
  real_t** vdev = map_upload(v_host, nrhs, vd, cv);
  real_t** wdev = map_upload(w_host, nrhs, wd, cw);

  // yp scratch: one Nin5*Nst_pad field per rhs, contiguous, with pointer array
  long yp_need = (long)nrhs * Nin5 * Nst_pad;
  if (yp_need > g_dw_yp_cap) {
    if (g_dw_yp_buf) cudaFree(g_dw_yp_buf);
    CHECK(cudaMalloc((void**)&g_dw_yp_buf, sizeof(real_t) * yp_need));
    g_dw_yp_cap = yp_need;
  }
  real_t* yp_host_ptrs[256];
  if (nrhs > 256) { printf("fineD_mrhs: nrhs > 256\n"); exit(1); }
  for (int r = 0; r < nrhs; ++r)
    yp_host_ptrs[r] = g_dw_yp_buf + (long)r * Nin5 * Nst_pad;
  static real_t** ypd = nullptr; static int cy = 0;
  if (nrhs > cy) { if (ypd) cudaFree(ypd);
    CHECK(cudaMalloc((void**)&ypd, sizeof(real_t*) * nrhs)); cy = nrhs; }
  // skip the H2D re-upload of the stable scratch pointers (see lu_scratch).
  static real_t* yp_lst = nullptr; static int yp_ln = 0, yp_li = 0, yp_lp = 0;
  if (g_dw_yp_buf != yp_lst || nrhs != yp_ln || Nin5 != yp_li || Nst_pad != yp_lp) {
    CHECK(cudaMemcpy(ypd, yp_host_ptrs, sizeof(real_t*) * nrhs, cudaMemcpyHostToDevice));
    yp_lst = g_dw_yp_buf; yp_ln = nrhs; yp_li = Nin5; yp_lp = Nst_pad;
  }

  // 1) 5dir: v = D5 w, yp = aux  (b/c read from __constant__ memory)
  int blockSize = VECTOR_LENGTH;
  int gridSize5 = (Nst_pad * NVC + blockSize - 1) / blockSize;
  mult_dw5din_5dir_mrhs_dev<<<gridSize5, blockSize>>>(
      vdev, ypd, wdev, nrhs, mq, M0, Ns, alpha, Nst_pad);

  // 2) hopb: v += U-hop(yp)   (flag=1 -> accumulate into v)
  int gridSizeH = (Nst_pad * Ns + blockSize - 1) / blockSize;
  mult_dw5din_hopb_mrhs_dev<<<gridSizeH, blockSize>>>(
      vdev, u_dev, ypd, nrhs, Ns,
      bc[0], bc[1], bc[2], bc[3],
      Nsize[0], Nsize[1], Nsize[2], Nsize[3],
      do_comm[0], do_comm[1], do_comm[2], do_comm[3],
      1, Nst_pad);

  afield_dd_kernel_sync();
}

//====================================================================
// MRHS 5d LU-inverse kernels (Prec = C^{-1} and its dagger).  Pure 5d-local
// serial recurrence over s, one thread per (site,ivc), NO gauge.  Reused
// verbatim from the production single-RHS kernels with an inner rhs loop.
// e/f/dpinv/dm are per-operator (mass-dependent) device coefficient arrays.
//====================================================================

// --- L^{-1} (forward substitution) ---
__global__ void mult_dw5din_Linv_mrhs_dev(
    real_t* const* __restrict__ vp_arr, real_t* const* __restrict__ wp_arr,
    int nrhs, int Ns, int Nin5,
    const real_t* __restrict__ e, const real_t* __restrict__ dpinv,
    const real_t* __restrict__ dm, int Nst_pad)
{
  const int ist = blockIdx.x * blockDim.x + threadIdx.x;
  const int GridSize = blockDim.x * gridDim.x;
  for (int idx = ist; idx < Nst_pad * NVC; idx += GridSize) {
    int idx2 = idx / NWP, idx_in = idx % NWP, ivc = idx2 % NVC, idx_out = idx2 / NVC;
    int site = idx_in + NWP * idx_out;
    for (int r = 0; r < nrhs; ++r) {
      real_t* __restrict__ vp = vp_arr[r];
      real_t* __restrict__ wp = wp_arr[r];
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
        vt1 += a * (xt1 + xt3); vt2 += a * (xt2 + xt4);
        vt3 += a * (xt3 + xt1); vt4 += a * (xt4 + xt2);
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
        vt1 += a * (xt1 + xt3); vt2 += a * (xt2 + xt4);
        vt3 += a * (xt3 + xt1); vt4 += a * (xt4 + xt2);
        vt1 += -0.5 * (yt1 - yt3); vt2 += -0.5 * (yt2 - yt4);
        vt3 += -0.5 * (yt3 - yt1); vt4 += -0.5 * (yt4 - yt2);
        vp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)] = vt1;
        vp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)] = vt2;
        vp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)] = vt3;
        vp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)] = vt4;
      }
    }
  }
}

// --- U^{-1} (backward substitution) ---
__global__ void mult_dw5din_Uinv_mrhs_dev(
    real_t* const* __restrict__ vp_arr, real_t* const* __restrict__ wp_arr,
    int nrhs, int Ns, int Nin5,
    const real_t* __restrict__ f, const real_t* __restrict__ dpinv,
    const real_t* __restrict__ dm, int Nst_pad)
{
  const int ist = blockIdx.x * blockDim.x + threadIdx.x;
  const int GridSize = blockDim.x * gridDim.x;
  for (int idx = ist; idx < Nst_pad * NVC; idx += GridSize) {
    int idx2 = idx / NWP, idx_in = idx % NWP, ivc = idx2 % NVC, idx_out = idx2 / NVC;
    int site = idx_in + NWP * idx_out;
    for (int r = 0; r < nrhs; ++r) {
      real_t* __restrict__ vp = vp_arr[r];
      real_t* __restrict__ wp = wp_arr[r];
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
      yt1 = 0.5 * (vt1 + vt3); yt2 = 0.5 * (vt2 + vt4);
      yt3 = 0.5 * (vt3 + vt1); yt4 = 0.5 * (vt4 + vt2);
      for (int is = Ns - 2; is >= 0; --is) {
        xt1 = vt1; xt2 = vt2; xt3 = vt3; xt4 = vt4;
        vt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)];
        vt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)];
        vt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)];
        vt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)];
        real_t a = real_t(0.5) * dm[is];
        vt1 += a * (xt1 - xt3); vt2 += a * (xt2 - xt4);
        vt3 += a * (xt3 - xt1); vt4 += a * (xt4 - xt2);
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
}

// --- Ldag^{-1} ---
__global__ void mult_dw5din_Ldaginv_mrhs_dev(
    real_t* const* __restrict__ vp_arr, real_t* const* __restrict__ wp_arr,
    int nrhs, int Ns, int Nin5,
    const real_t* __restrict__ e, const real_t* __restrict__ dpinv,
    const real_t* __restrict__ dm, int Nst_pad)
{
  const int ist = blockIdx.x * blockDim.x + threadIdx.x;
  const int GridSize = blockDim.x * gridDim.x;
  for (int idx = ist; idx < Nst_pad * NVC; idx += GridSize) {
    int idx2 = idx / NWP, idx_in = idx % NWP, ivc = idx2 % NVC, idx_out = idx2 / NVC;
    int site = idx_in + NWP * idx_out;
    for (int r = 0; r < nrhs; ++r) {
      real_t* __restrict__ vp = vp_arr[r];
      real_t* __restrict__ wp = wp_arr[r];
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
      yt1 = 0.5 * (vt1 - vt3); yt2 = 0.5 * (vt2 - vt4);
      yt3 = 0.5 * (vt3 - vt1); yt4 = 0.5 * (vt4 - vt2);
      for (int is = Ns - 2; is >= 0; --is) {
        xt1 = vt1; xt2 = vt2; xt3 = vt3; xt4 = vt4;
        vt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)];
        vt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)];
        vt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)];
        vt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)];
        real_t a = real_t(0.5) * dm[is + 1] * dpinv[is];
        vt1 += a * (xt1 + xt3); vt2 += a * (xt2 + xt4);
        vt3 += a * (xt3 + xt1); vt4 += a * (xt4 + xt2);
        real_t eis = e[is];
        vt1 += -eis * yt1; vt2 += -eis * yt2; vt3 += -eis * yt3; vt4 += -eis * yt4;
        vp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)] = vt1;
        vp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)] = vt2;
        vp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)] = vt3;
        vp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)] = vt4;
      }
    }
  }
}

// --- Udag^{-1} (forward substitution) ---
__global__ void mult_dw5din_Udaginv_mrhs_dev(
    real_t* const* __restrict__ vp_arr, real_t* const* __restrict__ wp_arr,
    int nrhs, int Ns, int Nin5,
    const real_t* __restrict__ f, const real_t* __restrict__ dpinv,
    const real_t* __restrict__ dm, int Nst_pad)
{
  const int ist = blockIdx.x * blockDim.x + threadIdx.x;
  const int GridSize = blockDim.x * gridDim.x;
  for (int idx = ist; idx < Nst_pad * NVC; idx += GridSize) {
    int idx2 = idx / NWP, idx_in = idx % NWP, ivc = idx2 % NVC, idx_out = idx2 / NVC;
    int site = idx_in + NWP * idx_out;
    for (int r = 0; r < nrhs; ++r) {
      real_t* __restrict__ vp = vp_arr[r];
      real_t* __restrict__ wp = wp_arr[r];
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
        vt1 += a * (xt1 - xt3); vt2 += a * (xt2 - xt4);
        vt3 += a * (xt3 - xt1); vt4 += a * (xt4 - xt2);
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
        vt1 += a * (xt1 - xt3); vt2 += a * (xt2 - xt4);
        vt3 += a * (xt3 - xt1); vt4 += a * (xt4 - xt2);
        vt1 += -0.5 * (yt1 + yt3); vt2 += -0.5 * (yt2 + yt4);
        vt3 += -0.5 * (yt3 + yt1); vt4 += -0.5 * (yt4 + yt2);
        real_t aa = dpinv[Ns - 1];
        vt1 *= aa; vt2 *= aa; vt3 *= aa; vt4 *= aa;
        vp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)] = vt1;
        vp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)] = vt2;
        vp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)] = vt3;
        vp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)] = vt4;
      }
    }
  }
}

// scratch for the 2-stage LU solves (one Nin5*Nst_pad field per rhs)
static real_t* g_dw_lu_buf = nullptr;
static long    g_dw_lu_cap = 0;
static real_t** g_dw_lu_ptrs = nullptr;
static int      g_dw_lu_pcap = 0;

static real_t** lu_scratch(int nrhs, int Nin5, int Nst_pad)
{
  long need = (long)nrhs * Nin5 * Nst_pad;
  if (need > g_dw_lu_cap) {
    if (g_dw_lu_buf) cudaFree(g_dw_lu_buf);
    CHECK(cudaMalloc((void**)&g_dw_lu_buf, sizeof(real_t) * need));
    g_dw_lu_cap = need;
  }
  real_t* hp[256];
  if (nrhs > 256) { printf("lu_scratch: nrhs > 256\n"); exit(1); }
  for (int r = 0; r < nrhs; ++r) hp[r] = g_dw_lu_buf + (long)r * Nin5 * Nst_pad;
  if (nrhs > g_dw_lu_pcap) { if (g_dw_lu_ptrs) cudaFree(g_dw_lu_ptrs);
    CHECK(cudaMalloc((void**)&g_dw_lu_ptrs, sizeof(real_t*) * nrhs)); g_dw_lu_pcap = nrhs; }
  // the scratch pointer values depend only on (buf, nrhs, Nin5, Nst_pad); skip
  // the H2D re-upload when none of those changed (stable across the solve).
  static real_t* lst_buf = nullptr; static int lst_n = 0, lst_in5 = 0, lst_pad = 0;
  if (g_dw_lu_buf != lst_buf || nrhs != lst_n || Nin5 != lst_in5 || Nst_pad != lst_pad) {
    CHECK(cudaMemcpy(g_dw_lu_ptrs, hp, sizeof(real_t*) * nrhs, cudaMemcpyHostToDevice));
    lst_buf = g_dw_lu_buf; lst_n = nrhs; lst_in5 = Nin5; lst_pad = Nst_pad;
  }
  return g_dw_lu_ptrs;
}

// Prec = C^{-1} = U^{-1} L^{-1}  (mass operator's LU; e/f/dpinv/dm host bases).
void finePrec_mrhs(real_t* const* v_host, real_t* const* w_host, int nrhs, int Ns,
                   real_t* e_host, real_t* f_host, real_t* dpinv_host, real_t* dm_host,
                   int* Nsize)
{
  int Nst = Nsize[0]*Nsize[1]*Nsize[2]*Nsize[3];
  int Nst_pad = ceil_nwp(Nst);
  int Nin5 = NVCD * Ns;
  static real_t **vd = nullptr, **wd = nullptr; static int cv = 0, cw = 0;
  real_t** vdev = map_upload(v_host, nrhs, vd, cv);
  real_t** wdev = map_upload(w_host, nrhs, wd, cw);
  real_t** sc   = lu_scratch(nrhs, Nin5, Nst_pad);
  real_t* e_dev = (real_t*)dev_ptr(e_host);
  real_t* f_dev = (real_t*)dev_ptr(f_host);
  real_t* dpinv_dev = (real_t*)dev_ptr(dpinv_host);
  real_t* dm_dev    = (real_t*)dev_ptr(dm_host);
  int blockSize = VECTOR_LENGTH;
  int gridSize  = (Nst_pad * NVC + blockSize - 1) / blockSize;
  // L^{-1} w -> sc ; U^{-1} sc -> v
  mult_dw5din_Linv_mrhs_dev<<<gridSize, blockSize>>>(sc, wdev, nrhs, Ns, Nin5, e_dev, dpinv_dev, dm_dev, Nst_pad);
  mult_dw5din_Uinv_mrhs_dev<<<gridSize, blockSize>>>(vdev, sc, nrhs, Ns, Nin5, f_dev, dpinv_dev, dm_dev, Nst_pad);
  afield_dd_kernel_sync();
}

// Precdag = (C^{-1})^dag = Ldag^{-1} Udag^{-1}.
void finePrecdag_mrhs(real_t* const* v_host, real_t* const* w_host, int nrhs, int Ns,
                      real_t* e_host, real_t* f_host, real_t* dpinv_host, real_t* dm_host,
                      int* Nsize)
{
  int Nst = Nsize[0]*Nsize[1]*Nsize[2]*Nsize[3];
  int Nst_pad = ceil_nwp(Nst);
  int Nin5 = NVCD * Ns;
  static real_t **vd = nullptr, **wd = nullptr; static int cv = 0, cw = 0;
  real_t** vdev = map_upload(v_host, nrhs, vd, cv);
  real_t** wdev = map_upload(w_host, nrhs, wd, cw);
  real_t** sc   = lu_scratch(nrhs, Nin5, Nst_pad);
  real_t* e_dev = (real_t*)dev_ptr(e_host);
  real_t* f_dev = (real_t*)dev_ptr(f_host);
  real_t* dpinv_dev = (real_t*)dev_ptr(dpinv_host);
  real_t* dm_dev    = (real_t*)dev_ptr(dm_host);
  int blockSize = VECTOR_LENGTH;
  int gridSize  = (Nst_pad * NVC + blockSize - 1) / blockSize;
  // Udag^{-1} w -> sc ; Ldag^{-1} sc -> v
  mult_dw5din_Udaginv_mrhs_dev<<<gridSize, blockSize>>>(sc, wdev, nrhs, Ns, Nin5, f_dev, dpinv_dev, dm_dev, Nst_pad);
  mult_dw5din_Ldaginv_mrhs_dev<<<gridSize, blockSize>>>(vdev, sc, nrhs, Ns, Nin5, e_dev, dpinv_dev, dm_dev, Nst_pad);
  afield_dd_kernel_sync();
}

//====================================================================
// MRHS 5dirdag (adjoint 5d block) + gm5, for the batched Ddag.
//====================================================================
__global__
void mult_dw5din_5dirdag_mrhs_dev(
    real_t* const* __restrict__ vp_arr, real_t* const* __restrict__ yp_arr,
    real_t* const* __restrict__ wp_arr,
    int nrhs, real_t mq, real_t M0, int Ns, real_t alpha, int Nst_pad)
{
  const real_t* b_arr = mrhs_cb;   // module __constant__ b/c (see 5dir kernel)
  const real_t* c_arr = mrhs_cc;
  const int Nin5 = NVCD * Ns;
  const int ist = blockIdx.x * blockDim.x + threadIdx.x;
  const int GridSize = blockDim.x * gridDim.x;
  for (int idx = ist; idx < Nst_pad * NVC; idx += GridSize) {
    int idx2 = idx / NWP, idx_in = idx % NWP, ivc = idx2 % NVC, idx_out = idx2 / NVC;
    int site = idx_in + NWP * idx_out;
    for (int r = 0; r < nrhs; ++r) {
      real_t* __restrict__ vp = vp_arr[r];
      real_t* __restrict__ yp = yp_arr[r];
      real_t* __restrict__ wp = wp_arr[r];
      for (int is = 0; is < Ns; ++is) {
        real_t B1 = b_arr[is] * (4.0 - M0) + 1.0;
        real_t a1 = -0.5 * b_arr[is];
        real_t wt1, wt2, wt3, wt4, yt1, yt2, yt3, yt4;
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
          real_t f1 = 0.5 * (1.0 + alpha), f2 = 0.5 * (-1.0 + alpha);
          vt1 = f1 * s1 + f2 * s3; vt2 = f1 * s2 + f2 * s4;
          vt3 = f1 * s3 + f2 * s1; vt4 = f1 * s4 + f2 * s2;
        } else if (is == Ns - 1) {
          real_t f1 = 0.5 * (1.0 + alpha), f2 = 0.5 * (1.0 - alpha);
          vt1 = f1 * s1 + f2 * s3; vt2 = f1 * s2 + f2 * s4;
          vt3 = f1 * s3 + f2 * s1; vt4 = f1 * s4 + f2 * s2;
        } else {
          vt1 = alpha * s1; vt2 = alpha * s2; vt3 = alpha * s3; vt4 = alpha * s4;
        }
        int is_up = (is + 1) % Ns;
        real_t C1 = c_arr[is_up] * (4.0 - M0) - 1.0;
        real_t aup = -0.5 * c_arr[is_up];
        wt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD * is_up), site)];
        wt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD * is_up), site)];
        wt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD * is_up), site)];
        wt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD * is_up), site)];
        yt1 = yp[IDX2(Nin5, (ID1 + ivc + NVCD * is_up), site)];
        yt2 = yp[IDX2(Nin5, (ID2 + ivc + NVCD * is_up), site)];
        yt3 = yp[IDX2(Nin5, (ID3 + ivc + NVCD * is_up), site)];
        yt4 = yp[IDX2(Nin5, (ID4 + ivc + NVCD * is_up), site)];
        real_t xu1 = C1 * wt1 + aup * yt3, xu2 = C1 * wt2 + aup * yt4;
        real_t xu3 = C1 * wt3 + aup * yt1, xu4 = C1 * wt4 + aup * yt2;
        real_t Fup = 0.5 * alpha;
        if (is == Ns - 1) Fup = -0.5 * mq;
        vt1 += Fup * (xu1 + xu3); vt2 += Fup * (xu2 + xu4);
        vt3 += Fup * (xu3 + xu1); vt4 += Fup * (xu4 + xu2);
        int is_dn = (is - 1 + Ns) % Ns;
        real_t C2 = c_arr[is_dn] * (4.0 - M0) - 1.0;
        real_t adn = -0.5 * c_arr[is_dn];
        wt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD * is_dn), site)];
        wt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD * is_dn), site)];
        wt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD * is_dn), site)];
        wt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD * is_dn), site)];
        yt1 = yp[IDX2(Nin5, (ID1 + ivc + NVCD * is_dn), site)];
        yt2 = yp[IDX2(Nin5, (ID2 + ivc + NVCD * is_dn), site)];
        yt3 = yp[IDX2(Nin5, (ID3 + ivc + NVCD * is_dn), site)];
        yt4 = yp[IDX2(Nin5, (ID4 + ivc + NVCD * is_dn), site)];
        real_t xd1 = C2 * wt1 + adn * yt3, xd2 = C2 * wt2 + adn * yt4;
        real_t xd3 = C2 * wt3 + adn * yt1, xd4 = C2 * wt4 + adn * yt2;
        real_t Fdn = 0.5 * alpha;
        if (is == 0) Fdn = -0.5 * mq;
        vt1 += Fdn * (xd1 - xd3); vt2 += Fdn * (xd2 - xd4);
        vt3 += Fdn * (xd3 - xd1); vt4 += Fdn * (xd4 - xd2);
        vp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)] = vt1;
        vp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)] = vt2;
        vp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)] = vt3;
        vp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)] = vt4;
      }
    }
  }
}

__global__
void mult_dw5din_gm5_mrhs_dev(real_t* const* __restrict__ vp_arr,
                              real_t* const* __restrict__ wp_arr,
                              int nrhs, int Ns, int Nst)
{
  const int Nin5 = NVCD * Ns;
  const int site = blockIdx.x * blockDim.x + threadIdx.x;
  if (site < Nst) {
    for (int r = 0; r < nrhs; ++r) {
      real_t* __restrict__ vp = vp_arr[r];
      real_t* __restrict__ wp = wp_arr[r];
      for (int is = 0; is < Ns; ++is) {
        real_t wt[NVCD];
        for (int ivcd = 0; ivcd < NVCD; ++ivcd)
          wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
        for (int ivc = 0; ivc < NVC; ++ivc) {
          vp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)] = wt[ID3 + ivc];
          vp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)] = wt[ID4 + ivc];
          vp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)] = wt[ID1 + ivc];
          vp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)] = wt[ID2 + ivc];
        }
      }
    }
  }
}

// second 5d scratch (Ddag needs two temporaries: gm5(w) and hopb(gm5(w)))
static real_t* g_dw_lu2_buf = nullptr;  static long g_dw_lu2_cap = 0;
static real_t** g_dw_lu2_ptrs = nullptr; static int g_dw_lu2_pcap = 0;
static real_t** lu_scratch2(int nrhs, int Nin5, int Nst_pad)
{
  long need = (long)nrhs * Nin5 * Nst_pad;
  if (need > g_dw_lu2_cap) { if (g_dw_lu2_buf) cudaFree(g_dw_lu2_buf);
    CHECK(cudaMalloc((void**)&g_dw_lu2_buf, sizeof(real_t) * need)); g_dw_lu2_cap = need; }
  real_t* hp[256];
  if (nrhs > 256) { printf("lu_scratch2: nrhs > 256\n"); exit(1); }
  for (int r = 0; r < nrhs; ++r) hp[r] = g_dw_lu2_buf + (long)r * Nin5 * Nst_pad;
  if (nrhs > g_dw_lu2_pcap) { if (g_dw_lu2_ptrs) cudaFree(g_dw_lu2_ptrs);
    CHECK(cudaMalloc((void**)&g_dw_lu2_ptrs, sizeof(real_t*) * nrhs)); g_dw_lu2_pcap = nrhs; }
  static real_t* l2_buf = nullptr; static int l2_n = 0, l2_in5 = 0, l2_pad = 0;
  if (g_dw_lu2_buf != l2_buf || nrhs != l2_n || Nin5 != l2_in5 || Nst_pad != l2_pad) {
    CHECK(cudaMemcpy(g_dw_lu2_ptrs, hp, sizeof(real_t*) * nrhs, cudaMemcpyHostToDevice));
    l2_buf = g_dw_lu2_buf; l2_n = nrhs; l2_in5 = Nin5; l2_pad = Nst_pad;
  }
  return g_dw_lu2_ptrs;
}

// Ddag w = 5dirdag( w, hopb( gm5 w ) ) : single-GPU no-comm path.
// b/c are read from __constant__ memory inside the kernels (no upload here).
void fineDdag_mrhs(real_t* const* v_host, real_t* const* w_host, real_t* u_field_host,
                   int nrhs, real_t mq, real_t M0, int Ns, real_t alpha,
                   int* Nsize, int* bc, int* do_comm)
{
  int Nst = Nsize[0]*Nsize[1]*Nsize[2]*Nsize[3];
  int Nst_pad = ceil_nwp(Nst);
  int Nin5 = NVCD * Ns;
  real_t* u_dev = (real_t*)dev_ptr(u_field_host);
  static real_t **vd = nullptr, **wd = nullptr; static int cv = 0, cw = 0;
  real_t** vdev = map_upload(v_host, nrhs, vd, cv);
  real_t** wdev = map_upload(w_host, nrhs, wd, cw);
  real_t** t1 = lu_scratch (nrhs, Nin5, Nst_pad);   // gm5(w)
  real_t** t2 = lu_scratch2(nrhs, Nin5, Nst_pad);   // hopb(gm5(w))

  int blockSize = VECTOR_LENGTH;
  int gridV = (Nst_pad * NVC + blockSize - 1) / blockSize;
  int gridS = (Nst + blockSize - 1) / blockSize;
  int gridH = (Nst_pad * Ns + blockSize - 1) / blockSize;

  // t1 = gm5 w
  mult_dw5din_gm5_mrhs_dev<<<gridS, blockSize>>>(t1, wdev, nrhs, Ns, Nst);
  // t2 = U-hopb(t1)  (flag=0 -> overwrite)
  mult_dw5din_hopb_mrhs_dev<<<gridH, blockSize>>>(
      t2, u_dev, t1, nrhs, Ns, bc[0], bc[1], bc[2], bc[3],
      Nsize[0], Nsize[1], Nsize[2], Nsize[3],
      do_comm[0], do_comm[1], do_comm[2], do_comm[3], 0, Nst_pad);
  // v = 5dirdag( w, t2 )  (b/c from __constant__)
  mult_dw5din_5dirdag_mrhs_dev<<<gridV, blockSize>>>(
      vdev, t2, wdev, nrhs, mq, M0, Ns, alpha, Nst_pad);

  afield_dd_kernel_sync();
}

//====================================================================
// On-device coarse-vector layout repack (replaces the V-cycle host repack).
//   MRHS coarse_dev:  mr = 2*(c + s*(j + nbasis*b)) interleaved re/im, j = 2i+ch
//   production coarse AField:  IDX2(Ncoarse_in=4*nvec, inr=2*(ch+2i), b),
//                              pr = (b%NWP)+NWP*(inr + Ncoarse_in*(b/NWP)), NWP=32
// One thread per (j, block); pure index permutation, no arithmetic, no host copy.
//====================================================================
__global__ void coarse_repack_m2p_dev(real_t* __restrict__ prod,
                                       const real_t* __restrict__ coarse,
                                       int c, int s, int nbasis, int nvec,
                                       long coarse_nvol)
{
  long idx = (long)blockIdx.x * blockDim.x + threadIdx.x;
  long total = (long)nbasis * coarse_nvol;
  if (idx >= total) return;
  int  j = (int)(idx / coarse_nvol);
  long b = idx % coarse_nvol;
  int  i = j / 2, ch = j % 2;
  const long Ncoarse_in = 4L * nvec;
  const int  MG_NWP = 32;
  long inr = 2L * (ch + 2 * i);
  long pr  = (b % MG_NWP) + MG_NWP * (inr     + Ncoarse_in * (b / MG_NWP));
  long pi  = (b % MG_NWP) + MG_NWP * (inr + 1 + Ncoarse_in * (b / MG_NWP));
  long mr  = 2L * (c + s * (j + (long)nbasis * b));
  prod[pr] = coarse[mr];
  prod[pi] = coarse[mr + 1];
}

__global__ void coarse_repack_p2m_dev(real_t* __restrict__ coarse,
                                      const real_t* __restrict__ prod,
                                      int c, int s, int nbasis, int nvec,
                                      long coarse_nvol)
{
  long idx = (long)blockIdx.x * blockDim.x + threadIdx.x;
  long total = (long)nbasis * coarse_nvol;
  if (idx >= total) return;
  int  j = (int)(idx / coarse_nvol);
  long b = idx % coarse_nvol;
  int  i = j / 2, ch = j % 2;
  const long Ncoarse_in = 4L * nvec;
  const int  MG_NWP = 32;
  long inr = 2L * (ch + 2 * i);
  long pr  = (b % MG_NWP) + MG_NWP * (inr     + Ncoarse_in * (b / MG_NWP));
  long pi  = (b % MG_NWP) + MG_NWP * (inr + 1 + Ncoarse_in * (b / MG_NWP));
  long mr  = 2L * (c + s * (j + (long)nbasis * b));
  coarse[mr]     = prod[pr];
  coarse[mr + 1] = prod[pi];
}

// coarse_dev (MRHS, column c) -> prod coarse AField device buffer, on device.
void coarse_repack_m2p(real_t* prod_field_host, const real_t* coarse_dev,
                       int c, int s, int nbasis, int nvec, long coarse_nvol)
{
  real_t* prod = (real_t*)dev_ptr(prod_field_host);
  long total = (long)nbasis * coarse_nvol;
  int bs = VECTOR_LENGTH, gs = (int)((total + bs - 1) / bs);
  coarse_repack_m2p_dev<<<gs, bs>>>(prod, coarse_dev, c, s, nbasis, nvec, coarse_nvol);
  afield_dd_kernel_sync();
}

// prod coarse AField (solution) -> coarse_dev (MRHS, column c), on device.
void coarse_repack_p2m(real_t* coarse_dev, real_t* prod_field_host,
                       int c, int s, int nbasis, int nvec, long coarse_nvol)
{
  real_t* prod = (real_t*)dev_ptr(prod_field_host);
  long total = (long)nbasis * coarse_nvol;
  int bs = VECTOR_LENGTH, gs = (int)((total + bs - 1) / bs);
  coarse_repack_p2m_dev<<<gs, bs>>>(coarse_dev, prod, c, s, nbasis, nvec, coarse_nvol);
  afield_dd_kernel_sync();
}

} // anonymous namespace
