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
#include <cstdlib>   // getenv (FP16_PROBE dynamic-range probe)
#include <cstring>   // memcpy (host bit-reinterpret)
#include <cmath>     // log10

namespace {  // anonymous: lives inside the BridgeACC TU (dev_ptr/macros visible)

//====================================================================
// FP16 dynamic-range probe (gated by env FP16_PROBE).  Scans one rhs of a
// fine field and reports max|.|, min nonzero |.|, and the decade span
// log10(max/min).  Decision metric for FP16-vs-BF16 storage of the
// preconditioner fields: FP16 normal range spans ~9.03 decades (65504 down to
// 6.10e-5), so after an optimal global rescale a field is underflow-free in
// FP16 iff its decade span < ~9; wider spans need BF16's FP32 exponent.
// Kept here (not perf path) so the same metric is read at every fine-apply
// boundary during a real solve.  No effect unless FP16_PROBE is set.
//====================================================================
__global__ void fp16_probe_kernel(const real_t* __restrict__ x, long n,
                                  unsigned* __restrict__ out)
{
  // out[0]=max bits, out[1]=min bits, out[2]=#(0<|x|<fp16min), out[3]=#nonzero
  long i = (long)blockIdx.x * blockDim.x + threadIdx.x;
  long stride = (long)blockDim.x * gridDim.x;
  for (; i < n; i += stride) {
    real_t a = fabsf(x[i]);
    if (a > 0.0f) {
      unsigned b = __float_as_uint(a);   // monotonic for non-negative floats
      atomicMax(&out[0], b);
      atomicMin(&out[1], b);
      atomicAdd(&out[3], 1u);
      if (a < 6.103515625e-05f) atomicAdd(&out[2], 1u);
    }
  }
}

static void fp16_probe(const char* tag, real_t* field_host, int Nin5, int Nst_pad)
{
  static int en = -1;
  static long cap = 0, seen = 0;
  if (en < 0) { const char* e = getenv("FP16_PROBE"); en = e ? 1 : 0;
                const char* c = getenv("FP16_PROBE_CAP"); cap = c ? atol(c) : 400; }
  if (!en || seen >= cap) return;
  ++seen;
  real_t* x = (real_t*)dev_ptr(field_host);
  long n = (long)Nin5 * Nst_pad;
  static unsigned* d = nullptr;
  if (!d) CHECK(cudaMalloc((void**)&d, 4 * sizeof(unsigned)));
  unsigned init[4] = { 0u, 0xFFFFFFFFu, 0u, 0u };
  CHECK(cudaMemcpy(d, init, 4 * sizeof(unsigned), cudaMemcpyHostToDevice));
  int bs = 256;
  int gs = (int)((n + bs - 1) / bs); if (gs > 1024) gs = 1024;
  fp16_probe_kernel<<<gs, bs>>>(x, n, d);
  unsigned h[4];
  CHECK(cudaDeviceSynchronize());
  CHECK(cudaMemcpy(h, d, 4 * sizeof(unsigned), cudaMemcpyDeviceToHost));
  float mx, mn; memcpy(&mx, &h[0], 4); memcpy(&mn, &h[1], 4);
  double dec = (h[3] && mn > 0.0f) ? log10((double)mx / (double)mn) : 0.0;
  printf("[FP16PROBE] %-16s max=%.3e minnz=%.3e decades=%.2f  uflow=%u/%u\n",
         tag, mx, mn, dec, h[2], h[3]);
  fflush(stdout);
}

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
// Gated per-wrapper kernel-time accounting (env FINE_TIME).  The batched mrhs
// kernels are invisible to nsys/ncu in this CUDA13 build, so time the four fine
// A-apply pieces with cudaEvents to get the V-cycle breakdown in a real solve.
// k: 0=fineD 1=fineDdag 2=finePrec 3=finePrecdag.  No effect unless FINE_TIME set.
//====================================================================
static double g_ft_ms[4] = {0,0,0,0};
static long   g_ft_n [4] = {0,0,0,0};
static cudaEvent_t g_ft_e0 = nullptr, g_ft_e1 = nullptr;
static int fine_time_on() { static int e = -1; if (e < 0) e = getenv("FINE_TIME") ? 1 : 0; return e; }
static void fine_time_begin() {
  if (!fine_time_on()) return;
  if (!g_ft_e0) { cudaEventCreate(&g_ft_e0); cudaEventCreate(&g_ft_e1); }
  cudaEventRecord(g_ft_e0);
}
static void fine_time_end(int k) {
  if (!fine_time_on()) return;
  cudaEventRecord(g_ft_e1); cudaEventSynchronize(g_ft_e1);
  float ms = 0; cudaEventElapsedTime(&ms, g_ft_e0, g_ft_e1);
  g_ft_ms[k] += ms; g_ft_n[k]++;
  long tot = g_ft_n[0] + g_ft_n[1] + g_ft_n[2] + g_ft_n[3];
  if (tot % 4000 == 0) {
    double s = g_ft_ms[0] + g_ft_ms[1] + g_ft_ms[2] + g_ft_ms[3];
    printf("[FINE_TIME] fineD=%.0fms(%ld) fineDdag=%.0f(%ld) finePrec=%.0f(%ld) "
           "finePrecdag=%.0f(%ld) | Prec share=%.0f%%\n",
           g_ft_ms[0], g_ft_n[0], g_ft_ms[1], g_ft_n[1], g_ft_ms[2], g_ft_n[2],
           g_ft_ms[3], g_ft_n[3], s>0 ? 100.0*(g_ft_ms[2]+g_ft_ms[3])/s : 0.0);
    fflush(stdout);
  }
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

  fp16_probe("fineD.in", w_host[0], Nin5, Nst_pad);
  fine_time_begin();

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
  fine_time_end(0);
  fp16_probe("fineD.out", v_host[0], Nin5, Nst_pad);
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

// Fused-LU buffer bound. The fused C^{-1}/C^{-dag} kernels keep the whole
// L^{-1}w (resp. Udag^{-1}w) Ls-vector for one (site,ivc,rhs) in thread-local
// storage between the forward and backward sweep, so the intermediate field
// never round-trips through DRAM. Buffer sized to this bound (covers all
// realistic Mobius Ls; our tests use Ns=8). Larger Ls falls back to the
// original two-kernel path in the drivers below.
#define MRHS_FUSE_LS_MAX 32

// --- C^{-1} = U^{-1} L^{-1}, fused (no global intermediate, one launch) ---
// gm5_arr (optional): if non-null, also store gamma5*(C^{-1}w) = the spin-swapped
// result (spin0<->2, spin1<->3). This pre-supplies the gm5(w) that fineDdag's hop
// needs, killing the separate gm5 kernel + its full-field read (the output is
// already in registers, so it is just an extra permuted store).
// Templated on compile-time NS (= Ls): the Ls buffer s[NS] is then promoted to
// REGISTERS (runtime-Ns forces it to local memory -> ~2.3x DRAM over-traffic,
// the dominant cost; ncu: local loads use 1 of 32 bytes/sector). Dispatched per
// Ns by the wrapper; uncommon Ns falls back to the non-fused Linv/Uinv path.
template<int NS>
__global__ void mult_dw5din_Cinv_mrhs_dev(
    real_t* const* __restrict__ vp_arr, real_t* const* __restrict__ wp_arr,
    real_t* const* __restrict__ gm5_arr,
    int nrhs, int Nin5,
    const real_t* __restrict__ e, const real_t* __restrict__ f,
    const real_t* __restrict__ dpinv, const real_t* __restrict__ dm, int Nst_pad)
{
  const int ist = blockIdx.x * blockDim.x + threadIdx.x;
  const int GridSize = blockDim.x * gridDim.x;
  for (int idx = ist; idx < Nst_pad * NVC; idx += GridSize) {
    int idx2 = idx / NWP, idx_in = idx % NWP, ivc = idx2 % NVC, idx_out = idx2 / NVC;
    int site = idx_in + NWP * idx_out;
    for (int r = 0; r < nrhs; ++r) {
      real_t* __restrict__ vp = vp_arr[r];
      real_t* __restrict__ wp = wp_arr[r];
      real_t* __restrict__ gp = gm5_arr ? gm5_arr[r] : nullptr;  // gamma5 copy out
      // thread-local L^{-1}w for this (site,ivc,rhs), all Ls -> registers (NS const)
      real_t s1[NS], s2[NS];
      real_t s3[NS], s4[NS];
      real_t vt1, vt2, vt3, vt4, yt1, yt2, yt3, yt4, xt1, xt2, xt3, xt4;
      // ---- forward sweep L^{-1}: w -> s (local) ----
      vt1 = wp[IDX2(Nin5, (ID1 + ivc), site)];
      vt2 = wp[IDX2(Nin5, (ID2 + ivc), site)];
      vt3 = wp[IDX2(Nin5, (ID3 + ivc), site)];
      vt4 = wp[IDX2(Nin5, (ID4 + ivc), site)];
      s1[0] = vt1; s2[0] = vt2; s3[0] = vt3; s4[0] = vt4;
      real_t e0 = e[0];
      yt1 = e0 * vt1; yt2 = e0 * vt2; yt3 = e0 * vt3; yt4 = e0 * vt4;
#pragma unroll
      for (int is = 1; is < NS - 1; ++is) {
        xt1 = vt1; xt2 = vt2; xt3 = vt3; xt4 = vt4;
        vt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)];
        vt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)];
        vt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)];
        vt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)];
        real_t a = real_t(0.5) * dm[is] * dpinv[is - 1];
        vt1 += a * (xt1 + xt3); vt2 += a * (xt2 + xt4);
        vt3 += a * (xt3 + xt1); vt4 += a * (xt4 + xt2);
        s1[is] = vt1; s2[is] = vt2; s3[is] = vt3; s4[is] = vt4;
        real_t eis = e[is];
        yt1 += eis * vt1; yt2 += eis * vt2; yt3 += eis * vt3; yt4 += eis * vt4;
      }
      {
        int is = NS - 1;
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
        s1[is] = vt1; s2[is] = vt2; s3[is] = vt3; s4[is] = vt4;
      }
      // ---- backward sweep U^{-1}: s (local) -> v ----
      int is0 = NS - 1;
      real_t a0 = dpinv[NS - 1];
      vt1 = a0 * s1[is0]; vt2 = a0 * s2[is0]; vt3 = a0 * s3[is0]; vt4 = a0 * s4[is0];
      vp[IDX2(Nin5, (ID1 + ivc + NVCD * is0), site)] = vt1;
      vp[IDX2(Nin5, (ID2 + ivc + NVCD * is0), site)] = vt2;
      vp[IDX2(Nin5, (ID3 + ivc + NVCD * is0), site)] = vt3;
      vp[IDX2(Nin5, (ID4 + ivc + NVCD * is0), site)] = vt4;
      if (gp) {  // gamma5: spin0<->2, spin1<->3
        gp[IDX2(Nin5, (ID1 + ivc + NVCD * is0), site)] = vt3;
        gp[IDX2(Nin5, (ID2 + ivc + NVCD * is0), site)] = vt4;
        gp[IDX2(Nin5, (ID3 + ivc + NVCD * is0), site)] = vt1;
        gp[IDX2(Nin5, (ID4 + ivc + NVCD * is0), site)] = vt2;
      }
      yt1 = 0.5 * (vt1 + vt3); yt2 = 0.5 * (vt2 + vt4);
      yt3 = 0.5 * (vt3 + vt1); yt4 = 0.5 * (vt4 + vt2);
#pragma unroll
      for (int is = NS - 2; is >= 0; --is) {
        xt1 = vt1; xt2 = vt2; xt3 = vt3; xt4 = vt4;
        vt1 = s1[is]; vt2 = s2[is]; vt3 = s3[is]; vt4 = s4[is];
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
        if (gp) {  // gamma5: spin0<->2, spin1<->3 (all Ls, not just is0)
          gp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)] = vt3;
          gp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)] = vt4;
          gp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)] = vt1;
          gp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)] = vt2;
        }
      }
    }
  }
}

// --- C^{-1} = U^{-1} L^{-1}, fused, BF16-STORAGE variant (Step 1 of the BF16
// smoother port).  Arithmetic is byte-for-byte the FP32 recurrence of
// mult_dw5din_Cinv_mrhs_dev above; the ONLY change is that the operand vectors
// (w in, v / gm5 out) live in global memory as __nv_bfloat16 (2 bytes).  This
// kernel is memory-bound, so halving the operand bytes is the ~2x lever on the
// fine-A apply (~65% of the V-cycle).  Accumulation stays in FP32 registers
// (ld bf16 -> promote -> fp32 compute -> demote -> st bf16); the LU coefficients
// e/f/dpinv/dm stay FP32.  The proven FP32 kernel is left untouched.
__device__ __forceinline__ float ld_bf16f(const __nv_bfloat16* __restrict__ p, long i)
{ return __bfloat162float(p[i]); }
__device__ __forceinline__ void  st_bf16f(__nv_bfloat16* __restrict__ p, long i, float v)
{ p[i] = __float2bfloat16(v); }

template<int NS>
__global__ void mult_dw5din_Cinv_mrhs_bf16_dev(
    __nv_bfloat16* const* __restrict__ vp_arr, __nv_bfloat16* const* __restrict__ wp_arr,
    __nv_bfloat16* const* __restrict__ gm5_arr,
    int nrhs, int Nin5,
    const real_t* __restrict__ e, const real_t* __restrict__ f,
    const real_t* __restrict__ dpinv, const real_t* __restrict__ dm, int Nst_pad)
{
  const int ist = blockIdx.x * blockDim.x + threadIdx.x;
  const int GridSize = blockDim.x * gridDim.x;
  for (int idx = ist; idx < Nst_pad * NVC; idx += GridSize) {
    int idx2 = idx / NWP, idx_in = idx % NWP, ivc = idx2 % NVC, idx_out = idx2 / NVC;
    int site = idx_in + NWP * idx_out;
    for (int r = 0; r < nrhs; ++r) {
      const __nv_bfloat16* __restrict__ wp = wp_arr[r];
      __nv_bfloat16* __restrict__ vp = vp_arr[r];
      __nv_bfloat16* __restrict__ gp = gm5_arr ? gm5_arr[r] : nullptr;
      real_t s1[NS], s2[NS];
      real_t s3[NS], s4[NS];
      real_t vt1, vt2, vt3, vt4, yt1, yt2, yt3, yt4, xt1, xt2, xt3, xt4;
      // ---- forward sweep L^{-1}: w -> s (local) ----
      vt1 = ld_bf16f(wp, IDX2(Nin5, (ID1 + ivc), site));
      vt2 = ld_bf16f(wp, IDX2(Nin5, (ID2 + ivc), site));
      vt3 = ld_bf16f(wp, IDX2(Nin5, (ID3 + ivc), site));
      vt4 = ld_bf16f(wp, IDX2(Nin5, (ID4 + ivc), site));
      s1[0] = vt1; s2[0] = vt2; s3[0] = vt3; s4[0] = vt4;
      real_t e0 = e[0];
      yt1 = e0 * vt1; yt2 = e0 * vt2; yt3 = e0 * vt3; yt4 = e0 * vt4;
#pragma unroll
      for (int is = 1; is < NS - 1; ++is) {
        xt1 = vt1; xt2 = vt2; xt3 = vt3; xt4 = vt4;
        vt1 = ld_bf16f(wp, IDX2(Nin5, (ID1 + ivc + NVCD * is), site));
        vt2 = ld_bf16f(wp, IDX2(Nin5, (ID2 + ivc + NVCD * is), site));
        vt3 = ld_bf16f(wp, IDX2(Nin5, (ID3 + ivc + NVCD * is), site));
        vt4 = ld_bf16f(wp, IDX2(Nin5, (ID4 + ivc + NVCD * is), site));
        real_t a = real_t(0.5) * dm[is] * dpinv[is - 1];
        vt1 += a * (xt1 + xt3); vt2 += a * (xt2 + xt4);
        vt3 += a * (xt3 + xt1); vt4 += a * (xt4 + xt2);
        s1[is] = vt1; s2[is] = vt2; s3[is] = vt3; s4[is] = vt4;
        real_t eis = e[is];
        yt1 += eis * vt1; yt2 += eis * vt2; yt3 += eis * vt3; yt4 += eis * vt4;
      }
      {
        int is = NS - 1;
        xt1 = vt1; xt2 = vt2; xt3 = vt3; xt4 = vt4;
        vt1 = ld_bf16f(wp, IDX2(Nin5, (ID1 + ivc + NVCD * is), site));
        vt2 = ld_bf16f(wp, IDX2(Nin5, (ID2 + ivc + NVCD * is), site));
        vt3 = ld_bf16f(wp, IDX2(Nin5, (ID3 + ivc + NVCD * is), site));
        vt4 = ld_bf16f(wp, IDX2(Nin5, (ID4 + ivc + NVCD * is), site));
        real_t a = real_t(0.5) * dm[is] * dpinv[is - 1];
        vt1 += a * (xt1 + xt3); vt2 += a * (xt2 + xt4);
        vt3 += a * (xt3 + xt1); vt4 += a * (xt4 + xt2);
        vt1 += -0.5 * (yt1 - yt3); vt2 += -0.5 * (yt2 - yt4);
        vt3 += -0.5 * (yt3 - yt1); vt4 += -0.5 * (yt4 - yt2);
        s1[is] = vt1; s2[is] = vt2; s3[is] = vt3; s4[is] = vt4;
      }
      // ---- backward sweep U^{-1}: s (local) -> v ----
      int is0 = NS - 1;
      real_t a0 = dpinv[NS - 1];
      vt1 = a0 * s1[is0]; vt2 = a0 * s2[is0]; vt3 = a0 * s3[is0]; vt4 = a0 * s4[is0];
      st_bf16f(vp, IDX2(Nin5, (ID1 + ivc + NVCD * is0), site), vt1);
      st_bf16f(vp, IDX2(Nin5, (ID2 + ivc + NVCD * is0), site), vt2);
      st_bf16f(vp, IDX2(Nin5, (ID3 + ivc + NVCD * is0), site), vt3);
      st_bf16f(vp, IDX2(Nin5, (ID4 + ivc + NVCD * is0), site), vt4);
      if (gp) {  // gamma5: spin0<->2, spin1<->3
        st_bf16f(gp, IDX2(Nin5, (ID1 + ivc + NVCD * is0), site), vt3);
        st_bf16f(gp, IDX2(Nin5, (ID2 + ivc + NVCD * is0), site), vt4);
        st_bf16f(gp, IDX2(Nin5, (ID3 + ivc + NVCD * is0), site), vt1);
        st_bf16f(gp, IDX2(Nin5, (ID4 + ivc + NVCD * is0), site), vt2);
      }
      yt1 = 0.5 * (vt1 + vt3); yt2 = 0.5 * (vt2 + vt4);
      yt3 = 0.5 * (vt3 + vt1); yt4 = 0.5 * (vt4 + vt2);
#pragma unroll
      for (int is = NS - 2; is >= 0; --is) {
        xt1 = vt1; xt2 = vt2; xt3 = vt3; xt4 = vt4;
        vt1 = s1[is]; vt2 = s2[is]; vt3 = s3[is]; vt4 = s4[is];
        real_t a = real_t(0.5) * dm[is];
        vt1 += a * (xt1 - xt3); vt2 += a * (xt2 - xt4);
        vt3 += a * (xt3 - xt1); vt4 += a * (xt4 - xt2);
        real_t fis = f[is];
        vt1 += -fis * yt1; vt2 += -fis * yt2; vt3 += -fis * yt3; vt4 += -fis * yt4;
        real_t aa = dpinv[is];
        vt1 *= aa; vt2 *= aa; vt3 *= aa; vt4 *= aa;
        st_bf16f(vp, IDX2(Nin5, (ID1 + ivc + NVCD * is), site), vt1);
        st_bf16f(vp, IDX2(Nin5, (ID2 + ivc + NVCD * is), site), vt2);
        st_bf16f(vp, IDX2(Nin5, (ID3 + ivc + NVCD * is), site), vt3);
        st_bf16f(vp, IDX2(Nin5, (ID4 + ivc + NVCD * is), site), vt4);
        if (gp) {  // gamma5: spin0<->2, spin1<->3 (all Ls)
          st_bf16f(gp, IDX2(Nin5, (ID1 + ivc + NVCD * is), site), vt3);
          st_bf16f(gp, IDX2(Nin5, (ID2 + ivc + NVCD * is), site), vt4);
          st_bf16f(gp, IDX2(Nin5, (ID3 + ivc + NVCD * is), site), vt1);
          st_bf16f(gp, IDX2(Nin5, (ID4 + ivc + NVCD * is), site), vt2);
        }
      }
    }
  }
}

// --- C^{-dag} = Ldag^{-1} Udag^{-1}, fused (no global intermediate) ---
// Templated on compile-time NS so the Ls buffer s[NS] lands in registers (see
// the Cinv note above); dispatched per Ns by the wrapper.
template<int NS>
__global__ void mult_dw5din_Cdaginv_mrhs_dev(
    real_t* const* __restrict__ vp_arr, real_t* const* __restrict__ wp_arr,
    int nrhs, int Nin5,
    const real_t* __restrict__ e, const real_t* __restrict__ f,
    const real_t* __restrict__ dpinv, const real_t* __restrict__ dm, int Nst_pad)
{
  const int ist = blockIdx.x * blockDim.x + threadIdx.x;
  const int GridSize = blockDim.x * gridDim.x;
  for (int idx = ist; idx < Nst_pad * NVC; idx += GridSize) {
    int idx2 = idx / NWP, idx_in = idx % NWP, ivc = idx2 % NVC, idx_out = idx2 / NVC;
    int site = idx_in + NWP * idx_out;
    for (int r = 0; r < nrhs; ++r) {
      real_t* __restrict__ vp = vp_arr[r];
      real_t* __restrict__ wp = wp_arr[r];
      real_t s1[NS], s2[NS];
      real_t s3[NS], s4[NS];
      real_t vt1, vt2, vt3, vt4, yt1, yt2, yt3, yt4, xt1, xt2, xt3, xt4;
      // ---- forward sweep Udag^{-1}: w -> s (local) ----
      real_t a0 = dpinv[0];
      vt1 = a0 * wp[IDX2(Nin5, (ID1 + ivc), site)];
      vt2 = a0 * wp[IDX2(Nin5, (ID2 + ivc), site)];
      vt3 = a0 * wp[IDX2(Nin5, (ID3 + ivc), site)];
      vt4 = a0 * wp[IDX2(Nin5, (ID4 + ivc), site)];
      s1[0] = vt1; s2[0] = vt2; s3[0] = vt3; s4[0] = vt4;
      real_t f0 = f[0];
      yt1 = f0 * vt1; yt2 = f0 * vt2; yt3 = f0 * vt3; yt4 = f0 * vt4;
#pragma unroll
      for (int is = 1; is < NS - 1; ++is) {
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
        s1[is] = vt1; s2[is] = vt2; s3[is] = vt3; s4[is] = vt4;
        real_t fis = f[is];
        yt1 += fis * vt1; yt2 += fis * vt2; yt3 += fis * vt3; yt4 += fis * vt4;
      }
      {
        int is = NS - 1;
        xt1 = vt1; xt2 = vt2; xt3 = vt3; xt4 = vt4;
        vt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)];
        vt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)];
        vt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)];
        vt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)];
        real_t a = real_t(0.5) * dm[NS - 2];
        vt1 += a * (xt1 - xt3); vt2 += a * (xt2 - xt4);
        vt3 += a * (xt3 - xt1); vt4 += a * (xt4 - xt2);
        vt1 += -0.5 * (yt1 + yt3); vt2 += -0.5 * (yt2 + yt4);
        vt3 += -0.5 * (yt3 + yt1); vt4 += -0.5 * (yt4 + yt2);
        real_t aa = dpinv[NS - 1];
        vt1 *= aa; vt2 *= aa; vt3 *= aa; vt4 *= aa;
        s1[is] = vt1; s2[is] = vt2; s3[is] = vt3; s4[is] = vt4;
      }
      // ---- backward sweep Ldag^{-1}: s (local) -> v ----
      int is0 = NS - 1;
      vt1 = s1[is0]; vt2 = s2[is0]; vt3 = s3[is0]; vt4 = s4[is0];
      vp[IDX2(Nin5, (ID1 + ivc + NVCD * is0), site)] = vt1;
      vp[IDX2(Nin5, (ID2 + ivc + NVCD * is0), site)] = vt2;
      vp[IDX2(Nin5, (ID3 + ivc + NVCD * is0), site)] = vt3;
      vp[IDX2(Nin5, (ID4 + ivc + NVCD * is0), site)] = vt4;
      yt1 = 0.5 * (vt1 - vt3); yt2 = 0.5 * (vt2 - vt4);
      yt3 = 0.5 * (vt3 - vt1); yt4 = 0.5 * (vt4 - vt2);
#pragma unroll
      for (int is = NS - 2; is >= 0; --is) {
        xt1 = vt1; xt2 = vt2; xt3 = vt3; xt4 = vt4;
        vt1 = s1[is]; vt2 = s2[is]; vt3 = s3[is]; vt4 = s4[is];
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

// --- C^{-dag} = Ldag^{-1} Udag^{-1}, fused, BF16-STORAGE variant (dag sibling
// of mult_dw5din_Cinv_mrhs_bf16_dev).  Arithmetic is byte-for-byte the FP32
// recurrence of mult_dw5din_Cdaginv_mrhs_dev above; only the operand vectors
// (w in, v out) live in global memory as __nv_bfloat16.  Accumulation stays in
// FP32 registers; the LU coefficients e/f/dpinv/dm stay FP32.  No gm5 output.
template<int NS>
__global__ void mult_dw5din_Cdaginv_mrhs_bf16_dev(
    __nv_bfloat16* const* __restrict__ vp_arr, __nv_bfloat16* const* __restrict__ wp_arr,
    int nrhs, int Nin5,
    const real_t* __restrict__ e, const real_t* __restrict__ f,
    const real_t* __restrict__ dpinv, const real_t* __restrict__ dm, int Nst_pad)
{
  const int ist = blockIdx.x * blockDim.x + threadIdx.x;
  const int GridSize = blockDim.x * gridDim.x;
  for (int idx = ist; idx < Nst_pad * NVC; idx += GridSize) {
    int idx2 = idx / NWP, idx_in = idx % NWP, ivc = idx2 % NVC, idx_out = idx2 / NVC;
    int site = idx_in + NWP * idx_out;
    for (int r = 0; r < nrhs; ++r) {
      __nv_bfloat16* __restrict__ vp = vp_arr[r];
      const __nv_bfloat16* __restrict__ wp = wp_arr[r];
      real_t s1[NS], s2[NS];
      real_t s3[NS], s4[NS];
      real_t vt1, vt2, vt3, vt4, yt1, yt2, yt3, yt4, xt1, xt2, xt3, xt4;
      // ---- forward sweep Udag^{-1}: w -> s (local) ----
      real_t a0 = dpinv[0];
      vt1 = a0 * ld_bf16f(wp, IDX2(Nin5, (ID1 + ivc), site));
      vt2 = a0 * ld_bf16f(wp, IDX2(Nin5, (ID2 + ivc), site));
      vt3 = a0 * ld_bf16f(wp, IDX2(Nin5, (ID3 + ivc), site));
      vt4 = a0 * ld_bf16f(wp, IDX2(Nin5, (ID4 + ivc), site));
      s1[0] = vt1; s2[0] = vt2; s3[0] = vt3; s4[0] = vt4;
      real_t f0 = f[0];
      yt1 = f0 * vt1; yt2 = f0 * vt2; yt3 = f0 * vt3; yt4 = f0 * vt4;
#pragma unroll
      for (int is = 1; is < NS - 1; ++is) {
        xt1 = vt1; xt2 = vt2; xt3 = vt3; xt4 = vt4;
        vt1 = ld_bf16f(wp, IDX2(Nin5, (ID1 + ivc + NVCD * is), site));
        vt2 = ld_bf16f(wp, IDX2(Nin5, (ID2 + ivc + NVCD * is), site));
        vt3 = ld_bf16f(wp, IDX2(Nin5, (ID3 + ivc + NVCD * is), site));
        vt4 = ld_bf16f(wp, IDX2(Nin5, (ID4 + ivc + NVCD * is), site));
        real_t a = real_t(0.5) * dm[is - 1];
        vt1 += a * (xt1 - xt3); vt2 += a * (xt2 - xt4);
        vt3 += a * (xt3 - xt1); vt4 += a * (xt4 - xt2);
        real_t aa = dpinv[is];
        vt1 *= aa; vt2 *= aa; vt3 *= aa; vt4 *= aa;
        s1[is] = vt1; s2[is] = vt2; s3[is] = vt3; s4[is] = vt4;
        real_t fis = f[is];
        yt1 += fis * vt1; yt2 += fis * vt2; yt3 += fis * vt3; yt4 += fis * vt4;
      }
      {
        int is = NS - 1;
        xt1 = vt1; xt2 = vt2; xt3 = vt3; xt4 = vt4;
        vt1 = ld_bf16f(wp, IDX2(Nin5, (ID1 + ivc + NVCD * is), site));
        vt2 = ld_bf16f(wp, IDX2(Nin5, (ID2 + ivc + NVCD * is), site));
        vt3 = ld_bf16f(wp, IDX2(Nin5, (ID3 + ivc + NVCD * is), site));
        vt4 = ld_bf16f(wp, IDX2(Nin5, (ID4 + ivc + NVCD * is), site));
        real_t a = real_t(0.5) * dm[NS - 2];
        vt1 += a * (xt1 - xt3); vt2 += a * (xt2 - xt4);
        vt3 += a * (xt3 - xt1); vt4 += a * (xt4 - xt2);
        vt1 += -0.5 * (yt1 + yt3); vt2 += -0.5 * (yt2 + yt4);
        vt3 += -0.5 * (yt3 + yt1); vt4 += -0.5 * (yt4 + yt2);
        real_t aa = dpinv[NS - 1];
        vt1 *= aa; vt2 *= aa; vt3 *= aa; vt4 *= aa;
        s1[is] = vt1; s2[is] = vt2; s3[is] = vt3; s4[is] = vt4;
      }
      // ---- backward sweep Ldag^{-1}: s (local) -> v ----
      int is0 = NS - 1;
      vt1 = s1[is0]; vt2 = s2[is0]; vt3 = s3[is0]; vt4 = s4[is0];
      st_bf16f(vp, IDX2(Nin5, (ID1 + ivc + NVCD * is0), site), vt1);
      st_bf16f(vp, IDX2(Nin5, (ID2 + ivc + NVCD * is0), site), vt2);
      st_bf16f(vp, IDX2(Nin5, (ID3 + ivc + NVCD * is0), site), vt3);
      st_bf16f(vp, IDX2(Nin5, (ID4 + ivc + NVCD * is0), site), vt4);
      yt1 = 0.5 * (vt1 - vt3); yt2 = 0.5 * (vt2 - vt4);
      yt3 = 0.5 * (vt3 - vt1); yt4 = 0.5 * (vt4 - vt2);
#pragma unroll
      for (int is = NS - 2; is >= 0; --is) {
        xt1 = vt1; xt2 = vt2; xt3 = vt3; xt4 = vt4;
        vt1 = s1[is]; vt2 = s2[is]; vt3 = s3[is]; vt4 = s4[is];
        real_t a = real_t(0.5) * dm[is + 1] * dpinv[is];
        vt1 += a * (xt1 + xt3); vt2 += a * (xt2 + xt4);
        vt3 += a * (xt3 + xt1); vt4 += a * (xt4 + xt2);
        real_t eis = e[is];
        vt1 += -eis * yt1; vt2 += -eis * yt2; vt3 += -eis * yt3; vt4 += -eis * yt4;
        st_bf16f(vp, IDX2(Nin5, (ID1 + ivc + NVCD * is), site), vt1);
        st_bf16f(vp, IDX2(Nin5, (ID2 + ivc + NVCD * is), site), vt2);
        st_bf16f(vp, IDX2(Nin5, (ID3 + ivc + NVCD * is), site), vt3);
        st_bf16f(vp, IDX2(Nin5, (ID4 + ivc + NVCD * is), site), vt4);
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

// Dispatch the fused LU over compile-time NS so the s[] Ls buffer lands in
// registers (runtime Ns -> local memory -> ~2.3x DRAM over-traffic). Common DWF
// Ls only; any other Ns falls through to the non-fused Linv/Uinv path. The macro
// names (vdev/wdev/.../gridSize) match the locals in both wrappers below.
#define CINV_CASE(N)    case N: mult_dw5din_Cinv_mrhs_dev<N><<<gridSize, blockSize>>>(vdev, wdev, gm5dev, nrhs, Nin5, e_dev, f_dev, dpinv_dev, dm_dev, Nst_pad); break
#define CDAGINV_CASE(N) case N: mult_dw5din_Cdaginv_mrhs_dev<N><<<gridSize, blockSize>>>(vdev, wdev, nrhs, Nin5, e_dev, f_dev, dpinv_dev, dm_dev, Nst_pad); break
#define FUSED_NS_LIST(CASE) CASE(8); CASE(12); CASE(16); CASE(24)

// Prec = C^{-1} = U^{-1} L^{-1}  (mass operator's LU; e/f/dpinv/dm host bases).
// gm5_host (optional): if non-null, ALSO write gamma5*(C^{-1}w) there -- the fused
// Cinv has the result in registers, so it is one extra permuted store, and it
// pre-supplies the gm5(w) that the following fineDdag needs (lets that call skip
// its separate gm5 full-field kernel). Only used for Ns in the FUSED_NS_LIST set.
void finePrec_mrhs(real_t* const* v_host, real_t* const* w_host, int nrhs, int Ns,
                   real_t* e_host, real_t* f_host, real_t* dpinv_host, real_t* dm_host,
                   int* Nsize, real_t* const* gm5_host)
{
  int Nst = Nsize[0]*Nsize[1]*Nsize[2]*Nsize[3];
  int Nst_pad = ceil_nwp(Nst);
  int Nin5 = NVCD * Ns;
  static real_t **vd = nullptr, **wd = nullptr; static int cv = 0, cw = 0;
  static real_t **gd = nullptr; static int cg = 0;
  real_t** vdev = map_upload(v_host, nrhs, vd, cv);
  real_t** wdev = map_upload(w_host, nrhs, wd, cw);
  real_t** gm5dev = gm5_host ? map_upload(gm5_host, nrhs, gd, cg) : nullptr;
  real_t** sc   = lu_scratch(nrhs, Nin5, Nst_pad);
  real_t* e_dev = (real_t*)dev_ptr(e_host);
  real_t* f_dev = (real_t*)dev_ptr(f_host);
  real_t* dpinv_dev = (real_t*)dev_ptr(dpinv_host);
  real_t* dm_dev    = (real_t*)dev_ptr(dm_host);
  int blockSize = VECTOR_LENGTH;
  int gridSize  = (Nst_pad * NVC + blockSize - 1) / blockSize;
  fp16_probe("finePrec.in", w_host[0], Nin5, Nst_pad);
  fine_time_begin();
  switch (Ns) {
    // fused C^{-1} with NS compile-time -> s[] in registers (gm5_arr=nullptr:
    // finePrec doesn't need the spin-swapped gm5 copy-out)
    FUSED_NS_LIST(CINV_CASE);
    default:
      // non-fused fallback (any Ns): L^{-1} w -> sc ; U^{-1} sc -> v.
      // The gm5 copy-out only exists on the fused path, so a gm5 request here
      // would leave the buffer stale -> fail loud rather than return wrong Ddag.
      if (gm5dev) { printf("finePrec_mrhs: gm5 fold unsupported for Ns=%d "
                           "(add it to FUSED_NS_LIST)\n", Ns); exit(1); }
      mult_dw5din_Linv_mrhs_dev<<<gridSize, blockSize>>>(sc, wdev, nrhs, Ns, Nin5, e_dev, dpinv_dev, dm_dev, Nst_pad);
      mult_dw5din_Uinv_mrhs_dev<<<gridSize, blockSize>>>(vdev, sc, nrhs, Ns, Nin5, f_dev, dpinv_dev, dm_dev, Nst_pad);
  }
  afield_dd_kernel_sync();
  fine_time_end(2);
  fp16_probe("finePrec.out", v_host[0], Nin5, Nst_pad);
}

// ---- BF16-storage scratch + convert kernels (Step-1 A/B harness) -----------
// Mirrors lu_scratch: one cached device bf16 buffer + per-rhs pointer array per
// role (0=w, 1=v, 2=gm5).
static __nv_bfloat16*  g_bf_buf[3]  = { nullptr, nullptr, nullptr };
static long            g_bf_cap[3]  = { 0, 0, 0 };
static __nv_bfloat16** g_bf_ptrs[3] = { nullptr, nullptr, nullptr };
static int             g_bf_pcap[3] = { 0, 0, 0 };
static __nv_bfloat16** bf16_scratch(int slot, int nrhs, long per)
{
  long need = (long)nrhs * per;
  if (need > g_bf_cap[slot]) {
    if (g_bf_buf[slot]) cudaFree(g_bf_buf[slot]);
    CHECK(cudaMalloc((void**)&g_bf_buf[slot], sizeof(__nv_bfloat16) * need));
    g_bf_cap[slot] = need;
  }
  __nv_bfloat16* hp[256];
  if (nrhs > 256) { printf("bf16_scratch: nrhs > 256\n"); exit(1); }
  for (int r = 0; r < nrhs; ++r) hp[r] = g_bf_buf[slot] + (long)r * per;
  if (nrhs > g_bf_pcap[slot]) {
    if (g_bf_ptrs[slot]) cudaFree(g_bf_ptrs[slot]);
    CHECK(cudaMalloc((void**)&g_bf_ptrs[slot], sizeof(__nv_bfloat16*) * nrhs));
    g_bf_pcap[slot] = nrhs;
  }
  CHECK(cudaMemcpy(g_bf_ptrs[slot], hp, sizeof(__nv_bfloat16*) * nrhs, cudaMemcpyHostToDevice));
  return g_bf_ptrs[slot];
}
__global__ void k_field_f2bf16(__nv_bfloat16* __restrict__ o, const real_t* __restrict__ in, long n)
{ long i = (long)blockIdx.x * blockDim.x + threadIdx.x; if (i < n) o[i] = __float2bfloat16(in[i]); }
__global__ void k_field_bf162f(real_t* __restrict__ o, const __nv_bfloat16* __restrict__ in, long n)
{ long i = (long)blockIdx.x * blockDim.x + threadIdx.x; if (i < n) o[i] = __bfloat162float(in[i]); }

// BF16-STORAGE finePrec (Step-1 harness): float in -> bf16 storage -> fused
// bf16 Cinv -> float out.  Matches finePrec_mrhs up to the bf16 operand rounding
// (~4e-3 rel), which the 2pt A/B showed sits inside the smoother budget.  The
// convert-in/out is for CORRECTNESS A/B only; the bandwidth win needs bf16-
// RESIDENT vectors across the V-cycle (Step-2 integration).
void finePrec_mrhs_bf16(real_t* const* v_host, real_t* const* w_host, int nrhs, int Ns,
                        real_t* e_host, real_t* f_host, real_t* dpinv_host, real_t* dm_host,
                        int* Nsize, real_t* const* gm5_host)
{
  int Nst = Nsize[0]*Nsize[1]*Nsize[2]*Nsize[3];
  int Nst_pad = ceil_nwp(Nst);
  int Nin5 = NVCD * Ns;
  long per = (long)Nin5 * Nst_pad;
  __nv_bfloat16** wbf = bf16_scratch(0, nrhs, per);
  __nv_bfloat16** vbf = bf16_scratch(1, nrhs, per);
  __nv_bfloat16** gbf = gm5_host ? bf16_scratch(2, nrhs, per) : nullptr;
  real_t* e_dev = (real_t*)dev_ptr(e_host);
  real_t* f_dev = (real_t*)dev_ptr(f_host);
  real_t* dpinv_dev = (real_t*)dev_ptr(dpinv_host);
  real_t* dm_dev    = (real_t*)dev_ptr(dm_host);
  const int bs = 256;
  int cg = (int)((per + bs - 1) / bs);
  for (int r = 0; r < nrhs; ++r)   // w: float -> bf16
    k_field_f2bf16<<<cg, bs>>>(g_bf_buf[0] + (long)r * per, (const real_t*)dev_ptr(w_host[r]), per);
  int blockSize = VECTOR_LENGTH;
  int gridSize  = (Nst_pad * NVC + blockSize - 1) / blockSize;
  fine_time_begin();
  switch (Ns) {
#define CINV_BF16_CASE(N) case N: mult_dw5din_Cinv_mrhs_bf16_dev<N><<<gridSize, blockSize>>>(vbf, wbf, gbf, nrhs, Nin5, e_dev, f_dev, dpinv_dev, dm_dev, Nst_pad); break
    FUSED_NS_LIST(CINV_BF16_CASE);
#undef CINV_BF16_CASE
    default:
      printf("finePrec_mrhs_bf16: Ns=%d not in FUSED_NS_LIST\n", Ns); exit(1);
  }
  fine_time_end(2);
  for (int r = 0; r < nrhs; ++r)   // v: bf16 -> float
    k_field_bf162f<<<cg, bs>>>((real_t*)dev_ptr(v_host[r]), g_bf_buf[1] + (long)r * per, per);
  if (gm5_host) for (int r = 0; r < nrhs; ++r)   // gm5: bf16 -> float
    k_field_bf162f<<<cg, bs>>>((real_t*)dev_ptr(gm5_host[r]), g_bf_buf[2] + (long)r * per, per);
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
  fp16_probe("finePrecdag.in", w_host[0], Nin5, Nst_pad);
  fine_time_begin();
  switch (Ns) {
    // fused C^{-dag} with NS compile-time -> s[] in registers
    FUSED_NS_LIST(CDAGINV_CASE);
    default:
      // non-fused fallback (any Ns): Udag^{-1} w -> sc ; Ldag^{-1} sc -> v
      mult_dw5din_Udaginv_mrhs_dev<<<gridSize, blockSize>>>(sc, wdev, nrhs, Ns, Nin5, f_dev, dpinv_dev, dm_dev, Nst_pad);
      mult_dw5din_Ldaginv_mrhs_dev<<<gridSize, blockSize>>>(vdev, sc, nrhs, Ns, Nin5, e_dev, dpinv_dev, dm_dev, Nst_pad);
  }
  afield_dd_kernel_sync();
  fine_time_end(3);
  fp16_probe("finePrecdag.out", v_host[0], Nin5, Nst_pad);
}

// BF16-STORAGE Precdag (dag sibling of finePrec_mrhs_bf16): float in -> bf16
// storage -> fused bf16 C^{-dag} -> float out.  Reuses bf16 scratch slots 0/1
// (the forward finePrec_bf16 has fully synced + converted out by the time the
// fineDdag/finePrecdag stage of apply_A_block runs).  No gm5.  FUSED_NS_LIST
// only (matches the FP32 fused dag path).
void finePrecdag_mrhs_bf16(real_t* const* v_host, real_t* const* w_host, int nrhs, int Ns,
                           real_t* e_host, real_t* f_host, real_t* dpinv_host, real_t* dm_host,
                           int* Nsize)
{
  int Nst = Nsize[0]*Nsize[1]*Nsize[2]*Nsize[3];
  int Nst_pad = ceil_nwp(Nst);
  int Nin5 = NVCD * Ns;
  long per = (long)Nin5 * Nst_pad;
  __nv_bfloat16** wbf = bf16_scratch(0, nrhs, per);
  __nv_bfloat16** vbf = bf16_scratch(1, nrhs, per);
  real_t* e_dev = (real_t*)dev_ptr(e_host);
  real_t* f_dev = (real_t*)dev_ptr(f_host);
  real_t* dpinv_dev = (real_t*)dev_ptr(dpinv_host);
  real_t* dm_dev    = (real_t*)dev_ptr(dm_host);
  const int bs = 256;
  int cg = (int)((per + bs - 1) / bs);
  for (int r = 0; r < nrhs; ++r)   // w: float -> bf16
    k_field_f2bf16<<<cg, bs>>>(g_bf_buf[0] + (long)r * per, (const real_t*)dev_ptr(w_host[r]), per);
  int blockSize = VECTOR_LENGTH;
  int gridSize  = (Nst_pad * NVC + blockSize - 1) / blockSize;
  fine_time_begin();
  switch (Ns) {
#define CDAGINV_BF16_CASE(N) case N: mult_dw5din_Cdaginv_mrhs_bf16_dev<N><<<gridSize, blockSize>>>(vbf, wbf, nrhs, Nin5, e_dev, f_dev, dpinv_dev, dm_dev, Nst_pad); break
    FUSED_NS_LIST(CDAGINV_BF16_CASE);
#undef CDAGINV_BF16_CASE
    default:
      printf("finePrecdag_mrhs_bf16: Ns=%d not in FUSED_NS_LIST\n", Ns); exit(1);
  }
  fine_time_end(3);
  for (int r = 0; r < nrhs; ++r)   // v: bf16 -> float
    k_field_bf162f<<<cg, bs>>>((real_t*)dev_ptr(v_host[r]), g_bf_buf[1] + (long)r * per, per);
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
// gm5_host (optional): precomputed gamma5*w (e.g. emitted by the preceding
// finePrec_mrhs). If supplied, use it directly as t1 and SKIP the separate gm5
// full-field kernel. Must equal gm5(w) for the given w, else results are wrong.
void fineDdag_mrhs(real_t* const* v_host, real_t* const* w_host, real_t* u_field_host,
                   int nrhs, real_t mq, real_t M0, int Ns, real_t alpha,
                   int* Nsize, int* bc, int* do_comm, real_t* const* gm5_host)
{
  int Nst = Nsize[0]*Nsize[1]*Nsize[2]*Nsize[3];
  int Nst_pad = ceil_nwp(Nst);
  int Nin5 = NVCD * Ns;
  real_t* u_dev = (real_t*)dev_ptr(u_field_host);
  static real_t **vd = nullptr, **wd = nullptr; static int cv = 0, cw = 0;
  static real_t **gd = nullptr; static int cg = 0;
  real_t** vdev = map_upload(v_host, nrhs, vd, cv);
  real_t** wdev = map_upload(w_host, nrhs, wd, cw);
  real_t** t2 = lu_scratch2(nrhs, Nin5, Nst_pad);   // hopb(gm5(w))
  // t1 = gm5(w): use the caller-supplied buffer if given, else compute it below.
  real_t** t1 = gm5_host ? map_upload(gm5_host, nrhs, gd, cg)
                         : lu_scratch(nrhs, Nin5, Nst_pad);

  int blockSize = VECTOR_LENGTH;
  int gridV = (Nst_pad * NVC + blockSize - 1) / blockSize;
  int gridS = (Nst + blockSize - 1) / blockSize;
  int gridH = (Nst_pad * Ns + blockSize - 1) / blockSize;

  fp16_probe("fineDdag.in", w_host[0], Nin5, Nst_pad);
  fine_time_begin();

  // t1 = gm5 w  (skipped when the caller pre-supplied it in gm5_host)
  if (!gm5_host)
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
  fine_time_end(1);
  fp16_probe("fineDdag.out", v_host[0], Nin5, Nst_pad);
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
