/*!
      @file    block_fgmres_dw.h
      @brief   Flexible block-GMRES (FBGMRES) for the domainwall AMG, dev8-local.
      @author  Wei-Lun Chen (wlchen)

      Solves  A X = B  with X,B = s columns (s right-hand sides) at once, where
      A is the float PVprec operator (m_afopr_A_F) and the (variable) right
      preconditioner is the float two-grid V-cycle (m_prec_mg->mult_single),
      applied per column.  The block-Arnoldi orthogonalisation is done with the
      VERIFIED MRHS GEMM kernels (mrhs_live::block_inner / block_update) and a
      Cholesky-QR block normalisation (block_inner -> host Cholesky -> block_update),
      so the expensive field-level block work runs as dense GEMMs (TF32 tensor
      cores when USE_TENSORCORE_MRHS).

      The small dense block-Hessenberg least-squares and the s x s Cholesky /
      triangular-inverse run on the host in std::complex<double>, matching the
      existing single-RHS FGMRES which also keeps its Hessenberg coefficients in
      double-complex even in the float solver.  No field/operator arithmetic is
      ever promoted to FP64 — only the O(m^2 s^2) host scalar reductions.

      AMG-local to dev8; nothing shared is touched.
*/
#ifndef BLOCK_FGMRES_DW_H
#define BLOCK_FGMRES_DW_H

#include <vector>
#include <complex>
#include <functional>
#include <cmath>

#include "mrhs_block_live.h"

template<typename AFIELD>
class BlockFGMRES_dw {
 public:
  using real_t    = typename AFIELD::real_t;        // float
  using complex_t = typename AFIELD::complex_t;     // std::complex<float>
  using cd        = std::complex<double>;

  // per-column operators (bound by the owner to A_F->mult and prec->mult_single)
  using Op = std::function<void(AFIELD&, const AFIELD&)>;

  BlockFGMRES_dw(int Nin, int Nvol, int Nex, const int Nsize[4])
    : m_Nin(Nin), m_Nvol(Nvol), m_Nex(Nex)
  { for (int i = 0; i < 4; ++i) m_Nsize[i] = Nsize[i]; }

  ~BlockFGMRES_dw() {
    if (m_G_dev)    mrhs_live::dev_free(m_G_dev);
    if (m_Hess_dev) mrhs_live::dev_free(m_Hess_dev);
    if (m_negLi_dev)mrhs_live::dev_free(m_negLi_dev);
    if (m_Y_dev)    mrhs_live::dev_free(m_Y_dev);
    if (m_L0_dev)   mrhs_live::dev_free(m_L0_dev);
  }

  // block preconditioner: applies M^{-1} to all s columns at once.
  using BlockOp = std::function<void(std::vector<AFIELD>&,
                                     const std::vector<AFIELD>&)>;

  void set_ops(Op A, Op Minv) { m_A = A; m_Minv = Minv; }
  // optional: batched (MRHS) preconditioner; if set, used instead of per-column Minv.
  void set_block_prec(BlockOp Mblk) { m_Mblk = Mblk; m_has_block_prec = true; }
  // optional: batched (MRHS) operator A applied to all s columns at once; if set,
  // used instead of the per-column m_A (the throughput path for the smoother).
  void set_block_A(BlockOp Ablk) { m_Ablk = Ablk; m_has_block_A = true; }
  void set_parameters(int restart_len, int max_restart, double stop_cond_sq) {
    m_m = restart_len; m_max_restart = max_restart; m_stop2 = stop_cond_sq;
  }

  // Solve A X = B (s columns each).  Returns the worst per-column relative
  // residual ‖A X_r - B_r‖/‖B_r‖ and fills relres[r] for each column.
  double solve(std::vector<AFIELD>& X, const std::vector<AFIELD>& B,
               std::vector<double>& relres, int& total_vcycles);

 private:
  int m_Nin, m_Nvol, m_Nex, m_Nsize[4];
  int m_m = 20, m_max_restart = 5;
  double m_stop2 = 1.0e-16;
  Op m_A, m_Minv;
  BlockOp m_Mblk, m_Ablk;
  bool m_has_block_prec = false;
  bool m_has_block_A = false;

  // persistent work buffers, reused across solve() calls (reset only when the
  // shape (s, m, Nin) changes) -- cuts the per-call AField/device allocations
  // that dominate the smoother hot loop.
  std::vector<std::vector<AFIELD> > m_V, m_Z;   // Krylov basis / preconditioned
  std::vector<AFIELD> m_R, m_AX, m_Wstep, m_QRtmp;
  AFIELD m_wtmp;
  // device-resident dense-LA buffers: the s x s Gram, the (m+1) x m block-
  // Hessenberg, the Cholesky inverse scratch, and the least-squares solution Y
  // all live on the GPU so NOTHING round-trips to the host per Krylov step.
  // Host only touches them ONCE per restart (download Hbar+L0, lstsq, upload Y).
  float *m_G_dev = nullptr; long m_G_cap = 0;     // s x s Gram scratch (blockQR)
  float *m_Hess_dev = nullptr; long m_Hess_cap = 0;// (m+1)*m blocks of 2*s*s
  float *m_negLi_dev = nullptr; long m_negLi_cap = 0; // 2*s*s Cholesky-inv scratch
  float *m_Y_dev = nullptr; long m_Y_cap = 0;     // m blocks of 2*s*s (-Y for update)
  float *m_L0_dev = nullptr; long m_L0_cap = 0;   // 2*s*s (R-factor of V[0])
  // persistent HOST scratch for the once-per-restart least-squares (resize() keeps
  // capacity, so no per-restart heap reallocation).
  std::vector<float> m_Hbuf, m_L0buf, m_negYbuf;
  std::vector<cd>    m_Hbar, m_Xi, m_Yv;
  int m_alloc_s = -1, m_alloc_m = -1, m_alloc_nin = -1;
  // device offset of Hessenberg block H[k][j] (k=0..m, j=0..m-1), each 2*s*s floats.
  long hess_off(int k, int j, int s) const { return (long)(k*m_m + j) * 2L * s * s; }

  // ---- small dense host helpers (std::complex<double>) ----
  // s x s Cholesky:  G = L^H L, L upper-triangular.  G,L row-major [i*s+j].
  static void choleskyU(const std::vector<cd>& G, std::vector<cd>& L, int s);
  // upper-triangular inverse (s x s).
  static void invU(const std::vector<cd>& L, std::vector<cd>& Li, int s);
  // least squares min_Y ‖Xi - H Y‖_F, H is M x N (M>=N), Xi is M x s, Y is N x s.
  static void lstsq(int M, int N, int s,
                    std::vector<cd> H, std::vector<cd> Xi, std::vector<cd>& Y);

  // ptr helper
  static void ptrs(std::vector<AFIELD>& blk, std::vector<real_t*>& p) {
    p.resize(blk.size());
    for (size_t i = 0; i < blk.size(); ++i) p[i] = blk[i].ptr(0);
  }

  // Cholesky-QR block normalise (FULLY ON DEVICE): R(s fields) <- orthonormal
  // block; the s x s upper R-factor L (G = L^H L) is written to the device buffer
  // L_dev_slot (2*s*s floats).  No host transfer.
  void blockQR(std::vector<AFIELD>& R, int s, float* L_dev_slot);
};

//====================================================================
template<typename AFIELD>
void BlockFGMRES_dw<AFIELD>::choleskyU(const std::vector<cd>& G,
                                       std::vector<cd>& L, int s)
{
  L.assign(s * s, cd(0.0, 0.0));
  for (int j = 0; j < s; ++j) {
    double d = G[j * s + j].real();
    for (int k = 0; k < j; ++k) d -= std::norm(L[k * s + j]);
    if (d < 1.0e-30) d = 1.0e-30;
    double Ljj = std::sqrt(d);
    L[j * s + j] = cd(Ljj, 0.0);
    for (int i = j + 1; i < s; ++i) {
      cd sum = G[j * s + i];
      for (int k = 0; k < j; ++k) sum -= std::conj(L[k * s + j]) * L[k * s + i];
      L[j * s + i] = sum / Ljj;
    }
  }
}

template<typename AFIELD>
void BlockFGMRES_dw<AFIELD>::invU(const std::vector<cd>& L,
                                  std::vector<cd>& Li, int s)
{
  Li.assign(s * s, cd(0.0, 0.0));
  for (int i = 0; i < s; ++i) Li[i * s + i] = cd(1.0, 0.0) / L[i * s + i];
  for (int j = 0; j < s; ++j) {
    for (int i = j - 1; i >= 0; --i) {
      cd sum(0.0, 0.0);
      for (int k = i + 1; k <= j; ++k) sum += L[i * s + k] * Li[k * s + j];
      Li[i * s + j] = -sum / L[i * s + i];
    }
  }
}

template<typename AFIELD>
void BlockFGMRES_dw<AFIELD>::lstsq(int M, int N, int s,
                                   std::vector<cd> H, std::vector<cd> Xi,
                                   std::vector<cd>& Y)
{
  // Householder QR of H (M x N), apply Q^H to Xi (M x s), back-substitute.
  for (int k = 0; k < N; ++k) {
    double xn = 0.0;
    for (int r = k; r < M; ++r) xn += std::norm(H[r * N + k]);
    xn = std::sqrt(xn);
    if (xn < 1.0e-300) continue;
    cd x0 = H[k * N + k];
    double a0 = std::abs(x0);
    cd phase = (a0 == 0.0) ? cd(1.0, 0.0) : x0 / a0;
    cd alpha = -phase * xn;
    std::vector<cd> v(M - k);
    v[0] = x0 - alpha;
    for (int r = k + 1; r < M; ++r) v[r - k] = H[r * N + k];
    double vn2 = 0.0;
    for (int r = 0; r < M - k; ++r) vn2 += std::norm(v[r]);
    if (vn2 < 1.0e-300) continue;
    for (int c = k; c < N; ++c) {
      cd dot(0.0, 0.0);
      for (int r = k; r < M; ++r) dot += std::conj(v[r - k]) * H[r * N + c];
      dot *= 2.0 / vn2;
      for (int r = k; r < M; ++r) H[r * N + c] -= v[r - k] * dot;
    }
    for (int c = 0; c < s; ++c) {
      cd dot(0.0, 0.0);
      for (int r = k; r < M; ++r) dot += std::conj(v[r - k]) * Xi[r * s + c];
      dot *= 2.0 / vn2;
      for (int r = k; r < M; ++r) Xi[r * s + c] -= v[r - k] * dot;
    }
    H[k * N + k] = alpha;
    for (int r = k + 1; r < M; ++r) H[r * N + k] = cd(0.0, 0.0);
  }
  Y.assign(N * s, cd(0.0, 0.0));
  for (int c = 0; c < s; ++c) {
    for (int i = N - 1; i >= 0; --i) {
      cd t = Xi[i * s + c];
      for (int j = i + 1; j < N; ++j) t -= H[i * N + j] * Y[j * s + c];
      Y[i * s + c] = t / H[i * N + i];
    }
  }
}

//====================================================================
template<typename AFIELD>
void BlockFGMRES_dw<AFIELD>::blockQR(std::vector<AFIELD>& R, int s,
                                     float* L_dev_slot)
{
  real_t* rp[256]; real_t* vp[256];                  // stack, no per-call heap
  for (int c = 0; c < s; ++c) rp[c] = R[c].ptr(0);

  // G = R^H R  (s x s) -> device Gram scratch.  NO host transfer.
  mrhs_live::block_inner(m_G_dev, rp, rp, s, m_Nin, m_Nex, m_Nsize);
  // device Cholesky+inverse:  L (R-factor) -> L_dev_slot,  -L^{-1} -> m_negLi_dev.
  mrhs_live::block_chol_inv(L_dev_slot, m_negLi_dev, m_G_dev, s);

  // Vout = R * L^{-1}  via block_update(Vout=0, R, -Li): Vout_j -= sum_i R_i (-Li)[i][j]
  for (int j = 0; j < s; ++j) { m_QRtmp[j].set(real_t(0.0)); vp[j] = m_QRtmp[j].ptr(0); }
  mrhs_live::block_update(vp, rp, m_negLi_dev, s, m_Nin, m_Nex, m_Nsize);

  mrhs_live::block_copy(rp, vp, s, (long)m_Nin*R[0].nvol_pad()*m_Nex);  // R <- orthonormal (batched)
}

//====================================================================
template<typename AFIELD>
double BlockFGMRES_dw<AFIELD>::solve(std::vector<AFIELD>& X,
                                     const std::vector<AFIELD>& B,
                                     std::vector<double>& relres,
                                     int& total_vcycles)
{
  const int s = (int)B.size();
  const int m = m_m;
  // floats per field = PADDED physical buffer (Nin * nvol_pad * Nex).  Padding is
  // zero-initialised and preserved by copy/axpy/scal (and the IDX2 operators only
  // touch logical sites), so flat batched BLAS-1 over the padded buffer is bit-
  // exact -- including norm2, whose padding terms are zero.
  const long nflt = (long)m_Nin * const_cast<AFIELD&>(B[0]).nvol_pad() * m_Nex;
  total_vcycles = 0;

  real_t* Bp[256]; real_t* Rp[256]; real_t* AXp[256];   // stack host-base ptr arrays
  for (int r = 0; r < s; ++r) Bp[r] = const_cast<AFIELD&>(B[r]).ptr(0);

  std::vector<double> bnorm(s);
  mrhs_live::block_norm2(bnorm.data(), Bp, s, nflt);    // s norms, one launch
  for (int r = 0; r < s; ++r) bnorm[r] = std::sqrt(bnorm[r]);

  for (int r = 0; r < s; ++r) { X[r].reset(m_Nin, m_Nvol, m_Nex); X[r].set(real_t(0.0)); }

  // Krylov / work storage: persistent across calls, (re)allocated only on shape
  // change.  V[0..m] blocks of s, Z[0..m-1] blocks of s, R/AX/Wstep blocks of s.
  if (s != m_alloc_s || m != m_alloc_m || m_Nin != m_alloc_nin) {
    // set(0) zeros the FULL padded buffer (m_nsize_pad), incl. the NWP padding
    // lanes.  The batched BLAS-1 sums over the padded buffer, so the padding MUST
    // be zero for norm2 to be correct on padded (coarse) fields -- operators write
    // only logical sites, and copy/axpy/scal preserve the zero padding thereafter.
    m_V.assign(m + 1, std::vector<AFIELD>());
    m_Z.assign(m,     std::vector<AFIELD>());
    for (int k = 0; k <= m; ++k) { m_V[k].resize(s); for (int c = 0; c < s; ++c) { m_V[k][c].reset(m_Nin, m_Nvol, m_Nex); m_V[k][c].set(real_t(0.0)); } }
    for (int k = 0; k <  m; ++k) { m_Z[k].resize(s); for (int c = 0; c < s; ++c) { m_Z[k][c].reset(m_Nin, m_Nvol, m_Nex); m_Z[k][c].set(real_t(0.0)); } }
    m_R.resize(s); m_AX.resize(s); m_Wstep.resize(s); m_QRtmp.resize(s);
    for (int c = 0; c < s; ++c) { m_R[c].reset(m_Nin, m_Nvol, m_Nex);     m_R[c].set(real_t(0.0));
                                  m_AX[c].reset(m_Nin, m_Nvol, m_Nex);    m_AX[c].set(real_t(0.0));
                                  m_Wstep[c].reset(m_Nin, m_Nvol, m_Nex); m_Wstep[c].set(real_t(0.0));
                                  m_QRtmp[c].reset(m_Nin, m_Nvol, m_Nex); m_QRtmp[c].set(real_t(0.0)); }
    m_wtmp.reset(m_Nin, m_Nvol, m_Nex); m_wtmp.set(real_t(0.0));
    // device dense-LA buffers (grow-only) -- everything stays on the GPU.
    long gss = 2L * s * s;
    if (gss > m_G_cap)      { if (m_G_dev)     mrhs_live::dev_free(m_G_dev);     m_G_dev     = mrhs_live::dev_alloc(gss); m_G_cap = gss; }
    if (gss > m_negLi_cap)  { if (m_negLi_dev) mrhs_live::dev_free(m_negLi_dev); m_negLi_dev = mrhs_live::dev_alloc(gss); m_negLi_cap = gss; }
    if (gss > m_L0_cap)     { if (m_L0_dev)    mrhs_live::dev_free(m_L0_dev);    m_L0_dev    = mrhs_live::dev_alloc(gss); m_L0_cap = gss; }
    long hcap = (long)(m + 1) * m * gss;
    if (hcap > m_Hess_cap)  { if (m_Hess_dev)  mrhs_live::dev_free(m_Hess_dev);  m_Hess_dev  = mrhs_live::dev_alloc(hcap); m_Hess_cap = hcap; }
    long ycap = (long)m * gss;
    if (ycap > m_Y_cap)     { if (m_Y_dev)     mrhs_live::dev_free(m_Y_dev);     m_Y_dev     = mrhs_live::dev_alloc(ycap); m_Y_cap = ycap; }
    m_alloc_s = s; m_alloc_m = m; m_alloc_nin = m_Nin;
  }
  std::vector<std::vector<AFIELD> >& V = m_V;
  std::vector<std::vector<AFIELD> >& Z = m_Z;
  std::vector<AFIELD>& R = m_R;
  AFIELD& wtmp = m_wtmp;

  double worst = 1.0;

  for (int r = 0; r < s; ++r) Rp[r] = R[r].ptr(0);
  std::vector<double> rn(s);

  for (int cyc = 0; cyc < m_max_restart; ++cyc) {
    // R0 = B - A X  (s columns; batched A + batched BLAS-1)
    if (m_has_block_A) {
      std::vector<AFIELD>& AX = m_AX;
      m_Ablk(AX, X);
      for (int r = 0; r < s; ++r) AXp[r] = AX[r].ptr(0);
      mrhs_live::block_copy(Rp, Bp, s, nflt);            // R = B
      mrhs_live::block_axpy(Rp, AXp, real_t(-1.0), s, nflt);  // R -= AX
    } else {
      for (int r = 0; r < s; ++r) {
        m_A(wtmp, X[r]);          // wtmp = A X_r
        copy(R[r], B[r]);
        axpy(R[r], real_t(-1.0), wtmp);
      }
    }
    // converged?  (s norms in one launch)
    mrhs_live::block_norm2(rn.data(), Rp, s, nflt);
    worst = 0.0;
    for (int r = 0; r < s; ++r) {
      double rr = std::sqrt(rn[r]) / (bnorm[r] > 0 ? bnorm[r] : 1.0);
      if (rr > worst) worst = rr;
    }
    if (worst * worst < m_stop2) break;

    // V[0] L0 = R0  -- L0 written to device buffer m_L0_dev (no host transfer)
    real_t* V0p[256];
    for (int c = 0; c < s; ++c) V0p[c] = V[0][c].ptr(0);
    mrhs_live::block_copy(V0p, Rp, s, nflt);             // V[0] = R (batched)
    blockQR(V[0], s, m_L0_dev);

    // block-Arnoldi: every Hessenberg block H[k][j] is computed and consumed
    // ENTIRELY ON DEVICE (block_inner writes the s x s block straight into its
    // slot in m_Hess_dev; block_update reads it back).  NO per-step host roundtrip.
    int jdone = 0;
    for (int j = 0; j < m; ++j) {
      // Z[j] = M^{-1} V[j]   (batched if a block prec is set, else per column)
      if (m_has_block_prec) m_Mblk(Z[j], V[j]);
      else for (int c = 0; c < s; ++c) m_Minv(Z[j][c], V[j][c]);
      // W = A Z[j]   (batched A if available, else per column)
      std::vector<AFIELD>& W = m_Wstep;
      real_t* wp[256]; real_t* vkp[256];          // stack, no per-step heap
      if (m_has_block_A) m_Ablk(W, Z[j]);
      else for (int c = 0; c < s; ++c) m_A(W[c], Z[j][c]);
      for (int c = 0; c < s; ++c) wp[c] = W[c].ptr(0);
      total_vcycles += s;

      // block modified Gram-Schmidt against V[0..j], all device-resident
      for (int k = 0; k <= j; ++k) {
        for (int c = 0; c < s; ++c) vkp[c] = V[k][c].ptr(0);
        float* Hslot = m_Hess_dev + hess_off(k, j, s);
        // H[k][j] = V[k]^H W  -> straight into the device Hessenberg slot
        mrhs_live::block_inner(Hslot, vkp, wp, s, m_Nin, m_Nex, m_Nsize);
        // W -= V[k] H[k][j]   (reads the same device slot)
        mrhs_live::block_update(wp, vkp, Hslot, s, m_Nin, m_Nex, m_Nsize);
      }
      // V[j+1] H[j+1][j] = W  -- subdiagonal R-factor written to its device slot
      real_t* Vj1p[256];
      for (int c = 0; c < s; ++c) Vj1p[c] = V[j + 1][c].ptr(0);
      mrhs_live::block_copy(Vj1p, wp, s, nflt);          // V[j+1] = W (batched)
      blockQR(V[j + 1], s, m_Hess_dev + hess_off(j + 1, j, s));
      jdone = j + 1;
    }

    // ---- ONE host transfer per RESTART (not per step): download the assembled
    // Hessenberg + L0, do the (M x N) least-squares on host, upload -Y. ----
    int Mr = (jdone + 1) * s, Nc = jdone * s;
    m_Hbuf.resize((long)(m + 1) * m * 2 * s * s);
    m_L0buf.resize(2 * s * s);
    mrhs_live::dev_to_host(m_Hbuf.data(),  m_Hess_dev, (long)(m + 1) * m * 2 * s * s);
    mrhs_live::dev_to_host(m_L0buf.data(), m_L0_dev,   2L * s * s);

    m_Hbar.assign(Mr * Nc, cd(0.0, 0.0));
    m_Xi.assign(Mr * s, cd(0.0, 0.0));
    std::vector<cd>& Hbar = m_Hbar; std::vector<cd>& Xi = m_Xi; std::vector<cd>& Y = m_Yv;
    for (int j = 0; j < jdone; ++j)
      for (int k = 0; k <= j + 1; ++k) {
        long off = hess_off(k, j, s);
        for (int a = 0; a < s; ++a)
          for (int b = 0; b < s; ++b)
            Hbar[(k * s + a) * Nc + (j * s + b)] =
                cd(m_Hbuf[off + 2 * (a * s + b)], m_Hbuf[off + 2 * (a * s + b) + 1]);
      }
    for (int a = 0; a < s; ++a)
      for (int b = 0; b < s; ++b)
        Xi[a * s + b] = cd(m_L0buf[2 * (a * s + b)], m_L0buf[2 * (a * s + b) + 1]);

    lstsq(Mr, Nc, s, Hbar, Xi, Y);

    // upload -Y once, then X += sum_j Z[j] Y[j] via device block_update:
    // X_r -= sum_b Z[j][b] (-Y[j])[b][r]  ==  X_r += sum_b Z[j][b] Y[j][b][r].
    m_negYbuf.resize((long)jdone * 2 * s * s);
    for (int j = 0; j < jdone; ++j)
      for (int b = 0; b < s; ++b)
        for (int r = 0; r < s; ++r) {
          cd y = Y[(j * s + b) * s + r];
          m_negYbuf[(long)j * 2 * s * s + 2 * (b * s + r)]     = -(float)y.real();
          m_negYbuf[(long)j * 2 * s * s + 2 * (b * s + r) + 1] = -(float)y.imag();
        }
    mrhs_live::host_to_dev(m_Y_dev, m_negYbuf.data(), (long)jdone * 2 * s * s);
    real_t* xp[256]; real_t* zjp[256];                 // stack, no per-restart heap
    for (int c = 0; c < s; ++c) xp[c] = X[c].ptr(0);
    for (int j = 0; j < jdone; ++j) {
      for (int c = 0; c < s; ++c) zjp[c] = Z[j][c].ptr(0);
      mrhs_live::block_update(xp, zjp, m_Y_dev + (long)j * 2 * s * s,
                              s, m_Nin, m_Nex, m_Nsize);
    }
  }

  // final true residual R = B - A X (batched A + batched BLAS-1); reuse m_AX/m_R.
  relres.assign(s, 0.0);
  worst = 0.0;
  if (m_has_block_A) {
    m_Ablk(m_AX, X);
    for (int r = 0; r < s; ++r) AXp[r] = m_AX[r].ptr(0);
    mrhs_live::block_copy(Rp, Bp, s, nflt);              // R = B
    mrhs_live::block_axpy(Rp, AXp, real_t(-1.0), s, nflt);    // R -= AX
    mrhs_live::block_norm2(rn.data(), Rp, s, nflt);
    for (int r = 0; r < s; ++r) {
      double v = std::sqrt(rn[r]) / (bnorm[r] > 0 ? bnorm[r] : 1.0);
      relres[r] = v; if (v > worst) worst = v;
    }
  } else {
    for (int r = 0; r < s; ++r) {
      m_A(wtmp, X[r]);
      AFIELD& rr = m_R[r];
      copy(rr, B[r]); axpy(rr, real_t(-1.0), wtmp);
      double v = std::sqrt(rr.norm2()) / (bnorm[r] > 0 ? bnorm[r] : 1.0);
      relres[r] = v; if (v > worst) worst = v;
    }
  }
  return worst;
}

#endif // BLOCK_FGMRES_DW_H
