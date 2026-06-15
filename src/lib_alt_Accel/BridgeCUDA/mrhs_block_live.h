/*!
      @file    mrhs_block_live.h
      @brief   C++-callable entry points for the MRHS / tensor-core AMG kernels.
      @author  Wei-Lun Chen (wlchen)

      Thin forwarders to the device host-drivers in
      src/lib_alt_Accel/BridgeCUDA/src/mrhs_block_tensorcore_cuda-inc.h, exposed
      so the dev8 block-FGMRES driver (plain C++/.cpp) can drive the verified
      MRHS restriction / prolongation / coarse-Dirac / block-Krylov GEMMs that
      live in libbridgecuda's BridgeACC namespace.

      Pointers are Bridge AField *host base pointers* (the same value the
      single-RHS coarse op feeds to dev_ptr); the device-side mapping is done
      inside each driver.  All arithmetic is FP32 (TF32 on the tensor-core path,
      -DUSE_TENSORCORE_MRHS); never FP64.  AMG-local to dev8: nothing shared is
      touched.
*/
#ifndef MRHS_BLOCK_LIVE_H
#define MRHS_BLOCK_LIVE_H

namespace mrhs_live {

// C[nbasis x nrhs] (per aggregate) = B^H R.  coarse_dev layout:
// 2*(r + nrhs*(j + nbasis*block)) interleaved (re,im).
void restrict_mrhs(float* coarse_dev,
                   float* const* basis_host, int nbasis,
                   float* const* rhs_host, int nrhs,
                   int Nin, int Nex, int* Nsize, int* block_size);

// fine[r] += fac * (B C)[r].
void prolong_mrhs(float* const* fine_host, int nrhs,
                  const float* coarse_dev,
                  float* const* basis_host, int nbasis, float fac,
                  int Nin, int Nex, int* Nsize, int* block_size);

// v2[r] = D_coarse v1[r] (8-direction coarse Dirac), r = 0..nrhs-1.
void coarse_mrhs(float* const* v2_host, float* const* v1_host,
                 const float* u_dev, int nrhs, int Ncol,
                 int* Nsize, int* bc);

// v2[r] = Clov^-1 v1[r] (site-local coarse clover-inverse, the "prec" source
// preconditioner), r = 0..nrhs-1.  ct_dev = Clov_inv raw device pointer.
void coarse_prec_mrhs(float* const* v2_host, float* const* v1_host,
                      const float* ct_dev, int nrhs, int Ncol, int* Nsize);

// v[r] = D_fine w[r] : batched fine Domainwall_5din D (5dir + hopb), single-GPU
// no-comm path, one launch each, fully on-device.  u_field_host = host base ptr
// of the fine gauge AField.  b/c = Moebius coeffs (length Ns).  do_comm all 0.
// upload the Moebius b/c coefficients (length Ns) into the dev8 module's
// __constant__ memory ONCE; the batched D/Ddag kernels read them from there
// (mass-independent, so a single call serves both the mq and PV operators).
void set_moebius_bc(const float* b, const float* c, int Ns);

// Runtime TF32(tensor-core) <-> FP32(CUDA-core reference) selector for the whole
// MRHS GEMM path (restrict / prolong / coarse / block_inner / block_update).
// Only has an effect when built with USE_TENSORCORE_MRHS (default true there);
// without the macro the path is FP32 reference and set_use_tc(true) is a no-op.
// Lets the live solve A/B TF32 vs FP32 transfer precision without a rebuild.
void set_use_tc(bool on);
bool get_use_tc();

// On-device coarse-vector layout repack (replaces the V-cycle host repack):
// MRHS coarse_dev[:,c] <-> production coarse AField (prod_field host base).
// j = 2i+ch, NWP=32, Ncoarse_in=4*nvec.  No host copy.
void coarse_repack_m2p(float* prod_field, const float* coarse_dev,
                       int c, int s, int nbasis, int nvec, long coarse_nvol);
void coarse_repack_p2m(float* coarse_dev, float* prod_field,
                       int c, int s, int nbasis, int nvec, long coarse_nvol);

void fineD_mrhs(float* const* v_host, float* const* w_host, float* u_field_host,
                int nrhs, float mq, float M0, int Ns, float alpha,
                int* Nsize, int* bc, int* do_comm);

// v[r] = C^{-1} w[r] = U^{-1} L^{-1} w[r] (Prec): 5d-local LU solve, batched,
// one launch per stage, on-device.  e/f/dpinv/dm = operator LU coeff host bases.
// gm5_host (optional): if non-null, also write gamma5*(C^{-1}w) there, to feed a
// following fineDdag_mrhs(gm5_host=...) so it can skip its own gm5 kernel.
void finePrec_mrhs(float* const* v_host, float* const* w_host, int nrhs, int Ns,
                   float* e, float* f, float* dpinv, float* dm, int* Nsize,
                   float* const* gm5_host = nullptr);

// v[r] = (C^{-1})^dag w[r] = Ldag^{-1} Udag^{-1} w[r] (Precdag).
void finePrecdag_mrhs(float* const* v_host, float* const* w_host, int nrhs, int Ns,
                      float* e, float* f, float* dpinv, float* dm, int* Nsize);

// v[r] = Ddag w[r] = 5dirdag( w, hopb( gm5 w ) ): batched fine Ddag, on-device.
// gm5_host (optional): precomputed gamma5*w; if given, skips the internal gm5 kernel.
void fineDdag_mrhs(float* const* v_host, float* const* w_host, float* u_field_host,
                   int nrhs, float mq, float M0, int Ns, float alpha,
                   int* Nsize, int* bc, int* do_comm, float* const* gm5_host = nullptr);

// G[i,j] = <V_i|W_j> (nv x nv complex Gram), reduced over the volume.
void block_inner(float* G_dev, float* const* V_host, float* const* W_host,
                 int nv, int Nin, int Nex, int* Nsize);

// W_j -= sum_i V_i G[i,j]  (block-Arnoldi orthogonalisation update).
void block_update(float* const* W_host, float* const* V_host, const float* G_dev,
                  int nv, int Nin, int Nex, int* Nsize);

// On-device s x s complex Cholesky+inverse: G (Hermitian PD, 2*s*s) -> upper L
// (G=L^H L) in Lout_dev and -L^{-1} in negLi_dev.  All RAW device buffers; keeps
// the block-Krylov dense LA on the GPU (NO per-step host roundtrip).
void block_chol_inv(float* Lout_dev, float* negLi_dev, const float* G_dev, int s);

// ---- batched BLAS-1 over s columns (nflt = floats per field = Nin*Nvol*Nex) --
// Collapse the per-column copy/axpy/scal/norm2 of the block-Krylov solver into a
// single kernel over all s columns (kills the per-column launch overhead that is
// the launch-bound floor on a small lattice).  y/x are AField host base ptrs.
void block_copy (float* const* y, float* const* x, int s, long nflt);  // y[c]=x[c]
void block_axpy (float* const* y, float* const* x, float a, int s, long nflt); // y[c]+=a*x[c]
void block_scal (float* const* y, float a, int s, long nflt);          // y[c]*=a
void block_norm2(double* out_host, float* const* x, int s, long nflt); // out[c]=||x[c]||^2

// ---- small device-memory helpers (keep all CUDA calls inside the .cu) ------
float* dev_alloc(long nfloat);                 // cudaMalloc nfloat*sizeof(float)
void   dev_free(float* p);                     // cudaFree
void   dev_to_host(float* host, const float* dev, long nfloat);   // raw device->host
void   host_to_dev(float* dev, const float* host, long nfloat);   // raw host->device
// map a registered Bridge AField host base pointer to its device pointer (the
// raw device address coarse_mrhs / fineD_mrhs want for the gauge link u_dev).
float* field_dev_ptr(float* field_hostbase);
// download a registered Bridge AField (host base pointer) via dev_ptr->memcpy.
void   field_to_host(float* host, float* field_hostbase, long nfloat);
// upload host buffer into a registered Bridge AField (host base pointer).
void   field_from_host(float* field_hostbase, const float* host, long nfloat);

} // namespace mrhs_live

#endif // MRHS_BLOCK_LIVE_H
