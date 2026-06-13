/*!
      @file    bridgeACC_mrhs_block_float.cu
      @brief   Compiles the MRHS / tensor-core AMG kernels into a BridgeCUDA
               library (libbridgecuda_mrhs) and exposes C++-callable forwarders
               (mrhs_live::*) for the AMG-on-A block-FGMRES solver.
      @author  Wei-Lun Chen (wlchen)

      Include recipe mirrors bridgeACC_Domainwall_coarse_float.cu (the proven
      BridgeACC translation unit) so dev_ptr / IDX2 / NWP / CHECK / real_t all
      resolve exactly as for the production single-RHS coarse operator.  The
      MRHS host-drivers live in an anonymous namespace inside the included
      header, so we bridge them out through small non-anonymous BridgeACC
      forwarders, which the global mrhs_live:: API then calls.

      Build: compiled as the separate static lib `bridgecuda_mrhs` (TF32, sm_80+,
      C++17) declared in BridgeCUDA/CMakeLists.txt -- additive, so the existing
      libbridgecuda and the shared fermion/Corr2pt path are NOT modified.  Tests
      that want the batched AMG kernels just link bridgecuda_mrhs.
*/

#include "inline/define_params.h"
#include "inline/define_index.h"

#include "bridgeACC.h"
#include "bridgeACC_AField.h"

#include "mrhs_block_live.h"

using namespace std;

typedef float real_t;

namespace BridgeACC {

// afield_dd_kernel_sync() is a (debug-only) no-op in afield_dd_cuda-inc.h; the
// MRHS header calls it after every launch.  Provide it locally so we do not pull
// the whole afield_dd translation unit into this small wrapper.
#ifndef AFIELD_DD_SYNC_EACH_KERNEL
#define AFIELD_DD_SYNC_EACH_KERNEL 0
#endif
static inline void afield_dd_kernel_sync()
{
#if AFIELD_DD_SYNC_EACH_KERNEL
  CHECK(cudaDeviceSynchronize());
#endif
}

#include "src/mrhs_block_tensorcore_cuda-inc.h"
#include "src/mult_Domainwall_5din_mrhs_cuda-inc.h"

// ---- non-anonymous bridges (callable from outside the anonymous namespace) --
void mrhs_fwd_restrict(real_t* coarse_dev, real_t* const* basis, int nbasis,
                       real_t* const* rhs, int nrhs,
                       int Nin, int Nex, int* Nsize, int* block_size)
{ block_restrict_mrhs(coarse_dev, basis, nbasis, rhs, nrhs,
                      Nin, Nex, Nsize, block_size); }

void mrhs_fwd_prolong(real_t* const* fine, int nrhs, const real_t* coarse_dev,
                      real_t* const* basis, int nbasis, real_t fac,
                      int Nin, int Nex, int* Nsize, int* block_size)
{ block_prolong_mrhs(fine, nrhs, coarse_dev, basis, nbasis, fac,
                     Nin, Nex, Nsize, block_size); }

void mrhs_fwd_coarse(real_t* const* v2, real_t* const* v1, const real_t* u_dev,
                     int nrhs, int Ncol, int* Nsize, int* bc)
{ mult_coarse_mrhs(v2, v1, u_dev, nrhs, Ncol, Nsize, bc); }

void mrhs_fwd_coarse_prec(real_t* const* v2, real_t* const* v1, const real_t* ct_dev,
                          int nrhs, int Ncol, int* Nsize)
{ mult_coarse_prec_mrhs(v2, v1, ct_dev, nrhs, Ncol, Nsize); }

void mrhs_fwd_inner(real_t* G_dev, real_t* const* V, real_t* const* W,
                    int nv, int Nin, int Nex, int* Nsize)
{ block_inner_mrhs(G_dev, V, W, nv, Nin, Nex, Nsize); }

void mrhs_fwd_update(real_t* const* W, real_t* const* V, const real_t* G_dev,
                     int nv, int Nin, int Nex, int* Nsize)
{ block_update_mrhs(W, V, G_dev, nv, Nin, Nex, Nsize); }

void mrhs_fwd_chol_inv(real_t* Lout_dev, real_t* negLi_dev, const real_t* G_dev, int s)
{ block_chol_inv(Lout_dev, negLi_dev, G_dev, s); }

void mrhs_fwd_set_moebius_bc(const real_t* b, const real_t* c, int Ns)
{ mrhs_set_moebius_bc_dev(b, c, Ns); }

// runtime TF32(tensor-core) <-> FP32(CUDA-core ref) selector for the GEMM path
void mrhs_fwd_set_use_tc(bool on) { mrhs_set_use_tc_dev(on); }
bool mrhs_fwd_get_use_tc()        { return mrhs_get_use_tc_dev(); }

void mrhs_fwd_coarse_repack_m2p(real_t* prod_field, const real_t* coarse_dev,
                                int c, int s, int nbasis, int nvec, long coarse_nvol)
{ coarse_repack_m2p(prod_field, coarse_dev, c, s, nbasis, nvec, coarse_nvol); }

void mrhs_fwd_coarse_repack_p2m(real_t* coarse_dev, real_t* prod_field,
                                int c, int s, int nbasis, int nvec, long coarse_nvol)
{ coarse_repack_p2m(coarse_dev, prod_field, c, s, nbasis, nvec, coarse_nvol); }

void mrhs_fwd_fineD(real_t* const* v, real_t* const* w, real_t* u_field_host,
                    int nrhs, real_t mq, real_t M0, int Ns, real_t alpha,
                    int* Nsize, int* bc, int* do_comm)
{ fineD_mrhs(v, w, u_field_host, nrhs, mq, M0, Ns, alpha, Nsize, bc, do_comm); }

void mrhs_fwd_finePrec(real_t* const* v, real_t* const* w, int nrhs, int Ns,
                       real_t* e, real_t* f, real_t* dpinv, real_t* dm, int* Nsize)
{ finePrec_mrhs(v, w, nrhs, Ns, e, f, dpinv, dm, Nsize); }

void mrhs_fwd_finePrecdag(real_t* const* v, real_t* const* w, int nrhs, int Ns,
                          real_t* e, real_t* f, real_t* dpinv, real_t* dm, int* Nsize)
{ finePrecdag_mrhs(v, w, nrhs, Ns, e, f, dpinv, dm, Nsize); }

void mrhs_fwd_fineDdag(real_t* const* v, real_t* const* w, real_t* u_field_host,
                       int nrhs, real_t mq, real_t M0, int Ns, real_t alpha,
                       int* Nsize, int* bc, int* do_comm)
{ fineDdag_mrhs(v, w, u_field_host, nrhs, mq, M0, Ns, alpha, Nsize, bc, do_comm); }

// device-memory helpers (dev_ptr is in this namespace)
float* mrhs_fwd_alloc(long n)
{ float* p = nullptr; CHECK(cudaMalloc((void**)&p, sizeof(float)*n)); return p; }
void mrhs_fwd_free(float* p) { if(p) cudaFree(p); }
void mrhs_fwd_d2h(float* host, const float* dev, long n)
{ CHECK(cudaMemcpy(host, dev, sizeof(float)*n, cudaMemcpyDeviceToHost)); }
void mrhs_fwd_h2d(float* dev, const float* host, long n)
{ CHECK(cudaMemcpy(dev, host, sizeof(float)*n, cudaMemcpyHostToDevice)); }
float* mrhs_fwd_field_devptr(float* field_hostbase)
{ return (float*)dev_ptr(field_hostbase); }
void mrhs_fwd_field2h(float* host, float* field_hostbase, long n)
{ float* dev = (float*)dev_ptr(field_hostbase);
  CHECK(cudaMemcpy(host, dev, sizeof(float)*n, cudaMemcpyDeviceToHost)); }
void mrhs_fwd_h2field(float* field_hostbase, const float* host, long n)
{ float* dev = (float*)dev_ptr(field_hostbase);
  CHECK(cudaMemcpy(dev, host, sizeof(float)*n, cudaMemcpyHostToDevice)); }

} // namespace BridgeACC

// ---- public C++ API (mrhs_live::) -----------------------------------------
namespace mrhs_live {

void restrict_mrhs(float* coarse_dev, float* const* basis_host, int nbasis,
                   float* const* rhs_host, int nrhs,
                   int Nin, int Nex, int* Nsize, int* block_size)
{ BridgeACC::mrhs_fwd_restrict(coarse_dev, basis_host, nbasis, rhs_host, nrhs,
                               Nin, Nex, Nsize, block_size); }

void prolong_mrhs(float* const* fine_host, int nrhs, const float* coarse_dev,
                  float* const* basis_host, int nbasis, float fac,
                  int Nin, int Nex, int* Nsize, int* block_size)
{ BridgeACC::mrhs_fwd_prolong(fine_host, nrhs, coarse_dev, basis_host, nbasis,
                              fac, Nin, Nex, Nsize, block_size); }

void coarse_mrhs(float* const* v2_host, float* const* v1_host,
                 const float* u_dev, int nrhs, int Ncol, int* Nsize, int* bc)
{ BridgeACC::mrhs_fwd_coarse(v2_host, v1_host, u_dev, nrhs, Ncol, Nsize, bc); }

void coarse_prec_mrhs(float* const* v2_host, float* const* v1_host,
                      const float* ct_dev, int nrhs, int Ncol, int* Nsize)
{ BridgeACC::mrhs_fwd_coarse_prec(v2_host, v1_host, ct_dev, nrhs, Ncol, Nsize); }

void block_inner(float* G_dev, float* const* V_host, float* const* W_host,
                 int nv, int Nin, int Nex, int* Nsize)
{ BridgeACC::mrhs_fwd_inner(G_dev, V_host, W_host, nv, Nin, Nex, Nsize); }

void block_update(float* const* W_host, float* const* V_host, const float* G_dev,
                  int nv, int Nin, int Nex, int* Nsize)
{ BridgeACC::mrhs_fwd_update(W_host, V_host, G_dev, nv, Nin, Nex, Nsize); }

void block_chol_inv(float* Lout_dev, float* negLi_dev, const float* G_dev, int s)
{ BridgeACC::mrhs_fwd_chol_inv(Lout_dev, negLi_dev, G_dev, s); }

void set_moebius_bc(const float* b, const float* c, int Ns)
{ BridgeACC::mrhs_fwd_set_moebius_bc(b, c, Ns); }

void set_use_tc(bool on) { BridgeACC::mrhs_fwd_set_use_tc(on); }
bool get_use_tc()        { return BridgeACC::mrhs_fwd_get_use_tc(); }

void coarse_repack_m2p(float* prod_field, const float* coarse_dev,
                       int c, int s, int nbasis, int nvec, long coarse_nvol)
{ BridgeACC::mrhs_fwd_coarse_repack_m2p(prod_field, coarse_dev, c, s, nbasis, nvec, coarse_nvol); }

void coarse_repack_p2m(float* coarse_dev, float* prod_field,
                       int c, int s, int nbasis, int nvec, long coarse_nvol)
{ BridgeACC::mrhs_fwd_coarse_repack_p2m(coarse_dev, prod_field, c, s, nbasis, nvec, coarse_nvol); }

void fineD_mrhs(float* const* v_host, float* const* w_host, float* u_field_host,
                int nrhs, float mq, float M0, int Ns, float alpha,
                int* Nsize, int* bc, int* do_comm)
{ BridgeACC::mrhs_fwd_fineD(v_host, w_host, u_field_host, nrhs, mq, M0, Ns, alpha,
                            Nsize, bc, do_comm); }

void finePrec_mrhs(float* const* v_host, float* const* w_host, int nrhs, int Ns,
                   float* e, float* f, float* dpinv, float* dm, int* Nsize)
{ BridgeACC::mrhs_fwd_finePrec(v_host, w_host, nrhs, Ns, e, f, dpinv, dm, Nsize); }

void finePrecdag_mrhs(float* const* v_host, float* const* w_host, int nrhs, int Ns,
                      float* e, float* f, float* dpinv, float* dm, int* Nsize)
{ BridgeACC::mrhs_fwd_finePrecdag(v_host, w_host, nrhs, Ns, e, f, dpinv, dm, Nsize); }

void fineDdag_mrhs(float* const* v_host, float* const* w_host, float* u_field_host,
                   int nrhs, float mq, float M0, int Ns, float alpha,
                   int* Nsize, int* bc, int* do_comm)
{ BridgeACC::mrhs_fwd_fineDdag(v_host, w_host, u_field_host, nrhs, mq, M0, Ns, alpha,
                               Nsize, bc, do_comm); }

float* dev_alloc(long nfloat) { return BridgeACC::mrhs_fwd_alloc(nfloat); }
void   dev_free(float* p) { BridgeACC::mrhs_fwd_free(p); }
void   dev_to_host(float* host, const float* dev, long nfloat)
{ BridgeACC::mrhs_fwd_d2h(host, dev, nfloat); }
void   host_to_dev(float* dev, const float* host, long nfloat)
{ BridgeACC::mrhs_fwd_h2d(dev, host, nfloat); }
float* field_dev_ptr(float* field_hostbase)
{ return BridgeACC::mrhs_fwd_field_devptr(field_hostbase); }
void   field_to_host(float* host, float* field_hostbase, long nfloat)
{ BridgeACC::mrhs_fwd_field2h(host, field_hostbase, nfloat); }
void   field_from_host(float* field_hostbase, const float* host, long nfloat)
{ BridgeACC::mrhs_fwd_h2field(field_hostbase, host, nfloat); }

} // namespace mrhs_live
