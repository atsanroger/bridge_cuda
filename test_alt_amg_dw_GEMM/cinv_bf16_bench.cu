// Self-contained microbench: AMG fine preconditioner C^{-1} (fused Moebius LU)
// with the FIELD STORAGE type templated (FP32 vs BF16), FP32 register math.
// Measures the per-kernel bandwidth win of halving DRAM storage and the
// accuracy cost of BF16 round-trip, on the real 12^4 x Ls=8 field layout.
// No MPI / no cmake:  nvcc -O3 -arch=sm_86 cinv_bf16_bench.cu -o cinv_bench
//
// The recurrence is copied verbatim from mult_dw5din_Cinv_mrhs_dev (the
// production fused C^{-1}); only the global load/store go through the storage
// type T.  Wei-Lun Chen.
#include <cuda_runtime.h>
#include <cuda_bf16.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>

#define NC    3
#define NVC   6
#define NVCD  24
#define NWP   32
#define ID1   0
#define ID2   6
#define ID3   12
#define ID4   18
#define IDX2(nin,in,ist)  (((ist)%NWP) + NWP*((in) + (nin)*((ist)/NWP)))
#ifndef MRHS_FUSE_LS_MAX
#define MRHS_FUSE_LS_MAX  32
#endif

#define CHECK(call) do { cudaError_t e=(call); if(e!=cudaSuccess){ \
  printf("CUDA error %s:%d %s\n",__FILE__,__LINE__,cudaGetErrorString(e)); exit(1);} } while(0)

template<typename T> __device__ __forceinline__ float ld(const T* p, long i){ return float(p[i]); }
template<typename T> __device__ __forceinline__ void  st(T* p, long i, float v){ p[i] = T(v); }

// Fused C^{-1} = U^{-1} L^{-1}, storage type T, math in float.  One field per
// rhs, contiguous with stride `strd`.  Grid = (site,ivc); loop rhs and is.
template<typename T>
__global__ void cinv_storage(T* vbase, const T* wbase, long strd,
                             int nrhs, int Ns, int Nin5,
                             const float* __restrict__ e, const float* __restrict__ f,
                             const float* __restrict__ dpinv, const float* __restrict__ dm,
                             int Nst_pad)
{
  const int ist = blockIdx.x * blockDim.x + threadIdx.x;
  const int GridSize = blockDim.x * gridDim.x;
  for (int idx = ist; idx < Nst_pad * NVC; idx += GridSize) {
    int idx2 = idx / NWP, idx_in = idx % NWP, ivc = idx2 % NVC, idx_out = idx2 / NVC;
    int site = idx_in + NWP * idx_out;
    for (int r = 0; r < nrhs; ++r) {
      const T* wp = wbase + (long)r * strd;
      T*       vp = vbase + (long)r * strd;
      float s1[MRHS_FUSE_LS_MAX], s2[MRHS_FUSE_LS_MAX], s3[MRHS_FUSE_LS_MAX], s4[MRHS_FUSE_LS_MAX];
      float vt1, vt2, vt3, vt4, yt1, yt2, yt3, yt4, xt1, xt2, xt3, xt4;
      // forward L^{-1}
      vt1 = ld(wp, IDX2(Nin5, (ID1 + ivc), site));
      vt2 = ld(wp, IDX2(Nin5, (ID2 + ivc), site));
      vt3 = ld(wp, IDX2(Nin5, (ID3 + ivc), site));
      vt4 = ld(wp, IDX2(Nin5, (ID4 + ivc), site));
      s1[0]=vt1; s2[0]=vt2; s3[0]=vt3; s4[0]=vt4;
      float e0 = e[0];
      yt1=e0*vt1; yt2=e0*vt2; yt3=e0*vt3; yt4=e0*vt4;
      for (int is = 1; is < Ns - 1; ++is) {
        xt1=vt1; xt2=vt2; xt3=vt3; xt4=vt4;
        vt1 = ld(wp, IDX2(Nin5, (ID1 + ivc + NVCD*is), site));
        vt2 = ld(wp, IDX2(Nin5, (ID2 + ivc + NVCD*is), site));
        vt3 = ld(wp, IDX2(Nin5, (ID3 + ivc + NVCD*is), site));
        vt4 = ld(wp, IDX2(Nin5, (ID4 + ivc + NVCD*is), site));
        float a = 0.5f * dm[is] * dpinv[is-1];
        vt1+=a*(xt1+xt3); vt2+=a*(xt2+xt4); vt3+=a*(xt3+xt1); vt4+=a*(xt4+xt2);
        s1[is]=vt1; s2[is]=vt2; s3[is]=vt3; s4[is]=vt4;
        float eis=e[is]; yt1+=eis*vt1; yt2+=eis*vt2; yt3+=eis*vt3; yt4+=eis*vt4;
      }
      {
        int is=Ns-1;
        xt1=vt1; xt2=vt2; xt3=vt3; xt4=vt4;
        vt1 = ld(wp, IDX2(Nin5, (ID1 + ivc + NVCD*is), site));
        vt2 = ld(wp, IDX2(Nin5, (ID2 + ivc + NVCD*is), site));
        vt3 = ld(wp, IDX2(Nin5, (ID3 + ivc + NVCD*is), site));
        vt4 = ld(wp, IDX2(Nin5, (ID4 + ivc + NVCD*is), site));
        float a = 0.5f * dm[is] * dpinv[is-1];
        vt1+=a*(xt1+xt3); vt2+=a*(xt2+xt4); vt3+=a*(xt3+xt1); vt4+=a*(xt4+xt2);
        vt1+=-0.5f*(yt1-yt3); vt2+=-0.5f*(yt2-yt4); vt3+=-0.5f*(yt3-yt1); vt4+=-0.5f*(yt4-yt2);
        s1[is]=vt1; s2[is]=vt2; s3[is]=vt3; s4[is]=vt4;
      }
      // backward U^{-1}
      int is0=Ns-1; float a0=dpinv[Ns-1];
      vt1=a0*s1[is0]; vt2=a0*s2[is0]; vt3=a0*s3[is0]; vt4=a0*s4[is0];
      st(vp, IDX2(Nin5, (ID1 + ivc + NVCD*is0), site), vt1);
      st(vp, IDX2(Nin5, (ID2 + ivc + NVCD*is0), site), vt2);
      st(vp, IDX2(Nin5, (ID3 + ivc + NVCD*is0), site), vt3);
      st(vp, IDX2(Nin5, (ID4 + ivc + NVCD*is0), site), vt4);
      yt1=0.5f*(vt1+vt3); yt2=0.5f*(vt2+vt4); yt3=0.5f*(vt3+vt1); yt4=0.5f*(vt4+vt2);
      for (int is = Ns-2; is >= 0; --is) {
        xt1=vt1; xt2=vt2; xt3=vt3; xt4=vt4;
        vt1=s1[is]; vt2=s2[is]; vt3=s3[is]; vt4=s4[is];
        float a=0.5f*dm[is];
        vt1+=a*(xt1-xt3); vt2+=a*(xt2-xt4); vt3+=a*(xt3-xt1); vt4+=a*(xt4-xt2);
        float fis=f[is];
        vt1+=-fis*yt1; vt2+=-fis*yt2; vt3+=-fis*yt3; vt4+=-fis*yt4;
        float aa=dpinv[is];
        vt1*=aa; vt2*=aa; vt3*=aa; vt4*=aa;
        st(vp, IDX2(Nin5, (ID1 + ivc + NVCD*is), site), vt1);
        st(vp, IDX2(Nin5, (ID2 + ivc + NVCD*is), site), vt2);
        st(vp, IDX2(Nin5, (ID3 + ivc + NVCD*is), site), vt3);
        st(vp, IDX2(Nin5, (ID4 + ivc + NVCD*is), site), vt4);
      }
    }
  }
}

// rhs-PARALLEL variant: one thread per (site,ivc,rhs).  Work index flattened as
// g = (site,ivc) + Nst_pad*NVC*r so a warp still spans 32 site-lanes (coalesced),
// rhs becomes the outer factor -> nrhs x more concurrent threads to hide latency.
template<typename T>
__global__ void cinv_storage_par(T* vbase, const T* wbase, long strd,
                                 int nrhs, int Ns, int Nin5,
                                 const float* __restrict__ e, const float* __restrict__ f,
                                 const float* __restrict__ dpinv, const float* __restrict__ dm,
                                 int Nst_pad)
{
  const long total = (long)Nst_pad * NVC * nrhs;
  const long g0 = (long)blockIdx.x * blockDim.x + threadIdx.x;
  const long step = (long)blockDim.x * gridDim.x;
  const int  sv = Nst_pad * NVC;
  for (long g = g0; g < total; g += step) {
    int r   = (int)(g / sv);
    int idx = (int)(g % sv);
    int idx2 = idx / NWP, idx_in = idx % NWP, ivc = idx2 % NVC, idx_out = idx2 / NVC;
    int site = idx_in + NWP * idx_out;
    const T* wp = wbase + (long)r * strd;
    T*       vp = vbase + (long)r * strd;
    float s1[MRHS_FUSE_LS_MAX], s2[MRHS_FUSE_LS_MAX], s3[MRHS_FUSE_LS_MAX], s4[MRHS_FUSE_LS_MAX];
    float vt1, vt2, vt3, vt4, yt1, yt2, yt3, yt4, xt1, xt2, xt3, xt4;
    vt1 = ld(wp, IDX2(Nin5, (ID1 + ivc), site));
    vt2 = ld(wp, IDX2(Nin5, (ID2 + ivc), site));
    vt3 = ld(wp, IDX2(Nin5, (ID3 + ivc), site));
    vt4 = ld(wp, IDX2(Nin5, (ID4 + ivc), site));
    s1[0]=vt1; s2[0]=vt2; s3[0]=vt3; s4[0]=vt4;
    float e0=e[0]; yt1=e0*vt1; yt2=e0*vt2; yt3=e0*vt3; yt4=e0*vt4;
    for (int is = 1; is < Ns - 1; ++is) {
      xt1=vt1; xt2=vt2; xt3=vt3; xt4=vt4;
      vt1 = ld(wp, IDX2(Nin5, (ID1 + ivc + NVCD*is), site));
      vt2 = ld(wp, IDX2(Nin5, (ID2 + ivc + NVCD*is), site));
      vt3 = ld(wp, IDX2(Nin5, (ID3 + ivc + NVCD*is), site));
      vt4 = ld(wp, IDX2(Nin5, (ID4 + ivc + NVCD*is), site));
      float a=0.5f*dm[is]*dpinv[is-1];
      vt1+=a*(xt1+xt3); vt2+=a*(xt2+xt4); vt3+=a*(xt3+xt1); vt4+=a*(xt4+xt2);
      s1[is]=vt1; s2[is]=vt2; s3[is]=vt3; s4[is]=vt4;
      float eis=e[is]; yt1+=eis*vt1; yt2+=eis*vt2; yt3+=eis*vt3; yt4+=eis*vt4;
    }
    { int is=Ns-1;
      xt1=vt1; xt2=vt2; xt3=vt3; xt4=vt4;
      vt1 = ld(wp, IDX2(Nin5, (ID1 + ivc + NVCD*is), site));
      vt2 = ld(wp, IDX2(Nin5, (ID2 + ivc + NVCD*is), site));
      vt3 = ld(wp, IDX2(Nin5, (ID3 + ivc + NVCD*is), site));
      vt4 = ld(wp, IDX2(Nin5, (ID4 + ivc + NVCD*is), site));
      float a=0.5f*dm[is]*dpinv[is-1];
      vt1+=a*(xt1+xt3); vt2+=a*(xt2+xt4); vt3+=a*(xt3+xt1); vt4+=a*(xt4+xt2);
      vt1+=-0.5f*(yt1-yt3); vt2+=-0.5f*(yt2-yt4); vt3+=-0.5f*(yt3-yt1); vt4+=-0.5f*(yt4-yt2);
      s1[is]=vt1; s2[is]=vt2; s3[is]=vt3; s4[is]=vt4; }
    int is0=Ns-1; float a0=dpinv[Ns-1];
    vt1=a0*s1[is0]; vt2=a0*s2[is0]; vt3=a0*s3[is0]; vt4=a0*s4[is0];
    st(vp, IDX2(Nin5, (ID1 + ivc + NVCD*is0), site), vt1);
    st(vp, IDX2(Nin5, (ID2 + ivc + NVCD*is0), site), vt2);
    st(vp, IDX2(Nin5, (ID3 + ivc + NVCD*is0), site), vt3);
    st(vp, IDX2(Nin5, (ID4 + ivc + NVCD*is0), site), vt4);
    yt1=0.5f*(vt1+vt3); yt2=0.5f*(vt2+vt4); yt3=0.5f*(vt3+vt1); yt4=0.5f*(vt4+vt2);
    for (int is = Ns-2; is >= 0; --is) {
      xt1=vt1; xt2=vt2; xt3=vt3; xt4=vt4;
      vt1=s1[is]; vt2=s2[is]; vt3=s3[is]; vt4=s4[is];
      float a=0.5f*dm[is];
      vt1+=a*(xt1-xt3); vt2+=a*(xt2-xt4); vt3+=a*(xt3-xt1); vt4+=a*(xt4-xt2);
      float fis=f[is]; vt1+=-fis*yt1; vt2+=-fis*yt2; vt3+=-fis*yt3; vt4+=-fis*yt4;
      float aa=dpinv[is]; vt1*=aa; vt2*=aa; vt3*=aa; vt4*=aa;
      st(vp, IDX2(Nin5, (ID1 + ivc + NVCD*is), site), vt1);
      st(vp, IDX2(Nin5, (ID2 + ivc + NVCD*is), site), vt2);
      st(vp, IDX2(Nin5, (ID3 + ivc + NVCD*is), site), vt3);
      st(vp, IDX2(Nin5, (ID4 + ivc + NVCD*is), site), vt4);
    }
  }
}

// Compile-time-Ns variant: arrays sized [NS] + unrolled loops -> the s[] Ls
// buffer is promoted to REGISTERS (no local memory, no hidden DRAM over-traffic).
// rhs-loop grid like the production kernel.
template<typename T, int NS>
__global__ void cinv_storage_ns(T* vbase, const T* wbase, long strd,
                                int nrhs, int Nin5,
                                const float* __restrict__ e, const float* __restrict__ f,
                                const float* __restrict__ dpinv, const float* __restrict__ dm,
                                int Nst_pad)
{
  const int ist = blockIdx.x * blockDim.x + threadIdx.x;
  const int GridSize = blockDim.x * gridDim.x;
  for (int idx = ist; idx < Nst_pad * NVC; idx += GridSize) {
    int idx2 = idx / NWP, idx_in = idx % NWP, ivc = idx2 % NVC, idx_out = idx2 / NVC;
    int site = idx_in + NWP * idx_out;
    for (int r = 0; r < nrhs; ++r) {
      const T* wp = wbase + (long)r * strd;
      T*       vp = vbase + (long)r * strd;
      float s1[NS], s2[NS], s3[NS], s4[NS];
      float vt1, vt2, vt3, vt4, yt1, yt2, yt3, yt4, xt1, xt2, xt3, xt4;
      vt1 = ld(wp, IDX2(Nin5, (ID1 + ivc), site));
      vt2 = ld(wp, IDX2(Nin5, (ID2 + ivc), site));
      vt3 = ld(wp, IDX2(Nin5, (ID3 + ivc), site));
      vt4 = ld(wp, IDX2(Nin5, (ID4 + ivc), site));
      s1[0]=vt1; s2[0]=vt2; s3[0]=vt3; s4[0]=vt4;
      float e0=e[0]; yt1=e0*vt1; yt2=e0*vt2; yt3=e0*vt3; yt4=e0*vt4;
#pragma unroll
      for (int is = 1; is < NS - 1; ++is) {
        xt1=vt1; xt2=vt2; xt3=vt3; xt4=vt4;
        vt1 = ld(wp, IDX2(Nin5, (ID1 + ivc + NVCD*is), site));
        vt2 = ld(wp, IDX2(Nin5, (ID2 + ivc + NVCD*is), site));
        vt3 = ld(wp, IDX2(Nin5, (ID3 + ivc + NVCD*is), site));
        vt4 = ld(wp, IDX2(Nin5, (ID4 + ivc + NVCD*is), site));
        float a=0.5f*dm[is]*dpinv[is-1];
        vt1+=a*(xt1+xt3); vt2+=a*(xt2+xt4); vt3+=a*(xt3+xt1); vt4+=a*(xt4+xt2);
        s1[is]=vt1; s2[is]=vt2; s3[is]=vt3; s4[is]=vt4;
        float eis=e[is]; yt1+=eis*vt1; yt2+=eis*vt2; yt3+=eis*vt3; yt4+=eis*vt4;
      }
      { int is=NS-1;
        xt1=vt1; xt2=vt2; xt3=vt3; xt4=vt4;
        vt1 = ld(wp, IDX2(Nin5, (ID1 + ivc + NVCD*is), site));
        vt2 = ld(wp, IDX2(Nin5, (ID2 + ivc + NVCD*is), site));
        vt3 = ld(wp, IDX2(Nin5, (ID3 + ivc + NVCD*is), site));
        vt4 = ld(wp, IDX2(Nin5, (ID4 + ivc + NVCD*is), site));
        float a=0.5f*dm[is]*dpinv[is-1];
        vt1+=a*(xt1+xt3); vt2+=a*(xt2+xt4); vt3+=a*(xt3+xt1); vt4+=a*(xt4+xt2);
        vt1+=-0.5f*(yt1-yt3); vt2+=-0.5f*(yt2-yt4); vt3+=-0.5f*(yt3-yt1); vt4+=-0.5f*(yt4-yt2);
        s1[is]=vt1; s2[is]=vt2; s3[is]=vt3; s4[is]=vt4; }
      int is0=NS-1; float a0=dpinv[NS-1];
      vt1=a0*s1[is0]; vt2=a0*s2[is0]; vt3=a0*s3[is0]; vt4=a0*s4[is0];
      st(vp, IDX2(Nin5, (ID1 + ivc + NVCD*is0), site), vt1);
      st(vp, IDX2(Nin5, (ID2 + ivc + NVCD*is0), site), vt2);
      st(vp, IDX2(Nin5, (ID3 + ivc + NVCD*is0), site), vt3);
      st(vp, IDX2(Nin5, (ID4 + ivc + NVCD*is0), site), vt4);
      yt1=0.5f*(vt1+vt3); yt2=0.5f*(vt2+vt4); yt3=0.5f*(vt3+vt1); yt4=0.5f*(vt4+vt2);
#pragma unroll
      for (int is = NS-2; is >= 0; --is) {
        xt1=vt1; xt2=vt2; xt3=vt3; xt4=vt4;
        vt1=s1[is]; vt2=s2[is]; vt3=s3[is]; vt4=s4[is];
        float a=0.5f*dm[is];
        vt1+=a*(xt1-xt3); vt2+=a*(xt2-xt4); vt3+=a*(xt3-xt1); vt4+=a*(xt4-xt2);
        float fis=f[is]; vt1+=-fis*yt1; vt2+=-fis*yt2; vt3+=-fis*yt3; vt4+=-fis*yt4;
        float aa=dpinv[is]; vt1*=aa; vt2*=aa; vt3*=aa; vt4*=aa;
        st(vp, IDX2(Nin5, (ID1 + ivc + NVCD*is), site), vt1);
        st(vp, IDX2(Nin5, (ID2 + ivc + NVCD*is), site), vt2);
        st(vp, IDX2(Nin5, (ID3 + ivc + NVCD*is), site), vt3);
        st(vp, IDX2(Nin5, (ID4 + ivc + NVCD*is), site), vt4);
      }
    }
  }
}

template<typename T>
float time_ns(T* v, const T* w, long strd, int nrhs, int Nin5,
              const float* e, const float* f, const float* dpinv, const float* dm,
              int Nst_pad, int iters)
{
  int bs=256; int gs=(Nst_pad*NVC + bs-1)/bs;
  for(int i=0;i<5;i++) cinv_storage_ns<T,8><<<gs,bs>>>(v,w,strd,nrhs,Nin5,e,f,dpinv,dm,Nst_pad);
  CHECK(cudaDeviceSynchronize());
  cudaEvent_t a,b; cudaEventCreate(&a); cudaEventCreate(&b);
  cudaEventRecord(a);
  for(int i=0;i<iters;i++) cinv_storage_ns<T,8><<<gs,bs>>>(v,w,strd,nrhs,Nin5,e,f,dpinv,dm,Nst_pad);
  cudaEventRecord(b); CHECK(cudaEventSynchronize(b));
  float ms=0; cudaEventElapsedTime(&ms,a,b); cudaEventDestroy(a); cudaEventDestroy(b);
  return ms/iters;
}

// Production-style gm5 kernel: ONE thread per site (loops is, ivc inside) -- the
// same poor parallelization (1 thread/site) as mult_dw5din_gm5_mrhs_dev. This is
// what the gm5-fold removes from fineDdag.
__global__ void gm5_kernel(float* vbase, const float* wbase, long strd,
                           int nrhs, int Ns, int Nin5, int Nst)
{
  int site = blockIdx.x * blockDim.x + threadIdx.x;
  if (site < Nst) {
    for (int r = 0; r < nrhs; ++r) {
      const float* wp = wbase + (long)r * strd; float* vp = vbase + (long)r * strd;
      for (int is = 0; is < Ns; ++is) {
        float wt[24];
        for (int i = 0; i < 24; ++i) wt[i] = wp[IDX2(Nin5, (i + NVCD*is), site)];
        for (int ivc = 0; ivc < NVC; ++ivc) {
          vp[IDX2(Nin5, (ID1 + ivc + NVCD*is), site)] = wt[ID3 + ivc];
          vp[IDX2(Nin5, (ID2 + ivc + NVCD*is), site)] = wt[ID4 + ivc];
          vp[IDX2(Nin5, (ID3 + ivc + NVCD*is), site)] = wt[ID1 + ivc];
          vp[IDX2(Nin5, (ID4 + ivc + NVCD*is), site)] = wt[ID2 + ivc];
        }
      }
    }
  }
}

static float time_gm5(float* v, const float* w, long strd, int nrhs, int Ns, int Nin5,
                      int Nst, int iters)
{
  int bs = 256; int gs = (Nst + bs - 1) / bs;          // 1 thread/site (as production)
  for (int i = 0; i < 5; i++) gm5_kernel<<<gs,bs>>>(v,w,strd,nrhs,Ns,Nin5,Nst);
  CHECK(cudaDeviceSynchronize());
  cudaEvent_t a,b; cudaEventCreate(&a); cudaEventCreate(&b);
  cudaEventRecord(a);
  for (int i = 0; i < iters; i++) gm5_kernel<<<gs,bs>>>(v,w,strd,nrhs,Ns,Nin5,Nst);
  cudaEventRecord(b); CHECK(cudaEventSynchronize(b));
  float ms=0; cudaEventElapsedTime(&ms,a,b); cudaEventDestroy(a); cudaEventDestroy(b);
  return ms/iters;
}

template<typename T> __global__ void cvt_from_f(T* o, const float* in, long n)
{ for(long i=blockIdx.x*(long)blockDim.x+threadIdx.x;i<n;i+=(long)blockDim.x*gridDim.x) o[i]=T(in[i]); }
template<typename T> __global__ void cvt_to_f(float* o, const T* in, long n)
{ for(long i=blockIdx.x*(long)blockDim.x+threadIdx.x;i<n;i+=(long)blockDim.x*gridDim.x) o[i]=float(in[i]); }

template<typename T>
float time_kernel(T* v, const T* w, long strd, int nrhs, int Ns, int Nin5,
                  const float* e, const float* f, const float* dpinv, const float* dm,
                  int Nst_pad, int iters)
{
  int bs=256; int gs=(Nst_pad*NVC + bs-1)/bs; if(gs>65535) gs=65535;
  for(int i=0;i<5;i++) cinv_storage<T><<<gs,bs>>>(v,w,strd,nrhs,Ns,Nin5,e,f,dpinv,dm,Nst_pad);
  CHECK(cudaDeviceSynchronize());
  cudaEvent_t a,b; cudaEventCreate(&a); cudaEventCreate(&b);
  cudaEventRecord(a);
  for(int i=0;i<iters;i++) cinv_storage<T><<<gs,bs>>>(v,w,strd,nrhs,Ns,Nin5,e,f,dpinv,dm,Nst_pad);
  cudaEventRecord(b); CHECK(cudaEventSynchronize(b));
  float ms=0; cudaEventElapsedTime(&ms,a,b);
  cudaEventDestroy(a); cudaEventDestroy(b);
  return ms/iters;
}

template<typename T>
float time_par(T* v, const T* w, long strd, int nrhs, int Ns, int Nin5,
               const float* e, const float* f, const float* dpinv, const float* dm,
               int Nst_pad, int iters)
{
  int bs=256; long total=(long)Nst_pad*NVC*nrhs;
  int gs=(int)((total + bs-1)/bs); if(gs>2147483000) gs=2147483000;
  for(int i=0;i<5;i++) cinv_storage_par<T><<<gs,bs>>>(v,w,strd,nrhs,Ns,Nin5,e,f,dpinv,dm,Nst_pad);
  CHECK(cudaDeviceSynchronize());
  cudaEvent_t a,b; cudaEventCreate(&a); cudaEventCreate(&b);
  cudaEventRecord(a);
  for(int i=0;i<iters;i++) cinv_storage_par<T><<<gs,bs>>>(v,w,strd,nrhs,Ns,Nin5,e,f,dpinv,dm,Nst_pad);
  cudaEventRecord(b); CHECK(cudaEventSynchronize(b));
  float ms=0; cudaEventElapsedTime(&ms,a,b);
  cudaEventDestroy(a); cudaEventDestroy(b);
  return ms/iters;
}

int main(int argc, char** argv)
{
  int Ns=8;
  int nrhs = (argc>1)? atoi(argv[1]) : 12;
  int iters= (argc>2)? atoi(argv[2]) : 200;
  int L    = (argc>3)? atoi(argv[3]) : 12;   // cubic L^4 lattice (size-dep check)
  int Nx=L,Ny=L,Nz=L,Nt=L;
  int Nst = Nx*Ny*Nz*Nt;
  int Nst_pad = ((Nst + NWP-1)/NWP)*NWP;
  int Nin5 = NVCD*Ns;
  long strd = (long)Nin5 * Nst_pad;          // floats per rhs
  long n    = strd * nrhs;
  double GB = 1.0e9;

  cudaDeviceProp prop; CHECK(cudaGetDeviceProperties(&prop,0));
  int memClkKHz=0, busBits=0;
  cudaDeviceGetAttribute(&memClkKHz, cudaDevAttrMemoryClockRate, 0);
  cudaDeviceGetAttribute(&busBits,  cudaDevAttrGlobalMemoryBusWidth, 0);
  double peakBW = 2.0*memClkKHz*1e3*(busBits/8.0)/GB; // GB/s (GDDR DDR factor 2)
  printf("GPU %s  peak BW ~ %.0f GB/s   |  12^4 Ls=%d  nrhs=%d  iters=%d\n",
         prop.name, Ns, nrhs, iters);
  printf("field/rhs = %ld floats (%.1f MB fp32)   total fp32 in+out = %.0f MB\n",
         strd, strd*4.0/1e6, 2.0*n*4.0/1e6);

  // host random input spanning a few decades (mimic real fine field range)
  std::vector<float> hw(n);
  srand(1234);
  for(long i=0;i<n;i++){ double u=(rand()/(double)RAND_MAX); double mag=pow(10.0,-6.0*u);
                         hw[i]=(float)((rand()&1?1:-1)*mag); }
  // coefficient arrays (representative Moebius-LU magnitudes; identical for both
  // precisions, so they don't bias the FP32-vs-BF16 comparison)
  std::vector<float> he(Ns),hf(Ns),hdp(Ns),hdm(Ns);
  for(int s=0;s<Ns;s++){ hdp[s]=0.8f; hdm[s]=0.3f; he[s]=powf(0.5f,s+1); hf[s]=powf(0.5f,Ns-s); }

  float *dw_f,*dv_f,*de,*df,*ddp,*ddm,*dref;
  CHECK(cudaMalloc(&dw_f,n*4)); CHECK(cudaMalloc(&dv_f,n*4)); CHECK(cudaMalloc(&dref,n*4));
  CHECK(cudaMalloc(&de,Ns*4)); CHECK(cudaMalloc(&df,Ns*4)); CHECK(cudaMalloc(&ddp,Ns*4)); CHECK(cudaMalloc(&ddm,Ns*4));
  CHECK(cudaMemcpy(dw_f,hw.data(),n*4,cudaMemcpyHostToDevice));
  CHECK(cudaMemcpy(de,he.data(),Ns*4,cudaMemcpyHostToDevice));
  CHECK(cudaMemcpy(df,hf.data(),Ns*4,cudaMemcpyHostToDevice));
  CHECK(cudaMemcpy(ddp,hdp.data(),Ns*4,cudaMemcpyHostToDevice));
  CHECK(cudaMemcpy(ddm,hdm.data(),Ns*4,cudaMemcpyHostToDevice));

  // ---- FP32 reference + timing ----
  float t_f = time_kernel<float>(dv_f,dw_f,strd,nrhs,Ns,Nin5,de,df,ddp,ddm,Nst_pad,iters);
  CHECK(cudaMemcpy(dref,dv_f,n*4,cudaMemcpyDeviceToDevice));
  double bytes_f = 2.0*n*4.0;
  printf("\nFP32  : %.4f ms/apply   %.0f GB/s  (%.0f%% peak)\n",
         t_f, bytes_f/(t_f*1e-3)/GB, 100.0*bytes_f/(t_f*1e-3)/GB/peakBW);

  // ---- BF16 storage + timing ----
  __nv_bfloat16 *dw_b,*dv_b;
  CHECK(cudaMalloc(&dw_b,n*2)); CHECK(cudaMalloc(&dv_b,n*2));
  { int bs=256; long gs=(n+bs-1)/bs; if(gs>65535)gs=65535;
    cvt_from_f<__nv_bfloat16><<<gs,bs>>>(dw_b,dw_f,n); CHECK(cudaDeviceSynchronize()); }
  float t_b = time_kernel<__nv_bfloat16>(dv_b,dw_b,strd,nrhs,Ns,Nin5,de,df,ddp,ddm,Nst_pad,iters);
  double bytes_b = 2.0*n*2.0;
  printf("BF16  : %.4f ms/apply   %.0f GB/s  (%.0f%% peak)   speedup x%.2f\n",
         t_b, bytes_b/(t_b*1e-3)/GB, 100.0*bytes_b/(t_b*1e-3)/GB/peakBW, t_f/t_b);

  // ---- accuracy: BF16 output vs FP32 reference ----
  float* dvb_f; CHECK(cudaMalloc(&dvb_f,n*4));
  { int bs=256; long gs=(n+bs-1)/bs; if(gs>65535)gs=65535;
    cvt_to_f<__nv_bfloat16><<<gs,bs>>>(dvb_f,dv_b,n); CHECK(cudaDeviceSynchronize()); }
  std::vector<float> vref(n), vbf(n);
  CHECK(cudaMemcpy(vref.data(),dref,n*4,cudaMemcpyDeviceToHost));
  CHECK(cudaMemcpy(vbf.data(),dvb_f,n*4,cudaMemcpyDeviceToHost));
  double num=0,den=0,maxrel=0;
  for(long i=0;i<n;i++){ double d=vbf[i]-vref[i]; num+=d*d; den+=(double)vref[i]*vref[i];
    double a=fabs(vref[i]); if(a>1e-6){ double rel=fabs(d)/a; if(rel>maxrel)maxrel=rel; } }
  printf("\nBF16 vs FP32 output:  rel L2 = %.2e   max elt rel (|ref|>1e-6) = %.2e\n",
         sqrt(num/den), maxrel);
  printf("(BF16 has 7 mantissa bits ~ 2^-8 = 3.9e-3 floor; preconditioner budget ~1e-3..1e-2)\n");

  // ============ rhs-PARALLEL variant (latency fix) ============
  printf("\n----- rhs-PARALLEL (one thread per site,ivc,rhs) -----\n");
  float tp_f = time_par<float>(dv_f,dw_f,strd,nrhs,Ns,Nin5,de,df,ddp,ddm,Nst_pad,iters);
  printf("FP32  : %.4f ms/apply   %.0f GB/s  (%.0f%% peak)   vs rhs-loop x%.2f\n",
         tp_f, bytes_f/(tp_f*1e-3)/GB, 100.0*bytes_f/(tp_f*1e-3)/GB/peakBW, t_f/tp_f);
  float tp_b = time_par<__nv_bfloat16>(dv_b,dw_b,strd,nrhs,Ns,Nin5,de,df,ddp,ddm,Nst_pad,iters);
  printf("BF16  : %.4f ms/apply   %.0f GB/s  (%.0f%% peak)   BF16/FP32 x%.2f   total vs FP32-rhsloop x%.2f\n",
         tp_b, bytes_b/(tp_b*1e-3)/GB, 100.0*bytes_b/(tp_b*1e-3)/GB/peakBW, tp_f/tp_b, t_f/tp_b);

  // ============ compile-time Ns=8 (s[] -> registers) ============
  printf("\n----- compile-time NS=8 (s[] promoted to registers) -----\n");
  float tn_f = time_ns<float>(dv_f,dw_f,strd,nrhs,Nin5,de,df,ddp,ddm,Nst_pad,iters);
  printf("FP32  : %.4f ms/apply   %.0f GB/s(useful)   vs FP32-rhsloop x%.2f\n",
         tn_f, bytes_f/(tn_f*1e-3)/GB, t_f/tn_f);
  { std::vector<float> vp(n); CHECK(cudaMemcpy(vp.data(),dv_f,n*4,cudaMemcpyDeviceToHost));
    double e2=0,d2=0; for(long i=0;i<n;i++){double d=vp[i]-vref[i];e2+=d*d;d2+=(double)vref[i]*vref[i];}
    printf("  NS=8 FP32 vs rhs-loop FP32: rel L2 = %.2e (should be ~0)\n", sqrt(e2/d2)); }
  float tn_b = time_ns<__nv_bfloat16>(dv_b,dw_b,strd,nrhs,Nin5,de,df,ddp,ddm,Nst_pad,iters);
  printf("BF16  : %.4f ms/apply   %.0f GB/s(useful)   BF16/FP32 x%.2f   total vs FP32-rhsloop x%.2f\n",
         tn_b, bytes_b/(tn_b*1e-3)/GB, tn_f/tn_b, t_f/tn_b);

  // ============ standalone gm5 kernel cost (what the fold removes) ============
  // fineDdag normally launches mult_dw5din_gm5_mrhs_dev (1 thread/site) before
  // hopb.  The finePrec->fineDdag gm5-fold lets Cinv write gm5(Cinv w) directly,
  // so this whole launch disappears (Cinv pays only ~1 extra write pass).
  printf("\n----- standalone gm5 kernel (removed by the finePrec->fineDdag fold) -----\n");
  float t_gm5 = time_gm5(dv_f,dw_f,strd,nrhs,Ns,Nin5,Nst_pad,iters);
  printf("gm5   : %.4f ms/apply   gm5/Cinv(NS=8 FP32) = %.2f   (fold removes ~this much)\n",
         t_gm5, t_gm5/tn_f);

  // verify rhs-parallel FP32 matches the rhs-loop reference (same math)
  { std::vector<float> vpar(n); CHECK(cudaMemcpy(vpar.data(),dv_f,n*4,cudaMemcpyDeviceToHost));
    double e2=0,d2=0; for(long i=0;i<n;i++){ double d=vpar[i]-vref[i]; e2+=d*d; d2+=(double)vref[i]*vref[i]; }
    printf("rhs-parallel FP32 vs rhs-loop FP32: rel L2 = %.2e (should be ~0)\n", sqrt(e2/d2)); }
  return 0;
}
