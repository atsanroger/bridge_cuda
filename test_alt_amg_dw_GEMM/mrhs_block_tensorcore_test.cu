#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <cuda_runtime.h>

typedef float real_t;
#define NWP 32
#define VECTOR_LENGTH 64
#define IDX2(nin, in, ist)  (((ist)%NWP) + NWP*((in) + (nin)*((ist)/NWP)))
#define CHECK(call) do { cudaError_t e=(call); if(e!=cudaSuccess){ \
  printf("CUDA err %s:%d %s\n",__FILE__,__LINE__,cudaGetErrorString(e)); exit(1);} } while(0)
void* dev_ptr(void* p){ return p; }
namespace { static inline void afield_dd_kernel_sync(){ CHECK(cudaDeviceSynchronize()); } }

#define MRHS_VARIANT_GUARD
#include "/home/wlchenhome/Bridge/with_alpha/src/lib_alt_Accel/BridgeCUDA/src/mrhs_block_tensorcore_cuda-inc.h"

int main(){
  int Ns[4]={4,4,4,4}, Bs[4]={2,2,2,2};
  int Nst=4*4*4*4;             // 256, multiple of NWP
  int Nin=8, Nex=1;            // 4 complex inner dof
  int nbasis=6, nrhs=5;
  int Nblock=(4/2)*(4/2)*(4/2)*(4/2);
  size_t fsz = (size_t)Nin*Nst*Nex;

  // host fields, deterministic pseudo-random
  std::vector<std::vector<real_t>> bh(nbasis, std::vector<real_t>(fsz));
  std::vector<std::vector<real_t>> rh(nrhs,   std::vector<real_t>(fsz));
  auto rnd=[](int s){ return (float)((s*2654435761u)%10007)/10007.f - 0.5f; };
  for(int j=0;j<nbasis;j++) for(size_t i=0;i<fsz;i++) bh[j][i]=rnd(j*7919+i);
  for(int r=0;r<nrhs;r++)   for(size_t i=0;i<fsz;i++) rh[r][i]=rnd(r*104729+i+13);

  std::vector<real_t*> bd(nbasis), rd(nrhs);
  for(int j=0;j<nbasis;j++){ CHECK(cudaMalloc(&bd[j],fsz*sizeof(real_t)));
    CHECK(cudaMemcpy(bd[j],bh[j].data(),fsz*sizeof(real_t),cudaMemcpyHostToDevice)); }
  for(int r=0;r<nrhs;r++){ CHECK(cudaMalloc(&rd[r],fsz*sizeof(real_t)));
    CHECK(cudaMemcpy(rd[r],rh[r].data(),fsz*sizeof(real_t),cudaMemcpyHostToDevice)); }

  size_t csz=2*(size_t)nrhs*nbasis*Nblock;
  real_t *cdev; CHECK(cudaMalloc(&cdev,csz*sizeof(real_t)));

  block_restrict_mrhs(cdev, bd.data(), nbasis, rd.data(), nrhs, Nin, Nex, Ns, Bs);
  std::vector<real_t> cgot(csz); CHECK(cudaMemcpy(cgot.data(),cdev,csz*sizeof(real_t),cudaMemcpyDeviceToHost));

  // CPU reference: C[j,r,block]=sum conj(b).r
  std::vector<real_t> cref(csz,0.f);
  int Bsize=2*2*2*2, NBx=2,NBy=2,NBz=2;
  for(int block=0;block<Nblock;block++){
    int ibx=block%NBx, iby=(block/NBx)%NBy, ibz=(block/(NBx*NBy))%NBz, ibt=block/(NBx*NBy*NBz);
    for(int j=0;j<nbasis;j++) for(int r=0;r<nrhs;r++){
      double cr=0,ci=0;
      for(int bsite=0;bsite<Bsize;bsite++){
        int kx=bsite%2,ix=kx+2*ibx; int ky=(bsite/2)%2,iy=ky+2*iby;
        int kz=(bsite/4)%2,iz=kz+2*ibz; int kt=(bsite/8),it=kt+2*ibt;
        int site=ix+4*(iy+4*(iz+4*it));
        for(int in2=0;in2<Nin/2;in2++){
          float br=bh[j][IDX2(Nin,2*in2,site)],   bi=bh[j][IDX2(Nin,2*in2+1,site)];
          float xr=rh[r][IDX2(Nin,2*in2,site)],   xi=rh[r][IDX2(Nin,2*in2+1,site)];
          cr+=br*xr+bi*xi; ci+=br*xi-bi*xr;
        }
      }
      int o=2*(r+nrhs*(j+nbasis*block)); cref[o]=cr; cref[o+1]=ci;
    }
  }
  double maxabs=0,maxrel=0,nrm=0;
  for(size_t i=0;i<csz;i++){ double d=fabs(cgot[i]-cref[i]); maxabs=fmax(maxabs,d);
    nrm=fmax(nrm,fabs(cref[i])); }
  for(size_t i=0;i<csz;i++){ double den=fabs(cref[i])>1e-4?fabs(cref[i]):1e-4;
    maxrel=fmax(maxrel,fabs(cgot[i]-cref[i])/den); }
  printf("coarse entries=%zu  max|ref|=%.4f  max abs diff=%.3e  max rel diff=%.3e\n",
         csz,nrm,maxabs,maxrel);
  double norm_err = maxabs/nrm;
  printf("normalized error (max abs diff / max|ref|) = %.3e\n", norm_err);
  printf(norm_err<1e-3 ? "PASS: TF32 tensor-core restriction matches FP32 oracle\n"
                       : "FAIL: dataflow/index error\n");
  // ---------- prolongation check: fine += fac * sum_j basis[j]*C[j,r] ----------
  float fac = 0.5f;
  std::vector<real_t*> fd(nrhs);
  std::vector<std::vector<real_t>> f0(nrhs, std::vector<real_t>(fsz));
  for(int r=0;r<nrhs;r++){ for(size_t i=0;i<fsz;i++) f0[r][i]=rnd(r*999+i+5);
    CHECK(cudaMalloc(&fd[r],fsz*sizeof(real_t)));
    CHECK(cudaMemcpy(fd[r],f0[r].data(),fsz*sizeof(real_t),cudaMemcpyHostToDevice)); }
  // use the FP32 reference coarse (cref) as input so we test prolong in isolation
  real_t* cin; CHECK(cudaMalloc(&cin,csz*sizeof(real_t)));
  CHECK(cudaMemcpy(cin,cref.data(),csz*sizeof(real_t),cudaMemcpyHostToDevice));
  block_prolong_mrhs(fd.data(), nrhs, cin, bd.data(), nbasis, fac, Nin, Nex, Ns, Bs);
  std::vector<std::vector<real_t>> fgot(nrhs, std::vector<real_t>(fsz));
  for(int r=0;r<nrhs;r++) CHECK(cudaMemcpy(fgot[r].data(),fd[r],fsz*sizeof(real_t),cudaMemcpyDeviceToHost));
  // CPU reference prolong on top of f0
  std::vector<std::vector<double>> fcpu(nrhs, std::vector<double>(fsz));
  for(int r=0;r<nrhs;r++) for(size_t i=0;i<fsz;i++) fcpu[r][i]=f0[r][i];
  for(int block=0;block<Nblock;block++){
    int ibx=block%NBx, iby=(block/NBx)%NBy, ibz=(block/(NBx*NBy))%NBz, ibt=block/(NBx*NBy*NBz);
    for(int bsite=0;bsite<Bsize;bsite++){
      int kx=bsite%2,ix=kx+2*ibx; int ky=(bsite/2)%2,iy=ky+2*iby;
      int kz=(bsite/4)%2,iz=kz+2*ibz; int kt=(bsite/8),it=kt+2*ibt;
      int site=ix+4*(iy+4*(iz+4*it));
      for(int in2=0;in2<Nin/2;in2++){
        for(int r=0;r<nrhs;r++){
          double vr=0,vi=0;
          for(int j=0;j<nbasis;j++){
            double cr=cref[2*(r+nrhs*(j+nbasis*block))], ci=cref[2*(r+nrhs*(j+nbasis*block))+1];
            double br=bh[j][IDX2(Nin,2*in2,site)], bi=bh[j][IDX2(Nin,2*in2+1,site)];
            vr += br*cr - bi*ci; vi += br*ci + bi*cr;
          }
          fcpu[r][IDX2(Nin,2*in2,site)]   += fac*vr;
          fcpu[r][IDX2(Nin,2*in2+1,site)] += fac*vi;
        }
      }
    }
  }
  double pmax=0,pnrm=0;
  for(int r=0;r<nrhs;r++) for(size_t i=0;i<fsz;i++){
    pmax=fmax(pmax,fabs(fgot[r][i]-fcpu[r][i])); pnrm=fmax(pnrm,fabs(fcpu[r][i])); }
  double pne=pmax/pnrm;
  printf("\nPROLONG: max|ref|=%.4f  max abs diff=%.3e  normalized=%.3e\n",pnrm,pmax,pne);
  printf(pne<1e-3 ? "PASS: TF32 tensor-core prolongation matches FP32 oracle\n"
                  : "FAIL: dataflow/index error\n");

  // ---------- coarse Dirac operator check: V2 = V1 + sum_dir U V1_nei ----------
  {
  int Ncol=8;                  // coarse dof per site (even, <=16)
  int Ncol2=2*Ncol, Ndf=2*Ncol*Ncol;
  size_t vsz=(size_t)Ncol2*Nst;          // one coarse vector (reals)
  size_t usz=(size_t)Ndf*Nst*8;          // 8-direction coarse gauge
  int cnrhs=6;
  std::vector<real_t> uh(usz);
  for(size_t i=0;i<usz;i++) uh[i]=rnd((int)(i*3+1))*0.3f;
  real_t* ud; CHECK(cudaMalloc(&ud,usz*sizeof(real_t)));
  CHECK(cudaMemcpy(ud,uh.data(),usz*sizeof(real_t),cudaMemcpyHostToDevice));
  std::vector<std::vector<real_t>> w1(cnrhs,std::vector<real_t>(vsz));
  for(int r=0;r<cnrhs;r++) for(size_t i=0;i<vsz;i++) w1[r][i]=rnd((int)(r*555+i+9));
  std::vector<real_t*> w1d(cnrhs), w2d(cnrhs);
  for(int r=0;r<cnrhs;r++){ CHECK(cudaMalloc(&w1d[r],vsz*sizeof(real_t)));
    CHECK(cudaMemcpy(w1d[r],w1[r].data(),vsz*sizeof(real_t),cudaMemcpyHostToDevice));
    CHECK(cudaMalloc(&w2d[r],vsz*sizeof(real_t))); }
  int bc4[4]={1,1,1,-1};
  mult_coarse_mrhs(w2d.data(), w1d.data(), ud, cnrhs, Ncol, Ns, bc4);
  std::vector<std::vector<real_t>> w2(cnrhs,std::vector<real_t>(vsz));
  for(int r=0;r<cnrhs;r++) CHECK(cudaMemcpy(w2[r].data(),w2d[r],vsz*sizeof(real_t),cudaMemcpyDeviceToHost));
  // CPU oracle replicating accum_mult_u
  int Ngst1=Ndf*Nst;
  auto idx2=[&](int ld,int in,int st){ return (st%NWP)+NWP*(in+ld*(st/NWP)); };
  double cmax=0,cnrm=0;
  for(int site=0;site<Nst;site++){
    int ix=site%4,iy=(site/4)%4,iz=(site/16)%4,it=site/64;
    for(int a=0;a<Ncol;a++) for(int r=0;r<cnrhs;r++){
      double ar=w1[r][idx2(Ncol2,2*a,site)], ai=w1[r][idx2(Ncol2,2*a+1,site)];
      for(int d=0;d<8;d++){
        int uoff,nei; double bc2=1;
        switch(d){
          case 0: uoff=Ngst1*0; nei=(ix==3)?site-3:site+1; if(ix==3)bc2=bc4[0]; break;
          case 1: uoff=Ngst1*4; nei=(ix==0)?site+3:site-1; if(ix==0)bc2=bc4[0]; break;
          case 2: uoff=Ngst1*1; nei=(iy==3)?site-4*3:site+4; if(iy==3)bc2=bc4[1]; break;
          case 3: uoff=Ngst1*5; nei=(iy==0)?site+4*3:site-4; if(iy==0)bc2=bc4[1]; break;
          case 4: uoff=Ngst1*2; nei=(iz==3)?site-16*3:site+16; if(iz==3)bc2=bc4[2]; break;
          case 5: uoff=Ngst1*6; nei=(iz==0)?site+16*3:site-16; if(iz==0)bc2=bc4[2]; break;
          case 6: uoff=Ngst1*3; nei=(it==3)?site-64*3:site+64; if(it==3)bc2=bc4[3]; break;
          default:uoff=Ngst1*7; nei=(it==0)?site+64*3:site-64; if(it==0)bc2=bc4[3]; break;
        }
        for(int b=0;b<Ncol;b++){
          double mr=uh[uoff+idx2(Ndf,Ncol2*b+2*a,site)], mi=uh[uoff+idx2(Ndf,Ncol2*b+2*a+1,site)];
          double vr=w1[r][idx2(Ncol2,2*b,nei)], vi=w1[r][idx2(Ncol2,2*b+1,nei)];
          ar+=bc2*(mr*vr-mi*vi); ai+=bc2*(mr*vi+mi*vr);
        }
      }
      double dr=fabs(w2[r][idx2(Ncol2,2*a,site)]-ar), di=fabs(w2[r][idx2(Ncol2,2*a+1,site)]-ai);
      cmax=fmax(cmax,fmax(dr,di)); cnrm=fmax(cnrm,fmax(fabs(ar),fabs(ai)));
    }
  }
  double cne=cmax/cnrm;
  printf("\nCOARSE OP: max|ref|=%.4f  max abs diff=%.3e  normalized=%.3e\n",cnrm,cmax,cne);
  printf(cne<1e-3 ? "PASS: TF32 tensor-core coarse operator matches FP32 oracle\n"
                  : "FAIL: dataflow/index error\n");
  }

  // ---------- (C) block-Krylov: inner product + orthogonalization ----------
  {
  int nv2=8;                         // block size (N_rhs)
  std::vector<std::vector<real_t>> Vv(nv2,std::vector<real_t>(fsz)), Ww(nv2,std::vector<real_t>(fsz));
  for(int i=0;i<nv2;i++) for(size_t q=0;q<fsz;q++){ Vv[i][q]=rnd((int)(i*271+q+3)); Ww[i][q]=rnd((int)(i*733+q+91)); }
  // orthonormalize V on host (modified GS) so the orthogonality test is meaningful
  auto dot=[&](std::vector<real_t>&a,std::vector<real_t>&b){ double r=0,im=0;
    for(int in2=0;in2<Nin/2;in2++) for(int st=0;st<Nst;st++){
      double ar=a[IDX2(Nin,2*in2,st)],ai=a[IDX2(Nin,2*in2+1,st)],br=b[IDX2(Nin,2*in2,st)],bi=b[IDX2(Nin,2*in2+1,st)];
      r+=ar*br+ai*bi; im+=ar*bi-ai*br; } return std::make_pair(r,im); };
  for(int i=0;i<nv2;i++){
    for(int k=0;k<i;k++){ auto g=dot(Vv[k],Vv[i]);
      for(int in2=0;in2<Nin/2;in2++) for(int st=0;st<Nst;st++){
        double vr=Vv[k][IDX2(Nin,2*in2,st)],vi=Vv[k][IDX2(Nin,2*in2+1,st)];
        Vv[i][IDX2(Nin,2*in2,st)]-=g.first*vr-g.second*vi;
        Vv[i][IDX2(Nin,2*in2+1,st)]-=g.first*vi+g.second*vr; } }
    auto n=dot(Vv[i],Vv[i]); double s=1.0/sqrt(n.first);
    for(size_t q=0;q<fsz;q++) Vv[i][q]*=s;
  }
  std::vector<real_t*> Vd(nv2),Wd(nv2);
  for(int i=0;i<nv2;i++){ CHECK(cudaMalloc(&Vd[i],fsz*sizeof(real_t))); CHECK(cudaMalloc(&Wd[i],fsz*sizeof(real_t)));
    CHECK(cudaMemcpy(Vd[i],Vv[i].data(),fsz*sizeof(real_t),cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(Wd[i],Ww[i].data(),fsz*sizeof(real_t),cudaMemcpyHostToDevice)); }
  real_t* Gd; CHECK(cudaMalloc(&Gd,2*nv2*nv2*sizeof(real_t)));
  // 1) inner product correctness  G = V^H W
  block_inner_mrhs(Gd, Vd.data(), Wd.data(), nv2, Nin, Nex, Ns);
  std::vector<real_t> Gh(2*nv2*nv2); CHECK(cudaMemcpy(Gh.data(),Gd,2*nv2*nv2*sizeof(real_t),cudaMemcpyDeviceToHost));
  double gmax=0,gnrm=0;
  for(int i=0;i<nv2;i++)for(int j=0;j<nv2;j++){ auto g=dot(Vv[i],Ww[j]);
    gmax=fmax(gmax,fmax(fabs(Gh[2*(j+nv2*i)]-g.first),fabs(Gh[2*(j+nv2*i)+1]-g.second)));
    gnrm=fmax(gnrm,fmax(fabs(g.first),fabs(g.second))); }
  printf("\nBLOCK INNER  G=V^H W : max|ref|=%.4f normalized=%.3e %s\n",gnrm,gmax/gnrm,
         gmax/gnrm<2e-3?"PASS":"FAIL");
  // 2) orthogonalization: W -= V (V^H W); then ||V^H W'|| should be ~0
  block_update_mrhs(Wd.data(), Vd.data(), Gd, nv2, Nin, Nex, Ns);
  real_t* G2; CHECK(cudaMalloc(&G2,2*nv2*nv2*sizeof(real_t)));
  block_inner_mrhs(G2, Vd.data(), Wd.data(), nv2, Nin, Nex, Ns);
  std::vector<real_t> G2h(2*nv2*nv2); CHECK(cudaMemcpy(G2h.data(),G2,2*nv2*nv2*sizeof(real_t),cudaMemcpyDeviceToHost));
  double resid=0; for(size_t q=0;q<G2h.size();q++) resid=fmax(resid,fabs(G2h[q]));
  printf("BLOCK ARNOLDI step: ||V^H (W - V V^H W)||_max = %.3e %s\n",resid,
         resid<5e-3?"PASS (W now orthogonal to V)":"FAIL");
  }

  return 0;
}
