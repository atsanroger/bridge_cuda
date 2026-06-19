/*!
      @file    mult_Domainwall_coarse_cuda-inc.h
      @brief
      @author  Wei-Lun Chen (wlchen)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2024-04-01 12:13:09 #$
      @version $LastChangedRevision: 2595 $
*/

#ifndef MULT_DOMAINWALL_COARSE_CUDA_INCLUDED
#define MULT_DOMAINWALL_COARSE_CUDA_INCLUDED

namespace{

__device__
inline void accum_mult_u(real_t *v2, real_t *u2, real_t *v1,
                         int ivec, int site, int nei,
                         real_t bc2, int Ndf, int Ncol2, int NcolH)
{
  int i1 = 2 * ivec;     // chirality -
  int i2 = 2 * ivec + 1; // chirality +
  real_t vt1r = 0.0;
  real_t vt1i = 0.0;
  real_t vt2r = 0.0;
  real_t vt2i = 0.0;

  for(int j = 0; j < NcolH; ++j){
    int j1 = 2 * j;
    int j2 = 2 * j + 1;

    real_t wt1r = v1[IDX2(Ncol2, 2*j1,   nei)];
    real_t wt1i = v1[IDX2(Ncol2, 2*j1+1, nei)];
    real_t wt2r = v1[IDX2(Ncol2, 2*j2,   nei)];
    real_t wt2i = v1[IDX2(Ncol2, 2*j2+1, nei)];

    real_t ut11r = u2[IDX2(Ndf, Ncol2 * j1 + 2*i1,   site)];
    real_t ut11i = u2[IDX2(Ndf, Ncol2 * j1 + 2*i1+1, site)];
    real_t ut21r = u2[IDX2(Ndf, Ncol2 * j1 + 2*i2,   site)];
    real_t ut21i = u2[IDX2(Ndf, Ncol2 * j1 + 2*i2+1, site)];
    vt1r += ut11r * wt1r - ut11i * wt1i;
    vt1i += ut11r * wt1i + ut11i * wt1r;
    vt2r += ut21r * wt1r - ut21i * wt1i;
    vt2i += ut21r * wt1i + ut21i * wt1r;

    real_t ut12r = u2[IDX2(Ndf, Ncol2 * j2 + 2*i1,   site)];
    real_t ut12i = u2[IDX2(Ndf, Ncol2 * j2 + 2*i1+1, site)];
    real_t ut22r = u2[IDX2(Ndf, Ncol2 * j2 + 2*i2,   site)];
    real_t ut22i = u2[IDX2(Ndf, Ncol2 * j2 + 2*i2+1, site)];
    vt1r += ut12r * wt2r - ut12i * wt2i;
    vt1i += ut12r * wt2i + ut12i * wt2r;
    vt2r += ut22r * wt2r - ut22i * wt2i;
    vt2i += ut22r * wt2i + ut22i * wt2r;
  }

  v2[IDX2(Ncol2, 2*i1,   site)] += bc2 * vt1r;
  v2[IDX2(Ncol2, 2*i1+1, site)] += bc2 * vt1i;
  v2[IDX2(Ncol2, 2*i2,   site)] += bc2 * vt2r;
  v2[IDX2(Ncol2, 2*i2+1, site)] += bc2 * vt2i;
}

} // namespace

//====================================================================
__global__ void mult_domainwall_coarse_mult_gm5_dev(
                                     real_t* v2, real_t* v1,
                                     int Ncol, 
                                     int Nx, int Ny, int Nz, int Nt)
{

  int Nst   = Nx * Ny * Nz * Nt;
  int Ncst  = Ncol * Nst;
  int Nvc   = 2 * Ncol;
  int icst  = blockIdx.x * blockDim.x + threadIdx.x;

  if( icst < Ncst ){
    int ic   = icst % Ncol;
    int site = icst / Ncol;
    int ch = ic % 2;
    if(ch == 0){  // negative chirality
      v2[IDX2(Nvc, 2*ic,   site)] = -v1[IDX2(Nvc, 2*ic,   site)];
      v2[IDX2(Nvc, 2*ic+1, site)] = -v1[IDX2(Nvc, 2*ic+1, site)];
    }else{        // positive chirality
      v2[IDX2(Nvc, 2*ic,   site)] =  v1[IDX2(Nvc, 2*ic,   site)];
      v2[IDX2(Nvc, 2*ic+1, site)] =  v1[IDX2(Nvc, 2*ic+1, site)];
    }
  }
}

void mult_domainwall_coarse_mult_gm5(real_t* v2,
                                     real_t* v1,
                                     int Ncol, int *Nsize)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst  = Nx * Ny * Nz * Nt;
  int Ncst = Ncol * Nst;

  real_t * v2_dev = (real_t*)dev_ptr(v2);
  real_t * v1_dev = (real_t*)dev_ptr(v1);

  int blockSize  = VECTOR_LENGTH;
  int gridSize   = (Ncst + blockSize - 1)/ blockSize;

  mult_domainwall_coarse_mult_gm5_dev<<<gridSize, blockSize>>>(
                                      v2_dev, v1_dev,
                                      Ncol, 
                                      Nx,  Ny,  Nz,  Nt);
  CHECK(cudaDeviceSynchronize());
}

//====================================================================
__global__ void mult_domainwall_coarse_bulk_dev(
                             real_t* v2, real_t* u, real_t* v1,
                             int Ncol, 
                             int Nx, int Ny, int Nz, int Nt, 
                             int bc_x, int bc_y, int bc_z, int bc_t)
{
  int Nst  = Nx * Ny * Nz * Nt;
  //int Nvst = 2 * Ncol * Nst;
  //int Ngst = 2 * Ncol * Ncol * Nst * 8;

  int Ncol2 = 2 * Ncol;
  int NcolH = Ncol/2;
  int Ndf   = 2 * Ncol * Ncol;
  int Ngst1 = 2 * Ncol * Ncol * Nst;
  int Nxy   = Nx * Ny;
  int Nxyz  = Nx * Ny * Nz;

  int Nloop = NcolH * Nst;
  int iloop = blockIdx.x * blockDim.x + threadIdx.x;
  //extern __shared__ real_t sharedMemory[];
  //real_t* u = sharedMemory;

  if( iloop < Nloop ){
      int ivec = iloop % NcolH;
      int site = iloop / NcolH;
      int ix = site % Nx;
      int iy = (site/Nx) % Ny;
      int iz = (site/Nxy) % Nz;
      int it = site/Nxyz;

      //u[iloop] = u_dev[Ngst1 * 8];

      for(int k = 0; k < 4; ++k){
        v2[IDX2(Ncol2, 4*ivec+k, site)] = v1[IDX2(Ncol2, 4*ivec+k, site)];
      }

      {
      int idir = 0;
      real_t *u2 = &u[Ngst1 * idir];

      int nei = site + 1;
      real_t bc2 = 1.0;
      if(ix == Nx-1){
        nei = site - (Nx-1);
        bc2 = bc_x;
      }
      accum_mult_u(v2, u2, v1, ivec, site, nei, bc2, Ndf, Ncol2, NcolH);

      u2 = &u[Ngst1 * (idir + 4)];
      nei = site - 1;
      bc2 = 1.0;
      if(ix == 0){
        nei = site + Nx-1;
        bc2 = bc_x;
      }
      accum_mult_u(v2, u2, v1, ivec, site, nei, bc2, Ndf, Ncol2, NcolH);
      }

     {
      int idir = 1;
      real_t *u2 = &u[Ngst1 * idir];

      int nei = site + Nx;
      real_t bc2 = 1.0;
      if(iy == Ny-1){
        nei = site - Nx*(Ny-1);
        bc2 = bc_y;
      }
      accum_mult_u(v2, u2, v1, ivec, site, nei, bc2, Ndf, Ncol2, NcolH);

      u2 = &u[Ngst1 * (idir + 4)];
      nei = site - Nx;
      bc2 = 1.0;
      if(iy == 0){
        nei = site + Nx*(Ny-1);
        bc2 = bc_y;
      }
      accum_mult_u(v2, u2, v1, ivec, site, nei, bc2, Ndf, Ncol2, NcolH);
     }

     {
      int idir = 2;
      real_t *u2 = &u[Ngst1 * idir];

      int nei = site + Nxy;
      real_t bc2 = 1.0;
      if(iz == Nz-1){
        nei = site - Nxy*(Nz-1);
        bc2 = bc_z;
      }
      accum_mult_u(v2, u2, v1, ivec, site, nei, bc2, Ndf, Ncol2, NcolH);

      u2 = &u[Ngst1 * (idir + 4)];
      nei = site - Nxy;
      bc2 = 1.0;
      if(iz == 0){
        nei = site + Nxy*(Nz-1);
        bc2 = bc_z;
      }
      accum_mult_u(v2, u2, v1, ivec, site, nei, bc2, Ndf, Ncol2, NcolH);
     }

     {
      int idir = 3;
      real_t *u2 = &u[Ngst1 * idir];

      int nei = site + Nxyz;
      real_t bc2 = 1.0;
      if(it == Nt-1){
        nei = site - Nxyz*(Nt-1);
        bc2 = bc_t;
      }

      accum_mult_u(v2, u2, v1, ivec, site, nei, bc2, Ndf, Ncol2, NcolH);

      u2 = &u[Ngst1 * (idir + 4)];
      nei = site - Nxyz;
      bc2 = 1.0;
      if(it == 0){
        nei = site + Nxyz*(Nt-1);
        bc2 = bc_t;
      }
      accum_mult_u(v2, u2, v1, ivec, site, nei, bc2, Ndf, Ncol2, NcolH);
     }
    }
}

void mult_domainwall_coarse_bulk(
                             real_t* v2, real_t* u,
                             real_t* v1,
                             int Ncol, int *Nsize, int *bc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];

  int Nst   = Nx * Ny * Nz * Nt;
  int Nloop = Nst * (Ncol/2);

  int bc_x = bc[0];
  int bc_y = bc[1];
  int bc_z = bc[2];
  int bc_t = bc[3];

  real_t * v2_dev = (real_t*)dev_ptr(v2);
  real_t * u_dev  = (real_t*)dev_ptr(u);
  real_t * v1_dev = (real_t*)dev_ptr(v1);

  int blockSize = VECTOR_LENGTH;
  int gridSize  = (Nloop + blockSize - 1)/ blockSize;
  //int sharedsize = (4 * Ncol * Ncol) * blockSize * sizeof(real_t);

  mult_domainwall_coarse_bulk_dev<<<gridSize, blockSize>>>(
                            v2_dev, u_dev, v1_dev,
                            Ncol, 
                            Nx, Ny, Nz, Nt, 
                            bc_x, bc_y, bc_z, bc_t);

  CHECK(cudaDeviceSynchronize());
}

//====================================================================
__global__ void mult_domainwall_coarse_prec_dev(
                             real_t* v2,
                             real_t* ct, real_t* v1,
                             int Ncol, 
                             int Nx, int Ny, int Nz, int Nt)
{
  int Nst   = Nx * Ny * Nz * Nt;
  //int Nvst  = 2 * Ncol * Nst;
  //int Nctst = 2 * Ncol * Ncol * Nst;

  int Ncol2 = 2 * Ncol;
  int NcolH = Ncol/2;
  int Ndf   = 2 * Ncol * Ncol;
  int Nloop = NcolH * Nst;
  int iloop = blockIdx.x * blockDim.x + threadIdx.x;

    if(iloop < Nloop){
      int ivec = iloop % NcolH;
      int site = iloop / NcolH;

      v2[IDX2(Ncol2, 4*ivec  , site)] = 0.0;
      v2[IDX2(Ncol2, 4*ivec+1, site)] = 0.0;
      v2[IDX2(Ncol2, 4*ivec+2, site)] = 0.0;
      v2[IDX2(Ncol2, 4*ivec+3, site)] = 0.0;

      accum_mult_u(v2, ct, v1, ivec, site, site, real_t(1.0),
                   Ndf, Ncol2, NcolH);

    }
}
void mult_domainwall_coarse_prec(
                             real_t* v2,
                             real_t* ct, real_t* v1,
                             int Ncol, int *Nsize)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];

  int Nst   = Nx * Ny * Nz * Nt;
  int Nloop = Nst * (Ncol / 2);

  real_t * v2_dev = (real_t*)dev_ptr(v2);
  real_t * ct_dev = (real_t*)dev_ptr(ct);
  real_t * v1_dev = (real_t*)dev_ptr(v1);

  int blockSize = VECTOR_LENGTH;
  int gridSize  = (Nloop + blockSize - 1)/ blockSize;

  mult_domainwall_coarse_prec_dev<<<gridSize, blockSize>>>(
                            v2_dev, ct_dev, v1_dev,
                            Ncol, 
                            Nx, Ny, Nz, Nt);

  CHECK(cudaDeviceSynchronize());

}
//====================================================================

__global__
void mult_domainwall_coarse_xpb_1_dev(
                            real_t* buf_xp, real_t* v1, 
                            int Ncol, int Nx, int Ny, int Nz, int Nt)
{
  int Nyzt  = Ny * Nz * Nt;
  int Ncol2 = 2 * Ncol; 

  int iyzt = blockIdx.x * blockDim.x + threadIdx.x;
  if(iyzt < Nyzt){
    int ix = 0;
    int ist  = ix + Nx * iyzt;
    for(int ic2 = 0; ic2 < Ncol2; ++ic2){
      buf_xp[IDX2(Ncol2, ic2, iyzt)] = v1[IDX2(Ncol2, ic2, ist)];
    }
  }
}

__global__
void mult_domainwall_coarse_xmb_1_dev(
                            real_t* buf_xm, real_t* v1,
                            int Ncol, int Nx, int Ny, int Nz, int Nt)
{
  int Nyzt  = Ny * Nz * Nt;
  int Ncol2 = 2 * Ncol; 

  int iyzt = blockIdx.x * blockDim.x + threadIdx.x;
  if(iyzt < Nyzt){
    int ix = Nx-1;
    int ist  = ix + Nx * iyzt;
    for(int ic2 = 0; ic2 < Ncol2; ++ic2){
      buf_xm[IDX2(Ncol2, ic2, iyzt)] = v1[IDX2(Ncol2, ic2, ist)];
    }
  }

}

__global__
void mult_domainwall_coarse_ypb_1_dev(
                            real_t* buf_yp, real_t* v1, 
                            int Ncol, int Nx, int Ny, int Nz, int Nt)
{
  int Nxzt  = Nx * Nz * Nt;
  int Ncol2 = 2 * Ncol; 

  int ixzt = blockIdx.x * blockDim.x + threadIdx.x;

    if(ixzt < Nxzt){
      int ix  = ixzt % Nx;
      int izt = ixzt / Nx;
      int iy  = 0;
      int ist = ix + Nx * (iy + Ny * izt);
      for(int ic2 = 0; ic2 < Ncol2; ++ic2){
        buf_yp[IDX2(Ncol2, ic2, ixzt)] = v1[IDX2(Ncol2, ic2, ist)];
      }
    }
}

__global__
void mult_domainwall_coarse_ymb_1_dev(
                            real_t* buf_ym, real_t* v1, 
                            int Ncol, int Nx, int Ny, int Nz, int Nt)
{
  int Nxzt  = Nx * Nz * Nt;
  int Ncol2 = 2 * Ncol; 

  int ixzt = blockIdx.x * blockDim.x + threadIdx.x;

  if(ixzt < Nxzt){
    int ix  = ixzt % Nx;
    int izt = ixzt / Nx;
    int iy  = Ny-1;
    int ist = ix + Nx * (iy + Ny * izt);
    for(int ic2 = 0; ic2 < Ncol2; ++ic2){
      buf_ym[IDX2(Ncol2, ic2, ixzt)] = v1[IDX2(Ncol2, ic2, ist)];
    }
  }
}

__global__
void mult_domainwall_coarse_zpb_1_dev(
                            real_t* buf_zp, real_t* v1, 
                            int Ncol, int Nx, int Ny, int Nz, int Nt)
{
  int Nxy   = Nx  * Ny;
  int Nxyt  = Nxy * Nt;
  int Ncol2 = 2 * Ncol; 

  int ixyt  = blockIdx.x * blockDim.x + threadIdx.x;

  if( ixyt < Nxyt ){
    int iz    = 0;
    int ixy = ixyt % Nxy;
    int it  = ixyt / Nxy;
    int ist  = ixy + Nxy *(iz + Nz*it);
    for(int ic2 = 0; ic2 < Ncol2; ++ic2){
      buf_zp[IDX2(Ncol2, ic2, ixyt)] = v1[IDX2(Ncol2, ic2, ist)];
    }
  }
}

__global__
void mult_domainwall_coarse_zmb_1_dev(
                            real_t* buf_zm, real_t* v1, 
                            int Ncol, int Nx, int Ny, int Nz, int Nt)
{
  int Nxy   = Nx  * Ny;
  int Nxyt  = Nxy * Nt;
  int Ncol2 = 2 * Ncol; 

  int ixyt  = blockIdx.x * blockDim.x + threadIdx.x;

  if( ixyt < Nxyt ){
    int iz    = Nz-1;
    int ixy = ixyt % Nxy;
    int it  = ixyt / Nxy;
    int ist  = ixy + Nxy *(iz + Nz*it);
    for(int ic2 = 0; ic2 < Ncol2; ++ic2){
      buf_zm[IDX2(Ncol2, ic2, ixyt)] = v1[IDX2(Ncol2, ic2, ist)];
    }
  }
}

__global__
void mult_domainwall_coarse_tpb_1_dev(
                            real_t* buf_tp, real_t* v1, 
                            int Ncol, int Nx, int Ny, int Nz, int Nt)
{
  int Nxyz  = Nx * Ny * Nz;
  int Ncol2 = 2 * Ncol; 

  int ixyz  = blockIdx.x * blockDim.x + threadIdx.x;

  if( ixyz < Nxyz){
    int it   = 0;
    int ist  = ixyz + Nxyz * it;
    for(int ic2 = 0; ic2 < Ncol2; ++ic2){
      buf_tp[IDX2(Ncol2, ic2, ixyz)] = v1[IDX2(Ncol2, ic2, ist)];
    }
  }
}

__global__
void mult_domainwall_coarse_tmb_1_dev(
                            real_t* buf_tm, real_t* v1, 
                            int Ncol, int Nx, int Ny, int Nz, int Nt)
{
  int Nxyz  = Nx * Ny * Nz;
  int Ncol2 = 2 * Ncol; 

  int ixyz  = blockIdx.x * blockDim.x + threadIdx.x;

  if( ixyz < Nxyz){
    int it   = Nt-1;
    int ist  = ixyz + Nxyz * it;
    for(int ic2 = 0; ic2 < Ncol2; ++ic2){
      buf_tm[IDX2(Ncol2, ic2, ixyz)] = v1[IDX2(Ncol2, ic2, ist)];
    }
  }
}

void mult_domainwall_coarse_1(real_t* buf_xp, real_t* buf_xm,
                              real_t* buf_yp, real_t* buf_ym,
                              real_t* buf_zp, real_t* buf_zm,
                              real_t* buf_tp, real_t* buf_tm,
                              real_t* v1,
                              int Ncol, int *Nsize, int *do_comm)
{
  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];

  int Nxyz  = Nx * Ny * Nz;
  int Nyzt  = Ny * Nz * Nt;
  int Nxyt  = Nx * Ny * Nt;
  int Nxzt  = Nx * Nz * Nt;

  int do_comm_x = do_comm[0];
  int do_comm_y = do_comm[1];
  int do_comm_z = do_comm[2];
  int do_comm_t = do_comm[3];

  int size_bx = 2 * Ncol * Ny * Nz * Nt;
  int size_by = 2 * Ncol * Nx * Nz * Nt;
  int size_bz = 2 * Ncol * Nx * Ny * Nt;
  int size_bt = 2 * Ncol * Nx * Ny * Nz;

  real_t * v1_dev     = (real_t*)dev_ptr(v1);

  int blockSize  = VECTOR_LENGTH;
  //int sharedSize = NDF * blockSize * sizeof(real_t); 

  if(do_comm_x > 0) {
    int gridSize   = (Nyzt + blockSize - 1)/ blockSize;
    mult_domainwall_coarse_xpb_1_dev<<<gridSize, blockSize>>>
                          (buf_xp, v1_dev, Ncol, Nx, Ny, Nz, Nt);

    mult_domainwall_coarse_xmb_1_dev<<<gridSize, blockSize>>>
                          (buf_xm, v1_dev, Ncol, Nx, Ny, Nz, Nt);
  }

  if(do_comm_y > 0) {
    int gridSize   = (Nxzt + blockSize - 1)/ blockSize;
    mult_domainwall_coarse_ypb_1_dev<<<gridSize, blockSize>>>
                          (buf_yp, v1_dev, Ncol, Nx, Ny, Nz, Nt);

    mult_domainwall_coarse_ymb_1_dev<<<gridSize, blockSize>>>
                          (buf_ym, v1_dev, Ncol, Nx, Ny, Nz, Nt);
  }

  if(do_comm_z > 0) {
    int gridSize   = (Nxyt + blockSize - 1)/ blockSize;
    mult_domainwall_coarse_zpb_1_dev<<<gridSize, blockSize>>>
                          (buf_zp, v1_dev, Ncol, Nx, Ny, Nz, Nt);

    mult_domainwall_coarse_zmb_1_dev<<<gridSize, blockSize>>>
                          (buf_zm, v1_dev, Ncol, Nx, Ny, Nz, Nt);
  }

  if(do_comm_t > 0) {
    int gridSize   = (Nxyz + blockSize - 1)/ blockSize;
     mult_domainwall_coarse_tpb_1_dev<<<gridSize, blockSize>>>
                          (buf_tp, v1_dev, Ncol, Nx, Ny, Nz, Nt);

    mult_domainwall_coarse_tmb_1_dev<<<gridSize, blockSize>>>
                          (buf_tm, v1_dev, Ncol, Nx, Ny, Nz, Nt);
  }

  CHECK(cudaDeviceSynchronize());
}

//====================================================================
__global__
void mult_domainwall_coarse_2_dev(real_t* v2, real_t* u, 
                                  real_t* buf_xp, real_t* buf_xm,
                                  real_t* buf_yp, real_t* buf_ym,
                                  real_t* buf_zp, real_t* buf_zm,
                                  real_t* buf_tp, real_t* buf_tm,
                                  int Ncol,
                                  int Nx, int Ny, int Nz, int Nt,
                                  int do_comm_x, int do_comm_y,
                                  int do_comm_z, int do_comm_t)
{

  int Nst   = Nx * Ny * Nz * Nt;
  int Ncol2 = 2 * Ncol;
  int NcolH = Ncol/2;
  int Ndf   = 2 * Ncol * Ncol;

  int Nloop = Nst * NcolH;
  int iloop = blockIdx.x * blockDim.x + threadIdx.x;

  if( iloop >= Nloop ) return;

  int ivec = iloop % NcolH;
  int site = iloop / NcolH;

  int ix = site % Nx;
  int iy = (site/Nx) % Ny;
  int iz = (site/(Nx*Ny)) % Nz;
  int it = site/(Nx*Ny*Nz);

  int idir = 0;
  if(do_comm_x > 0){
    if(ix == Nx-1){
      real_t *u2 = &u[Ndf * Nst * idir];
      int iyzt = site/Nx;
      accum_mult_u(v2, u2, buf_xp, ivec, site, iyzt, real_t(1.0),
                   Ndf, Ncol2, NcolH);
    }
    if(ix == 0){
      real_t *u2 = &u[Ndf * Nst * (idir+4)];
      int iyzt = site/Nx;
      accum_mult_u(v2, u2, buf_xm, ivec, site, iyzt, real_t(1.0),
                   Ndf, Ncol2, NcolH);
    }
  }

  idir = 1;
  if(do_comm_y > 0){
    if(iy == Ny-1){
      real_t *u2 = &u[Ndf * Nst * idir];
      int ixzt = ix + Nx * (iz + Nz*it);
      accum_mult_u(v2, u2, buf_yp, ivec, site, ixzt, real_t(1.0),
                   Ndf, Ncol2, NcolH);
    }
    if(iy == 0){
      real_t *u2 = &u[Ndf * Nst * (idir+4)];
      int ixzt = ix + Nx * (iz + Nz*it);
      accum_mult_u(v2, u2, buf_ym, ivec, site, ixzt, real_t(1.0),
                   Ndf, Ncol2, NcolH);
    }
  }

  idir = 2;
  if(do_comm_z > 0){
    if(iz == Nz-1){
      real_t *u2 = &u[Ndf * Nst * idir];
      int ixyt = ix + Nx * (iy + Ny * it);
      accum_mult_u(v2, u2, buf_zp, ivec, site, ixyt, real_t(1.0),
                   Ndf, Ncol2, NcolH);
    }
    if(iz == 0){
      real_t *u2 = &u[Ndf * Nst * (idir+4)];
      int ixyt = ix + Nx * (iy + Ny * it);
      accum_mult_u(v2, u2, buf_zm, ivec, site, ixyt, real_t(1.0),
                   Ndf, Ncol2, NcolH);
    }
  }

  idir = 3;
  if(do_comm_t > 0){
    if(it == Nt-1){
      real_t *u2 = &u[Ndf * Nst * idir];
      int ixyz = ix + Nx * (iy + Ny * iz);
      accum_mult_u(v2, u2, buf_tp, ivec, site, ixyz, real_t(1.0),
                   Ndf, Ncol2, NcolH);
    }
    if(it == 0){
      real_t *u2 = &u[Ndf * Nst * (idir+4)];
      int ixyz = ix + Nx * (iy + Ny * iz);
      accum_mult_u(v2, u2, buf_tm, ivec, site, ixyz, real_t(1.0),
                   Ndf, Ncol2, NcolH);
    }
  }

}

void mult_domainwall_coarse_2(real_t* v2, real_t* u, 
                              real_t* buf_xp, real_t* buf_xm,
                              real_t* buf_yp, real_t* buf_ym,
                              real_t* buf_zp, real_t* buf_zm,
                              real_t* buf_tp, real_t* buf_tm,
                              int Ncol, int *Nsize, int *do_comm)
{
  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nst  = Nx * Ny * Nz * Nt;

  int NcolH = Ncol/2;

  int do_comm_x = do_comm[0];
  int do_comm_y = do_comm[1];
  int do_comm_z = do_comm[2];
  int do_comm_t = do_comm[3];

  int size_bx = 2 * Ncol * Ny * Nz * Nt;
  int size_by = 2 * Ncol * Nx * Nz * Nt;
  int size_bz = 2 * Ncol * Nx * Ny * Nt;
  int size_bt = 2 * Ncol * Nx * Ny * Nz;

  real_t * v2_dev     = (real_t*)dev_ptr(v2);
  real_t * u_dev      = (real_t*)dev_ptr(u);

  int blockSize = VECTOR_LENGTH;
  int gridSize = (NcolH * Nst + blockSize - 1)/ blockSize;
  mult_domainwall_coarse_2_dev<<<gridSize, blockSize>>>(
                    v2_dev,  u_dev, 
                    buf_xp, buf_xm,
                    buf_yp, buf_ym,
                    buf_zp, buf_zm,
                    buf_tp, buf_tm,
                    Ncol,
                    Nx, Ny, Nz, Nt,
                    do_comm_x, do_comm_y, do_comm_z, do_comm_t);

}

#endif
//============================================================END=====
