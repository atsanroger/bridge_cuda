/*!
      @file    mult_Clover_coarse_cuda-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-08-20 14:25:12 #$
      @version $LastChangedRevision: 2535 $
*/

namespace{

__device__
inline void accum_mult_u(real_t *v2, real_t *u2, real_t *v1,
                         int ivec, int site, int nei, real_t bc2,
                         int Ndf, int Ncol2, int NcolH)
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

__device__
inline void accum_mult_udag(real_t *v2, real_t *u2, real_t *v1,
                            int ivec, int site, int nei, real_t bc2,
                            int Ndf, int Ncol2, int NcolH)
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

    real_t ut11r = u2[IDX2(Ndf, Ncol2 * i1 + 2*j1,   nei)];
    real_t ut11i = u2[IDX2(Ndf, Ncol2 * i1 + 2*j1+1, nei)];
    real_t ut21r = u2[IDX2(Ndf, Ncol2 * i2 + 2*j1,   nei)];
    real_t ut21i = u2[IDX2(Ndf, Ncol2 * i2 + 2*j1+1, nei)];
    vt1r += ut11r * wt1r + ut11i * wt1i;
    vt1i += ut11r * wt1i - ut11i * wt1r;
    vt2r -= ut21r * wt1r + ut21i * wt1i;
    vt2i -= ut21r * wt1i - ut21i * wt1r;

    real_t ut12r = u2[IDX2(Ndf, Ncol2 * i1 + 2*j2,   nei)];
    real_t ut12i = u2[IDX2(Ndf, Ncol2 * i1 + 2*j2+1, nei)];
    real_t ut22r = u2[IDX2(Ndf, Ncol2 * i2 + 2*j2,   nei)];
    real_t ut22i = u2[IDX2(Ndf, Ncol2 * i2 + 2*j2+1, nei)];
    vt1r -= ut12r * wt2r + ut12i * wt2i;
    vt1i -= ut12r * wt2i - ut12i * wt2r;
    vt2r += ut22r * wt2r + ut22i * wt2i;
    vt2i += ut22r * wt2i - ut22i * wt2r;
  }

  v2[IDX2(Ncol2, 2*i1,   site)] += bc2 * vt1r;
  v2[IDX2(Ncol2, 2*i1+1, site)] += bc2 * vt1i;
  v2[IDX2(Ncol2, 2*i2,   site)] += bc2 * vt2r;
  v2[IDX2(Ncol2, 2*i2+1, site)] += bc2 * vt2i;

}

} // namespace

//====================================================================
__global__
void mult_clover_coarse_mult_gm5_kernel(real_t* v2, real_t* v1,
                           int Ncol, int Nx, int Ny, int Nz, int Nt)
{
  int Nvc = 2 * Ncol;

  int icst = blockIdx.x * blockDim.x + threadIdx.x;
  //  for(int icst = 0; icst < Ncst; ++icst){

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

//====================================================================
void mult_clover_coarse_mult_gm5(real_t* v2, real_t* v1,
                                 int Ncol, int *Nsize)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst  = Nx * Ny * Nz * Nt;

  int Ncst = Ncol * Nst;

  real_t* v1_dev  = (real_t*)dev_ptr(v1);
  real_t* v2_dev  = (real_t*)dev_ptr(v2);

  int nth = VECTOR_LENGTH;
  int nbl = Ncst/nth;

  mult_clover_coarse_mult_gm5_kernel<<<nbl,nth>>>(
                                             v2_dev, v1_dev,
                                             Ncol, Nx, Ny, Nz, Nt);

  CHECK(cudaDeviceSynchronize());

}

//====================================================================
__global__
void mult_clover_coarse_bulk_kernel(
                            real_t* v2, real_t* u,
                            real_t* ct, real_t* v1,
                            int Ncol, int Nx, int Ny, int Nz, int Nt,
                            int bc_x, int bc_y, int bc_z, int bc_t)
{
  int Nst  = Nx * Ny * Nz * Nt;

  int Ncol2 = 2 * Ncol;
  int NcolH = Ncol/2;
  int Ndf   = 2 * Ncol * Ncol;
  int Ngst1 = 2 *  Ncol * Ncol * Nst;
  int Nxy   = Nx * Ny;
  int Nxyz  = Nx * Ny * Nz;

  int iloop = blockIdx.x * blockDim.x + threadIdx.x;
  //  int Nloop = NcolH * Nst;
  // for(int iloop = 0; iloop < Nloop; ++iloop){

  int ivec = iloop % NcolH;
  int site = iloop/NcolH;
  int ix = site % Nx;
  int iy = (site/Nx) % Ny;
  int iz = (site/Nxy) % Nz;
  int it = site/Nxyz;

  v2[IDX2(Ncol2, 4*ivec  , site)] = 0.0;
  v2[IDX2(Ncol2, 4*ivec+1, site)] = 0.0;
  v2[IDX2(Ncol2, 4*ivec+2, site)] = 0.0;
  v2[IDX2(Ncol2, 4*ivec+3, site)] = 0.0;

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

  nei = site - 1;
  bc2 = 1.0;
  if(ix == 0){
    nei = site + Nx-1;
    bc2 = bc_x;
  }
  accum_mult_udag(v2, u2, v1, ivec, site, nei, bc2, Ndf, Ncol2, NcolH);
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

  nei = site - Nx;
  bc2 = 1.0;
  if(iy == 0){
    nei = site + Nx*(Ny-1);
    bc2 = bc_y;
  }
  accum_mult_udag(v2, u2, v1, ivec, site, nei, bc2, Ndf, Ncol2, NcolH);
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

  nei = site - Nxy;
  bc2 = 1.0;
  if(iz == 0){
    nei = site + Nxy*(Nz-1);
    bc2 = bc_z;
  }
  accum_mult_udag(v2, u2, v1, ivec, site, nei, bc2, Ndf, Ncol2, NcolH);
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

  nei = site - Nxyz;
  bc2 = 1.0;
  if(it == 0){
    nei = site + Nxyz*(Nt-1);
    bc2 = bc_t;
  }
  accum_mult_udag(v2, u2, v1, ivec, site, nei, bc2, Ndf, Ncol2, NcolH);
 }

  accum_mult_u(v2, ct, v1, ivec, site, site, real_t(1.0),
               Ndf, Ncol2, NcolH);

}

//====================================================================
void mult_clover_coarse_bulk(real_t* v2, real_t* u,
                             real_t* ct, real_t* v1,
                             int Ncol, int *Nsize, int *bc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst  = Nx * Ny * Nz * Nt;

  int bc_x = bc[0];
  int bc_y = bc[1];
  int bc_z = bc[2];
  int bc_t = bc[3];

  real_t* v1_dev  = (real_t*)dev_ptr(v1);
  real_t* v2_dev  = (real_t*)dev_ptr(v2);
  real_t* u_dev   = (real_t*)dev_ptr(u);
  real_t* ct_dev  = (real_t*)dev_ptr(ct);

  int NcolH = Ncol/2;
  int Nloop = NcolH * Nst;

  int nth = VECTOR_LENGTH;
  int nbl = Nloop/nth;

  mult_clover_coarse_bulk_kernel<<<nbl,nth>>>(v2_dev, u_dev,
                                              ct_dev, v1_dev,
                                              Ncol, Nx, Ny, Nz, Nt,
					      bc_x, bc_y, bc_z, bc_t);

  CHECK(cudaDeviceSynchronize());

}

//====================================================================
__global__
void mult_clover_coarse_1x_kernel(
                            real_t* buf_xp, real_t* buf_xm,
                            real_t* u, real_t* v1,
                            int Ncol, int Nx, int Ny, int Nz, int Nt)
{
  int Nst  = Nx * Ny * Nz * Nt;

  int Ncol2 = 2 * Ncol; 
  int NcolH = Ncol/2;
  int Ndf   = 2 * Ncol * Ncol;

  int idir = 0;

  int iyzt = blockIdx.x * blockDim.x + threadIdx.x;

 {
  int ix = 0;
  int ist  = ix + Nx * iyzt;
  for(int ic2 = 0; ic2 < Ncol2; ++ic2){
    buf_xp[IDX2(Ncol2, ic2, iyzt)] = v1[IDX2(Ncol2, ic2, ist)];
  }
 }

 {
  int ix = Nx-1;
  real_t *u2 = &u[Ndf * Nst * idir];

  for(int ivec = 0; ivec < NcolH; ++ivec){
    int ist  = ix + Nx * iyzt;

    buf_xm[IDX2(Ncol2, 4*ivec  , iyzt)] = 0.0;
    buf_xm[IDX2(Ncol2, 4*ivec+1, iyzt)] = 0.0;
    buf_xm[IDX2(Ncol2, 4*ivec+2, iyzt)] = 0.0;
    buf_xm[IDX2(Ncol2, 4*ivec+3, iyzt)] = 0.0;
    accum_mult_udag(buf_xm, u2, v1, ivec, iyzt, ist, real_t(1.0),
                    Ndf, Ncol2, NcolH);
  }
 }

}

//====================================================================
__global__
void mult_clover_coarse_1y_kernel(
                            real_t* buf_yp, real_t* buf_ym,
                            real_t* u, real_t* v1,
                            int Ncol, int Nx, int Ny, int Nz, int Nt)
{
  int Nst  = Nx * Ny * Nz * Nt;
  int Nxzt = Nx * Nz * Nt;

  int Ncol2 = 2 * Ncol; 
  int NcolH = Ncol/2;
  int Ndf   = 2 * Ncol * Ncol;

  int idir = 1;

  int ixzt = blockIdx.x * blockDim.x + threadIdx.x;

 {
  int iy = 0;
  for(int ixzt = 0; ixzt < Nxzt; ++ixzt){
    int ix  = ixzt % Nx;
    int izt = ixzt / Nx;
    int ist  = ix + Nx * (iy + Ny * izt);
    for(int ic2 = 0; ic2 < Ncol2; ++ic2){
      buf_yp[IDX2(Ncol2, ic2, ixzt)] = v1[IDX2(Ncol2, ic2, ist)];
    }
  }
 }

 {
  int iy = Ny-1;

  real_t *u2 = &u[Ndf * Nst * idir];

  for(int ivec = 0; ivec < NcolH; ++ivec){
    int ix  = ixzt % Nx;
    int izt = ixzt / Nx;
    int ist  = ix + Nx * (iy + Ny * izt);

    buf_ym[IDX2(Ncol2, 4*ivec  , ixzt)] = 0.0;
    buf_ym[IDX2(Ncol2, 4*ivec+1, ixzt)] = 0.0;
    buf_ym[IDX2(Ncol2, 4*ivec+2, ixzt)] = 0.0;
    buf_ym[IDX2(Ncol2, 4*ivec+3, ixzt)] = 0.0;
    accum_mult_udag(buf_ym, u2, v1, ivec, ixzt, ist, real_t(1.0),
                    Ndf, Ncol2, NcolH);
  }
 }

}

//====================================================================
__global__
void mult_clover_coarse_1z_kernel(
                            real_t* buf_zp, real_t* buf_zm,
                            real_t* u, real_t* v1,
                            int Ncol, int Nx, int Ny, int Nz, int Nt)
{
  int Nst  = Nx * Ny * Nz * Nt;
  int Nxy  = Nx * Ny;

  int Ncol2 = 2 * Ncol; 
  int NcolH = Ncol/2;
  int Ndf   = 2 * Ncol * Ncol;

  int ixyt = blockIdx.x * blockDim.x + threadIdx.x;

  int  idir = 2;

 {
  int iz = 0;
  int ixy = ixyt % Nxy;
  int it  = ixyt/Nxy;
  int ist  = ixy + Nxy *(iz + Nz*it);
  for(int ic2 = 0; ic2 < Ncol2; ++ic2){
    buf_zp[IDX2(Ncol2, ic2, ixyt)] = v1[IDX2(Ncol2, ic2, ist)];
  }
 }

 {
  int iz = Nz-1;

  real_t *u2 = &u[Ndf * Nst * idir];

  for(int ivec = 0; ivec < NcolH; ++ivec){
    int ixy = ixyt % Nxy;
    int it  = ixyt/Nxy;
    int ist  = ixy + Nxy *(iz + Nz*it);

    buf_zm[IDX2(Ncol2, 4*ivec  , ixyt)] = 0.0;
    buf_zm[IDX2(Ncol2, 4*ivec+1, ixyt)] = 0.0;
    buf_zm[IDX2(Ncol2, 4*ivec+2, ixyt)] = 0.0;
    buf_zm[IDX2(Ncol2, 4*ivec+3, ixyt)] = 0.0;
    accum_mult_udag(buf_zm, u2, v1, ivec, ixyt, ist, real_t(1.0),
                    Ndf, Ncol2, NcolH);
  }
 }

}

//====================================================================
__global__
void mult_clover_coarse_1t_kernel(
                            real_t* buf_tp, real_t* buf_tm,
                            real_t* u, real_t* v1,
                            int Ncol, int Nx, int Ny, int Nz, int Nt)
{
  int Nst  = Nx * Ny * Nz * Nt;
  int Nxyz = Nx * Ny * Nz;

  int Ncol2 = 2 * Ncol; 
  int NcolH = Ncol/2;
  int Ndf   = 2 * Ncol * Ncol;

  int ixyz = blockIdx.x * blockDim.x + threadIdx.x;

  int idir = 3;

 {
  int it = 0;

  int ist  = ixyz + Nxyz * it;
  for(int ic2 = 0; ic2 < Ncol2; ++ic2){
    buf_tp[IDX2(Ncol2, ic2, ixyz)] = v1[IDX2(Ncol2, ic2, ist)];
  }
 }

 {
  int it = Nt-1;
  real_t *u2 = &u[Ndf * Nst * idir];

  for(int ivec = 0; ivec < NcolH; ++ivec){
    int ist  = ixyz + Nxyz * it;

    buf_tm[IDX2(Ncol2, 4*ivec  , ixyz)] = 0.0;
    buf_tm[IDX2(Ncol2, 4*ivec+1, ixyz)] = 0.0;
    buf_tm[IDX2(Ncol2, 4*ivec+2, ixyz)] = 0.0;
    buf_tm[IDX2(Ncol2, 4*ivec+3, ixyz)] = 0.0;
    accum_mult_udag(buf_tm, u2, v1, ivec, ixyz, ist, real_t(1.0),
                    Ndf, Ncol2, NcolH);
  }
 }

}

//====================================================================
void mult_clover_coarse_1(real_t* buf_xp, real_t* buf_xm,
                          real_t* buf_yp, real_t* buf_ym,
                          real_t* buf_zp, real_t* buf_zm,
                          real_t* buf_tp, real_t* buf_tm,
                          real_t* u, real_t* v1,
                          int Ncol,
                          int *Nsize, int *do_comm)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];

  real_t* u_dev  = (real_t*)dev_ptr(u);
  real_t* v1_dev = (real_t*)dev_ptr(v1);

  real_t* buf_xp_dev = (real_t*)dev_ptr(buf_xp);
  real_t* buf_xm_dev = (real_t*)dev_ptr(buf_xm);
  real_t* buf_yp_dev = (real_t*)dev_ptr(buf_yp);
  real_t* buf_ym_dev = (real_t*)dev_ptr(buf_ym);
  real_t* buf_zp_dev = (real_t*)dev_ptr(buf_zp);
  real_t* buf_zm_dev = (real_t*)dev_ptr(buf_zm);
  real_t* buf_tp_dev = (real_t*)dev_ptr(buf_tp);
  real_t* buf_tm_dev = (real_t*)dev_ptr(buf_tm);

  if(do_comm[0] > 0){
    int size = Ny * Nz * Nt;
    int nth = VECTOR_LENGTH;
    int nbl = size/nth;
    mult_clover_coarse_1x_kernel<<<nbl,nth>>>(
                                       buf_xp_dev, buf_xm_dev,
                                       u_dev, v1_dev,
                                       Ncol, Nx, Ny, Nz, Nt);
  }

  if(do_comm[1] > 0){
    int size = Nx * Nz * Nt;
    int nth = VECTOR_LENGTH;
    int nbl = size/nth;
    mult_clover_coarse_1y_kernel<<<nbl,nth>>>(
                                       buf_yp_dev, buf_ym_dev,
                                       u_dev, v1_dev,
                                       Ncol, Nx, Ny, Nz, Nt);
  }

  if(do_comm[2] > 0){
    int size = Nx * Ny * Nt;
    int nth = VECTOR_LENGTH;
    int nbl = size/nth;
    mult_clover_coarse_1z_kernel<<<nbl,nth>>>(
                                       buf_zp_dev, buf_zm_dev,
                                       u_dev, v1_dev,
                                       Ncol, Nx, Ny, Nz, Nt);
  }

  if(do_comm[3] > 0){
    int size = Nx * Ny * Nz;
    int nth = VECTOR_LENGTH;
    int nbl = size/nth;
    mult_clover_coarse_1t_kernel<<<nbl,nth>>>(
                                       buf_tp_dev, buf_tm_dev,
                                       u_dev, v1_dev,
                                       Ncol, Nx, Ny, Nz, Nt);
  }

  CHECK(cudaDeviceSynchronize());

  if(do_comm[0] > 0){
    int size_bx = Ncol * 2 * Ny * Nz * Nt;
    update_host(buf_xp, 0, size_bx);
    update_host(buf_xm, 0, size_bx);
  }

  if(do_comm[1] > 0){
    int size_by = Ncol * 2 * Nx * Nz * Nt;
    update_host(buf_yp, 0, size_by);
    update_host(buf_ym, 0, size_by);
  }

  if(do_comm[2] > 0){
    int size_bz = Ncol * 2 * Nx * Ny * Nt;
    update_host(buf_zp, 0, size_bz);
    update_host(buf_zm, 0, size_bz);
  }

  if(do_comm[3] > 0){
    int size_bt = Ncol * 2 * Nx * Ny * Nz;
    update_host(buf_tp, 0, size_bt);
    update_host(buf_tm, 0, size_bt);
  }

}

//====================================================================
__global__
void mult_clover_coarse_2_kernel(real_t* v2, real_t* u, 
                                 real_t* buf_xp, real_t* buf_xm,
                                 real_t* buf_yp, real_t* buf_ym,
                                 real_t* buf_zp, real_t* buf_zm,
                                 real_t* buf_tp, real_t* buf_tm,
                                 int Ncol,
                                 int Nx, int Ny, int Nz, int Nt,
                                 int do_comm_x, int do_comm_y,
                                 int do_comm_z, int do_comm_t)
{
  int Nst  = Nx * Ny * Nz * Nt;

  int Ncol2 = 2 * Ncol;
  int NcolH = Ncol/2;
  int Ndf   = 2 * Ncol * Ncol;

  int iloop = blockIdx.x * blockDim.x + threadIdx.x;
  //  for(int iloop = 0; iloop < Nloop; ++iloop){

  int ivec = iloop % NcolH;
  int site = iloop / NcolH;

  int ix = site % Nx;
  int iy = (site/Nx) % Ny;
  int iz = (site/(Nx*Ny)) % Nz;
  int it = site/(Nx*Ny*Nz);

  int idir = 0;
  if(do_comm_x > 0){

    real_t *u2 = &u[Ndf * Nst * idir];

    if(ix == Nx-1){
      int iyzt = site/Nx;
      accum_mult_u(v2, u2, buf_xp, ivec, site, iyzt, real_t(1.0),
                   Ndf, Ncol2, NcolH);
    }

    if(ix == 0){
      int iyzt = site/Nx;
      for(int i = 0; i < 4; ++i){
        v2[IDX2(Ncol2, 4*ivec + i, site)]
                        += buf_xm[IDX2(Ncol2, 4*ivec + i, iyzt)];
      }
    }

  }

  idir = 1;
  if(do_comm_y > 0){

    real_t *u2 = &u[Ndf * Nst * idir];

    if(iy == Ny-1){
      int ixzt = ix + Nx * (iz + Nz*it);
      accum_mult_u(v2, u2, buf_yp, ivec, site, ixzt, real_t(1.0),
                   Ndf, Ncol2, NcolH);
    }

    if(iy == 0){
      int ixzt = ix + Nx * (iz + Nz*it);
      for(int i = 0; i < 4; ++i){
        v2[IDX2(Ncol2, 4*ivec + i, site)]
                        += buf_ym[IDX2(Ncol2, 4*ivec + i, ixzt)];
      }
    }

  }

  idir = 2;
  if(do_comm_z > 0){

    real_t *u2 = &u[Ndf * Nst * idir];

    if(iz == Nz-1){
      int ixyt = ix + Nx * (iy + Ny * it);
      accum_mult_u(v2, u2, buf_zp, ivec, site, ixyt, real_t(1.0),
                   Ndf, Ncol2, NcolH);
    }

    if(iz == 0){
      int ixyt = ix + Nx * (iy + Ny * it);
      for(int i = 0; i < 4; ++i){
        v2[IDX2(Ncol2, 4*ivec + i, site)]
                        += buf_zm[IDX2(Ncol2, 4*ivec + i, ixyt)];
      }
    }

  }

  idir = 3;
  if(do_comm_t > 0){

    real_t *u2 = &u[Ndf * Nst * idir];

    if(it == Nt-1){
      int ixyz = ix + Nx * (iy + Ny * iz);
      accum_mult_u(v2, u2, buf_tp, ivec, site, ixyz, real_t(1.0),
                   Ndf, Ncol2, NcolH);
    }

    if(it == 0){
      int ixyz = ix + Nx * (iy + Ny * iz);
      for(int i = 0; i < 4; ++i){
        v2[IDX2(Ncol2, 4*ivec + i, site)]
                        += buf_tm[IDX2(Ncol2, 4*ivec + i, ixyz)];
      }
    }

  }
 
}

//====================================================================
void mult_clover_coarse_2(real_t* v2, real_t* u, 
                          real_t* buf_xp, real_t* buf_xm,
                          real_t* buf_yp, real_t* buf_ym,
                          real_t* buf_zp, real_t* buf_zm,
                          real_t* buf_tp, real_t* buf_tm,
                          int Ncol, int *Nsize, int *do_comm)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  if(do_comm[0] > 0){
    int size_bx = Ncol * 2 * Ny * Nz * Nt;
    update_device(buf_xp, 0, size_bx);
    update_device(buf_xm, 0, size_bx);
  }

  if(do_comm[1] > 0){
    int size_by = Ncol * 2 * Nx * Nz * Nt;
    update_device(buf_yp, 0, size_by);
    update_device(buf_ym, 0, size_by);
  }

  if(do_comm[2] > 0){
    int size_bz = Ncol * 2 * Nx * Ny * Nt;
    update_device(buf_zp, 0, size_bz);
    update_device(buf_zm, 0, size_bz);
  }

  if(do_comm[3] > 0){
    int size_bt = Ncol * 2 * Nx * Ny * Nz;
    update_device(buf_tp, 0, size_bt);
    update_device(buf_tm, 0, size_bt);
  }

  real_t* v2_dev = (real_t*)dev_ptr(v2);
  real_t* u_dev  = (real_t*)dev_ptr(u);

  real_t* buf_xp_dev = (real_t*)dev_ptr(buf_xp);
  real_t* buf_xm_dev = (real_t*)dev_ptr(buf_xm);
  real_t* buf_yp_dev = (real_t*)dev_ptr(buf_yp);
  real_t* buf_ym_dev = (real_t*)dev_ptr(buf_ym);
  real_t* buf_zp_dev = (real_t*)dev_ptr(buf_zp);
  real_t* buf_zm_dev = (real_t*)dev_ptr(buf_zm);
  real_t* buf_tp_dev = (real_t*)dev_ptr(buf_tp);
  real_t* buf_tm_dev = (real_t*)dev_ptr(buf_tm);

  int do_comm_x = do_comm[0];
  int do_comm_y = do_comm[1];
  int do_comm_z = do_comm[2];
  int do_comm_t = do_comm[3];

  int NcolH = Ncol/2;
  int Nloop = Nst * NcolH;

  int nth = VECTOR_LENGTH;
  int nbl = Nloop/nth;

  mult_clover_coarse_2_kernel<<<nbl,nth>>>(
                                    v2_dev, u_dev,
                                    buf_xp_dev, buf_xm_dev,
                                    buf_yp_dev, buf_ym_dev,
                                    buf_zp_dev, buf_zm_dev,
                                    buf_tp_dev, buf_tm_dev,
                                    Ncol, Nx, Ny, Nz, Nt,
                                    do_comm_x, do_comm_y,
                                    do_comm_z, do_comm_t);

  CHECK(cudaDeviceSynchronize());

}

//============================================================END=====
