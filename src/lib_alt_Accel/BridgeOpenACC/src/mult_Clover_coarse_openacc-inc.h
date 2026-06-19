/*!
      @file    mult_Clover_coarse_openacc-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2021-05-17 15:25:49 #$
      @version $LastChangedRevision: 2258 $
*/

#include "lib_alt_Accel/BridgeACC/define_index.h"

namespace{

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
void mult_clover_coarse_mult_gm5(real_t *restrict v2,
                                 real_t *restrict v1,
                                 int Ncol, int *Nsize)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst  = Nx * Ny * Nz * Nt;
  int Nvst = 2 * Ncol * Nst;

#pragma acc data present(v1[0:Nvst], v2[0:Nvst]) copyin(Ncol, Nst)

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
  int Ncst = Ncol * Nst;
  int Nvc  = 2 * Ncol;

#pragma acc loop
  for(int icst = 0; icst < Ncst; ++icst){
    int ic   = icst % Ncol;
    int site = icst/Ncol;
    int ch = ic % 2;
    if(ch == 0){  // negative chirality
      v2[IDX2(Nvc, 2*ic,   site)] = -v1[IDX2(Nvc, 2*ic,   site)];
      v2[IDX2(Nvc, 2*ic+1, site)] = -v1[IDX2(Nvc, 2*ic+1, site)];
    }else{        // positive chirality
      v2[IDX2(Nvc, 2*ic,   site)] =  v1[IDX2(Nvc, 2*ic,   site)];
      v2[IDX2(Nvc, 2*ic+1, site)] =  v1[IDX2(Nvc, 2*ic+1, site)];
    }
  }

 } // acc parallel

}

//====================================================================
void mult_clover_coarse_bulk(
                             real_t *restrict v2, real_t *restrict u,
                             real_t *restrict ct, real_t *restrict v1,
                             int Ncol, int *Nsize, int *bc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst  = Nx * Ny * Nz * Nt;
  int Nvst = 2 * Ncol * Nst;
  int Ngst = 2 * Ncol * Ncol * Nst * 4;
  int Nctst = 2 * Ncol * Ncol * Nst;

#pragma acc data present(v1[0:Nvst], u[0:Ngst], ct[0:Nctst], v2[0:Nvst]) \
                 copyin(bc[0:4], Ncol, Nx, Ny, Nz, Nt, Nst)

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
   int Ncol2 = 2 * Ncol;
   int NcolH = Ncol/2;
   int Ndf   = 2 * Ncol * Ncol;
   int Ngst1 = 2 *  Ncol * Ncol * Nst;
   int Nxy   = Nx * Ny;
   int Nxyz  = Nx * Ny * Nz;

   int Nloop = NcolH * Nst;

#pragma acc loop
    for(int iloop = 0; iloop < Nloop; ++iloop){
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
        bc2 = bc[idir];
      }
      accum_mult_u(v2, u2, v1, ivec, site, nei, bc2, Ndf, Ncol2, NcolH);

      nei = site - 1;
      bc2 = 1.0;
      if(ix == 0){
        nei = site + Nx-1;
        bc2 = bc[idir];
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
        bc2 = bc[idir];
      }
      accum_mult_u(v2, u2, v1, ivec, site, nei, bc2, Ndf, Ncol2, NcolH);

      nei = site - Nx;
      bc2 = 1.0;
      if(iy == 0){
        nei = site + Nx*(Ny-1);
        bc2 = bc[idir];
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
        bc2 = bc[idir];
      }
      accum_mult_u(v2, u2, v1, ivec, site, nei, bc2, Ndf, Ncol2, NcolH);

      nei = site - Nxy;
      bc2 = 1.0;
      if(iz == 0){
        nei = site + Nxy*(Nz-1);
        bc2 = bc[idir];
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
        bc2 = bc[idir];
      }

      accum_mult_u(v2, u2, v1, ivec, site, nei, bc2, Ndf, Ncol2, NcolH);

      nei = site - Nxyz;
      bc2 = 1.0;
      if(it == 0){
        nei = site + Nxyz*(Nt-1);
        bc2 = bc[idir];
      }
      accum_mult_udag(v2, u2, v1, ivec, site, nei, bc2, Ndf, Ncol2, NcolH);
     }

     accum_mult_u(v2, ct, v1, ivec, site, site, real_t(1.0),
                  Ndf, Ncol2, NcolH);

    }

 } // acc parallel

}

//====================================================================
void mult_clover_coarse_1(
                    real_t *restrict buf_xp, real_t *restrict buf_xm,
                    real_t *restrict buf_yp, real_t *restrict buf_ym,
                    real_t *restrict buf_zp, real_t *restrict buf_zm,
                    real_t *restrict buf_tp, real_t *restrict buf_tm,
                    real_t *restrict u, real_t *restrict v1,
                    int Ncol, int *Nsize, int *do_comm)
{
  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nst  = Nx * Ny * Nz * Nt;

  int size    = 2 * Ncol * Nst;
  int size_u  = 2 * Ncol * Ncol * Nst * 4;
  int size_bx = 2 * Ncol * Ny * Nz * Nt;
  int size_by = 2 * Ncol * Nx * Nz * Nt;
  int size_bz = 2 * Ncol * Nx * Ny * Nt;
  int size_bt = 2 * Ncol * Nx * Ny * Nz;

#pragma acc data present(buf_xp[0:size_bx], buf_xm[0:size_bx], \
                         buf_yp[0:size_by], buf_ym[0:size_by], \
                         buf_zp[0:size_bz], buf_zm[0:size_bz], \
                         buf_tp[0:size_bt], buf_tm[0:size_bt], \
                         u[0:size_u], v1[0:size]) \
                 copyin(do_comm[0:4], Nx, Ny, Nz, Nt, Ncol)
{

  int Ncol2 = 2 * Ncol; 
  int NcolH = Ncol/2;
  int Ndf   = 2 * Ncol * Ncol;

  int idir = 0;
  if(do_comm[idir] > 0){

#pragma acc parallel async \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
    int Nyzt = Ny * Nz * Nt;

    int ix = 0;
#pragma acc loop
    for(int iyzt = 0; iyzt < Nyzt; ++iyzt){
      int ist  = ix + Nx * iyzt;
      for(int ic2 = 0; ic2 < Ncol2; ++ic2){
        buf_xp[IDX2(Ncol2, ic2, iyzt)] = v1[IDX2(Ncol2, ic2, ist)];
      }
    }

   }

#pragma acc parallel async \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
    int Nyzt = Ny * Nz * Nt;
    real_t *u2 = &u[Ndf * Nst * idir];

    int ix = Nx-1;
    int Nloop = Nyzt * NcolH;
#pragma acc loop
    for(int iloop = 0; iloop < Nloop; ++iloop){
      int ivec = iloop % NcolH;
      int iyzt = iloop/NcolH;
      int ist  = ix + Nx * iyzt;
      buf_xm[IDX2(Ncol2, 4*ivec  , iyzt)] = 0.0;
      buf_xm[IDX2(Ncol2, 4*ivec+1, iyzt)] = 0.0;
      buf_xm[IDX2(Ncol2, 4*ivec+2, iyzt)] = 0.0;
      buf_xm[IDX2(Ncol2, 4*ivec+3, iyzt)] = 0.0;
      accum_mult_udag(buf_xm, u2, v1, ivec, iyzt, ist, real_t(1.0),
                      Ndf, Ncol2, NcolH);
    }
 }

  } // do_comm[0]

  idir = 1;
  if(do_comm[idir] > 0){

#pragma acc parallel async \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
    int Nxzt = Nx * Nz * Nt;

    int iy = 0;
#pragma acc loop
    for(int ixzt = 0; ixzt < Nxzt; ++ixzt){
      int ix  = ixzt % Nx;
      int izt = ixzt/Nx;
      int ist  = ix + Nx * (iy + Ny * izt);
      for(int ic2 = 0; ic2 < Ncol2; ++ic2){
        buf_yp[IDX2(Ncol2, ic2, ixzt)] = v1[IDX2(Ncol2, ic2, ist)];
      }
    }

 }

#pragma acc parallel async \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
    int Nxzt = Nx * Nz * Nt;
    real_t *u2 = &u[Ndf * Nst * idir];

    int iy = Ny-1;
    int Nloop = Nxzt * NcolH;
#pragma acc loop
    for(int iloop = 0; iloop < Nloop; ++iloop){
      int ivec = iloop % NcolH;
      int ixzt = iloop/NcolH;
      int ix  = ixzt % Nx;
      int izt = ixzt/Nx;
      int ist  = ix + Nx * (iy + Ny * izt);
      buf_ym[IDX2(Ncol2, 4*ivec  , ixzt)] = 0.0;
      buf_ym[IDX2(Ncol2, 4*ivec+1, ixzt)] = 0.0;
      buf_ym[IDX2(Ncol2, 4*ivec+2, ixzt)] = 0.0;
      buf_ym[IDX2(Ncol2, 4*ivec+3, ixzt)] = 0.0;
      accum_mult_udag(buf_ym, u2, v1, ivec, ixzt, ist, real_t(1.0),
                      Ndf, Ncol2, NcolH);
    }
 }

  } // do_comm[1]

   idir = 2;
   if(do_comm[idir] > 0){

#pragma acc parallel async \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
    int Nxy  = Nx * Ny;
    int Nxyt = Nxy * Nt;

    int iz = 0;
#pragma acc loop
    for(int ixyt = 0; ixyt < Nxyt; ++ixyt){
      int ixy = ixyt % Nxy;
      int it  = ixyt/Nxy;
      int ist  = ixy + Nxy *(iz + Nz*it);
      for(int ic2 = 0; ic2 < Ncol2; ++ic2){
        buf_zp[IDX2(Ncol2, ic2, ixyt)] = v1[IDX2(Ncol2, ic2, ist)];
      }
    }

 }

#pragma acc parallel async \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
    int Nxy  = Nx * Ny;
    int Nxyt = Nxy * Nt;
    real_t *u2 = &u[Ndf * Nst * idir];

    int iz = Nz-1;
    int Nloop = Nxyt * NcolH;
#pragma acc loop
    for(int iloop = 0; iloop < Nloop; ++iloop){
      int ivec = iloop % NcolH;
      int ixyt = iloop/NcolH;
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

 } // do_comm[2]

  idir = 3;
  if(do_comm[idir] > 0){

#pragma acc parallel async \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
    int Nxyz = Nx * Ny * Nz;

    int it = 0;
#pragma acc loop
    for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
      int ist  = ixyz + Nxyz * it;
      for(int ic2 = 0; ic2 < Ncol2; ++ic2){
        buf_tp[IDX2(Ncol2, ic2, ixyz)] = v1[IDX2(Ncol2, ic2, ist)];
      }
    }

 }

#pragma acc parallel async \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
    int Nxyz = Nx * Ny * Nz;
    real_t *u2 = &u[Ndf * Nst * idir];

    int it = Nt-1;
    int Nloop = Nxyz * NcolH;
#pragma acc loop
    for(int iloop = 0; iloop < Nloop; ++iloop){
      int ivec = iloop % NcolH;
      int ixyz = iloop/NcolH;
      int ist  = ixyz + Nxyz * it;
      buf_tm[IDX2(Ncol2, 4*ivec  , ixyz)] = 0.0;
      buf_tm[IDX2(Ncol2, 4*ivec+1, ixyz)] = 0.0;
      buf_tm[IDX2(Ncol2, 4*ivec+2, ixyz)] = 0.0;
      buf_tm[IDX2(Ncol2, 4*ivec+3, ixyz)] = 0.0;
      accum_mult_udag(buf_tm, u2, v1, ivec, ixyz, ist, real_t(1.0),
                      Ndf, Ncol2, NcolH);
    }

 }

 } // do_comm[3]

 #pragma acc wait

 } // acc data

 if(do_comm[0] > 0){
   #pragma acc update async host (buf_xp[0:size_bx])
   #pragma acc update async host (buf_xm[0:size_bx])
 }
 if(do_comm[1] > 0){
   #pragma acc update async host (buf_yp[0:size_by])
   #pragma acc update async host (buf_ym[0:size_by])
 }
 if(do_comm[2] > 0){
   #pragma acc update async host (buf_zp[0:size_bz])
   #pragma acc update async host (buf_zm[0:size_bz])
 }
 if(do_comm[3] > 0){
   #pragma acc update async host (buf_tp[0:size_bt])
   #pragma acc update async host (buf_tm[0:size_bt])
 }

 #pragma acc wait

}

//====================================================================
void mult_clover_coarse_2(
                    real_t *restrict v2, real_t *restrict u, 
                    real_t *restrict buf_xp, real_t *restrict buf_xm,
                    real_t *restrict buf_yp, real_t *restrict buf_ym,
                    real_t *restrict buf_zp, real_t *restrict buf_zm,
                    real_t *restrict buf_tp, real_t *restrict buf_tm,
                    int Ncol, int *Nsize, int *do_comm)
{
  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nst  = Nx * Ny * Nz * Nt;

  int size    = 2 * Ncol * Nst;
  int size_u  = 2 * Ncol * Ncol * Nst * 4;
  int size_bx = 2 * Ncol * Ny * Nz * Nt;
  int size_by = 2 * Ncol * Nx * Nz * Nt;
  int size_bz = 2 * Ncol * Nx * Ny * Nt;
  int size_bt = 2 * Ncol * Nx * Ny * Nz;

 if(do_comm[0] > 0){
   #pragma acc update async device (buf_xp[0:size_bx])
   #pragma acc update async device (buf_xm[0:size_bx])
 }
 if(do_comm[1] > 0){
   #pragma acc update async device (buf_yp[0:size_by])
   #pragma acc update async device (buf_ym[0:size_by])
 }
 if(do_comm[2] > 0){
   #pragma acc update async device (buf_zp[0:size_bz])
   #pragma acc update async device (buf_zm[0:size_bz])
 }
 if(do_comm[3] > 0){
   #pragma acc update async device (buf_tp[0:size_bt])
   #pragma acc update async device (buf_tm[0:size_bt])
 }

 #pragma acc wait

#pragma acc data present(buf_xp[0:size_bx], buf_xm[0:size_bx], \
                         buf_yp[0:size_by], buf_ym[0:size_by], \
                         buf_zp[0:size_bz], buf_zm[0:size_bz], \
                         buf_tp[0:size_bt], buf_tm[0:size_bt], \
			 v2[0:size], u[0:size_u]) \
                 copyin(do_comm[0:4], Ncol, Nx, Ny, Nz, Nt)

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
  int Nst   = Nx * Ny * Nz * Nt;
  int Ncol2 = 2 * Ncol;
  int NcolH = Ncol/2;
  int Ndf   = 2 * Ncol * Ncol;

  int Nloop = Nst * NcolH;

#pragma acc loop
  for(int iloop = 0; iloop < Nloop; ++iloop){
    int ivec = iloop % NcolH;
    int site = iloop/NcolH;

    int ix = site % Nx;
    int iy = (site/Nx) % Ny;
    int iz = (site/(Nx*Ny)) % Nz;
    int it = site/(Nx*Ny*Nz);

    int idir = 0;
    if(do_comm[idir] > 0){

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
    if(do_comm[idir] > 0){

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
    if(do_comm[idir] > 0){

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
    if(do_comm[idir] > 0){

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

 }

}

//============================================================END=====
