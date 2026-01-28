/*!
      @file    mult_Staggered_eo_openacc2R-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/


#define MULT_UV_R(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)   (u0*v0-u1*v1 + u2*v2-u3*v3 + u4*v4-u5*v5)
#define MULT_UV_I(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)   (u0*v1+u1*v0 + u2*v3+u3*v2 + u4*v5+u5*v4)

//====================================================================
void mult_staggered_Meo(real_t *restrict v2,
                        real_t *restrict u_up, real_t *restrict u_dn,
                        real_t *restrict v1, real_t *restrict x1, 
                        real_t mq2, int *bc, int *Nsize,
                        const int ieo, const int jeo, const int iflag)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst  = Nx * Ny * Nz * Nt;
  int nvst = NVC * Nst;
  int ngst = NDF * Nst * 4;

#pragma acc data present(v2[0:nvst], u_up[0:ngst], , u_dn[0:ngst], \
                         v1[0:nvst], x1[0:nvst]) \
                 copyin(bc[0:4], mq2, ieo, jeo, iflag, Nx, Ny, Nz, Nt, Nst)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
    for(int ist = 0; ist < Nst; ++ist){
      int ix = ist % Nx;
      int iy = (ist/Nx) % Ny;
      int iz = (ist/(Nx*Ny)) % Nz;
      int it = ist/(Nx*Ny*Nz);
      int keo  = (jeo + iy + iz + it) % 2;

      real_t bc2;
      int nei, istg, neig;

      real_t vt0, vt1, vt2, vt3, vt4, vt5;
      real_t ut0, ut1, ut2, ut3, ut4, ut5;
      real_t xt0, xt1, xt2, xt3, xt4, xt5;
      real_t wtr, wti;

      xt0 = 0.0;
      xt1 = 0.0;
      xt2 = 0.0;
      xt3 = 0.0;
      xt4 = 0.0;
      xt5 = 0.0;


      int idir = 0;

      // mult_xp
      int iyzt = iy + Ny * (iz + Nz * it);
      nei  = ((ix + keo) % Nx) + Nx * iyzt;
      istg = ist + Nst * (ieo + 2*idir);

      bc2 = 1.0;
      if(ix == Nx-1 && keo == 1) bc2 = bc[0];

#include "inc/mult_Staggered_openacc2R_uvup-inc.h"

      // mult_xm
      int ix2 = (ix - 1 + keo + Nx) % Nx;
      nei  = ix2 + Nx * iyzt;
      neig = nei + Nst * (1-ieo + 2*idir);

      bc2 = 1.0;
      if(ix == 0 && keo == 0) bc2 = bc[0];

#include "inc/mult_Staggered_openacc2R_uvdn-inc.h"


      idir = 1;

      // mult_yp
      int izt = ist/(Nx*Ny);
      int nn = (iy + 1) % Ny;

      nei  = ix + nn*Nx + izt*Nx*Ny;
      istg = ist + Nst * (ieo + 2*idir);

      bc2 = 1.0;
      if(iy == Ny-1) bc2 = bc[1];

#include "inc/mult_Staggered_openacc2R_uvup-inc.h"

      //mult_ym
      nn = (iy + Ny - 1) % Ny;

      nei  = ix + nn*Nx + izt*Nx*Ny;
      neig = nei + Nst * (1-ieo + 2*idir);

      bc2 = 1.0;
      if(iy == 0) bc2 = bc[1];

#include "inc/mult_Staggered_openacc2R_uvdn-inc.h"


      idir = 2;

      //mult_zp
      int ixy = ist % (Nx*Ny);

      nn = (iz + 1) % Nz;

      nei  = ixy + nn*Nx*Ny + it*Nx*Ny*Nz;
      istg = ist + Nst * (ieo + 2*idir);

      bc2 = 1.0;
      if(iz == Nz-1) bc2 = bc[2];

#include "inc/mult_Staggered_openacc2R_uvup-inc.h"

      //mult_zm
      nn = (iz + Nz - 1) % Nz;

      nei  = ixy + nn*Nx*Ny + it*Nx*Ny*Nz;
      neig = nei + Nst * (1-ieo + 2*idir);

      bc2 = 1.0;
      if(iz == 0) bc2 = bc[2];

#include "inc/mult_Staggered_openacc2R_uvdn-inc.h"

      idir = 3;

      //mult_tp
      nn = (it + 1) % Nt;

      int ixyz = ist % (Nx*Ny*Nz);

      nei  = ixyz + nn*Nx*Ny*Nz;
      istg = ist + Nst*(ieo + 2*idir);

      bc2 = 1.0;
      if(it == Nt-1) bc2 = bc[3];

#include "inc/mult_Staggered_openacc2R_uvup-inc.h"

      //mult_tm
      nn = (it + Nt - 1) % Nt;

      nei  = ixyz + nn*Nx*Ny*Nz;
      neig = nei + Nst*(1-ieo + 2*idir);

      bc2 = 1.0;
      if(it == 0) bc2 = bc[3];

#include "inc/mult_Staggered_openacc2R_uvdn-inc.h"

      // real_t fac = 0.5/mq;
      real_t fac = 0.5;
      if(iflag == 0){
        v2[IDX2_1SP_R(0, ist)] = fac * xt0;
        v2[IDX2_1SP_I(0, ist)] = fac * xt1;
        v2[IDX2_1SP_R(1, ist)] = fac * xt2;
        v2[IDX2_1SP_I(1, ist)] = fac * xt3;
        v2[IDX2_1SP_R(2, ist)] = fac * xt4;
        v2[IDX2_1SP_I(2, ist)] = fac * xt5;
      }else{
	/*
        v2[IDX2_1SP_R(0, ist)] = x1[IDX2_1SP_R(0, ist)] - fac * xt0;
        v2[IDX2_1SP_I(0, ist)] = x1[IDX2_1SP_I(0, ist)] - fac * xt1;
        v2[IDX2_1SP_R(1, ist)] = x1[IDX2_1SP_R(1, ist)] - fac * xt2;
        v2[IDX2_1SP_I(1, ist)] = x1[IDX2_1SP_I(1, ist)] - fac * xt3;
        v2[IDX2_1SP_R(2, ist)] = x1[IDX2_1SP_R(2, ist)] - fac * xt4;
        v2[IDX2_1SP_I(2, ist)] = x1[IDX2_1SP_I(2, ist)] - fac * xt5;
	*/
        v2[IDX2_1SP_R(0, ist)] = mq2 * x1[IDX2_1SP_R(0, ist)] - fac * xt0;
        v2[IDX2_1SP_I(0, ist)] = mq2 * x1[IDX2_1SP_I(0, ist)] - fac * xt1;
        v2[IDX2_1SP_R(1, ist)] = mq2 * x1[IDX2_1SP_R(1, ist)] - fac * xt2;
        v2[IDX2_1SP_I(1, ist)] = mq2 * x1[IDX2_1SP_I(1, ist)] - fac * xt3;
        v2[IDX2_1SP_R(2, ist)] = mq2 * x1[IDX2_1SP_R(2, ist)] - fac * xt4;
        v2[IDX2_1SP_I(2, ist)] = mq2 * x1[IDX2_1SP_I(2, ist)] - fac * xt5;
      }

    }

 }

}

//====================================================================
void mult_staggered_1eo(
                    real_t *restrict buf_xp, real_t *restrict buf_xm,
                    real_t *restrict buf_yp, real_t *restrict buf_ym,
                    real_t *restrict buf_zp, real_t *restrict buf_zm,
                    real_t *restrict buf_tp, real_t *restrict buf_tm,
                    real_t *restrict u_dn, real_t *restrict v1,
                    int *bc, int *Nsize, int *do_comm,
                    const int ieo, const int jeo)
{
  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nst  = Nx * Ny * Nz * Nt;

  int size    = NVC * Nst;
  int size_u  = NDF * Nst * 4;
  int size_bx = NVC * ((Ny * Nz * Nt + 1)/2);
  int size_by = NVC * Nx * Nz * Nt;
  int size_bz = NVC * Nx * Ny * Nt;
  int size_bt = NVC * Nx * Ny * Nz;

#pragma acc data present(buf_xp[0:size_bx], buf_xm[0:size_bx], \
                         buf_yp[0:size_by], buf_ym[0:size_by], \
                         buf_zp[0:size_bz], buf_zm[0:size_bz], \
                         buf_tp[0:size_bt], buf_tm[0:size_bt], \
                         u_dn[0:size_u], v1[0:size]) \
                 copyin(bc[0:4], do_comm[0:4], Nx, Ny, Nz, Nt)
{

 // int idir = 0;
 if(do_comm[0] > 0){

#pragma acc parallel async \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
   int Nyzt = Ny * Nz * Nt;

#pragma acc loop
  for(int iyzt = 0; iyzt < Nyzt; ++iyzt){
    int ix  = 0;
    int iy = iyzt % Ny;
    int iz = (iyzt/Ny) % Nz;
    int it = iyzt/(Ny * Nz);
    int ist = ix + Nx * iyzt;
    int keo = (jeo + iy + iz + it) % 2;

    if(keo == 1){
      int iyzt2 = iyzt/2;
      real_t bc2 = bc[0];
      real_t *wt = buf_xp;
      for(int ic = 0; ic < NC; ++ic){
        wt[IDX2_1SP_R(ic, iyzt2)] = bc2 * v1[ IDX2_1SP_R(ic, ist) ];
        wt[IDX2_1SP_I(ic, iyzt2)] = bc2 * v1[ IDX2_1SP_I(ic, ist) ];
      }
    }

  }
 }

#pragma acc parallel async \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
   int Nyzt = Ny * Nz * Nt;

#pragma acc loop
  for(int iyzt = 0; iyzt < Nyzt; ++iyzt){
    int ix  = Nx-1;
    int iy = iyzt % Ny;
    int iz = (iyzt/Ny) % Nz;
    int it = iyzt/(Ny * Nz);
    int ist = ix + Nx * iyzt;
    int istu = ist + Nst * (1-ieo + 2*0);
    int keo = (ieo + iy + iz + it) % 2;

    if(keo == 0){
      int ibf = iyzt/2;
      real_t *wt = buf_xm;

#include "inc/mult_Staggered_openacc2R_uvdn1-inc.h"

    }

  }

 }

 } // do_comm[0]

 // idir = 1;
 if(do_comm[1] > 0){

#pragma acc parallel async \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
   int Nzt = Nz * Nt;

#pragma acc loop collapse(2)
  for(int izt = 0; izt < Nzt; ++izt){
   for(int ix = 0; ix < Nx; ++ix){
     int iy   = 0;
     int ist  = ix + Nx * (iy + Ny * izt);
     int ixzt = ix + Nx * izt;
     real_t bc2 = bc[1];
     real_t *wt = buf_yp;
     for(int ic = 0; ic < NC; ++ic){
       wt[IDX2_1SP_R(ic, ixzt)] = bc2 * v1[ IDX2_1SP_R(ic, ist) ];
       wt[IDX2_1SP_I(ic, ixzt)] = bc2 * v1[ IDX2_1SP_I(ic, ist) ];
     }
   }
  }

 }

#pragma acc parallel async \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
   int Nzt = Nz * Nt;
   int Nst = Nx * Ny * Nzt;

#pragma acc loop collapse(2)
  for(int izt = 0; izt < Nzt; ++izt){
   for(int ix = 0; ix < Nx; ++ix){
     int iy   = Ny-1;
     int ist  = ix + Nx * (iy + Ny*izt);
     int istu = ist + Nst * (1-ieo + 2*1);
     int ibf = ix + Nx * izt;
     real_t *wt = buf_ym;

#include "inc/mult_Staggered_openacc2R_uvdn1-inc.h"

   }
  }

 }

 } // do_comm[1]

 //  idir = 2;
 if(do_comm[2] > 0){

#pragma acc parallel async \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
   int Nxy = Nx * Ny;

#pragma acc loop collapse(2)
  for(int it = 0; it < Nt; ++it){
   for(int ixy = 0; ixy < Nxy; ++ixy){
     int iz   = 0;
     int ist  = ixy + Nxy * (iz + Nz * it);
     int ixyt = ixy + Nxy * it;
     real_t bc2 = bc[2];
     real_t *wt = buf_zp;
     for(int ic = 0; ic < NC; ++ic){
       wt[IDX2_1SP_R(ic, ixyt)] = bc2 * v1[ IDX2_1SP_R(ic, ist) ];
       wt[IDX2_1SP_I(ic, ixyt)] = bc2 * v1[ IDX2_1SP_I(ic, ist) ];
     }
   }
  }

 }

#pragma acc parallel async \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
   int Nxy = Nx * Ny;
   int Nst = Nxy * Nz * Nt;

#pragma acc loop collapse(2)
  for(int it = 0; it < Nt; ++it){
   for(int ixy = 0; ixy < Nxy; ++ixy){
     int iz = Nz-1;
     int ist  = ixy + Nxy * (iz + Nz * it);
     int istu = ist + Nst * (1-ieo + 2*2);
     int ibf = ixy + Nxy * it;
     real_t *wt = buf_zm;

#include "inc/mult_Staggered_openacc2R_uvdn1-inc.h"

   }
  }

 }

 } // do_comm[2]

 //  idir = 3;
 if(do_comm[3] > 0){

#pragma acc parallel async \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
   int Nxyz = Nx * Ny * Nz;

#pragma acc loop
  for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
    int it = 0;
    int ist  = ixyz + Nxyz * it;
    real_t bc2 = bc[3];
    real_t *wt = buf_tp;
    for(int ic = 0; ic < NC; ++ic){
      wt[IDX2_1SP_R(ic, ixyz)] = bc2 * v1[ IDX2_1SP_R(ic, ist) ];
      wt[IDX2_1SP_I(ic, ixyz)] = bc2 * v1[ IDX2_1SP_I(ic, ist) ];
    }
  }

 }

#pragma acc parallel async \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
   int Nxyz = Nx * Ny * Nz;
   int Nst  = Nxyz * Nt;

#pragma acc loop
  for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
    int it   = Nt-1;
    int ist  = ixyz + Nxyz * it;
    int istu = ist + Nst * (1-ieo + 2*3);
    int ibf = ixyz;
    real_t *wt = buf_tm;

#include "inc/mult_Staggered_openacc2R_uvdn1-inc.h"

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
void mult_staggered_2eo(real_t *restrict v2, real_t *restrict u_up, 
                        real_t *restrict buf_xp, real_t *restrict buf_xm,
                        real_t *restrict buf_yp, real_t *restrict buf_ym,
                        real_t *restrict buf_zp, real_t *restrict buf_zm,
                        real_t *restrict buf_tp, real_t *restrict buf_tm,
                        real_t mq,int *bc,  int *Nsize, int *do_comm,
                        const int ieo, const int jeo, const int iflag)
{
  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nst  = Nx * Ny * Nz * Nt;

  int size    = NVC * Nst;
  int size_u  = NDF * Nst * 4;
  int size_bx = NVC * ((Ny * Nz * Nt + 1)/2);
  int size_by = NVC * Nx * Nz * Nt;
  int size_bz = NVC * Nx * Ny * Nt;
  int size_bt = NVC * Nx * Ny * Nz;

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
			 v2[0:size], u_up[0:size_u]) \
                 copyin(bc[0:4], do_comm[0:4], ieo, jeo, mq, Nx, Ny, Nz, Nt)
{

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
   int Nst = Nx * Ny * Nz * Nt;

#pragma acc loop
  for(int ist = 0; ist < Nst; ++ist){
    int ix = ist % Nx;
    int iy = (ist/Nx) % Ny;
    int iz = (ist/(Nx*Ny)) % Nz;
    int it = ist/(Nx*Ny*Nz);
    int keo  = (jeo + iy + iz + it) % 2;

    real_t xt0, xt1, xt2, xt3, xt4, xt5;
    xt0 = 0.0;
    xt1 = 0.0;
    xt2 = 0.0;
    xt3 = 0.0;
    xt4 = 0.0;
    xt5 = 0.0;

    // idir = 0
    if(do_comm[0] > 0){

      if(ix == Nx-1 && keo == 1){
        int iyzt = iy + Ny * (iz + Nz * it);
        int ibuf = iyzt/2;
        int istu = ist + Nst * (ieo + 2*0);
        real_t *buf = buf_xp;

#include "inc/mult_Staggered_openacc2R_uvup2-inc.h"

      }

      if(ix == 0 && keo == 0){
        int iyzt   = iy + Ny * (iz + Nz * it);
        int ibf = iyzt/2;
        real_t bc2 = bc[0];
        xt0 -= bc2 * buf_xm[IDX2_1SP_R(0, ibf)];
        xt1 -= bc2 * buf_xm[IDX2_1SP_I(0, ibf)];
        xt2 -= bc2 * buf_xm[IDX2_1SP_R(1, ibf)];
        xt3 -= bc2 * buf_xm[IDX2_1SP_I(1, ibf)];
        xt4 -= bc2 * buf_xm[IDX2_1SP_R(2, ibf)];
        xt5 -= bc2 * buf_xm[IDX2_1SP_I(2, ibf)];
      }

    }

    // idir = 1
    if(do_comm[1] > 0){

      if(iy == Ny-1){
        int ibuf = ix + Nx * (iz + Nz * it);
        int istu = ist + Nst * (ieo + 2*1);
        real_t *buf = buf_yp;

#include "inc/mult_Staggered_openacc2R_uvup2-inc.h"

      }

      if(iy == 0){
	int ibf = ix + Nx * (iz + Nz * it);
        real_t bc2 = bc[1];
        xt0 -= bc2 * buf_ym[IDX2_1SP_R(0, ibf)];
        xt1 -= bc2 * buf_ym[IDX2_1SP_I(0, ibf)];
        xt2 -= bc2 * buf_ym[IDX2_1SP_R(1, ibf)];
        xt3 -= bc2 * buf_ym[IDX2_1SP_I(1, ibf)];
        xt4 -= bc2 * buf_ym[IDX2_1SP_R(2, ibf)];
        xt5 -= bc2 * buf_ym[IDX2_1SP_I(2, ibf)];

      }

    }

    // idir = 2
    if(do_comm[2] > 0){

      if(iz == Nz-1){
        int ibuf = ix + Nx * (iy + Ny * it);
        int istu = ist + Nst * (ieo + 2*2);
        real_t *buf = buf_zp;

#include "inc/mult_Staggered_openacc2R_uvup2-inc.h"

      }

      if(iz == 0){
        int ibf = ix + Nx * (iy + Ny * it);
        real_t bc2 = bc[2];
        xt0 -= bc2 * buf_zm[IDX2_1SP_R(0, ibf)];
        xt1 -= bc2 * buf_zm[IDX2_1SP_I(0, ibf)];
        xt2 -= bc2 * buf_zm[IDX2_1SP_R(1, ibf)];
        xt3 -= bc2 * buf_zm[IDX2_1SP_I(1, ibf)];
        xt4 -= bc2 * buf_zm[IDX2_1SP_R(2, ibf)];
        xt5 -= bc2 * buf_zm[IDX2_1SP_I(2, ibf)];
      }

    }

    // idir = 3
    if(do_comm[3] > 0){

      if(it == Nt-1){
        int ibuf = ix + Nx * (iy + Ny * iz);
        int istu = ist + Nst * (ieo + 2*3);
        real_t *buf = buf_tp;

#include "inc/mult_Staggered_openacc2R_uvup2-inc.h"

      }

      if(it == 0){
        int ibf = ix + Nx * (iy + Ny * iz);
        real_t bc2 = bc[3];
        xt0 -= bc2 * buf_tm[IDX2_1SP_R(0, ibf)];
        xt1 -= bc2 * buf_tm[IDX2_1SP_I(0, ibf)];
        xt2 -= bc2 * buf_tm[IDX2_1SP_R(1, ibf)];
        xt3 -= bc2 * buf_tm[IDX2_1SP_I(1, ibf)];
        xt4 -= bc2 * buf_tm[IDX2_1SP_R(2, ibf)];
        xt5 -= bc2 * buf_tm[IDX2_1SP_I(2, ibf)];
      }

    }

    // real_t fac = 0.5/mq;
    real_t fac = 0.5;
    if(iflag != 0) fac = -fac;
    v2[IDX2_1SP_R(0, ist)] += fac * xt0;
    v2[IDX2_1SP_I(0, ist)] += fac * xt1;
    v2[IDX2_1SP_R(1, ist)] += fac * xt2;
    v2[IDX2_1SP_I(1, ist)] += fac * xt3;
    v2[IDX2_1SP_R(2, ist)] += fac * xt4;
    v2[IDX2_1SP_I(2, ist)] += fac * xt5;

  }

 }

 }

}

//============================================================END=====
