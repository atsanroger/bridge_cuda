/*!
      @file    afopr_Wilson_eo_openacc2-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2553 $
*/


#define MULT_GXr(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)   (u0*v0-u1*v1 + u2*v2-u3*v3 + u4*v4-u5*v5)
#define MULT_GXi(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)   (u0*v1+u1*v0 + u2*v3+u3*v2 + u4*v5+u5*v4)
#define MULT_GDXr(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)  (u0*v0+u1*v1 + u2*v2+u3*v3 + u4*v4+u5*v5)
#define MULT_GDXi(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)  (u0*v1-u1*v0 + u2*v3-u3*v2 + u4*v5-u5*v4)

#define EXT_IMG_R(v1r,v1i,v2r,v2i,w1r,w1i,w2r,w2i)  (v1r*w2r - v1i*w2i - v2r*w1r + v2i*w1i)
#define EXT_IMG_I(v1r,v1i,v2r,v2i,w1r,w1i,w2r,w2i)  (- v1r*w2i - v1i*w2r + v2r*w1i + v2i*w1r)


//====================================================================
void mult_wilson_Meo_dirac(real_t *restrict v2, real_t *restrict u,
                           real_t *restrict v1, real_t *restrict x1, 
                           const int ieo, const int jeo, const int iflag,
                           int *Nsize, int *bc, real_t kappa)
{
  int nx = Nsize[0];
  int ny = Nsize[1];
  int nz = Nsize[2];
  int nt = Nsize[3];
  int nst  = nx * ny * nz * nt;
  int nvst = NVC * ND * nst;
  int ngst = NDF * nst * 4;

  real_t *restrict u_up = u;
  real_t *restrict u_dn = u;

#pragma acc data present(v2[0:nvst], u_up[0:ngst], u_dn[0:ngst], \
                         v1[0:nvst], x1[0:nvst]) \
                 copyin(bc[0:4], kappa, ieo, jeo, iflag, nx, ny, nz, nt, nst)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
    for(int site = 0; site < nst; ++site){
      int ix = site % nx;
      int iy = (site/nx) % ny;
      int iz = (site/(nx*ny)) % nz;
      int it = site/(nx*ny*nz);
      int keo  = (jeo + iy + iz + it) % 2;

      real_t bc2;

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

#include "inc/mult_Wilson_eo_openacc2_xyz-inc.h"

#include "inc/mult_Wilson_eo_openacc2_t_dirac-inc.h"

      // aypx and write back to global memory
      if(iflag == 0){
        v2[IDX2_SP_R(0,0,site)] = - kappa * v2_01;
        v2[IDX2_SP_I(0,0,site)] = - kappa * v2_11;
        v2[IDX2_SP_R(1,0,site)] = - kappa * v2_21;
        v2[IDX2_SP_I(1,0,site)] = - kappa * v2_31;
        v2[IDX2_SP_R(2,0,site)] = - kappa * v2_41;
        v2[IDX2_SP_I(2,0,site)] = - kappa * v2_51;

        v2[IDX2_SP_R(0,1,site)] = - kappa * v2_02;
        v2[IDX2_SP_I(0,1,site)] = - kappa * v2_12;
        v2[IDX2_SP_R(1,1,site)] = - kappa * v2_22;
        v2[IDX2_SP_I(1,1,site)] = - kappa * v2_32;
        v2[IDX2_SP_R(2,1,site)] = - kappa * v2_42;
        v2[IDX2_SP_I(2,1,site)] = - kappa * v2_52;

        v2[IDX2_SP_R(0,2,site)] = - kappa * v2_03;
        v2[IDX2_SP_I(0,2,site)] = - kappa * v2_13;
        v2[IDX2_SP_R(1,2,site)] = - kappa * v2_23;
        v2[IDX2_SP_I(1,2,site)] = - kappa * v2_33;
        v2[IDX2_SP_R(2,2,site)] = - kappa * v2_43;
        v2[IDX2_SP_I(2,2,site)] = - kappa * v2_53;

        v2[IDX2_SP_R(0,3,site)] = - kappa * v2_04;
        v2[IDX2_SP_I(0,3,site)] = - kappa * v2_14;
        v2[IDX2_SP_R(1,3,site)] = - kappa * v2_24;
        v2[IDX2_SP_I(1,3,site)] = - kappa * v2_34;
        v2[IDX2_SP_R(2,3,site)] = - kappa * v2_44;
        v2[IDX2_SP_I(2,3,site)] = - kappa * v2_54;
      }else{
        v2[IDX2_SP_R(0,0,site)] = x1[IDX2_SP_R(0,0,site)] + kappa * v2_01;
        v2[IDX2_SP_I(0,0,site)] = x1[IDX2_SP_I(0,0,site)] + kappa * v2_11;
        v2[IDX2_SP_R(1,0,site)] = x1[IDX2_SP_R(1,0,site)] + kappa * v2_21;
        v2[IDX2_SP_I(1,0,site)] = x1[IDX2_SP_I(1,0,site)] + kappa * v2_31;
        v2[IDX2_SP_R(2,0,site)] = x1[IDX2_SP_R(2,0,site)] + kappa * v2_41;
        v2[IDX2_SP_I(2,0,site)] = x1[IDX2_SP_I(2,0,site)] + kappa * v2_51;

        v2[IDX2_SP_R(0,1,site)] = x1[IDX2_SP_R(0,1,site)] + kappa * v2_02;
        v2[IDX2_SP_I(0,1,site)] = x1[IDX2_SP_I(0,1,site)] + kappa * v2_12;
        v2[IDX2_SP_R(1,1,site)] = x1[IDX2_SP_R(1,1,site)] + kappa * v2_22;
        v2[IDX2_SP_I(1,1,site)] = x1[IDX2_SP_I(1,1,site)] + kappa * v2_32;
        v2[IDX2_SP_R(2,1,site)] = x1[IDX2_SP_R(2,1,site)] + kappa * v2_42;
        v2[IDX2_SP_I(2,1,site)] = x1[IDX2_SP_I(2,1,site)] + kappa * v2_52;

        v2[IDX2_SP_R(0,2,site)] = x1[IDX2_SP_R(0,2,site)] + kappa * v2_03;
        v2[IDX2_SP_I(0,2,site)] = x1[IDX2_SP_I(0,2,site)] + kappa * v2_13;
        v2[IDX2_SP_R(1,2,site)] = x1[IDX2_SP_R(1,2,site)] + kappa * v2_23;
        v2[IDX2_SP_I(1,2,site)] = x1[IDX2_SP_I(1,2,site)] + kappa * v2_33;
        v2[IDX2_SP_R(2,2,site)] = x1[IDX2_SP_R(2,2,site)] + kappa * v2_43;
        v2[IDX2_SP_I(2,2,site)] = x1[IDX2_SP_I(2,2,site)] + kappa * v2_53;

        v2[IDX2_SP_R(0,3,site)] = x1[IDX2_SP_R(0,3,site)] + kappa * v2_04;
        v2[IDX2_SP_I(0,3,site)] = x1[IDX2_SP_I(0,3,site)] + kappa * v2_14;
        v2[IDX2_SP_R(1,3,site)] = x1[IDX2_SP_R(1,3,site)] + kappa * v2_24;
        v2[IDX2_SP_I(1,3,site)] = x1[IDX2_SP_I(1,3,site)] + kappa * v2_34;
        v2[IDX2_SP_R(2,3,site)] = x1[IDX2_SP_R(2,3,site)] + kappa * v2_44;
        v2[IDX2_SP_I(2,3,site)] = x1[IDX2_SP_I(2,3,site)] + kappa * v2_54;
      }

    }

 }

}

//====================================================================
void mult_wilson_Meo_chiral(real_t *restrict v2, real_t *restrict u,
                            real_t *restrict v1, real_t *restrict x1, 
                            const int ieo, const int jeo, const int iflag,
                            int *Nsize, int *bc, real_t kappa)
{
  int nx = Nsize[0];
  int ny = Nsize[1];
  int nz = Nsize[2];
  int nt = Nsize[3];
  int nst  = nx*ny*nz*nt;
  int nvst = NVC * ND * nst;
  int ngst = NDF * nst * 4;

  real_t *restrict u_up = u;
  real_t *restrict u_dn = u;

#pragma acc data present(v2[0:nvst], u_up[0:ngst], u_dn[0:ngst], \
                         v1[0:nvst], x1[0:nvst]) \
                 copyin(bc[0:4], kappa, ieo, jeo, iflag, nx, ny, nz, nt, nst)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
    for(int site = 0; site < nst; ++site){
      int ix = site % nx;
      int iy = (site/nx) % ny;
      int iz = (site/(nx*ny)) % nz;
      int it = site/(nx*ny*nz);
      int keo  = (jeo + iy + iz + it) % 2;

      real_t bc2;

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

#include "inc/mult_Wilson_eo_openacc2_xyz-inc.h"

#include "inc/mult_Wilson_eo_openacc2_t_chiral-inc.h"

      // aypx and write back to global memory
      if(iflag == 0){
        v2[IDX2_SP_R(0,0,site)] = - kappa * v2_01;
        v2[IDX2_SP_I(0,0,site)] = - kappa * v2_11;
        v2[IDX2_SP_R(1,0,site)] = - kappa * v2_21;
        v2[IDX2_SP_I(1,0,site)] = - kappa * v2_31;
        v2[IDX2_SP_R(2,0,site)] = - kappa * v2_41;
        v2[IDX2_SP_I(2,0,site)] = - kappa * v2_51;

        v2[IDX2_SP_R(0,1,site)] = - kappa * v2_02;
        v2[IDX2_SP_I(0,1,site)] = - kappa * v2_12;
        v2[IDX2_SP_R(1,1,site)] = - kappa * v2_22;
        v2[IDX2_SP_I(1,1,site)] = - kappa * v2_32;
        v2[IDX2_SP_R(2,1,site)] = - kappa * v2_42;
        v2[IDX2_SP_I(2,1,site)] = - kappa * v2_52;

        v2[IDX2_SP_R(0,2,site)] = - kappa * v2_03;
        v2[IDX2_SP_I(0,2,site)] = - kappa * v2_13;
        v2[IDX2_SP_R(1,2,site)] = - kappa * v2_23;
        v2[IDX2_SP_I(1,2,site)] = - kappa * v2_33;
        v2[IDX2_SP_R(2,2,site)] = - kappa * v2_43;
        v2[IDX2_SP_I(2,2,site)] = - kappa * v2_53;

        v2[IDX2_SP_R(0,3,site)] = - kappa * v2_04;
        v2[IDX2_SP_I(0,3,site)] = - kappa * v2_14;
        v2[IDX2_SP_R(1,3,site)] = - kappa * v2_24;
        v2[IDX2_SP_I(1,3,site)] = - kappa * v2_34;
        v2[IDX2_SP_R(2,3,site)] = - kappa * v2_44;
        v2[IDX2_SP_I(2,3,site)] = - kappa * v2_54;
      }else{
        v2[IDX2_SP_R(0,0,site)] = x1[IDX2_SP_R(0,0,site)] + kappa * v2_01;
        v2[IDX2_SP_I(0,0,site)] = x1[IDX2_SP_I(0,0,site)] + kappa * v2_11;
        v2[IDX2_SP_R(1,0,site)] = x1[IDX2_SP_R(1,0,site)] + kappa * v2_21;
        v2[IDX2_SP_I(1,0,site)] = x1[IDX2_SP_I(1,0,site)] + kappa * v2_31;
        v2[IDX2_SP_R(2,0,site)] = x1[IDX2_SP_R(2,0,site)] + kappa * v2_41;
        v2[IDX2_SP_I(2,0,site)] = x1[IDX2_SP_I(2,0,site)] + kappa * v2_51;

        v2[IDX2_SP_R(0,1,site)] = x1[IDX2_SP_R(0,1,site)] + kappa * v2_02;
        v2[IDX2_SP_I(0,1,site)] = x1[IDX2_SP_I(0,1,site)] + kappa * v2_12;
        v2[IDX2_SP_R(1,1,site)] = x1[IDX2_SP_R(1,1,site)] + kappa * v2_22;
        v2[IDX2_SP_I(1,1,site)] = x1[IDX2_SP_I(1,1,site)] + kappa * v2_32;
        v2[IDX2_SP_R(2,1,site)] = x1[IDX2_SP_R(2,1,site)] + kappa * v2_42;
        v2[IDX2_SP_I(2,1,site)] = x1[IDX2_SP_I(2,1,site)] + kappa * v2_52;

        v2[IDX2_SP_R(0,2,site)] = x1[IDX2_SP_R(0,2,site)] + kappa * v2_03;
        v2[IDX2_SP_I(0,2,site)] = x1[IDX2_SP_I(0,2,site)] + kappa * v2_13;
        v2[IDX2_SP_R(1,2,site)] = x1[IDX2_SP_R(1,2,site)] + kappa * v2_23;
        v2[IDX2_SP_I(1,2,site)] = x1[IDX2_SP_I(1,2,site)] + kappa * v2_33;
        v2[IDX2_SP_R(2,2,site)] = x1[IDX2_SP_R(2,2,site)] + kappa * v2_43;
        v2[IDX2_SP_I(2,2,site)] = x1[IDX2_SP_I(2,2,site)] + kappa * v2_53;

        v2[IDX2_SP_R(0,3,site)] = x1[IDX2_SP_R(0,3,site)] + kappa * v2_04;
        v2[IDX2_SP_I(0,3,site)] = x1[IDX2_SP_I(0,3,site)] + kappa * v2_14;
        v2[IDX2_SP_R(1,3,site)] = x1[IDX2_SP_R(1,3,site)] + kappa * v2_24;
        v2[IDX2_SP_I(1,3,site)] = x1[IDX2_SP_I(1,3,site)] + kappa * v2_34;
        v2[IDX2_SP_R(2,3,site)] = x1[IDX2_SP_R(2,3,site)] + kappa * v2_44;
        v2[IDX2_SP_I(2,3,site)] = x1[IDX2_SP_I(2,3,site)] + kappa * v2_54;
      }

    }
   
 }

}

//====================================================================
void mult_wilson_eo_1_dirac(
                    real_t *restrict buf_xp, real_t *restrict buf_xm,
                    real_t *restrict buf_yp, real_t *restrict buf_ym,
                    real_t *restrict buf_zp, real_t *restrict buf_zm,
                    real_t *restrict buf_tp, real_t *restrict buf_tm,
                    real_t *restrict u, real_t *restrict v1,
                    const int ieo, const int jeo,
                    int *Nsize, int *bc, int *do_comm, int Nc)
{
  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nst  = Nx * Ny * Nz * Nt;

  int size    = NVC * ND * Nst;
  int size_u  = NDF * Nst * 4;
  int size_bx = NVC * 2 * ((Ny * Nz * Nt + 1)/2);
  int size_by = NVC * 2 * Nx * Nz * Nt;
  int size_bz = NVC * 2 * Nx * Ny * Nt;
  int size_bt = NVC * 2 * Nx * Ny * Nz;

#pragma acc data present(buf_xp[0:size_bx], buf_xm[0:size_bx], \
                         buf_yp[0:size_by], buf_ym[0:size_by], \
                         buf_zp[0:size_bz], buf_zm[0:size_bz], \
                         buf_tp[0:size_bt], buf_tm[0:size_bt], \
                         u[0:size_u], v1[0:size]) \
                 copyin(bc[0:4], do_comm[0:4], Nx, Ny, Nz, Nt, ieo, jeo)
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
      for(int ic = 0; ic < NC; ++ic){
        real_t vt1[2], vt2[2];
        vt1[0] = v1[IDX2_SP_R(ic, 0, ist)] - v1[IDX2_SP_I(ic, 3, ist)];
        vt1[1] = v1[IDX2_SP_I(ic, 0, ist)] + v1[IDX2_SP_R(ic, 3, ist)];
        vt2[0] = v1[IDX2_SP_R(ic, 1, ist)] - v1[IDX2_SP_I(ic, 2, ist)];
        vt2[1] = v1[IDX2_SP_I(ic, 1, ist)] + v1[IDX2_SP_R(ic, 2, ist)];
        buf_xp[IDXBF_R(ic, 0, iyzt2)] = bc2 * vt1[0];
        buf_xp[IDXBF_I(ic, 0, iyzt2)] = bc2 * vt1[1];
        buf_xp[IDXBF_R(ic, 1, iyzt2)] = bc2 * vt2[0];
        buf_xp[IDXBF_I(ic, 1, iyzt2)] = bc2 * vt2[1];
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
      int iyzt2 = iyzt/2;

      real_t vt1[NVC], vt2[NVC], ut[NDF];
      for(int ic = 0; ic < NC; ++ic){
        int icr = 2*ic;
        int ici = 2*ic + 1;
        vt1[icr] = v1[IDX2_SP_R(ic, 0, ist)] + v1[IDX2_SP_I(ic, 3, ist)];
        vt1[ici] = v1[IDX2_SP_I(ic, 0, ist)] - v1[IDX2_SP_R(ic, 3, ist)];
        vt2[icr] = v1[IDX2_SP_R(ic, 1, ist)] + v1[IDX2_SP_I(ic, 2, ist)];
        vt2[ici] = v1[IDX2_SP_I(ic, 1, ist)] - v1[IDX2_SP_R(ic, 2, ist)];
      }

      real_t wt1[2], wt2[2];

      for(int ic = 0; ic < NC; ++ic){

        for(int ic2 = 0; ic2 < NC; ++ic2){
          ut[2*ic2  ] =   u[IDX2_G_R(ic, ic2, istu)];
          ut[2*ic2+1] = - u[IDX2_G_I(ic, ic2, istu)];
        }

        wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                          vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
        wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                          vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
        wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                          vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
        wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                          vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

        buf_xm[IDXBF_R(ic, 0, iyzt2)] = wt1[0];
        buf_xm[IDXBF_I(ic, 0, iyzt2)] = wt1[1];
        buf_xm[IDXBF_R(ic, 1, iyzt2)] = wt2[0];
        buf_xm[IDXBF_I(ic, 1, iyzt2)] = wt2[1];
      }

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

    for(int ic = 0; ic < NC; ++ic){
      real_t vt1[2], vt2[2];
      vt1[0] = v1[IDX2_SP_R(ic, 0, ist)] + v1[IDX2_SP_R(ic, 3, ist)];
      vt1[1] = v1[IDX2_SP_I(ic, 0, ist)] + v1[IDX2_SP_I(ic, 3, ist)];
      vt2[0] = v1[IDX2_SP_R(ic, 1, ist)] - v1[IDX2_SP_R(ic, 2, ist)];
      vt2[1] = v1[IDX2_SP_I(ic, 1, ist)] - v1[IDX2_SP_I(ic, 2, ist)];
      buf_yp[IDXBF_R(ic, 0, ixzt)] = bc2 * vt1[0];
      buf_yp[IDXBF_I(ic, 0, ixzt)] = bc2 * vt1[1];
      buf_yp[IDXBF_R(ic, 1, ixzt)] = bc2 * vt2[0];
      buf_yp[IDXBF_I(ic, 1, ixzt)] = bc2 * vt2[1];
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
     int ixzt = ix + Nx * izt;

     real_t vt1[NVC], vt2[NVC];
     for(int ic = 0; ic < NC; ++ic){
       int icr = 2*ic;
       int ici = 2*ic + 1;
       vt1[icr] = v1[IDX2_SP_R(ic, 0, ist)] - v1[IDX2_SP_R(ic, 3, ist)];
       vt1[ici] = v1[IDX2_SP_I(ic, 0, ist)] - v1[IDX2_SP_I(ic, 3, ist)];
       vt2[icr] = v1[IDX2_SP_R(ic, 1, ist)] + v1[IDX2_SP_R(ic, 2, ist)];
       vt2[ici] = v1[IDX2_SP_I(ic, 1, ist)] + v1[IDX2_SP_I(ic, 2, ist)];
     }

     for(int ic = 0; ic < NC; ++ic){

       real_t ut[NVC];
       for(int ic2 = 0; ic2 < NC; ++ic2){
         ut[2*ic2  ] =   u[IDX2_G_R(ic, ic2, istu)];
         ut[2*ic2+1] = - u[IDX2_G_I(ic, ic2, istu)];
       }

       real_t wt1[2], wt2[2];
       wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
       wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
       wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
       wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

       buf_ym[IDXBF_R(ic, 0, ixzt)] = wt1[0];
       buf_ym[IDXBF_I(ic, 0, ixzt)] = wt1[1];
       buf_ym[IDXBF_R(ic, 1, ixzt)] = wt2[0];
       buf_ym[IDXBF_I(ic, 1, ixzt)] = wt2[1];

     }

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

     for(int ic = 0; ic < NC; ++ic){
       real_t wt1[2], wt2[2];
       wt1[0] = v1[IDX2_SP_R(ic, 0, ist)] - v1[IDX2_SP_I(ic, 2, ist)];
       wt1[1] = v1[IDX2_SP_I(ic, 0, ist)] + v1[IDX2_SP_R(ic, 2, ist)];
       wt2[0] = v1[IDX2_SP_R(ic, 1, ist)] + v1[IDX2_SP_I(ic, 3, ist)];
       wt2[1] = v1[IDX2_SP_I(ic, 1, ist)] - v1[IDX2_SP_R(ic, 3, ist)];
       buf_zp[IDXBF_R(ic, 0, ixyt)] = bc2 * wt1[0];
       buf_zp[IDXBF_I(ic, 0, ixyt)] = bc2 * wt1[1];
       buf_zp[IDXBF_R(ic, 1, ixyt)] = bc2 * wt2[0];
       buf_zp[IDXBF_I(ic, 1, ixyt)] = bc2 * wt2[1];
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
     int ixyt = ixy + Nxy * it;

     real_t vt1[NVC], vt2[NVC];
     for(int ic = 0; ic < NC; ++ic){
       int icr = 2*ic;
       int ici = 2*ic + 1;
       vt1[icr] = v1[IDX2_SP_R(ic, 0, ist)] + v1[IDX2_SP_I(ic, 2, ist)];
       vt1[ici] = v1[IDX2_SP_I(ic, 0, ist)] - v1[IDX2_SP_R(ic, 2, ist)];
       vt2[icr] = v1[IDX2_SP_R(ic, 1, ist)] - v1[IDX2_SP_I(ic, 3, ist)];
       vt2[ici] = v1[IDX2_SP_I(ic, 1, ist)] + v1[IDX2_SP_R(ic, 3, ist)];
     }

     for(int ic = 0; ic < NC; ++ic){

       real_t ut[NVC];
       for(int ic2 = 0; ic2 < NC; ++ic2){
         ut[2*ic2  ] =   u[IDX2_G_R(ic, ic2, istu)];
         ut[2*ic2+1] = - u[IDX2_G_I(ic, ic2, istu)];
       }

       real_t wt1[2], wt2[2];
       wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
       wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
       wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
       wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
       buf_zm[IDXBF_R(ic, 0, ixyt)] = wt1[0];
       buf_zm[IDXBF_I(ic, 0, ixyt)] = wt1[1];
       buf_zm[IDXBF_R(ic, 1, ixyt)] = wt2[0];
       buf_zm[IDXBF_I(ic, 1, ixyt)] = wt2[1];

    }

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
    real_t bc2 = 2.0 * bc[3];
    for(int ic = 0; ic < NC; ++ic){
      buf_tp[IDXBF_R(ic, 0, ixyz)] = bc2 * v1[IDX2_SP_R(ic, 2, ist)];
      buf_tp[IDXBF_I(ic, 0, ixyz)] = bc2 * v1[IDX2_SP_I(ic, 2, ist)];
      buf_tp[IDXBF_R(ic, 1, ixyz)] = bc2 * v1[IDX2_SP_R(ic, 3, ist)];
      buf_tp[IDXBF_I(ic, 1, ixyz)] = bc2 * v1[IDX2_SP_I(ic, 3, ist)];
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

    real_t vt1[NVC], vt2[NVC];
    for(int ic = 0; ic < NC; ++ic){
      int icr = 2*ic;
      int ici = 2*ic + 1;
      vt1[icr] = 2.0 * v1[IDX2_SP_R(ic, 0, ist)];
      vt1[ici] = 2.0 * v1[IDX2_SP_I(ic, 0, ist)];
      vt2[icr] = 2.0 * v1[IDX2_SP_R(ic, 1, ist)];
      vt2[ici] = 2.0 * v1[IDX2_SP_I(ic, 1, ist)];
    }

    for(int ic = 0; ic < NC; ++ic){

      real_t ut[NVC];
      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] =   u[IDX2_G_R(ic, ic2, istu)];
        ut[2*ic2+1] = - u[IDX2_G_I(ic, ic2, istu)];
      }

      real_t wt1[2], wt2[2];
      wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                        vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                        vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                        vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
      wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                        vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
      buf_tm[IDXBF_R(ic, 0, ixyz)] = wt1[0];
      buf_tm[IDXBF_I(ic, 0, ixyz)] = wt1[1];
      buf_tm[IDXBF_R(ic, 1, ixyz)] = wt2[0];
      buf_tm[IDXBF_I(ic, 1, ixyz)] = wt2[1];
    }

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
void mult_wilson_eo_1_chiral(
                    real_t *restrict buf_xp, real_t *restrict buf_xm,
                    real_t *restrict buf_yp, real_t *restrict buf_ym,
                    real_t *restrict buf_zp, real_t *restrict buf_zm,
                    real_t *restrict buf_tp, real_t *restrict buf_tm,
                    real_t *restrict u, real_t *restrict v1,
                    const int ieo, const int jeo,
                    int *Nsize, int *bc, int *do_comm, int Nc)
{
  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nst  = Nx * Ny * Nz * Nt;

  int size    = NVC * ND * Nst;
  int size_u  = NDF * Nst * 4;
  int size_bx = NVC * 2 * ((Ny * Nz * Nt + 1)/2);
  int size_by = NVC * 2 * Nx * Nz * Nt;
  int size_bz = NVC * 2 * Nx * Ny * Nt;
  int size_bt = NVC * 2 * Nx * Ny * Nz;

#pragma acc data present(buf_xp[0:size_bx], buf_xm[0:size_bx], \
                         buf_yp[0:size_by], buf_ym[0:size_by], \
                         buf_zp[0:size_bz], buf_zm[0:size_bz], \
                         buf_tp[0:size_bt], buf_tm[0:size_bt], \
                         u[0:size_u], v1[0:size]) \
                 copyin(bc[0:4], do_comm[0:4], ieo, jeo, Nx, Ny, Nz, Nt, Nst)
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
      for(int ic = 0; ic < NC; ++ic){
        real_t vt1[2], vt2[2];
        vt1[0] = v1[IDX2_SP_R(ic, 0, ist)] - v1[IDX2_SP_I(ic, 3, ist)];
        vt1[1] = v1[IDX2_SP_I(ic, 0, ist)] + v1[IDX2_SP_R(ic, 3, ist)];
        vt2[0] = v1[IDX2_SP_R(ic, 1, ist)] - v1[IDX2_SP_I(ic, 2, ist)];
        vt2[1] = v1[IDX2_SP_I(ic, 1, ist)] + v1[IDX2_SP_R(ic, 2, ist)];
        buf_xp[IDXBF_R(ic, 0, iyzt2)] = bc2 * vt1[0];
        buf_xp[IDXBF_I(ic, 0, iyzt2)] = bc2 * vt1[1];
        buf_xp[IDXBF_R(ic, 1, iyzt2)] = bc2 * vt2[0];
        buf_xp[IDXBF_I(ic, 1, iyzt2)] = bc2 * vt2[1];
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
      int iyzt2 = iyzt/2;

      real_t vt1[NVC], vt2[NVC], ut[NDF];
      for(int ic = 0; ic < NC; ++ic){
        int icr = 2*ic;
        int ici = 2*ic + 1;
        vt1[icr] = v1[IDX2_SP_R(ic, 0, ist)] + v1[IDX2_SP_I(ic, 3, ist)];
        vt1[ici] = v1[IDX2_SP_I(ic, 0, ist)] - v1[IDX2_SP_R(ic, 3, ist)];
        vt2[icr] = v1[IDX2_SP_R(ic, 1, ist)] + v1[IDX2_SP_I(ic, 2, ist)];
        vt2[ici] = v1[IDX2_SP_I(ic, 1, ist)] - v1[IDX2_SP_R(ic, 2, ist)];
      }

      real_t wt1[2], wt2[2];

      for(int ic = 0; ic < NC; ++ic){

        for(int ic2 = 0; ic2 < NC; ++ic2){
          ut[2*ic2  ] =   u[IDX2_G_R(ic, ic2, istu)];
          ut[2*ic2+1] = - u[IDX2_G_I(ic, ic2, istu)];
        }

        wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                          vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
        wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                          vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
        wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                          vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
        wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                          vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

        buf_xm[IDXBF_R(ic, 0, iyzt2)] = wt1[0];
        buf_xm[IDXBF_I(ic, 0, iyzt2)] = wt1[1];
        buf_xm[IDXBF_R(ic, 1, iyzt2)] = wt2[0];
        buf_xm[IDXBF_I(ic, 1, iyzt2)] = wt2[1];
      }

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

    for(int ic = 0; ic < NC; ++ic){
      real_t vt1[2], vt2[2];
      vt1[0] = v1[IDX2_SP_R(ic, 0, ist)] + v1[IDX2_SP_R(ic, 3, ist)];
      vt1[1] = v1[IDX2_SP_I(ic, 0, ist)] + v1[IDX2_SP_I(ic, 3, ist)];
      vt2[0] = v1[IDX2_SP_R(ic, 1, ist)] - v1[IDX2_SP_R(ic, 2, ist)];
      vt2[1] = v1[IDX2_SP_I(ic, 1, ist)] - v1[IDX2_SP_I(ic, 2, ist)];
      buf_yp[IDXBF_R(ic, 0, ixzt)] = bc2 * vt1[0];
      buf_yp[IDXBF_I(ic, 0, ixzt)] = bc2 * vt1[1];
      buf_yp[IDXBF_R(ic, 1, ixzt)] = bc2 * vt2[0];
      buf_yp[IDXBF_I(ic, 1, ixzt)] = bc2 * vt2[1];
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
     int ixzt = ix + Nx * izt;

     real_t vt1[NVC], vt2[NVC];
     for(int ic = 0; ic < NC; ++ic){
       int icr = 2*ic;
       int ici = 2*ic + 1;
       vt1[icr] = v1[IDX2_SP_R(ic, 0, ist)] - v1[IDX2_SP_R(ic, 3, ist)];
       vt1[ici] = v1[IDX2_SP_I(ic, 0, ist)] - v1[IDX2_SP_I(ic, 3, ist)];
       vt2[icr] = v1[IDX2_SP_R(ic, 1, ist)] + v1[IDX2_SP_R(ic, 2, ist)];
       vt2[ici] = v1[IDX2_SP_I(ic, 1, ist)] + v1[IDX2_SP_I(ic, 2, ist)];
     }

     for(int ic = 0; ic < NC; ++ic){

       real_t ut[NVC];
       for(int ic2 = 0; ic2 < NC; ++ic2){
         ut[2*ic2  ] =   u[IDX2_G_R(ic, ic2, istu)];
         ut[2*ic2+1] = - u[IDX2_G_I(ic, ic2, istu)];
       }

       real_t wt1[2], wt2[2];
       wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
       wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
       wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
       wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

       buf_ym[IDXBF_R(ic, 0, ixzt)] = wt1[0];
       buf_ym[IDXBF_I(ic, 0, ixzt)] = wt1[1];
       buf_ym[IDXBF_R(ic, 1, ixzt)] = wt2[0];
       buf_ym[IDXBF_I(ic, 1, ixzt)] = wt2[1];

     }

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

     for(int ic = 0; ic < NC; ++ic){
       real_t wt1[2], wt2[2];
       wt1[0] = v1[IDX2_SP_R(ic, 0, ist)] - v1[IDX2_SP_I(ic, 2, ist)];
       wt1[1] = v1[IDX2_SP_I(ic, 0, ist)] + v1[IDX2_SP_R(ic, 2, ist)];
       wt2[0] = v1[IDX2_SP_R(ic, 1, ist)] + v1[IDX2_SP_I(ic, 3, ist)];
       wt2[1] = v1[IDX2_SP_I(ic, 1, ist)] - v1[IDX2_SP_R(ic, 3, ist)];
       buf_zp[IDXBF_R(ic, 0, ixyt)] = bc2 * wt1[0];
       buf_zp[IDXBF_I(ic, 0, ixyt)] = bc2 * wt1[1];
       buf_zp[IDXBF_R(ic, 1, ixyt)] = bc2 * wt2[0];
       buf_zp[IDXBF_I(ic, 1, ixyt)] = bc2 * wt2[1];
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
     int ixyt = ixy + Nxy * it;

     real_t vt1[NVC], vt2[NVC];
     for(int ic = 0; ic < NC; ++ic){
       int icr = 2*ic;
       int ici = 2*ic + 1;
       vt1[icr] = v1[IDX2_SP_R(ic, 0, ist)] + v1[IDX2_SP_I(ic, 2, ist)];
       vt1[ici] = v1[IDX2_SP_I(ic, 0, ist)] - v1[IDX2_SP_R(ic, 2, ist)];
       vt2[icr] = v1[IDX2_SP_R(ic, 1, ist)] - v1[IDX2_SP_I(ic, 3, ist)];
       vt2[ici] = v1[IDX2_SP_I(ic, 1, ist)] + v1[IDX2_SP_R(ic, 3, ist)];
     }

     for(int ic = 0; ic < NC; ++ic){

       real_t ut[NVC];
       for(int ic2 = 0; ic2 < NC; ++ic2){
         ut[2*ic2  ] =   u[IDX2_G_R(ic, ic2, istu)];
         ut[2*ic2+1] = - u[IDX2_G_I(ic, ic2, istu)];
       }

       real_t wt1[2], wt2[2];
       wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
       wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
       wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
       wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
       buf_zm[IDXBF_R(ic, 0, ixyt)] = wt1[0];
       buf_zm[IDXBF_I(ic, 0, ixyt)] = wt1[1];
       buf_zm[IDXBF_R(ic, 1, ixyt)] = wt2[0];
       buf_zm[IDXBF_I(ic, 1, ixyt)] = wt2[1];

    }

   }
  }

 }

 } // do_comm[2]

 //  idir = 3;
 if(do_comm[3] > 0){

#pragma acc parallel async\
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
   int Nxyz = Nx * Ny * Nz;

#pragma acc loop
  for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
    int it = 0;
    int ist  = ixyz + Nxyz * it;
    real_t bc2 = bc[3];
    for(int ivc = 0; ivc < NVC; ++ivc){
      buf_tp[IDXBF(ivc, 0, ixyz)]
       = bc2 * (v1[IDX2_SP(ivc, 0, ist)] + v1[IDX2_SP(ivc, 2, ist)]);
      buf_tp[IDXBF(ivc, 1, ixyz)]
       = bc2 * (v1[IDX2_SP(ivc, 1, ist)] + v1[IDX2_SP(ivc, 3, ist)]);
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

    real_t vt1[NVC], vt2[NVC];
    for(int ivc = 0; ivc < NVC; ++ivc){
      vt1[ivc] = v1[IDX2_SP(ivc, 0, ist)] - v1[IDX2_SP(ivc, 2, ist)];
      vt2[ivc] = v1[IDX2_SP(ivc, 1, ist)] - v1[IDX2_SP(ivc, 3, ist)];
    }

    for(int ic = 0; ic < NC; ++ic){

      real_t ut[NVC];
      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] =   u[IDX2_G_R(ic, ic2, istu)];
        ut[2*ic2+1] = - u[IDX2_G_I(ic, ic2, istu)];
      }

      real_t wt1[2], wt2[2];
      wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                        vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                        vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                        vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
      wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                        vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
      buf_tm[IDXBF_R(ic, 0, ixyz)] = wt1[0];
      buf_tm[IDXBF_I(ic, 0, ixyz)] = wt1[1];
      buf_tm[IDXBF_R(ic, 1, ixyz)] = wt2[0];
      buf_tm[IDXBF_I(ic, 1, ixyz)] = wt2[1];

    }

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
void mult_wilson_eo_2_dirac(real_t *restrict v2, real_t *restrict u, 
                    real_t *restrict buf_xp, real_t *restrict buf_xm,
                    real_t *restrict buf_yp, real_t *restrict buf_ym,
                    real_t *restrict buf_zp, real_t *restrict buf_zm,
                    real_t *restrict buf_tp, real_t *restrict buf_tm,
                    real_t kappa, const int ieo, const int jeo,
                    const int iflag,
                    int *Nsize, int *bc, int *do_comm, int Nc)
{
  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nst  = Nx * Ny * Nz * Nt;

  int size    = NVC * ND * Nst;
  int size_u  = NDF * Nst * 4;
  int size_bx = NVC * 2 * ((Ny * Nz * Nt + 1)/2);
  int size_by = NVC * 2 * Nx * Nz * Nt;
  int size_bz = NVC * 2 * Nx * Ny * Nt;
  int size_bt = NVC * 2 * Nx * Ny * Nz;

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
                 copyin(bc[0:4], do_comm[0:4], ieo, jeo, kappa, Nx, Ny, Nz, Nt)
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

    real_t v2L[NVC*ND];
    for(int id = 0; id < ND; ++id){
      for(int ivc = 0; ivc < NVC; ++ivc){
        v2L[ivc + NVC * id] = 0.0;
      }
    }

#include "inc/mult_Wilson_eo_openacc2_2xyz-inc.h"

    // idir = 3
    if(do_comm[3] > 0){

      if(it == Nt-1){
        int ixyz = ix + Nx * (iy + Ny * iz);
        int istu = ist + Nst * (ieo + 2*3);

        real_t vt1[NVC], vt2[NVC];

        for(int ivc = 0; ivc < NVC; ++ivc){
          vt1[ivc] = buf_tp[IDXBF(ivc, 0, ixyz)];
          vt2[ivc] = buf_tp[IDXBF(ivc, 1, ixyz)];
        }

        for(int ic = 0; ic < NC; ++ic){

          real_t ut[NVC];
          for(int ic2 = 0; ic2 < NC; ++ic2){
            ut[2*ic2  ] = u[IDX2_G_R(ic2, ic, istu)];
            ut[2*ic2+1] = u[IDX2_G_I(ic2, ic, istu)];
          }

          real_t wt1[2], wt2[2];
          wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                            vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
          wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                            vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
          wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                            vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
          wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                            vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
          v2L[ic*2   + ID3] += wt1[0];
          v2L[ic*2+1 + ID3] += wt1[1];
          v2L[ic*2   + ID4] += wt2[0];
          v2L[ic*2+1 + ID4] += wt2[1];
	}

      }

      if(it == 0){
        int ixyz = ix + Nx * (iy + Ny * iz);
        real_t bc2 = bc[3];
        real_t wt1[2], wt2[2];

        for(int ic = 0; ic < NC; ++ic){
          wt1[0] = bc2 * buf_tm[IDXBF_R(ic, 0, ixyz)];
          wt1[1] = bc2 * buf_tm[IDXBF_I(ic, 0, ixyz)];
          wt2[0] = bc2 * buf_tm[IDXBF_R(ic, 1, ixyz)];
          wt2[1] = bc2 * buf_tm[IDXBF_I(ic, 1, ixyz)];
          v2L[ic*2   + ID1] += wt1[0];
          v2L[ic*2+1 + ID1] += wt1[1];
          v2L[ic*2   + ID2] += wt2[0];
          v2L[ic*2+1 + ID2] += wt2[1];
	}
      }

    }

    //axpy
    if(iflag == 0){
      for(int id = 0; id < ND; ++id){
        for(int ic = 0; ic < NC; ++ic){
          v2[IDX2_SP_R(ic, id, ist)] += -kappa * v2L[ic*2   + NVC*id];
          v2[IDX2_SP_I(ic, id, ist)] += -kappa * v2L[ic*2+1 + NVC*id];
        }
      }
    }else{
      for(int id = 0; id < ND; ++id){
        for(int ic = 0; ic < NC; ++ic){
          v2[IDX2_SP_R(ic, id, ist)] += kappa * v2L[ic*2   + NVC*id];
          v2[IDX2_SP_I(ic, id, ist)] += kappa * v2L[ic*2+1 + NVC*id];
        }
      }
    }

  }

 }

 }

}

//==================================================================== 
void mult_wilson_eo_2_chiral(real_t *restrict v2, real_t *restrict u,
                    real_t *restrict buf_xp, real_t *restrict buf_xm,
                    real_t *restrict buf_yp, real_t *restrict buf_ym,
                    real_t *restrict buf_zp, real_t *restrict buf_zm,
                    real_t *restrict buf_tp, real_t *restrict buf_tm,
                    real_t kappa, const int ieo, const int jeo,
                    const int iflag,
                    int *Nsize, int *bc, int *do_comm, int Nc)
{
  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nst  = Nx * Ny * Nz * Nt;

  int size    = NVC * ND * Nst;
  int size_u  = NDF * Nst * 4;
  int size_bx = NVC * 2 * ((Ny * Nz * Nt + 1)/2);
  int size_by = NVC * 2 * Nx * Nz * Nt;
  int size_bz = NVC * 2 * Nx * Ny * Nt;
  int size_bt = NVC * 2 * Nx * Ny * Nz;

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
                 copyin(bc[0:4], do_comm[0:4], ieo, jeo, kappa, Nx, Ny, Nz, Nt)
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

    real_t v2L[NVC*ND];
    for(int id = 0; id < ND; ++id){
      for(int ivc = 0; ivc < NVC; ++ivc){
        v2L[ivc + NVC * id] = 0.0;
      }
    }

#include "inc/mult_Wilson_eo_openacc2_2xyz-inc.h"

    // idir = 3
    if(do_comm[3] > 0){

      if(it == Nt-1){
        int ixyz = ix + Nx * (iy + Ny * iz);
        int istu = ist + Nst * (ieo + 2*3);

        real_t vt1[NVC], vt2[NVC], ut[NVC];
        real_t wt1[2], wt2[2];

        for(int ivc = 0; ivc < NVC; ++ivc){
          vt1[ivc] = buf_tp[IDXBF(ivc, 0, ixyz)];
          vt2[ivc] = buf_tp[IDXBF(ivc, 1, ixyz)];
        }

        for(int ic = 0; ic < NC; ++ic){

          for(int ic2 = 0; ic2 < NC; ++ic2){
            ut[2*ic2  ] = u[IDX2_G_R(ic2, ic, istu)];
            ut[2*ic2+1] = u[IDX2_G_I(ic2, ic, istu)];
          }

          wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                            vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
          wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                            vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
          wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                            vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
          wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                            vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

          v2L[ic*2   + ID1] += wt1[0];
          v2L[ic*2+1 + ID1] += wt1[1];
          v2L[ic*2   + ID2] += wt2[0];
          v2L[ic*2+1 + ID2] += wt2[1];
          v2L[ic*2   + ID3] += wt1[0];
          v2L[ic*2+1 + ID3] += wt1[1];
          v2L[ic*2   + ID4] += wt2[0];
          v2L[ic*2+1 + ID4] += wt2[1];
	}
      }

      if(it == 0){
        int ixyz = ix + Nx * (iy + Ny * iz);
        real_t bc2 = bc[3];
        real_t wt1[2], wt2[2];

        for(int ic = 0; ic < NC; ++ic){
          wt1[0] = bc2 * buf_tm[IDXBF_R(ic, 0, ixyz)];
          wt1[1] = bc2 * buf_tm[IDXBF_I(ic, 0, ixyz)];
          wt2[0] = bc2 * buf_tm[IDXBF_R(ic, 1, ixyz)];
          wt2[1] = bc2 * buf_tm[IDXBF_I(ic, 1, ixyz)];
          v2L[ic*2   + ID1] += wt1[0];
          v2L[ic*2+1 + ID1] += wt1[1];
          v2L[ic*2   + ID2] += wt2[0];
          v2L[ic*2+1 + ID2] += wt2[1];
          v2L[ic*2   + ID3] -= wt1[0];
          v2L[ic*2+1 + ID3] -= wt1[1];
          v2L[ic*2   + ID4] -= wt2[0];
          v2L[ic*2+1 + ID4] -= wt2[1];
	}
      }

    }

    //axpy
    if(iflag == 0){
      for(int id = 0; id < ND; ++id){
        for(int ic = 0; ic < NC; ++ic){
          v2[IDX2_SP_R(ic, id, ist)] += -kappa * v2L[ic*2   + NVC*id];
          v2[IDX2_SP_I(ic, id, ist)] += -kappa * v2L[ic*2+1 + NVC*id];
        }
      }
    }else{
      for(int id = 0; id < ND; ++id){
        for(int ic = 0; ic < NC; ++ic){
          v2[IDX2_SP_R(ic, id, ist)] += kappa * v2L[ic*2   + NVC*id];
          v2[IDX2_SP_I(ic, id, ist)] += kappa * v2L[ic*2+1 + NVC*id];
        }
      }
    }

  }

 }

 }

}

//============================================================END=====
