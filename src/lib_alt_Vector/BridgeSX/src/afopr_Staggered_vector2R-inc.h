/*!
      @file    afopr_Staggered_vector2-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2543 $
*/

#include "lib_alt_Vector/inline/define_params.h"
#include "lib_alt_Vector/inline/define_index.h"

//====================================================================
void mult_staggered_D(real_t *RESTRICT v2, real_t *RESTRICT u,
                      real_t *RESTRICT v1,
                      int *Nsize, int *bc, real_t mq, int jdag)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst  = Nx * Ny * Nz * Nt;

  real_t fac = real_t(jdag) * 0.5/mq;

  int nst2 = Nst/VLEN;

  real_t vt0[VLEN], vt1[VLEN], vt2[VLEN], vt3[VLEN], vt4[VLEN], vt5[VLEN];
  real_t ut0[VLEN], ut1[VLEN], ut2[VLEN], ut3[VLEN], ut4[VLEN], ut5[VLEN];
  real_t ur0[VLEN], ur1[VLEN], ur2[VLEN], ur3[VLEN], ur4[VLEN], ur5[VLEN];
  real_t us0[VLEN], us1[VLEN], us2[VLEN], us3[VLEN], us4[VLEN], us5[VLEN];
  real_t wt0[VLEN], wt1[VLEN], wt2[VLEN], wt3[VLEN], wt4[VLEN], wt5[VLEN];

  for(int ist2 = 0; ist2 < nst2; ++ist2){

#pragma _NEC ivdep
    for(int j = 0; j < VLEN; ++j){
      int ist = j + VLEN * ist2;
#pragma _NEC unroll(NVC)
      for(int ivc = 0; ivc < NVC; ++ivc){
        v2[IDXV_1SP(Nst, ivc, ist)] = v1[IDXV_1SP(Nst, ivc, ist)];
      }
    }

    int idir = 0;

#pragma _NEC ivdep
    for(int j = 0; j < VLEN; ++j){
      int ist = j + VLEN * ist2;
      int ix   = ist % Nx;
      int iyzt = ist/Nx;
      int nei  = ((ix+1) % Nx) + Nx * iyzt;
      vt0[j] = v1[IDXV_1SP(Nst, 0, nei)];
      vt1[j] = v1[IDXV_1SP(Nst, 1, nei)];
      vt2[j] = v1[IDXV_1SP(Nst, 2, nei)];
      vt3[j] = v1[IDXV_1SP(Nst, 3, nei)];
      vt4[j] = v1[IDXV_1SP(Nst, 4, nei)];
      vt5[j] = v1[IDXV_1SP(Nst, 5, nei)];
    }

    int ic = 0;
#pragma _NEC ivdep
    for(int j = 0; j < VLEN; ++j){
      int ist = j + VLEN * ist2;
      ut0[j] = u[IDXV_G_R(Nst, 0, ic, ist, idir)];
      ut1[j] = u[IDXV_G_I(Nst, 0, ic, ist, idir)];
      ut2[j] = u[IDXV_G_R(Nst, 1, ic, ist, idir)];
      ut3[j] = u[IDXV_G_I(Nst, 1, ic, ist, idir)];
      ut4[j] = u[IDXV_G_R(Nst, 2, ic, ist, idir)];
      ut5[j] = u[IDXV_G_I(Nst, 2, ic, ist, idir)];
    }
#pragma _NEC ivdep
    for(int j = 0; j < VLEN; ++j){
      wt0[j] =  ut0[j] * vt0[j] - ut1[j] * vt1[j]
              + ut2[j] * vt2[j] - ut3[j] * vt3[j]
              + ut4[j] * vt4[j] - ut5[j] * vt5[j];
      wt1[j] =  ut0[j] * vt1[j] + ut1[j] * vt0[j]
              + ut2[j] * vt3[j] + ut3[j] * vt2[j]
              + ut4[j] * vt5[j] + ut5[j] * vt4[j];
    }

    ic = 1;
#pragma _NEC ivdep
    for(int j = 0; j < VLEN; ++j){
      int ist = j + VLEN * ist2;
      ur0[j] = u[IDXV_G_R(Nst, 0, ic, ist, idir)];
      ur1[j] = u[IDXV_G_I(Nst, 0, ic, ist, idir)];
      ur2[j] = u[IDXV_G_R(Nst, 1, ic, ist, idir)];
      ur3[j] = u[IDXV_G_I(Nst, 1, ic, ist, idir)];
      ur4[j] = u[IDXV_G_R(Nst, 2, ic, ist, idir)];
      ur5[j] = u[IDXV_G_I(Nst, 2, ic, ist, idir)];
    }
#pragma _NEC ivdep
    for(int j = 0; j < VLEN; ++j){
      wt2[j] =  ur0[j] * vt0[j] - ur1[j] * vt1[j]
              + ur2[j] * vt2[j] - ur3[j] * vt3[j]
              + ur4[j] * vt4[j] - ur5[j] * vt5[j];
      wt3[j] =  ur0[j] * vt1[j] + ur1[j] * vt0[j]
              + ur2[j] * vt3[j] + ur3[j] * vt2[j]
              + ur4[j] * vt5[j] + ur5[j] * vt4[j];
    }

    ic = 2;
#pragma _NEC ivdep
    for(int j = 0; j < VLEN; ++j){
      int ist = j + VLEN * ist2;
      us0[j] = u[IDXV_G_R(Nst, 0, ic, ist, idir)];
      us1[j] = u[IDXV_G_I(Nst, 0, ic, ist, idir)];
      us2[j] = u[IDXV_G_R(Nst, 1, ic, ist, idir)];
      us3[j] = u[IDXV_G_I(Nst, 1, ic, ist, idir)];
      us4[j] = u[IDXV_G_R(Nst, 2, ic, ist, idir)];
      us5[j] = u[IDXV_G_I(Nst, 2, ic, ist, idir)];
    }
#pragma _NEC ivdep
    for(int j = 0; j < VLEN; ++j){
      wt4[j] =  us0[j] * vt0[j] - us1[j] * vt1[j]
              + us2[j] * vt2[j] - us3[j] * vt3[j]
              + us4[j] * vt4[j] - us5[j] * vt5[j];
      wt5[j] =  us0[j] * vt1[j] + us1[j] * vt0[j]
              + us2[j] * vt3[j] + us3[j] * vt2[j]
              + us4[j] * vt5[j] + us5[j] * vt4[j];
    }

#pragma _NEC ivdep
    for(int j = 0; j < VLEN; ++j){
      int ist = j + VLEN * ist2;
      int ix   = ist % Nx;
      real_t bc2 = fac * (ix < Nx-1 ? 1.0 : bc[idir]);
      v2[IDXV_1SP(Nst, 0, ist)] += bc2 * wt0[j];
      v2[IDXV_1SP(Nst, 1, ist)] += bc2 * wt1[j];
      v2[IDXV_1SP(Nst, 2, ist)] += bc2 * wt2[j];
      v2[IDXV_1SP(Nst, 3, ist)] += bc2 * wt3[j];
      v2[IDXV_1SP(Nst, 4, ist)] += bc2 * wt4[j];
      v2[IDXV_1SP(Nst, 5, ist)] += bc2 * wt5[j];
    }

    //qqqqqq
    /*
    for(int j = 0; j < VLEN; ++j){
      int ist = j + VLEN * ist2;
      int ix   = ist % Nx;
      int iyzt = ist/Nx;
      int nei = ((ix+1) % Nx) + Nx * iyzt;
      real_t bc2 = fac * (ix < Nx-1 ? 1.0 : bc[idir]);
      for(int ic = 0; ic < NC; ++ic){
        v2[IDXV_1SP_R(Nst, ic, ist)] += bc2 *
         (  u[IDXV_G_R(Nst, 0, ic, ist, idir)] * v1[IDXV_1SP_R(Nst, 0, nei)]
          - u[IDXV_G_I(Nst, 0, ic, ist, idir)] * v1[IDXV_1SP_I(Nst, 0, nei)]
          + u[IDXV_G_R(Nst, 1, ic, ist, idir)] * v1[IDXV_1SP_R(Nst, 1, nei)]
          - u[IDXV_G_I(Nst, 1, ic, ist, idir)] * v1[IDXV_1SP_I(Nst, 1, nei)]
          + u[IDXV_G_R(Nst, 2, ic, ist, idir)] * v1[IDXV_1SP_R(Nst, 2, nei)]
          - u[IDXV_G_I(Nst, 2, ic, ist, idir)] * v1[IDXV_1SP_I(Nst, 2, nei)]);
        v2[IDXV_1SP_I(Nst, ic, ist)] += bc2 *
         (  u[IDXV_G_R(Nst, 0, ic, ist, idir)] * v1[IDXV_1SP_I(Nst, 0, nei)]
          + u[IDXV_G_I(Nst, 0, ic, ist, idir)] * v1[IDXV_1SP_R(Nst, 0, nei)]
          + u[IDXV_G_R(Nst, 1, ic, ist, idir)] * v1[IDXV_1SP_I(Nst, 1, nei)]
          + u[IDXV_G_I(Nst, 1, ic, ist, idir)] * v1[IDXV_1SP_R(Nst, 1, nei)]
          + u[IDXV_G_R(Nst, 2, ic, ist, idir)] * v1[IDXV_1SP_I(Nst, 2, nei)]
          + u[IDXV_G_I(Nst, 2, ic, ist, idir)] * v1[IDXV_1SP_R(Nst, 2, nei)]);
      }

    }
    */

#pragma _NEC ivdep
    for(int j = 0; j < VLEN; ++j){
      int ist = j + VLEN * ist2;
      int ix   = ist % Nx;
      int iyzt = ist/Nx;
      int nei = ((ix-1+Nx) % Nx) + Nx * iyzt;
      real_t bc2 = fac * (ix > 0 ? 1.0 : bc[idir]);

      for(int ic = 0; ic < NC; ++ic){
        v2[IDXV_1SP_R(Nst, ic, ist)] -= bc2 *
         (  u[IDXV_G_R(Nst, ic, 0, nei, idir)] * v1[IDXV_1SP_R(Nst, 0, nei)]
          + u[IDXV_G_I(Nst, ic, 0, nei, idir)] * v1[IDXV_1SP_I(Nst, 0, nei)]
          + u[IDXV_G_R(Nst, ic, 1, nei, idir)] * v1[IDXV_1SP_R(Nst, 1, nei)]
          + u[IDXV_G_I(Nst, ic, 1, nei, idir)] * v1[IDXV_1SP_I(Nst, 1, nei)]
          + u[IDXV_G_R(Nst, ic, 2, nei, idir)] * v1[IDXV_1SP_R(Nst, 2, nei)]
          + u[IDXV_G_I(Nst, ic, 2, nei, idir)] * v1[IDXV_1SP_I(Nst, 2, nei)]);
        v2[IDXV_1SP_I(Nst, ic, ist)] -= bc2 *
         (  u[IDXV_G_R(Nst, ic, 0, nei, idir)] * v1[IDXV_1SP_I(Nst, 0, nei)]
          - u[IDXV_G_I(Nst, ic, 0, nei, idir)] * v1[IDXV_1SP_R(Nst, 0, nei)]
          + u[IDXV_G_R(Nst, ic, 1, nei, idir)] * v1[IDXV_1SP_I(Nst, 1, nei)]
          - u[IDXV_G_I(Nst, ic, 1, nei, idir)] * v1[IDXV_1SP_R(Nst, 1, nei)]
          + u[IDXV_G_R(Nst, ic, 2, nei, idir)] * v1[IDXV_1SP_I(Nst, 2, nei)]
          - u[IDXV_G_I(Nst, ic, 2, nei, idir)] * v1[IDXV_1SP_R(Nst, 2, nei)]);
      }

    }

    idir = 1;

#pragma _NEC ivdep
    for(int j = 0; j < VLEN; ++j){
      int ist = j + VLEN * ist2;

      int ix   = ist % Nx;
      int iyzt = ist/Nx;
      int iy  = iyzt % Ny;
      int izt = iyzt/Ny;
      int iy2 = (iy+1) % Ny;
      int nei = ix + Nx * (iy2 + Ny * izt);
      real_t bc2 = fac * (iy < Ny-1 ? 1.0 : bc[idir]);

      for(int ic = 0; ic < NC; ++ic){
        v2[IDXV_1SP_R(Nst, ic, ist)] += bc2 *
         (  u[IDXV_G_R(Nst, 0, ic, ist, idir)] * v1[IDXV_1SP_R(Nst, 0, nei)]
          - u[IDXV_G_I(Nst, 0, ic, ist, idir)] * v1[IDXV_1SP_I(Nst, 0, nei)]
          + u[IDXV_G_R(Nst, 1, ic, ist, idir)] * v1[IDXV_1SP_R(Nst, 1, nei)]
          - u[IDXV_G_I(Nst, 1, ic, ist, idir)] * v1[IDXV_1SP_I(Nst, 1, nei)]
          + u[IDXV_G_R(Nst, 2, ic, ist, idir)] * v1[IDXV_1SP_R(Nst, 2, nei)]
          - u[IDXV_G_I(Nst, 2, ic, ist, idir)] * v1[IDXV_1SP_I(Nst, 2, nei)]);
        v2[IDXV_1SP_I(Nst, ic, ist)] += bc2 *
         (  u[IDXV_G_R(Nst, 0, ic, ist, idir)] * v1[IDXV_1SP_I(Nst, 0, nei)]
          + u[IDXV_G_I(Nst, 0, ic, ist, idir)] * v1[IDXV_1SP_R(Nst, 0, nei)]
          + u[IDXV_G_R(Nst, 1, ic, ist, idir)] * v1[IDXV_1SP_I(Nst, 1, nei)]
          + u[IDXV_G_I(Nst, 1, ic, ist, idir)] * v1[IDXV_1SP_R(Nst, 1, nei)]
          + u[IDXV_G_R(Nst, 2, ic, ist, idir)] * v1[IDXV_1SP_I(Nst, 2, nei)]
          + u[IDXV_G_I(Nst, 2, ic, ist, idir)] * v1[IDXV_1SP_R(Nst, 2, nei)]);
      }

    }

#pragma _NEC ivdep
    for(int j = 0; j < VLEN; ++j){
      int ist = j + VLEN * ist2;

      int ix   = ist % Nx;
      int iyzt = ist/Nx;
      int iy  = iyzt % Ny;
      int izt = iyzt/Ny;

      int iy2 = (iy-1+Ny) % Ny;
      int nei = ix + Nx * (iy2 + Ny * izt);
      int bc2 = fac * (iy > 0 ? 1.0 : bc[idir]);

      for(int ic = 0; ic < NC; ++ic){
        v2[IDXV_1SP_R(Nst, ic, ist)] -= bc2 *
         (  u[IDXV_G_R(Nst, ic, 0, nei, idir)] * v1[IDXV_1SP_R(Nst, 0, nei)]
          + u[IDXV_G_I(Nst, ic, 0, nei, idir)] * v1[IDXV_1SP_I(Nst, 0, nei)]
          + u[IDXV_G_R(Nst, ic, 1, nei, idir)] * v1[IDXV_1SP_R(Nst, 1, nei)]
          + u[IDXV_G_I(Nst, ic, 1, nei, idir)] * v1[IDXV_1SP_I(Nst, 1, nei)]
          + u[IDXV_G_R(Nst, ic, 2, nei, idir)] * v1[IDXV_1SP_R(Nst, 2, nei)]
          + u[IDXV_G_I(Nst, ic, 2, nei, idir)] * v1[IDXV_1SP_I(Nst, 2, nei)]);
        v2[IDXV_1SP_I(Nst, ic, ist)] -= bc2 *
         (  u[IDXV_G_R(Nst, ic, 0, nei, idir)] * v1[IDXV_1SP_I(Nst, 0, nei)]
          - u[IDXV_G_I(Nst, ic, 0, nei, idir)] * v1[IDXV_1SP_R(Nst, 0, nei)]
          + u[IDXV_G_R(Nst, ic, 1, nei, idir)] * v1[IDXV_1SP_I(Nst, 1, nei)]
          - u[IDXV_G_I(Nst, ic, 1, nei, idir)] * v1[IDXV_1SP_R(Nst, 1, nei)]
          + u[IDXV_G_R(Nst, ic, 2, nei, idir)] * v1[IDXV_1SP_I(Nst, 2, nei)]
          - u[IDXV_G_I(Nst, ic, 2, nei, idir)] * v1[IDXV_1SP_R(Nst, 2, nei)]);
      }

    }

#pragma _NEC ivdep
    for(int j = 0; j < VLEN; ++j){
      int ist = j + VLEN * ist2;

      int ix   = ist % Nx;
      int iyzt = ist/Nx;
      int iy  = iyzt % Ny;
      int izt = iyzt/Ny;

      int idir = 2;
      int iz  = izt % Nz;
      int it  = izt/Nz;

    int iz2 = (iz+1) % Nz;

    int nei = ix + Nx * (iy + Ny * (iz2 + Nz * it));
    real_t bc2 = fac * (iz < Nz-1 ? 1.0 : bc[idir]);

    for(int ic = 0; ic < NC; ++ic){
      v2[IDXV_1SP_R(Nst, ic, ist)] += bc2 *
       (  u[IDXV_G_R(Nst, 0, ic, ist, idir)] * v1[IDXV_1SP_R(Nst, 0, nei)]
        - u[IDXV_G_I(Nst, 0, ic, ist, idir)] * v1[IDXV_1SP_I(Nst, 0, nei)]
        + u[IDXV_G_R(Nst, 1, ic, ist, idir)] * v1[IDXV_1SP_R(Nst, 1, nei)]
        - u[IDXV_G_I(Nst, 1, ic, ist, idir)] * v1[IDXV_1SP_I(Nst, 1, nei)]
        + u[IDXV_G_R(Nst, 2, ic, ist, idir)] * v1[IDXV_1SP_R(Nst, 2, nei)]
        - u[IDXV_G_I(Nst, 2, ic, ist, idir)] * v1[IDXV_1SP_I(Nst, 2, nei)]);
      v2[IDXV_1SP_I(Nst, ic, ist)] += bc2 *
       (  u[IDXV_G_R(Nst, 0, ic, ist, idir)] * v1[IDXV_1SP_I(Nst, 0, nei)]
        + u[IDXV_G_I(Nst, 0, ic, ist, idir)] * v1[IDXV_1SP_R(Nst, 0, nei)]
        + u[IDXV_G_R(Nst, 1, ic, ist, idir)] * v1[IDXV_1SP_I(Nst, 1, nei)]
        + u[IDXV_G_I(Nst, 1, ic, ist, idir)] * v1[IDXV_1SP_R(Nst, 1, nei)]
        + u[IDXV_G_R(Nst, 2, ic, ist, idir)] * v1[IDXV_1SP_I(Nst, 2, nei)]
        + u[IDXV_G_I(Nst, 2, ic, ist, idir)] * v1[IDXV_1SP_R(Nst, 2, nei)]);
    }

    }

#pragma _NEC ivdep
    for(int j = 0; j < VLEN; ++j){
      int ist = j + VLEN * ist2;

      int ix   = ist % Nx;
      int iyzt = ist/Nx;
      int iy  = iyzt % Ny;
      int izt = iyzt/Ny;

      int idir = 2;
      int iz  = izt % Nz;
      int it  = izt/Nz;

      int iz2 = (iz-1+Nz) % Nz;
      int nei = ix + Nx * (iy + Ny * (iz2 + Nz * it));
      real_t bc2 = fac * (iz > 0 ? 1.0 : bc[idir]);

    for(int ic = 0; ic < NC; ++ic){
      v2[IDXV_1SP_R(Nst, ic, ist)] -= bc2 *
       (  u[IDXV_G_R(Nst, ic, 0, nei, idir)] * v1[IDXV_1SP_R(Nst, 0, nei)]
        + u[IDXV_G_I(Nst, ic, 0, nei, idir)] * v1[IDXV_1SP_I(Nst, 0, nei)]
        + u[IDXV_G_R(Nst, ic, 1, nei, idir)] * v1[IDXV_1SP_R(Nst, 1, nei)]
        + u[IDXV_G_I(Nst, ic, 1, nei, idir)] * v1[IDXV_1SP_I(Nst, 1, nei)]
        + u[IDXV_G_R(Nst, ic, 2, nei, idir)] * v1[IDXV_1SP_R(Nst, 2, nei)]
        + u[IDXV_G_I(Nst, ic, 2, nei, idir)] * v1[IDXV_1SP_I(Nst, 2, nei)]);
      v2[IDXV_1SP_I(Nst, ic, ist)] -= bc2 *
       (  u[IDXV_G_R(Nst, ic, 0, nei, idir)] * v1[IDXV_1SP_I(Nst, 0, nei)]
        - u[IDXV_G_I(Nst, ic, 0, nei, idir)] * v1[IDXV_1SP_R(Nst, 0, nei)]
        + u[IDXV_G_R(Nst, ic, 1, nei, idir)] * v1[IDXV_1SP_I(Nst, 1, nei)]
        - u[IDXV_G_I(Nst, ic, 1, nei, idir)] * v1[IDXV_1SP_R(Nst, 1, nei)]
        + u[IDXV_G_R(Nst, ic, 2, nei, idir)] * v1[IDXV_1SP_I(Nst, 2, nei)]
        - u[IDXV_G_I(Nst, ic, 2, nei, idir)] * v1[IDXV_1SP_R(Nst, 2, nei)]);
    }

   }

#pragma _NEC ivdep
   for(int j = 0; j < VLEN; ++j){
     int ist = j + VLEN * ist2;

    int ix   = ist % Nx;
    int iyzt = ist/Nx;
    int iy  = iyzt % Ny;
    int izt = iyzt/Ny;

    int iz  = izt % Nz;
    int it  = izt/Nz;

    int idir = 3;
    int it2 = (it+1) % Nt;
    int nei = ix + Nx * (iy + Ny * (iz + Nz * it2));
    real_t bc2 = fac * (it < Nt-1 ? 1.0 : bc[idir]);

    for(int ic = 0; ic < NC; ++ic){
      v2[IDXV_1SP_R(Nst, ic, ist)] += bc2 *
       (  u[IDXV_G_R(Nst, 0, ic, ist, idir)] * v1[IDXV_1SP_R(Nst, 0, nei)]
        - u[IDXV_G_I(Nst, 0, ic, ist, idir)] * v1[IDXV_1SP_I(Nst, 0, nei)]
        + u[IDXV_G_R(Nst, 1, ic, ist, idir)] * v1[IDXV_1SP_R(Nst, 1, nei)]
        - u[IDXV_G_I(Nst, 1, ic, ist, idir)] * v1[IDXV_1SP_I(Nst, 1, nei)]
        + u[IDXV_G_R(Nst, 2, ic, ist, idir)] * v1[IDXV_1SP_R(Nst, 2, nei)]
        - u[IDXV_G_I(Nst, 2, ic, ist, idir)] * v1[IDXV_1SP_I(Nst, 2, nei)]);
      v2[IDXV_1SP_I(Nst, ic, ist)] += bc2 *
       (  u[IDXV_G_R(Nst, 0, ic, ist, idir)] * v1[IDXV_1SP_I(Nst, 0, nei)]
        + u[IDXV_G_I(Nst, 0, ic, ist, idir)] * v1[IDXV_1SP_R(Nst, 0, nei)]
        + u[IDXV_G_R(Nst, 1, ic, ist, idir)] * v1[IDXV_1SP_I(Nst, 1, nei)]
        + u[IDXV_G_I(Nst, 1, ic, ist, idir)] * v1[IDXV_1SP_R(Nst, 1, nei)]
        + u[IDXV_G_R(Nst, 2, ic, ist, idir)] * v1[IDXV_1SP_I(Nst, 2, nei)]
        + u[IDXV_G_I(Nst, 2, ic, ist, idir)] * v1[IDXV_1SP_R(Nst, 2, nei)]);
    }

    it2 = (it-1+Nt) % Nt;
    nei = ix + Nx * (iy + Ny * (iz + Nz * it2));
    bc2 = fac * (it > 0 ? 1.0 : bc[idir]);

    for(int ic = 0; ic < NC; ++ic){
      v2[IDXV_1SP_R(Nst, ic, ist)] -= bc2 *
       (  u[IDXV_G_R(Nst, ic, 0, nei, idir)] * v1[IDXV_1SP_R(Nst, 0, nei)]
        + u[IDXV_G_I(Nst, ic, 0, nei, idir)] * v1[IDXV_1SP_I(Nst, 0, nei)]
        + u[IDXV_G_R(Nst, ic, 1, nei, idir)] * v1[IDXV_1SP_R(Nst, 1, nei)]
        + u[IDXV_G_I(Nst, ic, 1, nei, idir)] * v1[IDXV_1SP_I(Nst, 1, nei)]
        + u[IDXV_G_R(Nst, ic, 2, nei, idir)] * v1[IDXV_1SP_R(Nst, 2, nei)]
        + u[IDXV_G_I(Nst, ic, 2, nei, idir)] * v1[IDXV_1SP_I(Nst, 2, nei)]);
      v2[IDXV_1SP_I(Nst, ic, ist)] -= bc2 *
       (  u[IDXV_G_R(Nst, ic, 0, nei, idir)] * v1[IDXV_1SP_I(Nst, 0, nei)]
        - u[IDXV_G_I(Nst, ic, 0, nei, idir)] * v1[IDXV_1SP_R(Nst, 0, nei)]
        + u[IDXV_G_R(Nst, ic, 1, nei, idir)] * v1[IDXV_1SP_I(Nst, 1, nei)]
        - u[IDXV_G_I(Nst, ic, 1, nei, idir)] * v1[IDXV_1SP_R(Nst, 1, nei)]
        + u[IDXV_G_R(Nst, ic, 2, nei, idir)] * v1[IDXV_1SP_I(Nst, 2, nei)]
        - u[IDXV_G_I(Nst, ic, 2, nei, idir)] * v1[IDXV_1SP_R(Nst, 2, nei)]);
    }

   }

  }

}

//====================================================================
void mult_staggered_1(
                  real_t *RESTRICT buf_xp, real_t *RESTRICT buf_xm,
                  real_t *RESTRICT buf_yp, real_t *RESTRICT buf_ym,
                  real_t *RESTRICT buf_zp, real_t *RESTRICT buf_zm,
                  real_t *RESTRICT buf_tp, real_t *RESTRICT buf_tm,
                  real_t *RESTRICT u, real_t *RESTRICT v1,
                  int *Nsize, int *bc, int *do_comm)
{
  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nst  = Nx * Ny * Nz * Nt;

  int Nyzt = Ny * Nz * Nt;
  int Nxzt = Nx * Nz * Nt;
  int Nxyt = Nx * Ny * Nt;
  int Nxyz = Nx * Ny * Nz;

  if(do_comm[0] > 0){

    int idir = 0;
    int Nyzt = Ny * Nz * Nt;

    for(int iyzt = 0; iyzt < Nyzt; ++iyzt){
      int ix  = 0;
      int ist = ix + Nx * iyzt;
      real_t bc2 = bc[0];

      real_t *wt = buf_xp;

#pragma _NEC unroll(NC)
      for(int ic = 0; ic < NC; ++ic){
        wt[IDXV_1SP_R(Nyzt, ic, iyzt)] = bc2 * v1[IDXV_1SP_R(Nst, ic, ist)];
        wt[IDXV_1SP_I(Nyzt, ic, iyzt)] = bc2 * v1[IDXV_1SP_I(Nst, ic, ist)];
      }
    }

    for(int iyzt = 0; iyzt < Nyzt; ++iyzt){
      int ix  = Nx-1;
      int ist  = ix + Nx * iyzt;
      int istu = ist;
      real_t *wt = buf_xm;
      int ibf = iyzt;
      int Nbf = Nyzt;

#pragma _NEC unroll(NC)
      for(int ic = 0; ic < NC; ++ic){
        wt[IDXV_1SP_R(Nbf, ic, ibf)] =
            u[IDXV_G_R(Nst, ic, 0, ist, idir)] * v1[IDXV_1SP_R(Nst, 0, ist)]
          + u[IDXV_G_I(Nst, ic, 0, ist, idir)] * v1[IDXV_1SP_I(Nst, 0, ist)]
          + u[IDXV_G_R(Nst, ic, 1, ist, idir)] * v1[IDXV_1SP_R(Nst, 1, ist)]
          + u[IDXV_G_I(Nst, ic, 1, ist, idir)] * v1[IDXV_1SP_I(Nst, 1, ist)]
          + u[IDXV_G_R(Nst, ic, 2, ist, idir)] * v1[IDXV_1SP_R(Nst, 2, ist)]
          + u[IDXV_G_I(Nst, ic, 2, ist, idir)] * v1[IDXV_1SP_I(Nst, 2, ist)];
       wt[IDXV_1SP_I(Nbf, ic, ibf)] =
            u[IDXV_G_R(Nst, ic, 0, ist, idir)] * v1[IDXV_1SP_I(Nst, 0, ist)]
          - u[IDXV_G_I(Nst, ic, 0, ist, idir)] * v1[IDXV_1SP_R(Nst, 0, ist)]
          + u[IDXV_G_R(Nst, ic, 1, ist, idir)] * v1[IDXV_1SP_I(Nst, 1, ist)]
          - u[IDXV_G_I(Nst, ic, 1, ist, idir)] * v1[IDXV_1SP_R(Nst, 1, ist)]
          + u[IDXV_G_R(Nst, ic, 2, ist, idir)] * v1[IDXV_1SP_I(Nst, 2, ist)]
          - u[IDXV_G_I(Nst, ic, 2, ist, idir)] * v1[IDXV_1SP_R(Nst, 2, ist)];
      }

    }

  } // do_comm[0]


  // idir = 1;
  if(do_comm[1] > 0){

    int Nzt  = Nz * Nt;
    int Nxzt = Nx * Nzt;

    int idir = 1;
    //  for(int ixzt = 0; ixzt < Nxzt; ++ixzt){
    //  int ix  = ixzt % Nx;
    //  int izt = ixzt/Nx;
    //  {
   for(int izt = 0; izt < Nzt; ++izt){
     for(int ix = 0; ix < Nx; ++ix){
       int iy   = 0;
       int ist  = ix + Nx * (iy + Ny * izt);
       int ixzt = ix + Nx * izt;
       real_t bc2 = bc[1];
       real_t *wt = buf_yp;
#pragma _NEC unroll(NC)
       for(int ic = 0; ic < NC; ++ic){
         wt[IDXV_1SP_R(Nxzt, ic, ixzt)] = bc2 * v1[IDXV_1SP_R(Nst, ic, ist)];
         wt[IDXV_1SP_I(Nxzt, ic, ixzt)] = bc2 * v1[IDXV_1SP_I(Nst, ic, ist)];
       }
     }
    }

    //for(int ixzt = 0; ixzt < Nxzt; ++ixzt){
    //  int ix  = ixzt % Nx;
    //  int izt = ixzt/Nx;
    //  {
    for(int izt = 0; izt < Nzt; ++izt){
     for(int ix = 0; ix < Nx; ++ix){
       int iy   = Ny-1;
       int ist  = ix + Nx * (iy + Ny*izt);
       int istu = ist + Nst;
       int ixzt = ix + Nx * izt;

       // int ibf = ix + Nx * izt;
       // int Nbf = Nxzt;
       real_t *wt = buf_ym;

      for(int ic = 0; ic < NC; ++ic){
        wt[IDXV_1SP_R(Nxzt, ic, ixzt)] =
            u[IDXV_G_R(Nst, ic, 0, ist, idir)] * v1[IDXV_1SP_R(Nst, 0, ist)]
          + u[IDXV_G_I(Nst, ic, 0, ist, idir)] * v1[IDXV_1SP_I(Nst, 0, ist)]
          + u[IDXV_G_R(Nst, ic, 1, ist, idir)] * v1[IDXV_1SP_R(Nst, 1, ist)]
          + u[IDXV_G_I(Nst, ic, 1, ist, idir)] * v1[IDXV_1SP_I(Nst, 1, ist)]
          + u[IDXV_G_R(Nst, ic, 2, ist, idir)] * v1[IDXV_1SP_R(Nst, 2, ist)]
          + u[IDXV_G_I(Nst, ic, 2, ist, idir)] * v1[IDXV_1SP_I(Nst, 2, ist)];
       wt[IDXV_1SP_I(Nxzt, ic, ixzt)] =
            u[IDXV_G_R(Nst, ic, 0, ist, idir)] * v1[IDXV_1SP_I(Nst, 0, ist)]
          - u[IDXV_G_I(Nst, ic, 0, ist, idir)] * v1[IDXV_1SP_R(Nst, 0, ist)]
          + u[IDXV_G_R(Nst, ic, 1, ist, idir)] * v1[IDXV_1SP_I(Nst, 1, ist)]
          - u[IDXV_G_I(Nst, ic, 1, ist, idir)] * v1[IDXV_1SP_R(Nst, 1, ist)]
          + u[IDXV_G_R(Nst, ic, 2, ist, idir)] * v1[IDXV_1SP_I(Nst, 2, ist)]
          - u[IDXV_G_I(Nst, ic, 2, ist, idir)] * v1[IDXV_1SP_R(Nst, 2, ist)];
      }

     }
    }

  } // do_comm[1]

  //  idir = 2;
  if(do_comm[2] > 0){

    int idir = 2;
    int Nxy = Nx * Ny;

#pragma _NEC ivdep
    // for(int ixyt = 0; ixyt < Nxyt; ++ixyt){
    //  int ixy = ixyt % Nxy;
    //  int it = ixyt/Nxy;
    //  {
    for(int it = 0; it < Nt; ++it){
     for(int ixy = 0; ixy < Nxy; ++ixy){
       int iz   = 0;
       int ist  = ixy + Nxy * (iz + Nz * it);
       int ixyt = ixy + Nxy * it;

       real_t bc2 = bc[2];

       real_t *wt = buf_zp;

       for(int ic = 0; ic < NC; ++ic){
         wt[IDXV_1SP_R(Nxyt, ic, ixyt)] = bc2 * v1[IDXV_1SP_R(Nst, ic, ist)];
         wt[IDXV_1SP_I(Nxyt, ic, ixyt)] = bc2 * v1[IDXV_1SP_I(Nst, ic, ist)];
       }

     }
    }

#pragma _NEC ivdep
    // for(int ixyt = 0; ixyt < Nxyt; ++ixyt){
    //  int ixy = ixyt % Nxy;
    //  int it = ixyt/Nxy;
    // {
    for(int it = 0; it < Nt; ++it){
     for(int ixy = 0; ixy < Nxy; ++ixy){
       int iz = Nz-1;
       int ist  = ixy + Nxy * (iz + Nz * it);
       int istu = ist + Nst * 2;

       real_t *wt = buf_zm;
       int ibf = ixy + Nxy * it;
       int Nbf = Nxyt;

      for(int ic = 0; ic < NC; ++ic){
        wt[IDXV_1SP_R(Nbf, ic, ibf)] =
            u[IDXV_G_R(Nst, ic, 0, ist, idir)] * v1[IDXV_1SP_R(Nst, 0, ist)]
          + u[IDXV_G_I(Nst, ic, 0, ist, idir)] * v1[IDXV_1SP_I(Nst, 0, ist)]
          + u[IDXV_G_R(Nst, ic, 1, ist, idir)] * v1[IDXV_1SP_R(Nst, 1, ist)]
          + u[IDXV_G_I(Nst, ic, 1, ist, idir)] * v1[IDXV_1SP_I(Nst, 1, ist)]
          + u[IDXV_G_R(Nst, ic, 2, ist, idir)] * v1[IDXV_1SP_R(Nst, 2, ist)]
          + u[IDXV_G_I(Nst, ic, 2, ist, idir)] * v1[IDXV_1SP_I(Nst, 2, ist)];
       wt[IDXV_1SP_I(Nbf, ic, ibf)] =
            u[IDXV_G_R(Nst, ic, 0, ist, idir)] * v1[IDXV_1SP_I(Nst, 0, ist)]
          - u[IDXV_G_I(Nst, ic, 0, ist, idir)] * v1[IDXV_1SP_R(Nst, 0, ist)]
          + u[IDXV_G_R(Nst, ic, 1, ist, idir)] * v1[IDXV_1SP_I(Nst, 1, ist)]
          - u[IDXV_G_I(Nst, ic, 1, ist, idir)] * v1[IDXV_1SP_R(Nst, 1, ist)]
          + u[IDXV_G_R(Nst, ic, 2, ist, idir)] * v1[IDXV_1SP_I(Nst, 2, ist)]
          - u[IDXV_G_I(Nst, ic, 2, ist, idir)] * v1[IDXV_1SP_R(Nst, 2, ist)];
      }

     }
    }

  } // do_comm[2]

   //  idir = 3;
  if(do_comm[3] > 0){

    int idir = 3;
    int Nxyz = Nx * Ny * Nz;

#pragma _NEC ivdep
    for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
      int it = 0;
      int ist  = ixyz + Nxyz * it;
      real_t bc2 = bc[3];

      real_t *wt = buf_tp;
      for(int ic = 0; ic < NC; ++ic){
        wt[IDXV_1SP_R(Nxyz, ic, ixyz)] = bc2 * v1[IDXV_1SP_R(Nst, ic, ist)];
        wt[IDXV_1SP_I(Nxyz, ic, ixyz)] = bc2 * v1[IDXV_1SP_I(Nst, ic, ist)];
      }

    }

#pragma _NEC ivdep
    for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
      int it   = Nt-1;
      int ist  = ixyz + Nxyz * it;
      int istu = ist + Nst * 3;

      real_t *wt = buf_tm;
      int ibf = ixyz;
      int Nbf = Nxyz;

      for(int ic = 0; ic < NC; ++ic){
        wt[IDXV_1SP_R(Nbf, ic, ibf)] =
            u[IDXV_G_R(Nst, ic, 0, ist, idir)] * v1[IDXV_1SP_R(Nst, 0, ist)]
          + u[IDXV_G_I(Nst, ic, 0, ist, idir)] * v1[IDXV_1SP_I(Nst, 0, ist)]
          + u[IDXV_G_R(Nst, ic, 1, ist, idir)] * v1[IDXV_1SP_R(Nst, 1, ist)]
          + u[IDXV_G_I(Nst, ic, 1, ist, idir)] * v1[IDXV_1SP_I(Nst, 1, ist)]
          + u[IDXV_G_R(Nst, ic, 2, ist, idir)] * v1[IDXV_1SP_R(Nst, 2, ist)]
          + u[IDXV_G_I(Nst, ic, 2, ist, idir)] * v1[IDXV_1SP_I(Nst, 2, ist)];
       wt[IDXV_1SP_I(Nbf, ic, ibf)] =
            u[IDXV_G_R(Nst, ic, 0, ist, idir)] * v1[IDXV_1SP_I(Nst, 0, ist)]
          - u[IDXV_G_I(Nst, ic, 0, ist, idir)] * v1[IDXV_1SP_R(Nst, 0, ist)]
          + u[IDXV_G_R(Nst, ic, 1, ist, idir)] * v1[IDXV_1SP_I(Nst, 1, ist)]
          - u[IDXV_G_I(Nst, ic, 1, ist, idir)] * v1[IDXV_1SP_R(Nst, 1, ist)]
          + u[IDXV_G_R(Nst, ic, 2, ist, idir)] * v1[IDXV_1SP_I(Nst, 2, ist)]
          - u[IDXV_G_I(Nst, ic, 2, ist, idir)] * v1[IDXV_1SP_R(Nst, 2, ist)];
      }


    }

  } // do_comm[3]

}

//====================================================================
void mult_staggered_2(real_t *RESTRICT v2, real_t *RESTRICT u,
		      real_t *RESTRICT buf_xp, real_t *RESTRICT buf_xm,
		      real_t *RESTRICT buf_yp, real_t *RESTRICT buf_ym,
		      real_t *RESTRICT buf_zp, real_t *RESTRICT buf_zm,
		      real_t *RESTRICT buf_tp, real_t *RESTRICT buf_tm,
                      int *Nsize, int *bc,  int *do_comm, real_t mq, int jdag)
{
  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nst  = Nx * Ny * Nz * Nt;

  int Nyzt = Ny * Nz * Nt;
  int Nxzt = Nx * Nz * Nt;
  int Nxyt = Nx * Ny * Nt;
  int Nxyz = Nx * Ny * Nz;

  real_t fac = real_t(jdag) * 0.5/mq;

  if(do_comm[0] > 0){

    int idir = 0;

#pragma _NEC ivdep
    for(int iyzt = 0; iyzt < Nyzt; ++iyzt){
      int ix = Nx-1;
      int ist = ix + Nx * iyzt;

      real_t *v1 = buf_xp;
      int Nbf = Nyzt;
      int ibf = iyzt;

      for(int ic = 0; ic < NC; ++ic){
        v2[IDXV_1SP_R(Nst, ic, ist)] += fac *
         (  u[IDXV_G_R(Nst, 0, ic, ist, idir)] * v1[IDXV_1SP_R(Nbf, 0, ibf)]
          - u[IDXV_G_I(Nst, 0, ic, ist, idir)] * v1[IDXV_1SP_I(Nbf, 0, ibf)]
          + u[IDXV_G_R(Nst, 1, ic, ist, idir)] * v1[IDXV_1SP_R(Nbf, 1, ibf)]
          - u[IDXV_G_I(Nst, 1, ic, ist, idir)] * v1[IDXV_1SP_I(Nbf, 1, ibf)]
          + u[IDXV_G_R(Nst, 2, ic, ist, idir)] * v1[IDXV_1SP_R(Nbf, 2, ibf)]
          - u[IDXV_G_I(Nst, 2, ic, ist, idir)] * v1[IDXV_1SP_I(Nbf, 2, ibf)]);
        v2[IDXV_1SP_I(Nst, ic, ist)] += fac *
         (  u[IDXV_G_R(Nst, 0, ic, ist, idir)] * v1[IDXV_1SP_I(Nbf, 0, ibf)]
          + u[IDXV_G_I(Nst, 0, ic, ist, idir)] * v1[IDXV_1SP_R(Nbf, 0, ibf)]
          + u[IDXV_G_R(Nst, 1, ic, ist, idir)] * v1[IDXV_1SP_I(Nbf, 1, ibf)]
          + u[IDXV_G_I(Nst, 1, ic, ist, idir)] * v1[IDXV_1SP_R(Nbf, 1, ibf)]
          + u[IDXV_G_R(Nst, 2, ic, ist, idir)] * v1[IDXV_1SP_I(Nbf, 2, ibf)]
          + u[IDXV_G_I(Nst, 2, ic, ist, idir)] * v1[IDXV_1SP_R(Nbf, 2, ibf)]);
      }

    }

    real_t bc2 = fac * bc[idir];

#pragma _NEC ivdep
    for(int iyzt = 0; iyzt < Nyzt; ++iyzt){
      int ix  = 0;
      int ist = ix + Nx * iyzt;
      int ibf = iyzt;

      for(int ivc = 0; ivc < NVC; ++ ivc){
        v2[IDXV_1SP(Nst, ivc, ist)] -= bc2 * buf_xm[IDXV_1SP(Nyzt, ivc, ibf)];
      }
    }

  } // do_comm[0]

  if(do_comm[1] > 0){

    int idir = 1;

#pragma _NEC ivdep
    for(int ixzt = 0; ixzt < Nxzt; ++ixzt){
      int iy  = Ny-1;
      int ix  = ixzt % Nx;
      int izt = ixzt/Nx;
      int ist = ix + Nx * (iy + Ny*izt);

      real_t *v1 = buf_yp;
      int Nbf = Nxzt;
      int ibf = ixzt;

      for(int ic = 0; ic < NC; ++ic){
        v2[IDXV_1SP_R(Nst, ic, ist)] += fac *
         (  u[IDXV_G_R(Nst, 0, ic, ist, idir)] * v1[IDXV_1SP_R(Nbf, 0, ibf)]
          - u[IDXV_G_I(Nst, 0, ic, ist, idir)] * v1[IDXV_1SP_I(Nbf, 0, ibf)]
          + u[IDXV_G_R(Nst, 1, ic, ist, idir)] * v1[IDXV_1SP_R(Nbf, 1, ibf)]
          - u[IDXV_G_I(Nst, 1, ic, ist, idir)] * v1[IDXV_1SP_I(Nbf, 1, ibf)]
          + u[IDXV_G_R(Nst, 2, ic, ist, idir)] * v1[IDXV_1SP_R(Nbf, 2, ibf)]
          - u[IDXV_G_I(Nst, 2, ic, ist, idir)] * v1[IDXV_1SP_I(Nbf, 2, ibf)]);
        v2[IDXV_1SP_I(Nst, ic, ist)] += fac *
         (  u[IDXV_G_R(Nst, 0, ic, ist, idir)] * v1[IDXV_1SP_I(Nbf, 0, ibf)]
          + u[IDXV_G_I(Nst, 0, ic, ist, idir)] * v1[IDXV_1SP_R(Nbf, 0, ibf)]
          + u[IDXV_G_R(Nst, 1, ic, ist, idir)] * v1[IDXV_1SP_I(Nbf, 1, ibf)]
          + u[IDXV_G_I(Nst, 1, ic, ist, idir)] * v1[IDXV_1SP_R(Nbf, 1, ibf)]
          + u[IDXV_G_R(Nst, 2, ic, ist, idir)] * v1[IDXV_1SP_I(Nbf, 2, ibf)]
          + u[IDXV_G_I(Nst, 2, ic, ist, idir)] * v1[IDXV_1SP_R(Nbf, 2, ibf)]);
      }

    }

    real_t bc2 = fac * bc[idir];

#pragma _NEC ivdep
    for(int ixzt = 0; ixzt < Nxzt; ++ixzt){
      int iy  = 0;
      int ix = ixzt % Nx;
      int izt = ixzt/Nx;
      int ist = ix + Nx * (iy + Ny*izt);

      for(int ivc = 0; ivc < NVC; ++ ivc){
        v2[IDXV_1SP(Nst, ivc, ist)] -= bc2 * buf_ym[IDXV_1SP(Nxzt, ivc, ixzt)];
      }
    }

  } // do_comm[1]

  if(do_comm[2] > 0){

    int idir = 2;
    int Nxy = Nx * Ny;

#pragma _NEC ivdep
    for(int ixyt = 0; ixyt < Nxyt; ++ixyt){
      int iz  = Nz-1;
      int ixy = ixyt % Nxy;
      int it  = ixyt/Nxy;
      int ist = ixy + Nxy * (iz + Nz * it);

      real_t *v1 = buf_zp;
      int Nbf = Nxyt;
      int ibf = ixyt;

      for(int ic = 0; ic < NC; ++ic){
        v2[IDXV_1SP_R(Nst, ic, ist)] += fac *
         (  u[IDXV_G_R(Nst, 0, ic, ist, idir)] * v1[IDXV_1SP_R(Nbf, 0, ibf)]
          - u[IDXV_G_I(Nst, 0, ic, ist, idir)] * v1[IDXV_1SP_I(Nbf, 0, ibf)]
          + u[IDXV_G_R(Nst, 1, ic, ist, idir)] * v1[IDXV_1SP_R(Nbf, 1, ibf)]
          - u[IDXV_G_I(Nst, 1, ic, ist, idir)] * v1[IDXV_1SP_I(Nbf, 1, ibf)]
          + u[IDXV_G_R(Nst, 2, ic, ist, idir)] * v1[IDXV_1SP_R(Nbf, 2, ibf)]
          - u[IDXV_G_I(Nst, 2, ic, ist, idir)] * v1[IDXV_1SP_I(Nbf, 2, ibf)]);
        v2[IDXV_1SP_I(Nst, ic, ist)] += fac *
         (  u[IDXV_G_R(Nst, 0, ic, ist, idir)] * v1[IDXV_1SP_I(Nbf, 0, ibf)]
          + u[IDXV_G_I(Nst, 0, ic, ist, idir)] * v1[IDXV_1SP_R(Nbf, 0, ibf)]
          + u[IDXV_G_R(Nst, 1, ic, ist, idir)] * v1[IDXV_1SP_I(Nbf, 1, ibf)]
          + u[IDXV_G_I(Nst, 1, ic, ist, idir)] * v1[IDXV_1SP_R(Nbf, 1, ibf)]
          + u[IDXV_G_R(Nst, 2, ic, ist, idir)] * v1[IDXV_1SP_I(Nbf, 2, ibf)]
          + u[IDXV_G_I(Nst, 2, ic, ist, idir)] * v1[IDXV_1SP_R(Nbf, 2, ibf)]);
      }

    }

    real_t bc2 = fac * bc[idir];

#pragma _NEC ivdep
    for(int ixyt = 0; ixyt < Nxyt; ++ixyt){
      int iz  = 0;
      int ixy = ixyt % Nxy;
      int it  = ixyt/Nxy;
      int ist = ixy + Nxy * (iz + Nz*it);

      for(int ivc = 0; ivc < NVC; ++ ivc){
        v2[IDXV_1SP(Nst, ivc, ist)] -= bc2 * buf_zm[IDXV_1SP(Nxyt, ivc, ixyt)];
      }
    }

  } // do_comm[2]

  if(do_comm[3] > 0){

    int idir = 3;

#pragma _NEC ivdep
    for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
      int it  = Nt-1;
      int ist = ixyz + Nxyz * it;

      real_t *v1 = buf_tp;
      int Nbf = Nxyz;
      int ibf = ixyz;

      for(int ic = 0; ic < NC; ++ic){
        v2[IDXV_1SP_R(Nst, ic, ist)] += fac *
         (  u[IDXV_G_R(Nst, 0, ic, ist, idir)] * v1[IDXV_1SP_R(Nbf, 0, ibf)]
          - u[IDXV_G_I(Nst, 0, ic, ist, idir)] * v1[IDXV_1SP_I(Nbf, 0, ibf)]
          + u[IDXV_G_R(Nst, 1, ic, ist, idir)] * v1[IDXV_1SP_R(Nbf, 1, ibf)]
          - u[IDXV_G_I(Nst, 1, ic, ist, idir)] * v1[IDXV_1SP_I(Nbf, 1, ibf)]
          + u[IDXV_G_R(Nst, 2, ic, ist, idir)] * v1[IDXV_1SP_R(Nbf, 2, ibf)]
          - u[IDXV_G_I(Nst, 2, ic, ist, idir)] * v1[IDXV_1SP_I(Nbf, 2, ibf)]);
        v2[IDXV_1SP_I(Nst, ic, ist)] += fac *
         (  u[IDXV_G_R(Nst, 0, ic, ist, idir)] * v1[IDXV_1SP_I(Nbf, 0, ibf)]
          + u[IDXV_G_I(Nst, 0, ic, ist, idir)] * v1[IDXV_1SP_R(Nbf, 0, ibf)]
          + u[IDXV_G_R(Nst, 1, ic, ist, idir)] * v1[IDXV_1SP_I(Nbf, 1, ibf)]
          + u[IDXV_G_I(Nst, 1, ic, ist, idir)] * v1[IDXV_1SP_R(Nbf, 1, ibf)]
          + u[IDXV_G_R(Nst, 2, ic, ist, idir)] * v1[IDXV_1SP_I(Nbf, 2, ibf)]
          + u[IDXV_G_I(Nst, 2, ic, ist, idir)] * v1[IDXV_1SP_R(Nbf, 2, ibf)]);
      }

    }

    real_t bc2 = fac * bc[idir];

#pragma _NEC ivdep
    for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
      int it  = 0;
      int ist = ixyz + Nxyz * it;

      for(int ivc = 0; ivc < NVC; ++ ivc){
        v2[IDXV_1SP(Nst, ivc, ist)] -= bc2 * buf_tm[IDXV_1SP(Nxyz, ivc, ixyz)];
      }
    }

  } // do_comm[3]

}


//====================================================================

