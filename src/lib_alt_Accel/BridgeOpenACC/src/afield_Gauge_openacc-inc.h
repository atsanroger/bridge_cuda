/*!
      @file    afield_Gauge_openacc-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2024-11-05 14:28:33 #$
      @version $LastChangedRevision: 2618 $
*/

#ifndef AFIELD_GAUGE_OPENACC_INC_INCLUDED
#define AFIELD_GAUGE_OPENACC_INC_INCLUDED

//#include <cstdlib> 

#include "inline/define_params.h"
#include "inline/define_index.h"

#define MULT_GXr(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)   (u0*v0-u1*v1 + u2*v2-u3*v3 + u4*v4-u5*v5)
#define MULT_GXi(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)   (u0*v1+u1*v0 + u2*v3+u3*v2 + u4*v5+u5*v4)

#define EXT_IMG_R(v1r,v1i,v2r,v2i,w1r,w1i,w2r,w2i)  (v1r*w2r - v1i*w2i - v2r*w1r + v2i*w1i)
#define EXT_IMG_I(v1r,v1i,v2r,v2i,w1r,w1i,w2r,w2i)  (- v1r*w2i - v1i*w2r + v2r*w1i + v2i*w1r)

#define MULT_R(v1r, v1i, v2r, v2i)  (v1r * v2r - v1i * v2i)
#define MULT_I(v1r, v1i, v2r, v2i)  (v1r * v2i + v1i * v2r)

//====================================================================
inline real_t mult_r(const real_t v1r, const real_t v1i,
                     const real_t v2r, const real_t v2i){
  return  v1r * v2r - v1i * v2i;
}

inline real_t mult_i(const real_t v1r, const real_t v1i,
                     const real_t v2r, const real_t v2i){
  return  v1r * v2i + v1i * v2r;
}

//====================================================================
void multadd_Gnn(real_t *restrict u, const int exu,
                 real_t *restrict v, const int exv,
                 real_t *restrict w, const int exw,
                 const real_t a, const int nst)
{
  int ngst = NDF * nst;
  int up = nst * exu;
  int vp = nst * exv;
  int wp = nst * exw;

#pragma acc data present(u[up:ngst], v[vp:ngst], w[wp:ngst]) \
                 copyin(a, nst, exu, exv, exw)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for(int ist = 0; ist < nst; ++ist){
    int igw = ist + nst * exw;
    int igv = ist + nst * exv;
    int igu = ist + nst * exu;

#ifdef SU3_3RD_ROW_RECONST
    real_t vt[NDF], wt[NVC];
    for(int ic2 = 0; ic2 < NC-1; ++ic2){
      for(int ic1 = 0; ic1 < NC; ++ic1){
        int icr = 2*ic1   + NVC * ic2;
        int ici = 2*ic1+1 + NVC * ic2;
        vt[icr] = v[IDX2_G_R(ic1,ic2,igv)];
        vt[ici] = v[IDX2_G_I(ic1,ic2,igv)];
      }
    }

    vt[12] = EXT_IMG_R(vt[ 2], vt[ 3], vt[ 4], vt[ 5],
                       vt[ 8], vt[ 9], vt[10], vt[11]);
    vt[13] = EXT_IMG_I(vt[ 2], vt[ 3], vt[ 4], vt[ 5],
                       vt[ 8], vt[ 9], vt[10], vt[11]);
    vt[14] = EXT_IMG_R(vt[ 4], vt[ 5], vt[ 0], vt[ 1],
                       vt[10], vt[11], vt[ 6], vt[ 7]);
    vt[15] = EXT_IMG_I(vt[ 4], vt[ 5], vt[ 0], vt[ 1],
                       vt[10], vt[11], vt[ 6], vt[ 7]);
    vt[16] = EXT_IMG_R(vt[ 0], vt[ 1], vt[ 2], vt[ 3],
                       vt[ 6], vt[ 7], vt[ 8], vt[ 9]);
    vt[17] = EXT_IMG_I(vt[ 0], vt[ 1], vt[ 2], vt[ 3],
                       vt[ 6], vt[ 7], vt[ 8], vt[ 9]);
#else
    real_t vt[NDF], wt[NVC];
    for(int ic2 = 0; ic2 < NC; ++ic2){
      for(int ic1 = 0; ic1 < NC; ++ic1){
        int icr = 2*ic1   + NVC * ic2;
        int ici = 2*ic1+1 + NVC * ic2;
        vt[icr] = v[IDX2_G_R(ic1,ic2,igv)];
        vt[ici] = v[IDX2_G_I(ic1,ic2,igv)];
      }
    }
#endif

    for(int ic2 = 0; ic2 < NC; ++ic2){

      for(int ic1 = 0; ic1 < NC; ++ic1){
        int icr = 2*ic1;
        int ici = 2*ic1+1;
        wt[icr] = w[IDX2_G_R(ic2,ic1,igw)];
        wt[ici] = w[IDX2_G_I(ic2,ic1,igw)];
      }

      for(int ic1 = 0; ic1 < NC; ++ic1){
        int j = NVC * ic1;
        real_t ur, ui;
        ur = MULT_GXr(vt[0+j], vt[1+j], vt[2+j], vt[3+j], vt[4+j], vt[5+j],
                      wt[0],   wt[1],   wt[2],   wt[3],   wt[4], wt[5]);
        ui = MULT_GXi(vt[0+j], vt[1+j], vt[2+j], vt[3+j], vt[4+j], vt[5+j],
                      wt[0],   wt[1],   wt[2],   wt[3],   wt[4], wt[5]);
        u[IDX2_G_R(ic2,ic1,igu)] += a * ur;
        u[IDX2_G_I(ic2,ic1,igu)] += a * ui;
      }

    }

  }

 } // acc parallel

}

//====================================================================
void mult_Gnn(real_t *restrict u, const int exu,
              real_t *restrict v, const int exv,
              real_t *restrict w, const int exw, const int nst)
{
  int ngst = NDF * nst;
  int up = nst * exu;
  int vp = nst * exv;
  int wp = nst * exw;

#pragma acc data present(u[up:ngst], v[vp:ngst], w[wp:ngst]) \
                 copyin(nst, exu, exv, exw)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for(int ist = 0; ist < nst; ++ist){
    int igw = ist + nst * exw;
    int igv = ist + nst * exv;
    int igu = ist + nst * exu;

#ifdef SU3_3RD_ROW_RECONST
    real_t vt[NDF], wt[NVC];
    for(int ic2 = 0; ic2 < NC-1; ++ic2){
      for(int ic1 = 0; ic1 < NC; ++ic1){
        int icr = 2*ic1   + NVC * ic2;
        int ici = 2*ic1+1 + NVC * ic2;
        vt[icr] = v[IDX2_G_R(ic1,ic2,igv)];
        vt[ici] = v[IDX2_G_I(ic1,ic2,igv)];
      }
    }

    vt[12] = EXT_IMG_R(vt[ 2], vt[ 3], vt[ 4], vt[ 5],
                       vt[ 8], vt[ 9], vt[10], vt[11]);
    vt[13] = EXT_IMG_I(vt[ 2], vt[ 3], vt[ 4], vt[ 5],
                       vt[ 8], vt[ 9], vt[10], vt[11]);
    vt[14] = EXT_IMG_R(vt[ 4], vt[ 5], vt[ 0], vt[ 1],
                       vt[10], vt[11], vt[ 6], vt[ 7]);
    vt[15] = EXT_IMG_I(vt[ 4], vt[ 5], vt[ 0], vt[ 1],
                       vt[10], vt[11], vt[ 6], vt[ 7]);
    vt[16] = EXT_IMG_R(vt[ 0], vt[ 1], vt[ 2], vt[ 3],
                       vt[ 6], vt[ 7], vt[ 8], vt[ 9]);
    vt[17] = EXT_IMG_I(vt[ 0], vt[ 1], vt[ 2], vt[ 3],
                       vt[ 6], vt[ 7], vt[ 8], vt[ 9]);
#else
    real_t vt[NDF], wt[NVC];
    for(int ic2 = 0; ic2 < NC; ++ic2){
      for(int ic1 = 0; ic1 < NC; ++ic1){
        int icr = 2*ic1   + NVC * ic2;
        int ici = 2*ic1+1 + NVC * ic2;
        vt[icr] = v[IDX2_G_R(ic1,ic2,igv)];
        vt[ici] = v[IDX2_G_I(ic1,ic2,igv)];
      }
    }
#endif

    for(int ic2 = 0; ic2 < NC; ++ic2){

      for(int ic1 = 0; ic1 < NC; ++ic1){
        int icr = 2*ic1;
        int ici = 2*ic1+1;
        wt[icr] = w[IDX2_G_R(ic2,ic1,igw)];
        wt[ici] = w[IDX2_G_I(ic2,ic1,igw)];
      }

      for(int ic1 = 0; ic1 < NC; ++ic1){
        int j = NVC * ic1;
        real_t ur, ui;
        ur = MULT_GXr(vt[0+j], vt[1+j], vt[2+j], vt[3+j], vt[4+j], vt[5+j],
                      wt[0],   wt[1],   wt[2],   wt[3],   wt[4],   wt[5]);
        ui = MULT_GXi(vt[0+j], vt[1+j], vt[2+j], vt[3+j], vt[4+j], vt[5+j],
                      wt[0],   wt[1],   wt[2],   wt[3],   wt[4],   wt[5]);
        u[IDX2_G_R(ic2,ic1,igu)] = ur;
        u[IDX2_G_I(ic2,ic1,igu)] = ui;
      }

    }

  }

 } // acc parallel

}

//====================================================================
void multadd_Gnd(real_t *restrict u, const int exu,
                 real_t *restrict v, const int exv,
                 real_t *restrict w, const int exw,
                 const real_t a, const int nst)
{
  int ngst = NDF * nst;
  int up = nst * exu;
  int vp = nst * exv;
  int wp = nst * exw;

#pragma acc data present(u[up:ngst], v[vp:ngst], w[wp:ngst]) \
                 copyin(a, nst, exu, exv, exw)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for(int ist = 0; ist < nst; ++ist){
    int igw = ist + nst * exw;
    int igv = ist + nst * exv;
    int igu = ist + nst * exu;

#ifdef SU3_3RD_ROW_RECONST
    real_t vt[NDF], wt[NVC];
    for(int ic2 = 0; ic2 < NC-1; ++ic2){
      for(int ic1 = 0; ic1 < NC; ++ic1){
        int icr = 2*ic1   + NVC * ic2;
        int ici = 2*ic1+1 + NVC * ic2;
        vt[icr] = v[IDX2_G_R(ic1,ic2,igv)];
        vt[ici] = v[IDX2_G_I(ic1,ic2,igv)];
      }
    }

    vt[12] = EXT_IMG_R(vt[ 2], vt[ 3], vt[ 4], vt[ 5],
                       vt[ 8], vt[ 9], vt[10], vt[11]);
    vt[13] = EXT_IMG_I(vt[ 2], vt[ 3], vt[ 4], vt[ 5],
                       vt[ 8], vt[ 9], vt[10], vt[11]);
    vt[14] = EXT_IMG_R(vt[ 4], vt[ 5], vt[ 0], vt[ 1],
                       vt[10], vt[11], vt[ 6], vt[ 7]);
    vt[15] = EXT_IMG_I(vt[ 4], vt[ 5], vt[ 0], vt[ 1],
                       vt[10], vt[11], vt[ 6], vt[ 7]);
    vt[16] = EXT_IMG_R(vt[ 0], vt[ 1], vt[ 2], vt[ 3],
                       vt[ 6], vt[ 7], vt[ 8], vt[ 9]);
    vt[17] = EXT_IMG_I(vt[ 0], vt[ 1], vt[ 2], vt[ 3],
                       vt[ 6], vt[ 7], vt[ 8], vt[ 9]);
#else
    real_t vt[NDF], wt[NVC];
    for(int ic2 = 0; ic2 < NC; ++ic2){
      for(int ic1 = 0; ic1 < NC; ++ic1){
        int icr = 2*ic1   + NVC * ic2;
        int ici = 2*ic1+1 + NVC * ic2;
        vt[icr] = v[IDX2_G_R(ic1,ic2,igv)];
        vt[ici] = v[IDX2_G_I(ic1,ic2,igv)];
      }
    }
#endif

    for(int ic2 = 0; ic2 < NC; ++ic2){

      for(int ic1 = 0; ic1 < NC; ++ic1){
        int icr = 2*ic1;
        int ici = 2*ic1+1;
        wt[icr] =  w[IDX2_G_R(ic1,ic2,igw)];
        wt[ici] = -w[IDX2_G_I(ic1,ic2,igw)];
      }

      for(int ic1 = 0; ic1 < NC; ++ic1){
        int j = NVC * ic1;
        real_t ur, ui;
        ur = MULT_GXr(vt[0+j], vt[1+j], vt[2+j], vt[3+j], vt[4+j], vt[5+j],
                      wt[0],   wt[1],   wt[2],   wt[3],   wt[4],   wt[5]);
        ui = MULT_GXi(vt[0+j], vt[1+j], vt[2+j], vt[3+j], vt[4+j], vt[5+j],
                      wt[0],   wt[1],   wt[2],   wt[3],   wt[4],   wt[5]);
        u[IDX2_G_R(ic2,ic1,igu)] += a * ur;
        u[IDX2_G_I(ic2,ic1,igu)] += a * ui;
      }

    }

  }

 } // acc parallel

}

//====================================================================
void mult_Gnd(real_t *restrict u, const int exu,
              real_t *restrict v, const int exv,
              real_t *restrict w, const int exw, const int nst)
{
  int ngst = NDF * nst;
  int up = nst * exu;
  int vp = nst * exv;
  int wp = nst * exw;

#pragma acc data present(u[up:ngst], v[vp:ngst], w[wp:ngst]) \
                 copyin(nst, exu, exv, exw)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for(int ist = 0; ist < nst; ++ist){
    int igw = ist + nst * exw;
    int igv = ist + nst * exv;
    int igu = ist + nst * exu;

#ifdef SU3_3RD_ROW_RECONST
    real_t vt[NDF], wt[NVC];
    for(int ic2 = 0; ic2 < NC-1; ++ic2){
      for(int ic1 = 0; ic1 < NC; ++ic1){
        int icr = 2*ic1   + NVC * ic2;
        int ici = 2*ic1+1 + NVC * ic2;
        vt[icr] = v[IDX2_G_R(ic1,ic2,igv)];
        vt[ici] = v[IDX2_G_I(ic1,ic2,igv)];
      }
    }

    vt[12] = EXT_IMG_R(vt[ 2], vt[ 3], vt[ 4], vt[ 5],
                       vt[ 8], vt[ 9], vt[10], vt[11]);
    vt[13] = EXT_IMG_I(vt[ 2], vt[ 3], vt[ 4], vt[ 5],
                       vt[ 8], vt[ 9], vt[10], vt[11]);
    vt[14] = EXT_IMG_R(vt[ 4], vt[ 5], vt[ 0], vt[ 1],
                       vt[10], vt[11], vt[ 6], vt[ 7]);
    vt[15] = EXT_IMG_I(vt[ 4], vt[ 5], vt[ 0], vt[ 1],
                       vt[10], vt[11], vt[ 6], vt[ 7]);
    vt[16] = EXT_IMG_R(vt[ 0], vt[ 1], vt[ 2], vt[ 3],
                       vt[ 6], vt[ 7], vt[ 8], vt[ 9]);
    vt[17] = EXT_IMG_I(vt[ 0], vt[ 1], vt[ 2], vt[ 3],
                       vt[ 6], vt[ 7], vt[ 8], vt[ 9]);
#else
    real_t vt[NDF], wt[NVC];
    for(int ic2 = 0; ic2 < NC; ++ic2){
      for(int ic1 = 0; ic1 < NC; ++ic1){
        int icr = 2*ic1   + NVC * ic2;
        int ici = 2*ic1+1 + NVC * ic2;
        vt[icr] = v[IDX2_G_R(ic1,ic2,igv)];
        vt[ici] = v[IDX2_G_I(ic1,ic2,igv)];
      }
    }
#endif

    for(int ic2 = 0; ic2 < NC; ++ic2){

      for(int ic1 = 0; ic1 < NC; ++ic1){
        int icr = 2*ic1;
        int ici = 2*ic1+1;
        wt[icr] =  w[IDX2_G_R(ic1,ic2,igw)];
        wt[ici] = -w[IDX2_G_I(ic1,ic2,igw)];
      }

      for(int ic1 = 0; ic1 < NC; ++ic1){
        int j = NVC * ic1;
        real_t ur, ui;
        ur = MULT_GXr(vt[0+j], vt[1+j], vt[2+j], vt[3+j], vt[4+j], vt[5+j],
                      wt[0],   wt[1],   wt[2],   wt[3],   wt[4],   wt[5]);
        ui = MULT_GXi(vt[0+j], vt[1+j], vt[2+j], vt[3+j], vt[4+j], vt[5+j],
                      wt[0],   wt[1],   wt[2],   wt[3],   wt[4],   wt[5]);
        u[IDX2_G_R(ic2,ic1,igu)] = ur;
        u[IDX2_G_I(ic2,ic1,igu)] = ui;
      }

    }

  }

 } // acc parallel

}

//====================================================================
void multadd_Gdn(real_t *restrict u, const int exu,
                 real_t *restrict v, const int exv,
                 real_t *restrict w, const int exw,
                  const real_t a, const int nst)
{
  int ngst = NDF * nst;
  int up = nst * exu;
  int vp = nst * exv;
  int wp = nst * exw;

#pragma acc data present(u[up:ngst], v[vp:ngst], w[wp:ngst]) \
                 copyin(a, nst, exu, exv, exw)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for(int ist = 0; ist < nst; ++ist){
    int igw = ist + nst * exw;
    int igv = ist + nst * exv;
    int igu = ist + nst * exu;

#ifdef SU3_3RD_ROW_RECONST
    real_t vt[NDF], wt[NVC];
    for(int ic2 = 0; ic2 < NC-1; ++ic2){
      for(int ic1 = 0; ic1 < NC; ++ic1){
        int icr = 2*ic1   + NVC * ic2;
        int ici = 2*ic1+1 + NVC * ic2;
        vt[icr] =  v[IDX2_G_R(ic2,ic1,igv)];
        vt[ici] = -v[IDX2_G_I(ic2,ic1,igv)];
      }
    }

    vt[12] = EXT_IMG_R(vt[ 2], vt[ 3], vt[ 4], vt[ 5],
                       vt[ 8], vt[ 9], vt[10], vt[11]);
    vt[13] = EXT_IMG_I(vt[ 2], vt[ 3], vt[ 4], vt[ 5],
                       vt[ 8], vt[ 9], vt[10], vt[11]);
    vt[14] = EXT_IMG_R(vt[ 4], vt[ 5], vt[ 0], vt[ 1],
                       vt[10], vt[11], vt[ 6], vt[ 7]);
    vt[15] = EXT_IMG_I(vt[ 4], vt[ 5], vt[ 0], vt[ 1],
                       vt[10], vt[11], vt[ 6], vt[ 7]);
    vt[16] = EXT_IMG_R(vt[ 0], vt[ 1], vt[ 2], vt[ 3],
                       vt[ 6], vt[ 7], vt[ 8], vt[ 9]);
    vt[17] = EXT_IMG_I(vt[ 0], vt[ 1], vt[ 2], vt[ 3],
                       vt[ 6], vt[ 7], vt[ 8], vt[ 9]);
#else
    real_t vt[NDF], wt[NVC];
    for(int ic2 = 0; ic2 < NC; ++ic2){
      for(int ic1 = 0; ic1 < NC; ++ic1){
        int icr = 2*ic1   + NVC * ic2;
        int ici = 2*ic1+1 + NVC * ic2;
        vt[icr] =  v[IDX2_G_R(ic2,ic1,igv)];
        vt[ici] = -v[IDX2_G_I(ic2,ic1,igv)];
      }
    }
#endif

    for(int ic2 = 0; ic2 < NC; ++ic2){

      for(int ic1 = 0; ic1 < NC; ++ic1){
        int icr = 2*ic1;
        int ici = 2*ic1+1;
        wt[icr] = w[IDX2_G_R(ic2,ic1,igw)];
        wt[ici] = w[IDX2_G_I(ic2,ic1,igw)];
      }

      for(int ic1 = 0; ic1 < NC; ++ic1){
        int j = NVC * ic1;
        real_t ur, ui;
        ur = MULT_GXr(vt[0+j], vt[1+j], vt[2+j], vt[3+j], vt[4+j], vt[5+j],
                      wt[0],   wt[1],   wt[2],   wt[3],   wt[4],   wt[5]);
        ui = MULT_GXi(vt[0+j], vt[1+j], vt[2+j], vt[3+j], vt[4+j], vt[5+j],
                      wt[0],   wt[1],   wt[2],   wt[3],   wt[4],   wt[5]);
        u[IDX2_G_R(ic2,ic1,igu)] += a * ur;
        u[IDX2_G_I(ic2,ic1,igu)] += a * ui;
      }

    }

  }

 } // acc parallel

}

//====================================================================
void mult_Gdn(real_t *restrict u, const int exu,
              real_t *restrict v, const int exv,
              real_t *restrict w, const int exw, const int nst)
{
  int ngst = NDF * nst;
  int up = nst * exu;
  int vp = nst * exv;
  int wp = nst * exw;

#pragma acc data present(u[up:ngst], v[vp:ngst], w[wp:ngst]) \
                 copyin(nst, exu, exv, exw)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for(int ist = 0; ist < nst; ++ist){
    int igw = ist + nst * exw;
    int igv = ist + nst * exv;
    int igu = ist + nst * exu;

#ifdef SU3_3RD_ROW_RECONST
    real_t vt[NDF], wt[NVC];
    for(int ic2 = 0; ic2 < NC-1; ++ic2){
      for(int ic1 = 0; ic1 < NC; ++ic1){
        int icr = 2*ic1   + NVC * ic2;
        int ici = 2*ic1+1 + NVC * ic2;
        vt[icr] =  v[IDX2_G_R(ic2,ic1,igv)];
        vt[ici] = -v[IDX2_G_I(ic2,ic1,igv)];
      }
    }

    vt[12] = EXT_IMG_R(vt[ 2], vt[ 3], vt[ 4], vt[ 5],
                       vt[ 8], vt[ 9], vt[10], vt[11]);
    vt[13] = EXT_IMG_I(vt[ 2], vt[ 3], vt[ 4], vt[ 5],
                       vt[ 8], vt[ 9], vt[10], vt[11]);
    vt[14] = EXT_IMG_R(vt[ 4], vt[ 5], vt[ 0], vt[ 1],
                       vt[10], vt[11], vt[ 6], vt[ 7]);
    vt[15] = EXT_IMG_I(vt[ 4], vt[ 5], vt[ 0], vt[ 1],
                       vt[10], vt[11], vt[ 6], vt[ 7]);
    vt[16] = EXT_IMG_R(vt[ 0], vt[ 1], vt[ 2], vt[ 3],
                       vt[ 6], vt[ 7], vt[ 8], vt[ 9]);
    vt[17] = EXT_IMG_I(vt[ 0], vt[ 1], vt[ 2], vt[ 3],
                       vt[ 6], vt[ 7], vt[ 8], vt[ 9]);
#else
    real_t vt[NDF], wt[NVC];
    for(int ic2 = 0; ic2 < NC; ++ic2){
      for(int ic1 = 0; ic1 < NC; ++ic1){
        int icr = 2*ic1   + NVC * ic2;
        int ici = 2*ic1+1 + NVC * ic2;
        vt[icr] =  v[IDX2_G_R(ic2,ic1,igv)];
        vt[ici] = -v[IDX2_G_I(ic2,ic1,igv)];
      }
    }
#endif

    for(int ic2 = 0; ic2 < NC; ++ic2){

      for(int ic1 = 0; ic1 < NC; ++ic1){
        int icr = 2*ic1;
        int ici = 2*ic1+1;
        wt[icr] = w[IDX2_G_R(ic2,ic1,igw)];
        wt[ici] = w[IDX2_G_I(ic2,ic1,igw)];
      }

      for(int ic1 = 0; ic1 < NC; ++ic1){
        int j = NVC * ic1;
        real_t ur, ui;
        ur = MULT_GXr(vt[0+j], vt[1+j], vt[2+j], vt[3+j], vt[4+j], vt[5+j],
                      wt[0],   wt[1],   wt[2],   wt[3],   wt[4],   wt[5]);
        ui = MULT_GXi(vt[0+j], vt[1+j], vt[2+j], vt[3+j], vt[4+j], vt[5+j],
                      wt[0],   wt[1],   wt[2],   wt[3],   wt[4],   wt[5]);
        u[IDX2_G_R(ic2,ic1,igu)] = ur;
        u[IDX2_G_I(ic2,ic1,igu)] = ui;
      }

    }

  }

 } // acc parallel

}

//====================================================================
void mult_Gdd(real_t *restrict u, const int exu,
              real_t *restrict v, const int exv,
              real_t *restrict w, const int exw, const int nst)
{
  int ngst = NDF * nst;
  int up = nst * exu;
  int vp = nst * exv;
  int wp = nst * exw;

#pragma acc data present(u[up:ngst], v[vp:ngst], w[wp:ngst]) \
                 copyin(nst, exu, exv, exw)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for(int ist = 0; ist < nst; ++ist){
    int igw = ist + nst * exw;
    int igv = ist + nst * exv;
    int igu = ist + nst * exu;

#ifdef SU3_3RD_ROW_RECONST
    real_t vt[NDF], wt[NVC];
    for(int ic2 = 0; ic2 < NC-1; ++ic2){
      for(int ic1 = 0; ic1 < NC; ++ic1){
        int icr = 2*ic1   + NVC * ic2;
        int ici = 2*ic1+1 + NVC * ic2;
        vt[icr] =  v[IDX2_G_R(ic2,ic1,igv)];
        vt[ici] = -v[IDX2_G_I(ic2,ic1,igv)];
      }
    }

    vt[12] = EXT_IMG_R(vt[ 2], vt[ 3], vt[ 4], vt[ 5],
                       vt[ 8], vt[ 9], vt[10], vt[11]);
    vt[13] = EXT_IMG_I(vt[ 2], vt[ 3], vt[ 4], vt[ 5],
                       vt[ 8], vt[ 9], vt[10], vt[11]);
    vt[14] = EXT_IMG_R(vt[ 4], vt[ 5], vt[ 0], vt[ 1],
                       vt[10], vt[11], vt[ 6], vt[ 7]);
    vt[15] = EXT_IMG_I(vt[ 4], vt[ 5], vt[ 0], vt[ 1],
                       vt[10], vt[11], vt[ 6], vt[ 7]);
    vt[16] = EXT_IMG_R(vt[ 0], vt[ 1], vt[ 2], vt[ 3],
                       vt[ 6], vt[ 7], vt[ 8], vt[ 9]);
    vt[17] = EXT_IMG_I(vt[ 0], vt[ 1], vt[ 2], vt[ 3],
                       vt[ 6], vt[ 7], vt[ 8], vt[ 9]);
#else
    real_t vt[NDF], wt[NVC];
    for(int ic2 = 0; ic2 < NC; ++ic2){
      for(int ic1 = 0; ic1 < NC; ++ic1){
        int icr = 2*ic1   + NVC * ic2;
        int ici = 2*ic1+1 + NVC * ic2;
        vt[icr] =  v[IDX2_G_R(ic2,ic1,igv)];
        vt[ici] = -v[IDX2_G_I(ic2,ic1,igv)];
      }
    }
#endif

    for(int ic2 = 0; ic2 < NC; ++ic2){

      for(int ic1 = 0; ic1 < NC; ++ic1){
        int icr = 2*ic1;
        int ici = 2*ic1+1;
        wt[icr] =  w[IDX2_G_R(ic1,ic2,igw)];
        wt[ici] = -w[IDX2_G_I(ic1,ic2,igw)];
      }

      for(int ic1 = 0; ic1 < NC; ++ic1){
        int j = NVC * ic1;
        real_t ur, ui;
        ur = MULT_GXr(vt[0+j], vt[1+j], vt[2+j], vt[3+j], vt[4+j], vt[5+j],
                      wt[0],   wt[1],   wt[2],   wt[3],   wt[4],   wt[5]);
        ui = MULT_GXi(vt[0+j], vt[1+j], vt[2+j], vt[3+j], vt[4+j], vt[5+j],
                      wt[0],   wt[1],   wt[2],   wt[3],   wt[4],   wt[5]);
        u[IDX2_G_R(ic2,ic1,igu)] = ur;
        u[IDX2_G_I(ic2,ic1,igu)] = ui;
      }

    }

  }

 } // acc parallel

}

//====================================================================
// anti-hermitian
void ah_G(real_t *u, const int ex, const int nst)
{
  int ngst = NDF * nst;
  int up   = NDF * nst * ex;

#pragma acc data present(u[up:ngst]) copyin(ex, nst)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for(int ist = 0; ist < nst; ++ist){
    int igu = ist + nst * ex;

    real_t ut[NDF];
    for(int ic2 = 0; ic2 < NC; ++ic2){
      for(int ic1 = 0; ic1 < NC; ++ic1){
        int icr = 2*ic1   + NVC * ic2;
        int ici = 2*ic1+1 + NVC * ic2;
        ut[icr] = u[IDX2_G_R(ic1,ic2,igu)];
        ut[ici] = u[IDX2_G_I(ic1,ic2,igu)];
      }
    }

    real_t vt[NDF];
    vt[ 0] = 0.0;
    vt[ 1] = ut[ 1];
    vt[ 2] = 0.5 *(ut[ 2] - ut[ 6]);
    vt[ 3] = 0.5 *(ut[ 3] + ut[ 7]);
    vt[ 4] = 0.5 *(ut[ 4] - ut[12]);
    vt[ 5] = 0.5 *(ut[ 5] + ut[13]);
    vt[ 6] = -vt[ 2];
    vt[ 7] =  vt[ 3];
    vt[ 8] = 0.0;
    vt[ 9] = ut[ 9];
    vt[10] = 0.5 *(ut[10] - ut[14]);
    vt[11] = 0.5 *(ut[11] + ut[15]);
    vt[12] = -vt[ 4];
    vt[13] =  vt[ 5];
    vt[14] = -vt[10];
    vt[15] =  vt[11];
    vt[16] = 0.0;
    vt[17] = ut[17];

    for(int ic2 = 0; ic2 < NC; ++ic2){
      for(int ic1 = 0; ic1 < NC; ++ic1){
        int icr = 2*ic1   + NVC * ic2;
        int ici = 2*ic1+1 + NVC * ic2;
        u[IDX2_G_R(ic1,ic2,igu)] = vt[icr];
        u[IDX2_G_I(ic1,ic2,igu)] = vt[ici];
      }
    }

  }

 } // acc parallel

}

//====================================================================
// anti-hermitian traceless
void at_G(real_t *u, const int ex, const int nst)
{
  int ngst = NDF * nst;
  int up   = NDF * nst * ex;

#pragma acc data present(u[up:ngst]) copyin(ex, nst)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for(int ist = 0; ist < nst; ++ist){
    int igu = ist + nst * ex;

    real_t ut[NDF];
    for(int ic2 = 0; ic2 < NC; ++ic2){
      for(int ic1 = 0; ic1 < NC; ++ic1){
        int icr = 2*ic1   + NVC * ic2;
        int ici = 2*ic1+1 + NVC * ic2;
        ut[icr] = u[IDX2_G_R(ic1,ic2,igu)];
        ut[ici] = u[IDX2_G_I(ic1,ic2,igu)];
      }
    }

    real_t vt[NDF];
    real_t tr = (ut[ 1] + ut[9] + ut[17])/3.0;
    vt[ 0] = 0.0;
    vt[ 1] = ut[ 1] - tr;
    vt[ 2] = 0.5 *(ut[ 2] - ut[ 6]);
    vt[ 3] = 0.5 *(ut[ 3] + ut[ 7]);
    vt[ 4] = 0.5 *(ut[ 4] - ut[12]);
    vt[ 5] = 0.5 *(ut[ 5] + ut[13]);
    vt[ 6] = -vt[ 2];
    vt[ 7] =  vt[ 3];
    vt[ 8] = 0.0;
    vt[ 9] = ut[ 9] - tr;
    vt[10] = 0.5 *(ut[10] - ut[14]);
    vt[11] = 0.5 *(ut[11] + ut[15]);
    vt[12] = -vt[ 4];
    vt[13] =  vt[ 5];
    vt[14] = -vt[10];
    vt[15] =  vt[11];
    vt[16] = 0.0;
    vt[17] = ut[17] - tr;

    for(int ic2 = 0; ic2 < NC; ++ic2){
      for(int ic1 = 0; ic1 < NC; ++ic1){
        int icr = 2*ic1   + NVC * ic2;
        int ici = 2*ic1+1 + NVC * ic2;
        u[IDX2_G_R(ic1,ic2,igu)] = vt[icr];
        u[IDX2_G_I(ic1,ic2,igu)] = vt[ici];
      }
    }

  }

 } // acc parallel

}

//====================================================================
void add_unit(real_t *u, const int ex, const real_t a, const int nst)
{         // u = u + a * I (I: unit matrix)

  int ngst = NDF * nst;
  int up   = NDF * nst * ex;

#pragma acc data present(u[up:ngst]) copyin(a, ex, nst)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for(int ist = 0; ist < nst; ++ist){
    int igu = ist + nst * ex;

    for(int ic = 0; ic < NC; ++ic){
      real_t ut = u[IDX2_G_R(ic,ic,igu)];
      ut += a;
      u[IDX2_G_R(ic,ic,igu)] = ut;
    }

  }

 } // acc parallel

}

//====================================================================
void inverse_dag(real_t *restrict u, const int ex1,
                 const real_t *restrict v, const int ex2,
                 const int nst)
{ // calculate u = ((v)^inv)^dag

  int ngst = NDF * nst;
  int up   = NDF * nst * ex1;
  int vp   = NDF * nst * ex2;

#pragma acc data present(v[vp:ngst], u[up:ngst]) copyin(ex1, ex2, nst)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for(int ist = 0; ist < nst; ++ist){

    int ist1 = ist + nst * ex1;
    int ist2 = ist + nst * ex2;

    real_t vt[18], wt[18], ut[18];

    vt[ 0] = v[IDX2_G_R(0, 0, ist2)];
    vt[ 1] = v[IDX2_G_I(0, 0, ist2)];
    vt[ 2] = v[IDX2_G_R(1, 0, ist2)];
    vt[ 3] = v[IDX2_G_I(1, 0, ist2)];
    vt[ 4] = v[IDX2_G_R(2, 0, ist2)];
    vt[ 5] = v[IDX2_G_I(2, 0, ist2)];

    vt[ 6] = v[IDX2_G_R(0, 1, ist2)];
    vt[ 7] = v[IDX2_G_I(0, 1, ist2)];
    vt[ 8] = v[IDX2_G_R(1, 1, ist2)];
    vt[ 9] = v[IDX2_G_I(1, 1, ist2)];
    vt[10] = v[IDX2_G_R(2, 1, ist2)];
    vt[11] = v[IDX2_G_I(2, 1, ist2)];

    vt[12] = v[IDX2_G_R(0, 2, ist2)];
    vt[13] = v[IDX2_G_I(0, 2, ist2)];
    vt[14] = v[IDX2_G_R(1, 2, ist2)];
    vt[15] = v[IDX2_G_I(1, 2, ist2)];
    vt[16] = v[IDX2_G_R(2, 2, ist2)];
    vt[17] = v[IDX2_G_I(2, 2, ist2)];

    wt[ 0] =  MULT_R(vt[ 8], vt[ 9], vt[16], vt[17])
            - MULT_R(vt[10], vt[11], vt[14], vt[15]);
    wt[ 1] =  MULT_I(vt[ 8], vt[ 9], vt[16], vt[17])
            - MULT_I(vt[10], vt[11], vt[14], vt[15]);

    wt[ 2] =  MULT_R(vt[10], vt[11], vt[12], vt[13])
            - MULT_R(vt[ 6], vt[ 7], vt[16], vt[17]);
    wt[ 3] =  MULT_I(vt[10], vt[11], vt[12], vt[13])
            - MULT_I(vt[ 6], vt[ 7], vt[16], vt[17]);

    wt[ 4] =  MULT_R(vt[ 6], vt[ 7], vt[14], vt[15])
            - MULT_R(vt[ 8], vt[ 9], vt[12], vt[13]);
    wt[ 5] =  MULT_I(vt[ 6], vt[ 7], vt[14], vt[15])
            - MULT_I(vt[ 8], vt[ 9], vt[12], vt[13]);

    wt[ 6] =  MULT_R(vt[14], vt[15], vt[ 4], vt[ 5])
            - MULT_R(vt[16], vt[17], vt[ 2], vt[ 3]);
    wt[ 7] =  MULT_I(vt[14], vt[15], vt[ 4], vt[ 5])
            - MULT_I(vt[16], vt[17], vt[ 2], vt[ 3]);

    wt[ 8] =  MULT_R(vt[16], vt[17], vt[ 0], vt[ 1])
            - MULT_R(vt[12], vt[13], vt[ 4], vt[ 5]);
    wt[ 9] =  MULT_I(vt[16], vt[17], vt[ 0], vt[ 1])
            - MULT_I(vt[12], vt[13], vt[ 4], vt[ 5]);

    wt[10] =  MULT_R(vt[12], vt[13], vt[ 2], vt[ 3])
            - MULT_R(vt[14], vt[15], vt[ 0], vt[ 1]);
    wt[11] =  MULT_I(vt[12], vt[13], vt[ 2], vt[ 3])
            - MULT_I(vt[14], vt[15], vt[ 0], vt[ 1]);

    wt[12] =  MULT_R(vt[ 2], vt[ 3], vt[10], vt[11])
            - MULT_R(vt[ 4], vt[ 5], vt[ 8], vt[ 9]);
    wt[13] =  MULT_I(vt[ 2], vt[ 3], vt[10], vt[11])
            - MULT_I(vt[ 4], vt[ 5], vt[ 8], vt[ 9]);

    wt[14] =  MULT_R(vt[ 4], vt[ 5], vt[ 6], vt[ 7])
            - MULT_R(vt[ 0], vt[ 1], vt[10], vt[11]);
    wt[15] =  MULT_I(vt[ 4], vt[ 5], vt[ 6], vt[ 7])
            - MULT_I(vt[ 0], vt[ 1], vt[10], vt[11]);

    wt[16] =  MULT_R(vt[ 0], vt[ 1], vt[ 8], vt[ 9])
            - MULT_R(vt[ 2], vt[ 3], vt[ 6], vt[ 7]);
    wt[17] =  MULT_I(vt[ 0], vt[ 1], vt[ 8], vt[ 9])
            - MULT_I(vt[ 2], vt[ 3], vt[ 6], vt[ 7]);

    real_t det_r =  MULT_R(vt[0], vt[1], wt[0], wt[1])
                  + MULT_R(vt[2], vt[3], wt[2], wt[3])
                  + MULT_R(vt[4], vt[5], wt[4], wt[5]);
    real_t det_i =  MULT_I(vt[0], vt[1], wt[0], wt[1])
                  + MULT_I(vt[2], vt[3], wt[2], wt[3])
                  + MULT_I(vt[4], vt[5], wt[4], wt[5]);

    real_t det2 = det_r * det_r + det_i * det_i;
    real_t detinv_r =   det_r/det2;
    real_t detinv_i = - det_i/det2;

    ut[ 0] = MULT_R(wt[ 0], wt[ 1], detinv_r, detinv_i);
    ut[ 1] = MULT_I(wt[ 0], wt[ 1], detinv_r, detinv_i);

    ut[ 2] = MULT_R(wt[ 6], wt[ 7], detinv_r, detinv_i);
    ut[ 3] = MULT_I(wt[ 6], wt[ 7], detinv_r, detinv_i);

    ut[ 4] = MULT_R(wt[12], wt[13], detinv_r, detinv_i);
    ut[ 5] = MULT_I(wt[12], wt[13], detinv_r, detinv_i);

    ut[ 6] = MULT_R(wt[ 2], wt[ 3], detinv_r, detinv_i);
    ut[ 7] = MULT_I(wt[ 2], wt[ 3], detinv_r, detinv_i);

    ut[ 8] = MULT_R(wt[ 8], wt[ 9], detinv_r, detinv_i);
    ut[ 9] = MULT_I(wt[ 8], wt[ 9], detinv_r, detinv_i);

    ut[10] = MULT_R(wt[14], wt[15], detinv_r, detinv_i);
    ut[11] = MULT_I(wt[14], wt[15], detinv_r, detinv_i);

    ut[12] = MULT_R(wt[ 4], wt[ 5], detinv_r, detinv_i);
    ut[13] = MULT_I(wt[ 4], wt[ 5], detinv_r, detinv_i);

    ut[14] = MULT_R(wt[10], wt[11], detinv_r, detinv_i);
    ut[15] = MULT_I(wt[10], wt[11], detinv_r, detinv_i);

    ut[16] = MULT_R(wt[16], wt[17], detinv_r, detinv_i);
    ut[17] = MULT_I(wt[16], wt[17], detinv_r, detinv_i);

    u[IDX2_G_R(0, 0, ist1)] =  ut[ 0];
    u[IDX2_G_I(0, 0, ist1)] = -ut[ 1];
    u[IDX2_G_R(1, 0, ist1)] =  ut[ 6];
    u[IDX2_G_I(1, 0, ist1)] = -ut[ 7];
    u[IDX2_G_R(2, 0, ist1)] =  ut[12];
    u[IDX2_G_I(2, 0, ist1)] = -ut[13];

    u[IDX2_G_R(0, 1, ist1)] =  ut[ 2];
    u[IDX2_G_I(0, 1, ist1)] = -ut[ 3];
    u[IDX2_G_R(1, 1, ist1)] =  ut[ 8];
    u[IDX2_G_I(1, 1, ist1)] = -ut[ 9];
    u[IDX2_G_R(2, 1, ist1)] =  ut[14];
    u[IDX2_G_I(2, 1, ist1)] = -ut[15];

    u[IDX2_G_R(0, 2, ist1)] =  ut[ 4];
    u[IDX2_G_I(0, 2, ist1)] = -ut[ 5];
    u[IDX2_G_R(1, 2, ist1)] =  ut[10];
    u[IDX2_G_I(1, 2, ist1)] = -ut[11];
    u[IDX2_G_R(2, 2, ist1)] =  ut[16];
    u[IDX2_G_I(2, 2, ist1)] = -ut[17];

  }

 } // acc parallel

}

//============================================================END=====
#endif
