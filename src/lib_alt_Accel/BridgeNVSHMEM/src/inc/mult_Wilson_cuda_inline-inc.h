/*!
      @file    mult_Wilson_openacc_inline-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2606 $
*/

#ifndef MULT_WILSON_OPENACC_INLINE_INCLUDED
#define MULT_WILSON_OPENACC_INLINE_INCLUDED

#define MULT_UV_R(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)   (u0*v0-u1*v1 + u2*v2-u3*v3 + u4*v4-u5*v5)
#define MULT_UV_I(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)   (u0*v1+u1*v0 + u2*v3+u3*v2 + u4*v5+u5*v4)

#define MULT_GXr(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)   (u0*v0-u1*v1 + u2*v2-u3*v3 + u4*v4-u5*v5)
#define MULT_GXi(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)   (u0*v1+u1*v0 + u2*v3+u3*v2 + u4*v5+u5*v4)
#define MULT_GDXr(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)  (u0*v0+u1*v1 + u2*v2+u3*v3 + u4*v4+u5*v5)
#define MULT_GDXi(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)  (u0*v1-u1*v0 + u2*v3-u3*v2 + u4*v5-u5*v4)

#define EXT_IMG_R(v1r,v1i,v2r,v2i,w1r,w1i,w2r,w2i)  (v1r*w2r - v1i*w2i - v2r*w1r + v2i*w1i)
#define EXT_IMG_I(v1r,v1i,v2r,v2i,w1r,w1i,w2r,w2i)  (- v1r*w2i - v1i*w2r + v2r*w1i + v2i*w1r)


//====================================================================
__device__ inline void load_u(real_t* ut, real_t* up, int site)
{
#ifdef SU3_3RD_ROW_RECONST
  for(int idf = 0; idf < 2*NVC; ++idf){
    ut[idf] = up[IDX2(NDF, idf, site)];
  }

  ut[12] = EXT_IMG_R(ut[2], ut[3], ut[4], ut[5], ut[ 8], ut[ 9], ut[10], ut[11]);
  ut[13] = EXT_IMG_I(ut[2], ut[3], ut[4], ut[5], ut[ 8], ut[ 9], ut[10], ut[11]);
  ut[14] = EXT_IMG_R(ut[4], ut[5], ut[0], ut[1], ut[10], ut[11], ut[ 6], ut[ 7]);
  ut[15] = EXT_IMG_I(ut[4], ut[5], ut[0], ut[1], ut[10], ut[11], ut[ 6], ut[ 7]);
  ut[16] = EXT_IMG_R(ut[0], ut[1], ut[2], ut[3], ut[ 6], ut[ 7], ut[ 8], ut[ 9]);
  ut[17] = EXT_IMG_I(ut[0], ut[1], ut[2], ut[3], ut[ 6], ut[ 7], ut[ 8], ut[ 9]);
#else
  for(int idf = 0; idf < NDF; ++idf){
    ut[idf] = up[IDX2(NDF, idf, site)];
  }
#endif

}

//====================================================================
__device__ inline void mult_wilson_xpb(real_t* vt, real_t* ut, real_t* wt)
{
  real_t vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5;
  real_t vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5;

  vt1_0 = wt[ID1 + 0] - wt[ID4 + 1];
  vt1_1 = wt[ID1 + 1] + wt[ID4 + 0];
  vt1_2 = wt[ID1 + 2] - wt[ID4 + 3];
  vt1_3 = wt[ID1 + 3] + wt[ID4 + 2];
  vt1_4 = wt[ID1 + 4] - wt[ID4 + 5];
  vt1_5 = wt[ID1 + 5] + wt[ID4 + 4];

  vt2_0 = wt[ID2 + 0] - wt[ID3 + 1];
  vt2_1 = wt[ID2 + 1] + wt[ID3 + 0];
  vt2_2 = wt[ID2 + 2] - wt[ID3 + 3];
  vt2_3 = wt[ID2 + 3] + wt[ID3 + 2];
  vt2_4 = wt[ID2 + 4] - wt[ID3 + 5];
  vt2_5 = wt[ID2 + 5] + wt[ID3 + 4];

  for(int ic = 0; ic < NC; ++ic){

    real_t wt1r, wt1i, wt2r, wt2i;
    real_t* up = &ut[NVC * ic];

    wt1r = MULT_GXr(up[0], up[1], up[2], up[3], up[4], up[5],
                    vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
    wt1i = MULT_GXi(up[0], up[1], up[2], up[3], up[4], up[5],
                    vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
    wt2r = MULT_GXr(up[0], up[1], up[2], up[3], up[4], up[5],
                    vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
    wt2i = MULT_GXi(up[0], up[1], up[2], up[3], up[4], up[5],
                    vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

    vt[ID1 + 2*ic  ] +=  wt1r;
    vt[ID1 + 2*ic+1] +=  wt1i;
    vt[ID2 + 2*ic  ] +=  wt2r;
    vt[ID2 + 2*ic+1] +=  wt2i;
    vt[ID3 + 2*ic  ] +=  wt2i;
    vt[ID3 + 2*ic+1] += -wt2r;
    vt[ID4 + 2*ic  ] +=  wt1i;
    vt[ID4 + 2*ic+1] += -wt1r;
  }

}

//====================================================================
__device__ inline void mult_wilson_xp1(real_t* vt, real_t* wt)
{
  vt[ 0] = wt[ID1 + 0] - wt[ID4 + 1];
  vt[ 1] = wt[ID1 + 1] + wt[ID4 + 0];
  vt[ 2] = wt[ID1 + 2] - wt[ID4 + 3];
  vt[ 3] = wt[ID1 + 3] + wt[ID4 + 2];
  vt[ 4] = wt[ID1 + 4] - wt[ID4 + 5];
  vt[ 5] = wt[ID1 + 5] + wt[ID4 + 4];

  vt[ 6] = wt[ID2 + 0] - wt[ID3 + 1];
  vt[ 7] = wt[ID2 + 1] + wt[ID3 + 0];
  vt[ 8] = wt[ID2 + 2] - wt[ID3 + 3];
  vt[ 9] = wt[ID2 + 3] + wt[ID3 + 2];
  vt[10] = wt[ID2 + 4] - wt[ID3 + 5];
  vt[11] = wt[ID2 + 5] + wt[ID3 + 4];
}

//====================================================================
__device__ inline void mult_wilson_xp2(real_t* vt, real_t* ut, real_t* wt)
{
  real_t* wt1 = &wt[0];
  real_t* wt2 = &wt[NVC];

  for(int ic = 0; ic < NC; ++ic){

    real_t wt1r, wt1i, wt2r, wt2i;
    real_t* up = &ut[NVC * ic];

    wt1r = MULT_GXr( up[0],  up[1],  up[2],  up[3],  up[4],  up[5],
                    wt1[0], wt1[1], wt1[2], wt1[3], wt1[4], wt1[5]);
    wt1i = MULT_GXi( up[0],  up[1],  up[2],  up[3],  up[4],  up[5],
                    wt1[0], wt1[1], wt1[2], wt1[3], wt1[4], wt1[5]);
    wt2r = MULT_GXr( up[0],  up[1],  up[2],  up[3],  up[4],  up[5],
                    wt2[0], wt2[1], wt2[2], wt2[3], wt2[4], wt2[5]);
    wt2i = MULT_GXi( up[0],  up[1],  up[2],  up[3],  up[4],  up[5],
                    wt2[0], wt2[1], wt2[2], wt2[3], wt2[4], wt2[5]);

    vt[ID1 + 2*ic  ] +=  wt1r;
    vt[ID1 + 2*ic+1] +=  wt1i;
    vt[ID2 + 2*ic  ] +=  wt2r;
    vt[ID2 + 2*ic+1] +=  wt2i;
    vt[ID3 + 2*ic  ] +=  wt2i;
    vt[ID3 + 2*ic+1] += -wt2r;
    vt[ID4 + 2*ic  ] +=  wt1i;
    vt[ID4 + 2*ic+1] += -wt1r;
  }

}

//====================================================================
__device__ inline void mult_wilson_xmb(real_t* vt, real_t* ut, real_t* wt)
{
  real_t vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5;
  real_t vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5;

  vt1_0 = wt[ID1 + 0] + wt[ID4 + 1];
  vt1_1 = wt[ID1 + 1] - wt[ID4 + 0];
  vt1_2 = wt[ID1 + 2] + wt[ID4 + 3];
  vt1_3 = wt[ID1 + 3] - wt[ID4 + 2];
  vt1_4 = wt[ID1 + 4] + wt[ID4 + 5];
  vt1_5 = wt[ID1 + 5] - wt[ID4 + 4];

  vt2_0 = wt[ID2 + 0] + wt[ID3 + 1];
  vt2_1 = wt[ID2 + 1] - wt[ID3 + 0];
  vt2_2 = wt[ID2 + 2] + wt[ID3 + 3];
  vt2_3 = wt[ID2 + 3] - wt[ID3 + 2];
  vt2_4 = wt[ID2 + 4] + wt[ID3 + 5];
  vt2_5 = wt[ID2 + 5] - wt[ID3 + 4];

  for(int ic = 0; ic < NC; ++ic){

    real_t wt1r, wt1i, wt2r, wt2i;
    real_t* up = &ut[2 * ic];

    wt1r = MULT_GDXr(up[0], up[1], up[6], up[7], up[12], up[13],
                     vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
    wt1i = MULT_GDXi(up[0], up[1], up[6], up[7], up[12], up[13],
                     vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
    wt2r = MULT_GDXr(up[0], up[1], up[6], up[7], up[12], up[13],
                     vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
    wt2i = MULT_GDXi(up[0], up[1], up[6], up[7], up[12], up[13],
                     vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

    vt[ID1 + 2*ic  ] +=  wt1r;
    vt[ID1 + 2*ic+1] +=  wt1i;
    vt[ID2 + 2*ic  ] +=  wt2r;
    vt[ID2 + 2*ic+1] +=  wt2i;
    vt[ID3 + 2*ic  ] += -wt2i;
    vt[ID3 + 2*ic+1] +=  wt2r;
    vt[ID4 + 2*ic  ] += -wt1i;
    vt[ID4 + 2*ic+1] +=  wt1r;
  }

}

//====================================================================
__device__ inline void mult_wilson_xm1(real_t* vt, real_t* ut, real_t* wt)
{
  real_t vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5;
  real_t vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5;

  vt1_0 = wt[ID1 + 0] + wt[ID4 + 1];
  vt1_1 = wt[ID1 + 1] - wt[ID4 + 0];
  vt1_2 = wt[ID1 + 2] + wt[ID4 + 3];
  vt1_3 = wt[ID1 + 3] - wt[ID4 + 2];
  vt1_4 = wt[ID1 + 4] + wt[ID4 + 5];
  vt1_5 = wt[ID1 + 5] - wt[ID4 + 4];

  vt2_0 = wt[ID2 + 0] + wt[ID3 + 1];
  vt2_1 = wt[ID2 + 1] - wt[ID3 + 0];
  vt2_2 = wt[ID2 + 2] + wt[ID3 + 3];
  vt2_3 = wt[ID2 + 3] - wt[ID3 + 2];
  vt2_4 = wt[ID2 + 4] + wt[ID3 + 5];
  vt2_5 = wt[ID2 + 5] - wt[ID3 + 4];

  for(int ic = 0; ic < NC; ++ic){
    real_t* up = &ut[2 * ic];

    vt[2*ic] =
           MULT_GDXr(up[0], up[1], up[6], up[7], up[12], up[13],
                     vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
    vt[2*ic+1] =
           MULT_GDXi(up[0], up[1], up[6], up[7], up[12], up[13],
                     vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
    vt[2*ic + NVC] =
           MULT_GDXr(up[0], up[1], up[6], up[7], up[12], up[13],
                     vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
    vt[2*ic+1 + NVC] =
            MULT_GDXi(up[0], up[1], up[6], up[7], up[12], up[13],
                     vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
  }

}

//====================================================================
__device__ inline void mult_wilson_xm2(real_t* vt, real_t* wt)
{
  real_t* wt1 = &wt[0];
  real_t* wt2 = &wt[NVC];
  for(int ic = 0; ic < NC; ++ic){
    vt[ID1 + 2*ic  ] +=  wt1[2*ic];
    vt[ID1 + 2*ic+1] +=  wt1[2*ic+1];
    vt[ID2 + 2*ic  ] +=  wt2[2*ic];
    vt[ID2 + 2*ic+1] +=  wt2[2*ic+1];
    vt[ID3 + 2*ic  ] += -wt2[2*ic+1];
    vt[ID3 + 2*ic+1] +=  wt2[2*ic];
    vt[ID4 + 2*ic  ] += -wt1[2*ic+1];
    vt[ID4 + 2*ic+1] +=  wt1[2*ic];
  }

}

//====================================================================
__device__ inline void mult_wilson_ypb(real_t* vt, real_t* ut, real_t* wt)
{
  real_t vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5;
  real_t vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5;

  vt1_0 = wt[ID1 + 0] + wt[ID4 + 0];
  vt1_1 = wt[ID1 + 1] + wt[ID4 + 1];
  vt1_2 = wt[ID1 + 2] + wt[ID4 + 2];
  vt1_3 = wt[ID1 + 3] + wt[ID4 + 3];
  vt1_4 = wt[ID1 + 4] + wt[ID4 + 4];
  vt1_5 = wt[ID1 + 5] + wt[ID4 + 5];

  vt2_0 = wt[ID2 + 0] - wt[ID3 + 0];
  vt2_1 = wt[ID2 + 1] - wt[ID3 + 1];
  vt2_2 = wt[ID2 + 2] - wt[ID3 + 2];
  vt2_3 = wt[ID2 + 3] - wt[ID3 + 3];
  vt2_4 = wt[ID2 + 4] - wt[ID3 + 4];
  vt2_5 = wt[ID2 + 5] - wt[ID3 + 5];

  for(int ic = 0; ic < NC; ++ic){

    real_t wt1r, wt1i, wt2r, wt2i;
    real_t* up = &ut[NVC * ic];

    wt1r = MULT_GXr(up[0], up[1], up[2], up[3], up[4], up[5],
                    vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
    wt1i = MULT_GXi(up[0], up[1], up[2], up[3], up[4], up[5],
                    vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
    wt2r = MULT_GXr(up[0], up[1], up[2], up[3], up[4], up[5],
                    vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
    wt2i = MULT_GXi(up[0], up[1], up[2], up[3], up[4], up[5],
                    vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

    vt[ID1 + 2*ic  ] +=  wt1r;
    vt[ID1 + 2*ic+1] +=  wt1i;
    vt[ID2 + 2*ic  ] +=  wt2r;
    vt[ID2 + 2*ic+1] +=  wt2i;
    vt[ID3 + 2*ic  ] += -wt2r;
    vt[ID3 + 2*ic+1] += -wt2i;
    vt[ID4 + 2*ic  ] +=  wt1r;
    vt[ID4 + 2*ic+1] +=  wt1i;
  }

}

//====================================================================
__device__ inline void mult_wilson_yp1(real_t* vt, real_t* wt)
{
  vt[ 0] = wt[ID1 + 0] + wt[ID4 + 0];
  vt[ 1] = wt[ID1 + 1] + wt[ID4 + 1];
  vt[ 2] = wt[ID1 + 2] + wt[ID4 + 2];
  vt[ 3] = wt[ID1 + 3] + wt[ID4 + 3];
  vt[ 4] = wt[ID1 + 4] + wt[ID4 + 4];
  vt[ 5] = wt[ID1 + 5] + wt[ID4 + 5];

  vt[ 6] = wt[ID2 + 0] - wt[ID3 + 0];
  vt[ 7] = wt[ID2 + 1] - wt[ID3 + 1];
  vt[ 8] = wt[ID2 + 2] - wt[ID3 + 2];
  vt[ 9] = wt[ID2 + 3] - wt[ID3 + 3];
  vt[10] = wt[ID2 + 4] - wt[ID3 + 4];
  vt[11] = wt[ID2 + 5] - wt[ID3 + 5];
}

//====================================================================
__device__ inline void mult_wilson_yp2(real_t* vt, real_t* ut, real_t* wt)
{
  real_t* wt1 = &wt[0];
  real_t* wt2 = &wt[NVC];

  for(int ic = 0; ic < NC; ++ic){

    real_t wt1r, wt1i, wt2r, wt2i;
    real_t* up = &ut[NVC * ic];

    wt1r = MULT_GXr( up[0],  up[1],  up[2],  up[3],  up[4],  up[5],
                    wt1[0], wt1[1], wt1[2], wt1[3], wt1[4], wt1[5]);
    wt1i = MULT_GXi( up[0],  up[1],  up[2],  up[3],  up[4],  up[5],
                    wt1[0], wt1[1], wt1[2], wt1[3], wt1[4], wt1[5]);
    wt2r = MULT_GXr( up[0],  up[1],  up[2],  up[3],  up[4],  up[5],
                    wt2[0], wt2[1], wt2[2], wt2[3], wt2[4], wt2[5]);
    wt2i = MULT_GXi( up[0],  up[1],  up[2],  up[3],  up[4],  up[5],
                    wt2[0], wt2[1], wt2[2], wt2[3], wt2[4], wt2[5]);

    vt[ID1 + 2*ic  ] +=  wt1r;
    vt[ID1 + 2*ic+1] +=  wt1i;
    vt[ID2 + 2*ic  ] +=  wt2r;
    vt[ID2 + 2*ic+1] +=  wt2i;
    vt[ID3 + 2*ic  ] += -wt2r;
    vt[ID3 + 2*ic+1] += -wt2i;
    vt[ID4 + 2*ic  ] +=  wt1r;
    vt[ID4 + 2*ic+1] +=  wt1i;
  }

}

//====================================================================
__device__ inline void mult_wilson_ymb(real_t* vt, real_t* ut, real_t* wt)
{
  real_t vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5;
  real_t vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5;

  vt1_0 = wt[ID1 + 0] - wt[ID4 + 0];
  vt1_1 = wt[ID1 + 1] - wt[ID4 + 1];
  vt1_2 = wt[ID1 + 2] - wt[ID4 + 2];
  vt1_3 = wt[ID1 + 3] - wt[ID4 + 3];
  vt1_4 = wt[ID1 + 4] - wt[ID4 + 4];
  vt1_5 = wt[ID1 + 5] - wt[ID4 + 5];

  vt2_0 = wt[ID2 + 0] + wt[ID3 + 0];
  vt2_1 = wt[ID2 + 1] + wt[ID3 + 1];
  vt2_2 = wt[ID2 + 2] + wt[ID3 + 2];
  vt2_3 = wt[ID2 + 3] + wt[ID3 + 3];
  vt2_4 = wt[ID2 + 4] + wt[ID3 + 4];
  vt2_5 = wt[ID2 + 5] + wt[ID3 + 5];

  for(int ic = 0; ic < NC; ++ic){

    real_t wt1r, wt1i, wt2r, wt2i;
    real_t* up = &ut[2 * ic];

    wt1r = MULT_GDXr(up[0], up[1], up[6], up[7], up[12], up[13],
                     vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
    wt1i = MULT_GDXi(up[0], up[1], up[6], up[7], up[12], up[13],
                     vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
    wt2r = MULT_GDXr(up[0], up[1], up[6], up[7], up[12], up[13],
                     vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
    wt2i = MULT_GDXi(up[0], up[1], up[6], up[7], up[12], up[13],
                     vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

    vt[ID1 + 2*ic  ] +=  wt1r;
    vt[ID1 + 2*ic+1] +=  wt1i;
    vt[ID2 + 2*ic  ] +=  wt2r;
    vt[ID2 + 2*ic+1] +=  wt2i;
    vt[ID3 + 2*ic  ] +=  wt2r;
    vt[ID3 + 2*ic+1] +=  wt2i;
    vt[ID4 + 2*ic  ] += -wt1r;
    vt[ID4 + 2*ic+1] += -wt1i;
  }

}

//====================================================================
__device__ inline void mult_wilson_ym1(real_t* vt, real_t* ut, real_t* wt)
{
  real_t vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5;
  real_t vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5;

  vt1_0 = wt[ID1 + 0] - wt[ID4 + 0];
  vt1_1 = wt[ID1 + 1] - wt[ID4 + 1];
  vt1_2 = wt[ID1 + 2] - wt[ID4 + 2];
  vt1_3 = wt[ID1 + 3] - wt[ID4 + 3];
  vt1_4 = wt[ID1 + 4] - wt[ID4 + 4];
  vt1_5 = wt[ID1 + 5] - wt[ID4 + 5];

  vt2_0 = wt[ID2 + 0] + wt[ID3 + 0];
  vt2_1 = wt[ID2 + 1] + wt[ID3 + 1];
  vt2_2 = wt[ID2 + 2] + wt[ID3 + 2];
  vt2_3 = wt[ID2 + 3] + wt[ID3 + 3];
  vt2_4 = wt[ID2 + 4] + wt[ID3 + 4];
  vt2_5 = wt[ID2 + 5] + wt[ID3 + 5];

  for(int ic = 0; ic < NC; ++ic){
    real_t* up = &ut[2 * ic];

    vt[2*ic] =
           MULT_GDXr(up[0], up[1], up[6], up[7], up[12], up[13],
                     vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
    vt[2*ic+1] =
           MULT_GDXi(up[0], up[1], up[6], up[7], up[12], up[13],
                     vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
    vt[2*ic + NVC] =
           MULT_GDXr(up[0], up[1], up[6], up[7], up[12], up[13],
                     vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
    vt[2*ic+1 + NVC] =
            MULT_GDXi(up[0], up[1], up[6], up[7], up[12], up[13],
                     vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
  }

}

//====================================================================
__device__ inline void mult_wilson_ym2(real_t* vt, real_t* wt)
{
  real_t* wt1 = &wt[0];
  real_t* wt2 = &wt[NVC];
  for(int ic = 0; ic < NC; ++ic){
    vt[ID1 + 2*ic  ] +=  wt1[2*ic];
    vt[ID1 + 2*ic+1] +=  wt1[2*ic+1];
    vt[ID2 + 2*ic  ] +=  wt2[2*ic];
    vt[ID2 + 2*ic+1] +=  wt2[2*ic+1];
    vt[ID3 + 2*ic  ] +=  wt2[2*ic];
    vt[ID3 + 2*ic+1] +=  wt2[2*ic+1];
    vt[ID4 + 2*ic  ] += -wt1[2*ic];
    vt[ID4 + 2*ic+1] += -wt1[2*ic+1];
  }

}

//====================================================================
__device__ inline void mult_wilson_zpb(real_t* vt, real_t* ut, real_t* wt)
{
  real_t vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5;
  real_t vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5;

  vt1_0 = wt[ID1 + 0] - wt[ID3 + 1];
  vt1_1 = wt[ID1 + 1] + wt[ID3 + 0];
  vt1_2 = wt[ID1 + 2] - wt[ID3 + 3];
  vt1_3 = wt[ID1 + 3] + wt[ID3 + 2];
  vt1_4 = wt[ID1 + 4] - wt[ID3 + 5];
  vt1_5 = wt[ID1 + 5] + wt[ID3 + 4];

  vt2_0 = wt[ID2 + 0] + wt[ID4 + 1];
  vt2_1 = wt[ID2 + 1] - wt[ID4 + 0];
  vt2_2 = wt[ID2 + 2] + wt[ID4 + 3];
  vt2_3 = wt[ID2 + 3] - wt[ID4 + 2];
  vt2_4 = wt[ID2 + 4] + wt[ID4 + 5];
  vt2_5 = wt[ID2 + 5] - wt[ID4 + 4];

  for(int ic = 0; ic < NC; ++ic){

    real_t wt1r, wt1i, wt2r, wt2i;
    real_t* up = &ut[NVC * ic];

    wt1r = MULT_GXr(up[0], up[1], up[2], up[3], up[4], up[5],
                    vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
    wt1i = MULT_GXi(up[0], up[1], up[2], up[3], up[4], up[5],
                    vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
    wt2r = MULT_GXr(up[0], up[1], up[2], up[3], up[4], up[5],
                    vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
    wt2i = MULT_GXi(up[0], up[1], up[2], up[3], up[4], up[5],
                    vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

    vt[ID1 + 2*ic  ] +=  wt1r;
    vt[ID1 + 2*ic+1] +=  wt1i;
    vt[ID2 + 2*ic  ] +=  wt2r;
    vt[ID2 + 2*ic+1] +=  wt2i;
    vt[ID3 + 2*ic  ] +=  wt1i;
    vt[ID3 + 2*ic+1] += -wt1r;
    vt[ID4 + 2*ic  ] += -wt2i;
    vt[ID4 + 2*ic+1] +=  wt2r;
  }

}

//====================================================================
__device__ inline void mult_wilson_zp1(real_t* vt, real_t* wt)
{
  vt[ 0] = wt[ID1 + 0] - wt[ID3 + 1];
  vt[ 1] = wt[ID1 + 1] + wt[ID3 + 0];
  vt[ 2] = wt[ID1 + 2] - wt[ID3 + 3];
  vt[ 3] = wt[ID1 + 3] + wt[ID3 + 2];
  vt[ 4] = wt[ID1 + 4] - wt[ID3 + 5];
  vt[ 5] = wt[ID1 + 5] + wt[ID3 + 4];

  vt[ 6] = wt[ID2 + 0] + wt[ID4 + 1];
  vt[ 7] = wt[ID2 + 1] - wt[ID4 + 0];
  vt[ 8] = wt[ID2 + 2] + wt[ID4 + 3];
  vt[ 9] = wt[ID2 + 3] - wt[ID4 + 2];
  vt[10] = wt[ID2 + 4] + wt[ID4 + 5];
  vt[11] = wt[ID2 + 5] - wt[ID4 + 4];
}

//====================================================================
__device__ inline void mult_wilson_zp2(real_t* vt, real_t* ut, real_t* wt)
{
  real_t* wt1 = &wt[0];
  real_t* wt2 = &wt[NVC];

  for(int ic = 0; ic < NC; ++ic){

    real_t wt1r, wt1i, wt2r, wt2i;
    real_t* up = &ut[NVC * ic];

    wt1r = MULT_GXr( up[0],  up[1],  up[2],  up[3],  up[4],  up[5],
                    wt1[0], wt1[1], wt1[2], wt1[3], wt1[4], wt1[5]);
    wt1i = MULT_GXi( up[0],  up[1],  up[2],  up[3],  up[4],  up[5],
                    wt1[0], wt1[1], wt1[2], wt1[3], wt1[4], wt1[5]);
    wt2r = MULT_GXr( up[0],  up[1],  up[2],  up[3],  up[4],  up[5],
                    wt2[0], wt2[1], wt2[2], wt2[3], wt2[4], wt2[5]);
    wt2i = MULT_GXi( up[0],  up[1],  up[2],  up[3],  up[4],  up[5],
                    wt2[0], wt2[1], wt2[2], wt2[3], wt2[4], wt2[5]);

    vt[ID1 + 2*ic  ] +=  wt1r;
    vt[ID1 + 2*ic+1] +=  wt1i;
    vt[ID2 + 2*ic  ] +=  wt2r;
    vt[ID2 + 2*ic+1] +=  wt2i;
    vt[ID3 + 2*ic  ] +=  wt1i;
    vt[ID3 + 2*ic+1] += -wt1r;
    vt[ID4 + 2*ic  ] += -wt2i;
    vt[ID4 + 2*ic+1] +=  wt2r;
  }

}

//====================================================================
__device__ inline void mult_wilson_zmb(real_t* vt, real_t* ut, real_t* wt)
{
  real_t vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5;
  real_t vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5;

  vt1_0 = wt[ID1 + 0] + wt[ID3 + 1];
  vt1_1 = wt[ID1 + 1] - wt[ID3 + 0];
  vt1_2 = wt[ID1 + 2] + wt[ID3 + 3];
  vt1_3 = wt[ID1 + 3] - wt[ID3 + 2];
  vt1_4 = wt[ID1 + 4] + wt[ID3 + 5];
  vt1_5 = wt[ID1 + 5] - wt[ID3 + 4];

  vt2_0 = wt[ID2 + 0] - wt[ID4 + 1];
  vt2_1 = wt[ID2 + 1] + wt[ID4 + 0];
  vt2_2 = wt[ID2 + 2] - wt[ID4 + 3];
  vt2_3 = wt[ID2 + 3] + wt[ID4 + 2];
  vt2_4 = wt[ID2 + 4] - wt[ID4 + 5];
  vt2_5 = wt[ID2 + 5] + wt[ID4 + 4];

  for(int ic = 0; ic < NC; ++ic){

    real_t wt1r, wt1i, wt2r, wt2i;
    real_t* up = &ut[2 * ic];

    wt1r = MULT_GDXr(up[0], up[1], up[6], up[7], up[12], up[13],
                     vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
    wt1i = MULT_GDXi(up[0], up[1], up[6], up[7], up[12], up[13],
                     vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
    wt2r = MULT_GDXr(up[0], up[1], up[6], up[7], up[12], up[13],
                     vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
    wt2i = MULT_GDXi(up[0], up[1], up[6], up[7], up[12], up[13],
                     vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

    vt[ID1 + 2*ic  ] +=  wt1r;
    vt[ID1 + 2*ic+1] +=  wt1i;
    vt[ID2 + 2*ic  ] +=  wt2r;
    vt[ID2 + 2*ic+1] +=  wt2i;
    vt[ID3 + 2*ic  ] += -wt1i;
    vt[ID3 + 2*ic+1] +=  wt1r;
    vt[ID4 + 2*ic  ] +=  wt2i;
    vt[ID4 + 2*ic+1] += -wt2r;
  }

}

//====================================================================
__device__ inline void mult_wilson_zm1(real_t* vt, real_t* ut, real_t* wt)
{
  real_t vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5;
  real_t vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5;

  vt1_0 = wt[ID1 + 0] + wt[ID3 + 1];
  vt1_1 = wt[ID1 + 1] - wt[ID3 + 0];
  vt1_2 = wt[ID1 + 2] + wt[ID3 + 3];
  vt1_3 = wt[ID1 + 3] - wt[ID3 + 2];
  vt1_4 = wt[ID1 + 4] + wt[ID3 + 5];
  vt1_5 = wt[ID1 + 5] - wt[ID3 + 4];

  vt2_0 = wt[ID2 + 0] - wt[ID4 + 1];
  vt2_1 = wt[ID2 + 1] + wt[ID4 + 0];
  vt2_2 = wt[ID2 + 2] - wt[ID4 + 3];
  vt2_3 = wt[ID2 + 3] + wt[ID4 + 2];
  vt2_4 = wt[ID2 + 4] - wt[ID4 + 5];
  vt2_5 = wt[ID2 + 5] + wt[ID4 + 4];

  for(int ic = 0; ic < NC; ++ic){
    real_t* up = &ut[2 * ic];

    vt[2*ic] =
           MULT_GDXr(up[0], up[1], up[6], up[7], up[12], up[13],
                     vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
    vt[2*ic+1] =
           MULT_GDXi(up[0], up[1], up[6], up[7], up[12], up[13],
                     vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
    vt[2*ic + NVC] =
           MULT_GDXr(up[0], up[1], up[6], up[7], up[12], up[13],
                     vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
    vt[2*ic+1 + NVC] =
            MULT_GDXi(up[0], up[1], up[6], up[7], up[12], up[13],
                     vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
  }

}

//====================================================================
__device__ inline void mult_wilson_zm2(real_t* vt, real_t* wt)
{
  real_t* wt1 = &wt[0];
  real_t* wt2 = &wt[NVC];
  for(int ic = 0; ic < NC; ++ic){
    vt[ID1 + 2*ic  ] +=  wt1[2*ic];
    vt[ID1 + 2*ic+1] +=  wt1[2*ic+1];
    vt[ID2 + 2*ic  ] +=  wt2[2*ic];
    vt[ID2 + 2*ic+1] +=  wt2[2*ic+1];
    vt[ID3 + 2*ic  ] += -wt1[2*ic+1];
    vt[ID3 + 2*ic+1] +=  wt1[2*ic];
    vt[ID4 + 2*ic  ] +=  wt2[2*ic+1];
    vt[ID4 + 2*ic+1] += -wt2[2*ic];
  }

}

//====================================================================
__device__ inline void mult_wilson_tpb_dirac(real_t* vt, real_t* ut, real_t* wt)
{
  real_t vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5;
  real_t vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5;

  vt1_0 = 2.0 * wt[ID3 + 0];
  vt1_1 = 2.0 * wt[ID3 + 1];
  vt1_2 = 2.0 * wt[ID3 + 2];
  vt1_3 = 2.0 * wt[ID3 + 3];
  vt1_4 = 2.0 * wt[ID3 + 4];
  vt1_5 = 2.0 * wt[ID3 + 5];

  vt2_0 = 2.0 * wt[ID4 + 0];
  vt2_1 = 2.0 * wt[ID4 + 1];
  vt2_2 = 2.0 * wt[ID4 + 2];
  vt2_3 = 2.0 * wt[ID4 + 3];
  vt2_4 = 2.0 * wt[ID4 + 4];
  vt2_5 = 2.0 * wt[ID4 + 5];

  for(int ic = 0; ic < NC; ++ic){

    real_t wt1r, wt1i, wt2r, wt2i;
    real_t* up = &ut[NVC * ic];

    wt1r = MULT_GXr(up[0], up[1], up[2], up[3], up[4], up[5],
                    vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
    wt1i = MULT_GXi(up[0], up[1], up[2], up[3], up[4], up[5],
                    vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
    wt2r = MULT_GXr(up[0], up[1], up[2], up[3], up[4], up[5],
                    vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
    wt2i = MULT_GXi(up[0], up[1], up[2], up[3], up[4], up[5],
                    vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

    vt[ID3 + 2*ic  ] += wt1r;
    vt[ID3 + 2*ic+1] += wt1i;
    vt[ID4 + 2*ic  ] += wt2r;
    vt[ID4 + 2*ic+1] += wt2i;
  }

}

//====================================================================
__device__ inline void mult_wilson_tp1_dirac(real_t* vt, real_t* wt)
{
  vt[ 0] = 2.0 * wt[ID3 + 0];
  vt[ 1] = 2.0 * wt[ID3 + 1];
  vt[ 2] = 2.0 * wt[ID3 + 2];
  vt[ 3] = 2.0 * wt[ID3 + 3];
  vt[ 4] = 2.0 * wt[ID3 + 4];
  vt[ 5] = 2.0 * wt[ID3 + 5];

  vt[ 6] = 2.0 * wt[ID4 + 0];
  vt[ 7] = 2.0 * wt[ID4 + 1];
  vt[ 8] = 2.0 * wt[ID4 + 2];
  vt[ 9] = 2.0 * wt[ID4 + 3];
  vt[10] = 2.0 * wt[ID4 + 4];
  vt[11] = 2.0 * wt[ID4 + 5];
}

//====================================================================
__device__ inline void mult_wilson_tp2_dirac(real_t* vt, real_t* ut, real_t* wt)
{
  real_t* wt1 = &wt[0];
  real_t* wt2 = &wt[NVC];

  for(int ic = 0; ic < NC; ++ic){

    real_t wt1r, wt1i, wt2r, wt2i;
    real_t* up = &ut[NVC * ic];

    wt1r = MULT_GXr( up[0],  up[1],  up[2],  up[3],  up[4],  up[5],
                    wt1[0], wt1[1], wt1[2], wt1[3], wt1[4], wt1[5]);
    wt1i = MULT_GXi( up[0],  up[1],  up[2],  up[3],  up[4],  up[5],
                    wt1[0], wt1[1], wt1[2], wt1[3], wt1[4], wt1[5]);
    wt2r = MULT_GXr( up[0],  up[1],  up[2],  up[3],  up[4],  up[5],
                    wt2[0], wt2[1], wt2[2], wt2[3], wt2[4], wt2[5]);
    wt2i = MULT_GXi( up[0],  up[1],  up[2],  up[3],  up[4],  up[5],
                    wt2[0], wt2[1], wt2[2], wt2[3], wt2[4], wt2[5]);

    vt[ID3 + 2*ic  ] += wt1r;
    vt[ID3 + 2*ic+1] += wt1i;
    vt[ID4 + 2*ic  ] += wt2r;
    vt[ID4 + 2*ic+1] += wt2i;
  }

}

//====================================================================
__device__ inline void mult_wilson_tmb_dirac(real_t* vt, real_t* ut, real_t* wt)
{
  real_t vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5;
  real_t vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5;

  vt1_0 = 2.0 * wt[ID1 + 0];
  vt1_1 = 2.0 * wt[ID1 + 1];
  vt1_2 = 2.0 * wt[ID1 + 2];
  vt1_3 = 2.0 * wt[ID1 + 3];
  vt1_4 = 2.0 * wt[ID1 + 4];
  vt1_5 = 2.0 * wt[ID1 + 5];

  vt2_0 = 2.0 * wt[ID2 + 0];
  vt2_1 = 2.0 * wt[ID2 + 1];
  vt2_2 = 2.0 * wt[ID2 + 2];
  vt2_3 = 2.0 * wt[ID2 + 3];
  vt2_4 = 2.0 * wt[ID2 + 4];
  vt2_5 = 2.0 * wt[ID2 + 5];

  for(int ic = 0; ic < NC; ++ic){

    real_t wt1r, wt1i, wt2r, wt2i;
    real_t* up = &ut[2 * ic];

    wt1r = MULT_GDXr(up[0], up[1], up[6], up[7], up[12], up[13],
                     vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
    wt1i = MULT_GDXi(up[0], up[1], up[6], up[7], up[12], up[13],
                     vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
    wt2r = MULT_GDXr(up[0], up[1], up[6], up[7], up[12], up[13],
                     vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
    wt2i = MULT_GDXi(up[0], up[1], up[6], up[7], up[12], up[13],
                     vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

    vt[ID1 + 2*ic  ] +=  wt1r;
    vt[ID1 + 2*ic+1] +=  wt1i;
    vt[ID2 + 2*ic  ] +=  wt2r;
    vt[ID2 + 2*ic+1] +=  wt2i;
  }

}

//====================================================================
__device__ inline void mult_wilson_tm1_dirac(real_t* vt, real_t* ut, real_t* wt)
{
  real_t vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5;
  real_t vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5;

  vt1_0 = 2.0 * wt[ID1 + 0];
  vt1_1 = 2.0 * wt[ID1 + 1];
  vt1_2 = 2.0 * wt[ID1 + 2];
  vt1_3 = 2.0 * wt[ID1 + 3];
  vt1_4 = 2.0 * wt[ID1 + 4];
  vt1_5 = 2.0 * wt[ID1 + 5];

  vt2_0 = 2.0 * wt[ID2 + 0];
  vt2_1 = 2.0 * wt[ID2 + 1];
  vt2_2 = 2.0 * wt[ID2 + 2];
  vt2_3 = 2.0 * wt[ID2 + 3];
  vt2_4 = 2.0 * wt[ID2 + 4];
  vt2_5 = 2.0 * wt[ID2 + 5];

  for(int ic = 0; ic < NC; ++ic){
    real_t* up = &ut[2 * ic];

    vt[2*ic] =
           MULT_GDXr(up[0], up[1], up[6], up[7], up[12], up[13],
                     vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
    vt[2*ic+1] =
           MULT_GDXi(up[0], up[1], up[6], up[7], up[12], up[13],
                     vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
    vt[2*ic + NVC] =
           MULT_GDXr(up[0], up[1], up[6], up[7], up[12], up[13],
                     vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
    vt[2*ic+1 + NVC] =
            MULT_GDXi(up[0], up[1], up[6], up[7], up[12], up[13],
                     vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
  }

}

//====================================================================
__device__ inline void mult_wilson_tm2_dirac(real_t* vt, real_t* wt)
{
  real_t* wt1 = &wt[0];
  real_t* wt2 = &wt[NVC];
  for(int ic = 0; ic < NC; ++ic){
    vt[ID1 + 2*ic  ] +=  wt1[2*ic];
    vt[ID1 + 2*ic+1] +=  wt1[2*ic+1];
    vt[ID2 + 2*ic  ] +=  wt2[2*ic];
    vt[ID2 + 2*ic+1] +=  wt2[2*ic+1];
  }

}

#endif
//============================================================END=====
