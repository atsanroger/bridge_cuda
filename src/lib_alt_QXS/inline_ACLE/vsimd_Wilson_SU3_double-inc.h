/*!
      @file    vsimd_Wilson_SU3_double-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2562 $
*/

#ifndef QXS_VSIMD_WILSON_SU3_INCLUDED
#define QXS_VSIMD_WILSON_SU3_INCLUDED

#define  ID1     0
#define  ID2     2
#define  ID3     4
#define  ID4     6

namespace {

  /*
void check_Nc()
  {
    vout.general(CommonParameters::Vlevel(),
                 "Fopr_Wilson_impl: implementation for SU(3).\n");
  }
  */

inline void load_u(svbool_t pg,
                   svreal_t& u0, svreal_t& u1,
                   svreal_t& u2, svreal_t& u3,
                   svreal_t& u4, svreal_t& u5,
                   real_t *u)
{
  load_vec(pg, u0, &u[VLEN * 0]);
  load_vec(pg, u1, &u[VLEN * 1]);
  load_vec(pg, u2, &u[VLEN * (0 + NVC)]);
  load_vec(pg, u3, &u[VLEN * (1 + NVC)]);
  load_vec(pg, u4, &u[VLEN * (0 + 2 * NVC)]);
  load_vec(pg, u5, &u[VLEN * (1 + 2 * NVC)]);
}

inline void load_udag(svbool_t pg,
                      svreal_t& u0, svreal_t& u1,
                      svreal_t& u2, svreal_t& u3,
                      svreal_t& u4, svreal_t& u5,
                      real_t *u)
{
  load_vec(pg, u0, &u[VLEN * 0]);
  load_vec(pg, u1, &u[VLEN * 1]);
  load_vec(pg, u2, &u[VLEN * 2]);
  load_vec(pg, u3, &u[VLEN * 3]);
  load_vec(pg, u4, &u[VLEN * 4]);
  load_vec(pg, u5, &u[VLEN * 5]);
}

/*
inline void load_udag(svbool_t pg,
                      svreal_t& u0, svreal_t& u1,
                      svreal_t& u2, svreal_t& u3,
                      svreal_t& u4, svreal_t& u5,
                      real_t *u, svint_t index)
{
  load_vec_gather(pg, u0, &u[VLEN * 0], index);
  load_vec_gather(pg, u1, &u[VLEN * 1], index);
  load_vec_gather(pg, u2, &u[VLEN * 2], index);
  load_vec_gather(pg, u3, &u[VLEN * 3], index);
  load_vec_gather(pg, u4, &u[VLEN * 4], index);
  load_vec_gather(pg, u5, &u[VLEN * 5], index);
}
*/

inline void mult_uv(svbool_t pg,
                    svreal_t& yr, svreal_t& yi,
                    svreal_t u0, svreal_t u1,
                    svreal_t u2, svreal_t u3,
                    svreal_t u4, svreal_t u5,
                    svreal_t v0, svreal_t v1,
                    svreal_t v2, svreal_t v3,
                    svreal_t v4, svreal_t v5)
{
  yr = svmul_m(pg,     u0, v0);
  yr = svmls_m(pg, yr, u1, v1);
  yr = svmla_m(pg, yr, u2, v2);
  yr = svmls_m(pg, yr, u3, v3);
  yr = svmla_m(pg, yr, u4, v4);
  yr = svmls_m(pg, yr, u5, v5);
  yi = svmul_m(pg,     u0, v1);
  yi = svmla_m(pg, yi, u1, v0);
  yi = svmla_m(pg, yi, u2, v3);
  yi = svmla_m(pg, yi, u3, v2);
  yi = svmla_m(pg, yi, u4, v5);
  yi = svmla_m(pg, yi, u5, v4);
}

inline void mult_udv(svbool_t pg,
                     svreal_t& yr, svreal_t& yi,
                     svreal_t u0, svreal_t u1,
                     svreal_t u2, svreal_t u3,
                     svreal_t u4, svreal_t u5,
                     svreal_t v0, svreal_t v1,
                     svreal_t v2, svreal_t v3,
                     svreal_t v4, svreal_t v5)
{
  yr = svmul_m(pg,     u0, v0);
  yr = svmla_m(pg, yr, u1, v1);
  yr = svmla_m(pg, yr, u2, v2);
  yr = svmla_m(pg, yr, u3, v3);
  yr = svmla_m(pg, yr, u4, v4);
  yr = svmla_m(pg, yr, u5, v5);
  yi = svmul_m(pg,     u0, v1);
  yi = svmls_m(pg, yi, u1, v0);
  yi = svmla_m(pg, yi, u2, v3);
  yi = svmls_m(pg, yi, u3, v2);
  yi = svmla_m(pg, yi, u4, v5);
  yi = svmls_m(pg, yi, u5, v4);
}

inline void mult_ctv(Vsimd_t *x, Vsimd_t *u, Vsimd_t *w, int Nc)
{
  for(int k = 0; k < VLEN; ++k){
    x[0].v[k] =
        u[ 0].v[k] * w[0].v[k]    - u[ 1].v[k] * w[1].v[k]
      + u[ 6].v[k] * w[2*ND].v[k] - u[ 7].v[k] * w[2*ND+1].v[k]
      + u[12].v[k] * w[4*ND].v[k] - u[13].v[k] * w[4*ND+1].v[k];
    x[1].v[k] =
        u[ 0].v[k] * w[1].v[k]      + u[ 1].v[k] * w[0].v[k]
      + u[ 6].v[k] * w[2*ND+1].v[k] + u[ 7].v[k] * w[2*ND].v[k]
      + u[12].v[k] * w[4*ND+1].v[k] + u[13].v[k] * w[4*ND].v[k];
  }
}

inline void mult_uv(Vsimd_t *x, Vsimd_t *u, Vsimd_t *w, int Nc)
{
  for(int k = 0; k < VLEN; ++k){
    x[0].v[k] =
        u[ 0].v[k] * w[0].v[k] - u[ 1].v[k] * w[1].v[k]
      + u[ 6].v[k] * w[2].v[k] - u[ 7].v[k] * w[3].v[k]
      + u[12].v[k] * w[4].v[k] - u[13].v[k] * w[5].v[k];
    x[1].v[k] =
        u[ 0].v[k] * w[1].v[k] + u[ 1].v[k] * w[0].v[k]
      + u[ 6].v[k] * w[3].v[k] + u[ 7].v[k] * w[2].v[k]
      + u[12].v[k] * w[5].v[k] + u[13].v[k] * w[4].v[k];
  }
}

inline void mult_udagv(Vsimd_t *x, Vsimd_t *u, Vsimd_t *w, int Nc)
{
  for(int k = 0; k < VLEN; ++k){
    x[0].v[k] =
        u[0].v[k] * w[0].v[k] + u[1].v[k] * w[1].v[k]
      + u[2].v[k] * w[2].v[k] + u[3].v[k] * w[3].v[k]
      + u[4].v[k] * w[4].v[k] + u[5].v[k] * w[5].v[k];
    x[1].v[k] =
        u[0].v[k] * w[1].v[k] - u[1].v[k] * w[0].v[k]
      + u[2].v[k] * w[3].v[k] - u[3].v[k] * w[2].v[k]
      + u[4].v[k] * w[5].v[k] - u[5].v[k] * w[4].v[k];
  }
}


template <typename REALTYPE>
inline void load_sp4(svbool_t pg, 
                     svreal_t &vt1r, svreal_t &vt1i,
                     svreal_t &vt2r, svreal_t &vt2i,
                     svreal_t &vt3r, svreal_t &vt3i,
                     svreal_t &vt4r, svreal_t &vt4i,
                     REALTYPE *v, int ic)
 {
   int ic2 = ND * 2 * ic;
   load_vec(pg, vt1r, &v[VLEN*(ic2 +   ID1)]);
   load_vec(pg, vt1i, &v[VLEN*(ic2 + 1+ID1)]);
   load_vec(pg, vt2r, &v[VLEN*(ic2 +   ID2)]);
   load_vec(pg, vt2i, &v[VLEN*(ic2 + 1+ID2)]);
   load_vec(pg, vt3r, &v[VLEN*(ic2 +   ID3)]);
   load_vec(pg, vt3i, &v[VLEN*(ic2 + 1+ID3)]);
   load_vec(pg, vt4r, &v[VLEN*(ic2 +   ID4)]);
   load_vec(pg, vt4i, &v[VLEN*(ic2 + 1+ID4)]);
 }

template <typename REALTYPE>
inline void set_sp4(svbool_t pg, REALTYPE *v,
                    svreal_t vt1r, svreal_t vt1i,
                    svreal_t vt2r, svreal_t vt2i,
                    svreal_t vt3r, svreal_t vt3i,
                    svreal_t vt4r, svreal_t vt4i,
                    int ic)
 {
   int ic2 = ND * 2 * ic;
   save_vec(pg, &v[VLEN*(ic2 +   ID1)], vt1r);
   save_vec(pg, &v[VLEN*(ic2 + 1+ID1)], vt1i);
   save_vec(pg, &v[VLEN*(ic2 +   ID2)], vt2r);
   save_vec(pg, &v[VLEN*(ic2 + 1+ID2)], vt2i);
   save_vec(pg, &v[VLEN*(ic2 +   ID3)], vt3r);
   save_vec(pg, &v[VLEN*(ic2 + 1+ID3)], vt3i);
   save_vec(pg, &v[VLEN*(ic2 +   ID4)], vt4r);
   save_vec(pg, &v[VLEN*(ic2 + 1+ID4)], vt4i);
 }

template <typename REALTYPE>
inline void set_sp2_xp(Vsimd_t *vt1, Vsimd_t *vt2, REALTYPE *w)
{
  for(int ic = 0; ic < NC; ++ic){
    int icr = ND * 2 * ic;
    int ici = ND * 2 * ic + 1;
    for(int k = 0; k < VLEN; ++k){
      vt1[2*ic  ].v[k] = w[k + VLEN*(icr + ID1)] + w[k + VLEN*(ici + ID4)];
      vt1[2*ic+1].v[k] = w[k + VLEN*(ici + ID1)] - w[k + VLEN*(icr + ID4)];
      vt2[2*ic  ].v[k] = w[k + VLEN*(icr + ID2)] + w[k + VLEN*(ici + ID3)];
      vt2[2*ic+1].v[k] = w[k + VLEN*(ici + ID2)] - w[k + VLEN*(icr + ID3)];
    }
  }
}

 inline void set_sp2_xp(svbool_t pg,
			svreal_t& vt1r, svreal_t& vt1i,
			svreal_t& vt2r, svreal_t& vt2i,
			real_t* v, int ic)
 {
   int icr = ND * 2 * ic;
   int ici = ND * 2 * ic + 1;
   svreal_t w1r, w1i, w2r, w2i, w3r, w3i, w4r, w4i;
   load_vec(pg, w1r, &v[VLEN * (icr + ID1)]);
   load_vec(pg, w1i, &v[VLEN * (ici + ID1)]);
   load_vec(pg, w2r, &v[VLEN * (icr + ID2)]);
   load_vec(pg, w2i, &v[VLEN * (ici + ID2)]);
   load_vec(pg, w3r, &v[VLEN * (icr + ID3)]);
   load_vec(pg, w3i, &v[VLEN * (ici + ID3)]);
   load_vec(pg, w4r, &v[VLEN * (icr + ID4)]);
   load_vec(pg, w4i, &v[VLEN * (ici + ID4)]);

   vt1r = svadd_m(pg, w1r, w4i);
   vt1i = svsub_m(pg, w1i, w4r);
   vt2r = svadd_m(pg, w2r, w3i);
   vt2i = svsub_m(pg, w2i, w3r);
 }

/*
 inline void set_sp2_xp(svbool_t pg,
			svreal_t& vt1r, svreal_t& vt1i,
			svreal_t& vt2r, svreal_t& vt2i,
			real_t* v, int ic, svint_t index)
 {
   int icr = ND * 2 * ic;
   int ici = ND * 2 * ic + 1;
   svreal_t w1r, w1i, w2r, w2i, w3r, w3i, w4r, w4i;
   load_vec_gather(pg, w1r, &v[VLEN * (icr + ID1)], index);
   load_vec_gather(pg, w1i, &v[VLEN * (ici + ID1)], index);
   load_vec_gather(pg, w2r, &v[VLEN * (icr + ID2)], index);
   load_vec_gather(pg, w2i, &v[VLEN * (ici + ID2)], index);
   load_vec_gather(pg, w3r, &v[VLEN * (icr + ID3)], index);
   load_vec_gather(pg, w3i, &v[VLEN * (ici + ID3)], index);
   load_vec_gather(pg, w4r, &v[VLEN * (icr + ID4)], index);
   load_vec_gather(pg, w4i, &v[VLEN * (ici + ID4)], index);

   vt1r = svadd_m(pg, w1r, w4i);
   vt1i = svsub_m(pg, w1i, w4r);
   vt2r = svadd_m(pg, w2r, w3i);
   vt2i = svsub_m(pg, w2i, w3r);
 }
*/

template <typename REALTYPE>
inline void set_sp2_xp1(REALTYPE *vt, REALTYPE *w)
{
  for(int ic = 0; ic < NC; ++ic){
    int icr = ND * 2 * ic;
    int ici = ND * 2 * ic + 1;
    for(int ky = 0; ky < VLENY; ++ky){
      vt[ky + VLENY * (2*ic)]         =  w[ky + VLENY*(icr + ID1)]
                                       + w[ky + VLENY*(ici + ID4)];
      vt[ky + VLENY * (2*ic+1)]       =  w[ky + VLENY*(ici + ID1)]
                                       - w[ky + VLENY*(icr + ID4)];
      vt[ky + VLENY * (2*ic   + NVC)] =  w[ky + VLENY*(icr + ID2)]
                                       + w[ky + VLENY*(ici + ID3)];
      vt[ky + VLENY * (2*ic+1 + NVC)] =  w[ky + VLENY*(ici + ID2)]
                                       - w[ky + VLENY*(icr + ID3)];
    }
  }
}

template <typename REALTYPE>
inline void set_sp2_xm(Vsimd_t *vt1, Vsimd_t *vt2, REALTYPE *w)
{
  for(int ic = 0; ic < NC; ++ic){
    int icr = ND * 2 * ic;
    int ici = ND * 2 * ic + 1;
    for(int k = 0; k < VLEN; ++k){
      vt1[2*ic  ].v[k] = w[k + VLEN*(icr + ID1)] - w[k + VLEN*(ici + ID4)];
      vt1[2*ic+1].v[k] = w[k + VLEN*(ici + ID1)] + w[k + VLEN*(icr + ID4)];
      vt2[2*ic  ].v[k] = w[k + VLEN*(icr + ID2)] - w[k + VLEN*(ici + ID3)];
      vt2[2*ic+1].v[k] = w[k + VLEN*(ici + ID2)] + w[k + VLEN*(icr + ID3)];
    }
  }
}

 inline void set_sp2_xm(svbool_t pg,
			svreal_t& vt1r, svreal_t& vt1i,
			svreal_t& vt2r, svreal_t& vt2i,
			real_t* v, int ic)
 {
   int icr = ND * 2 * ic;
   int ici = ND * 2 * ic + 1;
   svreal_t w1r, w1i, w2r, w2i, w3r, w3i, w4r, w4i;
   load_vec(pg, w1r, &v[VLEN * (icr + ID1)]);
   load_vec(pg, w1i, &v[VLEN * (ici + ID1)]);
   load_vec(pg, w2r, &v[VLEN * (icr + ID2)]);
   load_vec(pg, w2i, &v[VLEN * (ici + ID2)]);
   load_vec(pg, w3r, &v[VLEN * (icr + ID3)]);
   load_vec(pg, w3i, &v[VLEN * (ici + ID3)]);
   load_vec(pg, w4r, &v[VLEN * (icr + ID4)]);
   load_vec(pg, w4i, &v[VLEN * (ici + ID4)]);

   vt1r = svsub_m(pg, w1r, w4i);
   vt1i = svadd_m(pg, w1i, w4r);
   vt2r = svsub_m(pg, w2r, w3i);
   vt2i = svadd_m(pg, w2i, w3r);
 }

/*
 inline void set_sp2_xm(svbool_t pg,
			svreal_t& vt1r, svreal_t& vt1i,
			svreal_t& vt2r, svreal_t& vt2i,
			real_t* v, int ic, svint_t index)
 {
   int icr = ND * 2 * ic;
   int ici = ND * 2 * ic + 1;
   svreal_t w1r, w1i, w2r, w2i, w3r, w3i, w4r, w4i;
   load_vec_gather(pg, w1r, &v[VLEN * (icr + ID1)], index);
   load_vec_gather(pg, w1i, &v[VLEN * (ici + ID1)], index);
   load_vec_gather(pg, w2r, &v[VLEN * (icr + ID2)], index);
   load_vec_gather(pg, w2i, &v[VLEN * (ici + ID2)], index);
   load_vec_gather(pg, w3r, &v[VLEN * (icr + ID3)], index);
   load_vec_gather(pg, w3i, &v[VLEN * (ici + ID3)], index);
   load_vec_gather(pg, w4r, &v[VLEN * (icr + ID4)], index);
   load_vec_gather(pg, w4i, &v[VLEN * (ici + ID4)], index);

   vt1r = svsub_m(pg, w1r, w4i);
   vt1i = svadd_m(pg, w1i, w4r);
   vt2r = svsub_m(pg, w2r, w3i);
   vt2i = svadd_m(pg, w2i, w3r);
 }
*/

template <typename REALTYPE>
inline void set_sp2_yp(Vsimd_t *vt1, Vsimd_t *vt2, REALTYPE *w)
{
  for(int ic = 0; ic < NC; ++ic){
    int icr = ND * 2 * ic;
    int ici = ND * 2 * ic + 1;
    for(int k = 0; k < VLEN; ++k){
      vt1[2*ic  ].v[k] = w[k + VLEN*(icr + ID1)] - w[k + VLEN*(icr + ID4)];
      vt1[2*ic+1].v[k] = w[k + VLEN*(ici + ID1)] - w[k + VLEN*(ici + ID4)];
      vt2[2*ic  ].v[k] = w[k + VLEN*(icr + ID2)] + w[k + VLEN*(icr + ID3)];
      vt2[2*ic+1].v[k] = w[k + VLEN*(ici + ID2)] + w[k + VLEN*(ici + ID3)];
    }
  }
}

 inline void set_sp2_yp(svbool_t pg,
			svreal_t& vt1r, svreal_t& vt1i,
			svreal_t& vt2r, svreal_t& vt2i,
			real_t* v, int ic)
 {
   int icr = ND * 2 * ic;
   int ici = ND * 2 * ic + 1;
   svreal_t w1r, w1i, w2r, w2i, w3r, w3i, w4r, w4i;
   load_vec(pg, w1r, &v[VLEN * (icr + ID1)]);
   load_vec(pg, w1i, &v[VLEN * (ici + ID1)]);
   load_vec(pg, w2r, &v[VLEN * (icr + ID2)]);
   load_vec(pg, w2i, &v[VLEN * (ici + ID2)]);
   load_vec(pg, w3r, &v[VLEN * (icr + ID3)]);
   load_vec(pg, w3i, &v[VLEN * (ici + ID3)]);
   load_vec(pg, w4r, &v[VLEN * (icr + ID4)]);
   load_vec(pg, w4i, &v[VLEN * (ici + ID4)]);

   vt1r = svsub_m(pg, w1r, w4r);
   vt1i = svsub_m(pg, w1i, w4i);
   vt2r = svadd_m(pg, w2r, w3r);
   vt2i = svadd_m(pg, w2i, w3i);
 }

/*
 inline void set_sp2_yp(svbool_t pg,
			svreal_t& vt1r, svreal_t& vt1i,
			svreal_t& vt2r, svreal_t& vt2i,
			real_t* v, int ic, svint_t index)
 {
   int icr = ND * 2 * ic;
   int ici = ND * 2 * ic + 1;
   svreal_t w1r, w1i, w2r, w2i, w3r, w3i, w4r, w4i;
   load_vec_gather(pg, w1r, &v[VLEN * (icr + ID1)], index);
   load_vec_gather(pg, w1i, &v[VLEN * (ici + ID1)], index);
   load_vec_gather(pg, w2r, &v[VLEN * (icr + ID2)], index);
   load_vec_gather(pg, w2i, &v[VLEN * (ici + ID2)], index);
   load_vec_gather(pg, w3r, &v[VLEN * (icr + ID3)], index);
   load_vec_gather(pg, w3i, &v[VLEN * (ici + ID3)], index);
   load_vec_gather(pg, w4r, &v[VLEN * (icr + ID4)], index);
   load_vec_gather(pg, w4i, &v[VLEN * (ici + ID4)], index);

   vt1r = svsub_m(pg, w1r, w4r);
   vt1i = svsub_m(pg, w1i, w4i);
   vt2r = svadd_m(pg, w2r, w3r);
   vt2i = svadd_m(pg, w2i, w3i);
 }
*/

template <typename REALTYPE>
inline void set_sp2_yp1(REALTYPE *vt, REALTYPE *w)
{
  for(int ic = 0; ic < NC; ++ic){
    int icr = ND * 2 * ic;
    int ici = ND * 2 * ic + 1;
    for(int kx = 0; kx < VLENX; ++kx){
      vt[kx + VLENX * (2*ic)]         =  w[kx + VLENX*(icr + ID1)]
                                       - w[kx + VLENX*(icr + ID4)];
      vt[kx + VLENX * (2*ic+1)]       =  w[kx + VLENX*(ici + ID1)]
                                       - w[kx + VLENX*(ici + ID4)];
      vt[kx + VLENX * (2*ic   + NVC)] =  w[kx + VLENX*(icr + ID2)]
                                       + w[kx + VLENX*(icr + ID3)];
      vt[kx + VLENX * (2*ic+1 + NVC)] =  w[kx + VLENX*(ici + ID2)]
                                       + w[kx + VLENX*(ici + ID3)];
    }
  }
}

template <typename REALTYPE>
inline void set_sp2_ym(Vsimd_t *vt1, Vsimd_t *vt2, REALTYPE *w)
{
  for(int ic = 0; ic < NC; ++ic){
    int icr = ND * 2 * ic;
    int ici = ND * 2 * ic + 1;
    for(int k = 0; k < VLEN; ++k){
      vt1[2*ic  ].v[k] = w[k + VLEN*(icr + ID1)] + w[k + VLEN*(icr + ID4)];
      vt1[2*ic+1].v[k] = w[k + VLEN*(ici + ID1)] + w[k + VLEN*(ici + ID4)];
      vt2[2*ic  ].v[k] = w[k + VLEN*(icr + ID2)] - w[k + VLEN*(icr + ID3)];
      vt2[2*ic+1].v[k] = w[k + VLEN*(ici + ID2)] - w[k + VLEN*(ici + ID3)];
    }
  }
}

 inline void set_sp2_ym(svbool_t pg,
			svreal_t& vt1r, svreal_t& vt1i,
			svreal_t& vt2r, svreal_t& vt2i,
			real_t* v, int ic)
 {
   int icr = ND * 2 * ic;
   int ici = ND * 2 * ic + 1;
   svreal_t w1r, w1i, w2r, w2i, w3r, w3i, w4r, w4i;
   load_vec(pg, w1r, &v[VLEN * (icr + ID1)]);
   load_vec(pg, w1i, &v[VLEN * (ici + ID1)]);
   load_vec(pg, w2r, &v[VLEN * (icr + ID2)]);
   load_vec(pg, w2i, &v[VLEN * (ici + ID2)]);
   load_vec(pg, w3r, &v[VLEN * (icr + ID3)]);
   load_vec(pg, w3i, &v[VLEN * (ici + ID3)]);
   load_vec(pg, w4r, &v[VLEN * (icr + ID4)]);
   load_vec(pg, w4i, &v[VLEN * (ici + ID4)]);

   vt1r = svadd_m(pg, w1r, w4r);
   vt1i = svadd_m(pg, w1i, w4i);
   vt2r = svsub_m(pg, w2r, w3r);
   vt2i = svsub_m(pg, w2i, w3i);
 }

/*
 inline void set_sp2_ym(svbool_t pg,
			svreal_t& vt1r, svreal_t& vt1i,
			svreal_t& vt2r, svreal_t& vt2i,
			real_t* v, int ic, svint_t index)
 {
   int icr = ND * 2 * ic;
   int ici = ND * 2 * ic + 1;
   svreal_t w1r, w1i, w2r, w2i, w3r, w3i, w4r, w4i;
   load_vec_gather(pg, w1r, &v[VLEN * (icr + ID1)], index);
   load_vec_gather(pg, w1i, &v[VLEN * (ici + ID1)], index);
   load_vec_gather(pg, w2r, &v[VLEN * (icr + ID2)], index);
   load_vec_gather(pg, w2i, &v[VLEN * (ici + ID2)], index);
   load_vec_gather(pg, w3r, &v[VLEN * (icr + ID3)], index);
   load_vec_gather(pg, w3i, &v[VLEN * (ici + ID3)], index);
   load_vec_gather(pg, w4r, &v[VLEN * (icr + ID4)], index);
   load_vec_gather(pg, w4i, &v[VLEN * (ici + ID4)], index);

   vt1r = svadd_m(pg, w1r, w4r);
   vt1i = svadd_m(pg, w1i, w4i);
   vt2r = svsub_m(pg, w2r, w3r);
   vt2i = svsub_m(pg, w2i, w3i);
 }
*/

template <typename REALTYPE>
inline void set_sp2_zp(Vsimd_t *vt1, Vsimd_t *vt2, REALTYPE *w)
{
  for(int ic = 0; ic < NC; ++ic){
    int icr = ND * 2 * ic;
    int ici = ND * 2 * ic + 1;
    for(int k = 0; k < VLEN; ++k){
      vt1[2*ic  ].v[k] = w[k + VLEN*(icr + ID1)] + w[k + VLEN*(ici + ID3)];
      vt1[2*ic+1].v[k] = w[k + VLEN*(ici + ID1)] - w[k + VLEN*(icr + ID3)];
      vt2[2*ic  ].v[k] = w[k + VLEN*(icr + ID2)] - w[k + VLEN*(ici + ID4)];
      vt2[2*ic+1].v[k] = w[k + VLEN*(ici + ID2)] + w[k + VLEN*(icr + ID4)];
    }
  }
}

inline void set_sp2_zp(svbool_t pg,
                       svreal_t& vt1r, svreal_t& vt1i,
                       svreal_t& vt2r, svreal_t& vt2i,
                       real_t* v, int ic)
{
  int icr = ND * 2 * ic;
  int ici = ND * 2 * ic + 1;
  svreal_t w1r, w1i, w2r, w2i, w3r, w3i, w4r, w4i;
  load_vec(pg, w1r, &v[VLEN * (icr + ID1)]);
  load_vec(pg, w1i, &v[VLEN * (ici + ID1)]);
  load_vec(pg, w2r, &v[VLEN * (icr + ID2)]);
  load_vec(pg, w2i, &v[VLEN * (ici + ID2)]);
  load_vec(pg, w3r, &v[VLEN * (icr + ID3)]);
  load_vec(pg, w3i, &v[VLEN * (ici + ID3)]);
  load_vec(pg, w4r, &v[VLEN * (icr + ID4)]);
  load_vec(pg, w4i, &v[VLEN * (ici + ID4)]);

  vt1r = svadd_m(pg, w1r, w3i);
  vt1i = svsub_m(pg, w1i, w3r);
  vt2r = svsub_m(pg, w2r, w4i);
  vt2i = svadd_m(pg, w2i, w4r);
}

template <typename REALTYPE>
inline void set_sp2_zm(Vsimd_t *vt1, Vsimd_t *vt2, REALTYPE *w)
{
  for(int ic = 0; ic < NC; ++ic){
    int icr = ND * 2 * ic;
    int ici = ND * 2 * ic + 1;
    for(int k = 0; k < VLEN; ++k){
      vt1[2*ic  ].v[k] = w[k + VLEN*(icr + ID1)] - w[k + VLEN*(ici + ID3)];
      vt1[2*ic+1].v[k] = w[k + VLEN*(ici + ID1)] + w[k + VLEN*(icr + ID3)];
      vt2[2*ic  ].v[k] = w[k + VLEN*(icr + ID2)] + w[k + VLEN*(ici + ID4)];
      vt2[2*ic+1].v[k] = w[k + VLEN*(ici + ID2)] - w[k + VLEN*(icr + ID4)];
    }
  }
}

inline void set_sp2_zm(svbool_t pg,
                       svreal_t& vt1r, svreal_t& vt1i,
                       svreal_t& vt2r, svreal_t& vt2i,
                       real_t* v, int ic)
{
  int icr = ND * 2 * ic;
  int ici = ND * 2 * ic + 1;
  svreal_t w1r, w1i, w2r, w2i, w3r, w3i, w4r, w4i;
  load_vec(pg, w1r, &v[VLEN * (icr + ID1)]);
  load_vec(pg, w1i, &v[VLEN * (ici + ID1)]);
  load_vec(pg, w2r, &v[VLEN * (icr + ID2)]);
  load_vec(pg, w2i, &v[VLEN * (ici + ID2)]);
  load_vec(pg, w3r, &v[VLEN * (icr + ID3)]);
  load_vec(pg, w3i, &v[VLEN * (ici + ID3)]);
  load_vec(pg, w4r, &v[VLEN * (icr + ID4)]);
  load_vec(pg, w4i, &v[VLEN * (ici + ID4)]);

  vt1r = svsub_m(pg, w1r, w3i);
  vt1i = svadd_m(pg, w1i, w3r);
  vt2r = svadd_m(pg, w2r, w4i);
  vt2i = svsub_m(pg, w2i, w4r);
}

template <typename REALTYPE>
inline void set_sp2_tp_dirac(Vsimd_t *vt1, Vsimd_t *vt2, REALTYPE *w)
{
  for(int ic = 0; ic < NC; ++ic){
    int icr = ND * 2 * ic;
    int ici = ND * 2 * ic + 1;
    for(int k = 0; k < VLEN; ++k){
      vt1[2*ic  ].v[k] = 2.0 * w[k + VLEN*(icr + ID3)];
      vt1[2*ic+1].v[k] = 2.0 * w[k + VLEN*(ici + ID3)];
      vt2[2*ic  ].v[k] = 2.0 * w[k + VLEN*(icr + ID4)];
      vt2[2*ic+1].v[k] = 2.0 * w[k + VLEN*(ici + ID4)];
    }
  }
}

inline void set_sp2_tp_dirac(svbool_t pg,
                             svreal_t& vt1r, svreal_t& vt1i,
                             svreal_t& vt2r, svreal_t& vt2i,
                             real_t* v, int ic)
{
  int icr = ND * 2 * ic;
  int ici = ND * 2 * ic + 1;
  svreal_t w3r, w3i, w4r, w4i;
  load_vec(pg, w3r, &v[VLEN * (icr + ID3)]);
  load_vec(pg, w3i, &v[VLEN * (ici + ID3)]);
  load_vec(pg, w4r, &v[VLEN * (icr + ID4)]);
  load_vec(pg, w4i, &v[VLEN * (ici + ID4)]);

  vt1r = svmul_m(pg, w3r, real_t(2.0));
  vt1i = svmul_m(pg, w3i, real_t(2.0));
  vt2r = svmul_m(pg, w4r, real_t(2.0));
  vt2i = svmul_m(pg, w4i, real_t(2.0));
}

template <typename REALTYPE>
inline void set_sp2_tm_dirac(Vsimd_t *vt1, Vsimd_t *vt2, REALTYPE *w)
{
  for(int ic = 0; ic < NC; ++ic){
    int icr = ND * 2 * ic;
    int ici = ND * 2 * ic + 1;
    for(int k = 0; k < VLEN; ++k){
      vt1[2*ic  ].v[k] = 2.0 * w[k + VLEN*(icr + ID1)];
      vt1[2*ic+1].v[k] = 2.0 * w[k + VLEN*(ici + ID1)];
      vt2[2*ic  ].v[k] = 2.0 * w[k + VLEN*(icr + ID2)];
      vt2[2*ic+1].v[k] = 2.0 * w[k + VLEN*(ici + ID2)];
    }
  }
}

inline void set_sp2_tm_dirac(svbool_t pg,
                             svreal_t& vt1r, svreal_t& vt1i,
                             svreal_t& vt2r, svreal_t& vt2i,
                             real_t* v, int ic)
{
  int icr = ND * 2 * ic;
  int ici = ND * 2 * ic + 1;
  svreal_t w1r, w1i, w2r, w2i;
  load_vec(pg, w1r, &v[VLEN * (icr + ID1)]);
  load_vec(pg, w1i, &v[VLEN * (ici + ID1)]);
  load_vec(pg, w2r, &v[VLEN * (icr + ID2)]);
  load_vec(pg, w2i, &v[VLEN * (ici + ID2)]);

  vt1r = svmul_m(pg, w1r, real_t(2.0));
  vt1i = svmul_m(pg, w1i, real_t(2.0));
  vt2r = svmul_m(pg, w2r, real_t(2.0));
  vt2i = svmul_m(pg, w2i, real_t(2.0));
}

inline void set_sp4_xp(Vsimd_t *x, Vsimd_t *wt1, Vsimd_t *wt2)
{
    for(int k = 0; k < VLEN; ++k){
      x[  ID1].v[k] +=  wt1[0].v[k];
      x[1+ID1].v[k] +=  wt1[1].v[k];
      x[  ID2].v[k] +=  wt2[0].v[k];
      x[1+ID2].v[k] +=  wt2[1].v[k];
      x[  ID3].v[k] += -wt2[1].v[k];
      x[1+ID3].v[k] +=  wt2[0].v[k];
      x[  ID4].v[k] += -wt1[1].v[k];
      x[1+ID4].v[k] +=  wt1[0].v[k];
    }

 }

 inline void set_sp4_xp(svbool_t pg, Vsimd_t *v,
			svreal_t wt1r, svreal_t wt1i,
			svreal_t wt2r, svreal_t wt2i, int ic)
 {
   int ic2 = ND * 2 * ic;
   svreal_t vtr, vti;
   load_vec(pg, vtr, &v[ic2 +   ID1].v[0]);
   load_vec(pg, vti, &v[ic2 + 1+ID1].v[0]);
   vtr = svadd_m(pg, vtr, wt1r);
   vti = svadd_m(pg, vti, wt1i);
   save_vec(pg, &v[ic2 +   ID1].v[0], vtr);
   save_vec(pg, &v[ic2 + 1+ID1].v[0], vti);

   load_vec(pg, vtr, &v[ic2 +   ID2].v[0]);
   load_vec(pg, vti, &v[ic2 + 1+ID2].v[0]);
   vtr = svadd_m(pg, vtr, wt2r);
   vti = svadd_m(pg, vti, wt2i);
   save_vec(pg, &v[ic2 +   ID2].v[0], vtr);
   save_vec(pg, &v[ic2 + 1+ID2].v[0], vti);

   load_vec(pg, vtr, &v[ic2 +   ID3].v[0]);
   load_vec(pg, vti, &v[ic2 + 1+ID3].v[0]);
   vtr = svsub_m(pg, vtr, wt2i);
   vti = svadd_m(pg, vti, wt2r);
   save_vec(pg, &v[ic2 +   ID3].v[0], vtr);
   save_vec(pg, &v[ic2 + 1+ID3].v[0], vti);

   load_vec(pg, vtr, &v[ic2 +   ID4].v[0]);
   load_vec(pg, vti, &v[ic2 + 1+ID4].v[0]);
   vtr = svsub_m(pg, vtr, wt1i);
   vti = svadd_m(pg, vti, wt1r);
   save_vec(pg, &v[ic2 +   ID4].v[0], vtr);
   save_vec(pg, &v[ic2 + 1+ID4].v[0], vti);
 }

 inline void set_sp4_xp(svbool_t pg, real_t *v,
			svreal_t wt1r, svreal_t wt1i,
			svreal_t wt2r, svreal_t wt2i, int ic)
 {
   int ic2 = ND * 2 * ic;
   svreal_t vtr, vti;
   load_vec(pg, vtr, &v[VLEN*(ic2 +   ID1)]);
   load_vec(pg, vti, &v[VLEN*(ic2 + 1+ID1)]);
   vtr = svadd_m(pg, vtr, wt1r);
   vti = svadd_m(pg, vti, wt1i);
   save_vec(pg, &v[VLEN*(ic2 +   ID1)], vtr);
   save_vec(pg, &v[VLEN*(ic2 + 1+ID1)], vti);

   load_vec(pg, vtr, &v[VLEN*(ic2 +   ID2)]);
   load_vec(pg, vti, &v[VLEN*(ic2 + 1+ID2)]);
   vtr = svadd_m(pg, vtr, wt2r);
   vti = svadd_m(pg, vti, wt2i);
   save_vec(pg, &v[VLEN*(ic2 +   ID2)], vtr);
   save_vec(pg, &v[VLEN*(ic2 + 1+ID2)], vti);

   load_vec(pg, vtr, &v[VLEN*(ic2 +   ID3)]);
   load_vec(pg, vti, &v[VLEN*(ic2 + 1+ID3)]);
   vtr = svsub_m(pg, vtr, wt2i);
   vti = svadd_m(pg, vti, wt2r);
   save_vec(pg, &v[VLEN*(ic2 +   ID3)], vtr);
   save_vec(pg, &v[VLEN*(ic2 + 1+ID3)], vti);

   load_vec(pg, vtr, &v[VLEN*(ic2 +   ID4)]);
   load_vec(pg, vti, &v[VLEN*(ic2 + 1+ID4)]);
   vtr = svsub_m(pg, vtr, wt1i);
   vti = svadd_m(pg, vti, wt1r);
   save_vec(pg, &v[VLEN*(ic2 +   ID4)], vtr);
   save_vec(pg, &v[VLEN*(ic2 + 1+ID4)], vti);
 }

inline void set_sp4_xm(Vsimd_t *x, Vsimd_t *wt1, Vsimd_t *wt2)
{
    for(int k = 0; k < VLEN; ++k){
      x[  ID1].v[k] +=  wt1[0].v[k];
      x[1+ID1].v[k] +=  wt1[1].v[k];
      x[  ID2].v[k] +=  wt2[0].v[k];
      x[1+ID2].v[k] +=  wt2[1].v[k];
      x[  ID3].v[k] +=  wt2[1].v[k];
      x[1+ID3].v[k] += -wt2[0].v[k];
      x[  ID4].v[k] +=  wt1[1].v[k];
      x[1+ID4].v[k] += -wt1[0].v[k];
    }

 }

 inline void set_sp4_xm(svbool_t pg, Vsimd_t *v,
			svreal_t wt1r, svreal_t wt1i,
			svreal_t wt2r, svreal_t wt2i, int ic)
 {
   int ic2 = ND * 2 * ic;
   svreal_t vtr, vti;
   load_vec(pg, vtr, &v[ic2 +   ID1].v[0]);
   load_vec(pg, vti, &v[ic2 + 1+ID1].v[0]);
   vtr = svadd_m(pg, vtr, wt1r);
   vti = svadd_m(pg, vti, wt1i);
   save_vec(pg, &v[ic2 +   ID1].v[0], vtr);
   save_vec(pg, &v[ic2 + 1+ID1].v[0], vti);

   load_vec(pg, vtr, &v[ic2 +   ID2].v[0]);
   load_vec(pg, vti, &v[ic2 + 1+ID2].v[0]);
   vtr = svadd_m(pg, vtr, wt2r);
   vti = svadd_m(pg, vti, wt2i);
   save_vec(pg, &v[ic2 +   ID2].v[0], vtr);
   save_vec(pg, &v[ic2 + 1+ID2].v[0], vti);

   load_vec(pg, vtr, &v[ic2 +   ID3].v[0]);
   load_vec(pg, vti, &v[ic2 + 1+ID3].v[0]);
   vtr = svadd_m(pg, vtr, wt2i);
   vti = svsub_m(pg, vti, wt2r);
   save_vec(pg, &v[ic2 +   ID3].v[0], vtr);
   save_vec(pg, &v[ic2 + 1+ID3].v[0], vti);

   load_vec(pg, vtr, &v[ic2 +   ID4].v[0]);
   load_vec(pg, vti, &v[ic2 + 1+ID4].v[0]);
   vtr = svadd_m(pg, vtr, wt1i);
   vti = svsub_m(pg, vti, wt1r);
   save_vec(pg, &v[ic2 +   ID4].v[0], vtr);
   save_vec(pg, &v[ic2 + 1+ID4].v[0], vti);
 }

 inline void set_sp4_xm(svbool_t pg, real_t *v,
			svreal_t wt1r, svreal_t wt1i,
			svreal_t wt2r, svreal_t wt2i, int ic)
 {
   int ic2 = ND * 2 * ic;
   svreal_t vtr, vti;
   load_vec(pg, vtr, &v[VLEN*(ic2 +   ID1)]);
   load_vec(pg, vti, &v[VLEN*(ic2 + 1+ID1)]);
   vtr = svadd_m(pg, vtr, wt1r);
   vti = svadd_m(pg, vti, wt1i);
   save_vec(pg, &v[VLEN*(ic2 +   ID1)], vtr);
   save_vec(pg, &v[VLEN*(ic2 + 1+ID1)], vti);

   load_vec(pg, vtr, &v[VLEN*(ic2 +   ID2)]);
   load_vec(pg, vti, &v[VLEN*(ic2 + 1+ID2)]);
   vtr = svadd_m(pg, vtr, wt2r);
   vti = svadd_m(pg, vti, wt2i);
   save_vec(pg, &v[VLEN*(ic2 +   ID2)], vtr);
   save_vec(pg, &v[VLEN*(ic2 + 1+ID2)], vti);

   load_vec(pg, vtr, &v[VLEN*(ic2 +   ID3)]);
   load_vec(pg, vti, &v[VLEN*(ic2 + 1+ID3)]);
   vtr = svadd_m(pg, vtr, wt2i);
   vti = svsub_m(pg, vti, wt2r);
   save_vec(pg, &v[VLEN*(ic2 +   ID3)], vtr);
   save_vec(pg, &v[VLEN*(ic2 + 1+ID3)], vti);

   load_vec(pg, vtr, &v[VLEN*(ic2 +   ID4)]);
   load_vec(pg, vti, &v[VLEN*(ic2 + 1+ID4)]);
   vtr = svadd_m(pg, vtr, wt1i);
   vti = svsub_m(pg, vti, wt1r);
   save_vec(pg, &v[VLEN*(ic2 +   ID4)], vtr);
   save_vec(pg, &v[VLEN*(ic2 + 1+ID4)], vti);
 }


inline void set_sp4_yp(Vsimd_t *v, Vsimd_t *wt1, Vsimd_t *wt2)
{
  for(int k = 0; k < VLEN; ++k){
    v[  ID1].v[k] +=  wt1[0].v[k];
    v[1+ID1].v[k] +=  wt1[1].v[k];
    v[  ID2].v[k] +=  wt2[0].v[k];
    v[1+ID2].v[k] +=  wt2[1].v[k];
    v[  ID3].v[k] +=  wt2[0].v[k];
    v[1+ID3].v[k] +=  wt2[1].v[k];
    v[  ID4].v[k] += -wt1[0].v[k];
    v[1+ID4].v[k] += -wt1[1].v[k];
  }
}

 inline void set_sp4_yp(svbool_t pg, Vsimd_t *v,
			svreal_t wt1r, svreal_t wt1i,
			svreal_t wt2r, svreal_t wt2i, int ic)
 {
   int ic2 = ND * 2 * ic;
   svreal_t vtr, vti;
   load_vec(pg, vtr, &v[ic2 +   ID1].v[0]);
   load_vec(pg, vti, &v[ic2 + 1+ID1].v[0]);
   vtr = svadd_m(pg, vtr, wt1r);
   vti = svadd_m(pg, vti, wt1i);
   save_vec(pg, &v[ic2 +   ID1].v[0], vtr);
   save_vec(pg, &v[ic2 + 1+ID1].v[0], vti);

   load_vec(pg, vtr, &v[ic2 +   ID2].v[0]);
   load_vec(pg, vti, &v[ic2 + 1+ID2].v[0]);
   vtr = svadd_m(pg, vtr, wt2r);
   vti = svadd_m(pg, vti, wt2i);
   save_vec(pg, &v[ic2 +   ID2].v[0], vtr);
   save_vec(pg, &v[ic2 + 1+ID2].v[0], vti);

   load_vec(pg, vtr, &v[ic2 +   ID3].v[0]);
   load_vec(pg, vti, &v[ic2 + 1+ID3].v[0]);
   vtr = svadd_m(pg, vtr, wt2r);
   vti = svadd_m(pg, vti, wt2i);
   save_vec(pg, &v[ic2 +   ID3].v[0], vtr);
   save_vec(pg, &v[ic2 + 1+ID3].v[0], vti);

   load_vec(pg, vtr, &v[ic2 +   ID4].v[0]);
   load_vec(pg, vti, &v[ic2 + 1+ID4].v[0]);
   vtr = svsub_m(pg, vtr, wt1r);
   vti = svsub_m(pg, vti, wt1i);
   save_vec(pg, &v[ic2 +   ID4].v[0], vtr);
   save_vec(pg, &v[ic2 + 1+ID4].v[0], vti);
 }


 inline void set_sp4_yp(svbool_t pg, real_t *v,
			svreal_t wt1r, svreal_t wt1i,
			svreal_t wt2r, svreal_t wt2i, int ic)
 {
   int ic2 = ND * 2 * ic;
   svreal_t vtr, vti;

   load_vec(pg, vtr, &v[VLEN*(ic2 +   ID1)]);
   load_vec(pg, vti, &v[VLEN*(ic2 + 1+ID1)]);
   vtr = svadd_m(pg, vtr, wt1r);
   vti = svadd_m(pg, vti, wt1i);
   save_vec(pg, &v[VLEN*(ic2 +   ID1)], vtr);
   save_vec(pg, &v[VLEN*(ic2 + 1+ID1)], vti);

   load_vec(pg, vtr, &v[VLEN*(ic2 +   ID2)]);
   load_vec(pg, vti, &v[VLEN*(ic2 + 1+ID2)]);
   vtr = svadd_m(pg, vtr, wt2r);
   vti = svadd_m(pg, vti, wt2i);
   save_vec(pg, &v[VLEN*(ic2 +   ID2)], vtr);
   save_vec(pg, &v[VLEN*(ic2 + 1+ID2)], vti);

   load_vec(pg, vtr, &v[VLEN*(ic2 +   ID3)]);
   load_vec(pg, vti, &v[VLEN*(ic2 + 1+ID3)]);
   vtr = svadd_m(pg, vtr, wt2r);
   vti = svadd_m(pg, vti, wt2i);
   save_vec(pg, &v[VLEN*(ic2 +   ID3)], vtr);
   save_vec(pg, &v[VLEN*(ic2 + 1+ID3)], vti);

   load_vec(pg, vtr, &v[VLEN*(ic2 +   ID4)]);
   load_vec(pg, vti, &v[VLEN*(ic2 + 1+ID4)]);
   vtr = svsub_m(pg, vtr, wt1r);
   vti = svsub_m(pg, vti, wt1i);
   save_vec(pg, &v[VLEN*(ic2 +   ID4)], vtr);
   save_vec(pg, &v[VLEN*(ic2 + 1+ID4)], vti);
 }


inline void set_sp4_ym(Vsimd_t *v, Vsimd_t *wt1, Vsimd_t *wt2)
{
  for(int k = 0; k < VLEN; ++k){
    v[  ID1].v[k] +=  wt1[0].v[k];
    v[1+ID1].v[k] +=  wt1[1].v[k];
    v[  ID2].v[k] +=  wt2[0].v[k];
    v[1+ID2].v[k] +=  wt2[1].v[k];
    v[  ID3].v[k] += -wt2[0].v[k];
    v[1+ID3].v[k] += -wt2[1].v[k];
    v[  ID4].v[k] +=  wt1[0].v[k];
    v[1+ID4].v[k] +=  wt1[1].v[k];
  }
}

 inline void set_sp4_ym(svbool_t pg, Vsimd_t *v,
			svreal_t wt1r, svreal_t wt1i,
			svreal_t wt2r, svreal_t wt2i, int ic)
 {
   int ic2 = ND * 2 * ic;
   svreal_t vtr, vti;
   load_vec(pg, vtr, &v[ic2 +   ID1].v[0]);
   load_vec(pg, vti, &v[ic2 + 1+ID1].v[0]);
   vtr = svadd_m(pg, vtr, wt1r);
   vti = svadd_m(pg, vti, wt1i);
   save_vec(pg, &v[ic2 +   ID1].v[0], vtr);
   save_vec(pg, &v[ic2 + 1+ID1].v[0], vti);

   load_vec(pg, vtr, &v[ic2 +   ID2].v[0]);
   load_vec(pg, vti, &v[ic2 + 1+ID2].v[0]);
   vtr = svadd_m(pg, vtr, wt2r);
   vti = svadd_m(pg, vti, wt2i);
   save_vec(pg, &v[ic2 +   ID2].v[0], vtr);
   save_vec(pg, &v[ic2 + 1+ID2].v[0], vti);

   load_vec(pg, vtr, &v[ic2 +   ID3].v[0]);
   load_vec(pg, vti, &v[ic2 + 1+ID3].v[0]);
   vtr = svsub_m(pg, vtr, wt2r);
   vti = svsub_m(pg, vti, wt2i);
   save_vec(pg, &v[ic2 +   ID3].v[0], vtr);
   save_vec(pg, &v[ic2 + 1+ID3].v[0], vti);

   load_vec(pg, vtr, &v[ic2 +   ID4].v[0]);
   load_vec(pg, vti, &v[ic2 + 1+ID4].v[0]);
   vtr = svadd_m(pg, vtr, wt1r);
   vti = svadd_m(pg, vti, wt1i);
   save_vec(pg, &v[ic2 +   ID4].v[0], vtr);
   save_vec(pg, &v[ic2 + 1+ID4].v[0], vti);
 }

 inline void set_sp4_ym(svbool_t pg, real_t *v,
			svreal_t wt1r, svreal_t wt1i,
			svreal_t wt2r, svreal_t wt2i, int ic)
 {
   int ic2 = ND * 2 * ic;
   svreal_t vtr, vti;
   load_vec(pg, vtr, &v[VLEN*(ic2 +   ID1)]);
   load_vec(pg, vti, &v[VLEN*(ic2 + 1+ID1)]);
   vtr = svadd_m(pg, vtr, wt1r);
   vti = svadd_m(pg, vti, wt1i);
   save_vec(pg, &v[VLEN*(ic2 +   ID1)], vtr);
   save_vec(pg, &v[VLEN*(ic2 + 1+ID1)], vti);

   load_vec(pg, vtr, &v[VLEN*(ic2 +   ID2)]);
   load_vec(pg, vti, &v[VLEN*(ic2 + 1+ID2)]);
   vtr = svadd_m(pg, vtr, wt2r);
   vti = svadd_m(pg, vti, wt2i);
   save_vec(pg, &v[VLEN*(ic2 +   ID2)], vtr);
   save_vec(pg, &v[VLEN*(ic2 + 1+ID2)], vti);

   load_vec(pg, vtr, &v[VLEN*(ic2 +   ID3)]);
   load_vec(pg, vti, &v[VLEN*(ic2 + 1+ID3)]);
   vtr = svsub_m(pg, vtr, wt2r);
   vti = svsub_m(pg, vti, wt2i);
   save_vec(pg, &v[VLEN*(ic2 +   ID3)], vtr);
   save_vec(pg, &v[VLEN*(ic2 + 1+ID3)], vti);

   load_vec(pg, vtr, &v[VLEN*(ic2 +   ID4)]);
   load_vec(pg, vti, &v[VLEN*(ic2 + 1+ID4)]);
   vtr = svadd_m(pg, vtr, wt1r);
   vti = svadd_m(pg, vti, wt1i);
   save_vec(pg, &v[VLEN*(ic2 +   ID4)], vtr);
   save_vec(pg, &v[VLEN*(ic2 + 1+ID4)], vti);
 }

inline void set_sp4_zp(Vsimd_t *v, Vsimd_t *wt1, Vsimd_t *wt2)
{
  for(int k = 0; k < VLEN; ++k){
    v[  ID1].v[k] +=  wt1[0].v[k];
    v[1+ID1].v[k] +=  wt1[1].v[k];
    v[  ID2].v[k] +=  wt2[0].v[k];
    v[1+ID2].v[k] +=  wt2[1].v[k];
    v[  ID3].v[k] += -wt1[1].v[k];
    v[1+ID3].v[k] +=  wt1[0].v[k];
    v[  ID4].v[k] +=  wt2[1].v[k];
    v[1+ID4].v[k] += -wt2[0].v[k];
  }
}

inline void set_sp4_zp(svbool_t pg, Vsimd_t *v,
			svreal_t wt1r, svreal_t wt1i,
			svreal_t wt2r, svreal_t wt2i, int ic)
{
   int ic2 = ND * 2 * ic;
   svreal_t vtr, vti;
   load_vec(pg, vtr, &v[ic2 +   ID1].v[0]);
   load_vec(pg, vti, &v[ic2 + 1+ID1].v[0]);
   vtr = svadd_m(pg, vtr, wt1r);
   vti = svadd_m(pg, vti, wt1i);
   save_vec(pg, &v[ic2 +   ID1].v[0], vtr);
   save_vec(pg, &v[ic2 + 1+ID1].v[0], vti);

   load_vec(pg, vtr, &v[ic2 +   ID2].v[0]);
   load_vec(pg, vti, &v[ic2 + 1+ID2].v[0]);
   vtr = svadd_m(pg, vtr, wt2r);
   vti = svadd_m(pg, vti, wt2i);
   save_vec(pg, &v[ic2 +   ID2].v[0], vtr);
   save_vec(pg, &v[ic2 + 1+ID2].v[0], vti);

   load_vec(pg, vtr, &v[ic2 +   ID3].v[0]);
   load_vec(pg, vti, &v[ic2 + 1+ID3].v[0]);
   vtr = svsub_m(pg, vtr, wt1i);
   vti = svadd_m(pg, vti, wt1r);
   save_vec(pg, &v[ic2 +   ID3].v[0], vtr);
   save_vec(pg, &v[ic2 + 1+ID3].v[0], vti);

   load_vec(pg, vtr, &v[ic2 +   ID4].v[0]);
   load_vec(pg, vti, &v[ic2 + 1+ID4].v[0]);
   vtr = svadd_m(pg, vtr, wt2i);
   vti = svsub_m(pg, vti, wt2r);
   save_vec(pg, &v[ic2 +   ID4].v[0], vtr);
   save_vec(pg, &v[ic2 + 1+ID4].v[0], vti);
}

inline void set_sp4_zp(svbool_t pg, real_t *v,
			svreal_t wt1r, svreal_t wt1i,
			svreal_t wt2r, svreal_t wt2i, int ic)
{
   int ic2 = ND * 2 * ic;
   svreal_t vtr, vti;
   load_vec(pg, vtr, &v[VLEN*(ic2 +   ID1)]);
   load_vec(pg, vti, &v[VLEN*(ic2 + 1+ID1)]);
   vtr = svadd_m(pg, vtr, wt1r);
   vti = svadd_m(pg, vti, wt1i);
   save_vec(pg, &v[VLEN*(ic2 +   ID1)], vtr);
   save_vec(pg, &v[VLEN*(ic2 + 1+ID1)], vti);

   load_vec(pg, vtr, &v[VLEN*(ic2 +   ID2)]);
   load_vec(pg, vti, &v[VLEN*(ic2 + 1+ID2)]);
   vtr = svadd_m(pg, vtr, wt2r);
   vti = svadd_m(pg, vti, wt2i);
   save_vec(pg, &v[VLEN*(ic2 +   ID2)], vtr);
   save_vec(pg, &v[VLEN*(ic2 + 1+ID2)], vti);

   load_vec(pg, vtr, &v[VLEN*(ic2 +   ID3)]);
   load_vec(pg, vti, &v[VLEN*(ic2 + 1+ID3)]);
   vtr = svsub_m(pg, vtr, wt1i);
   vti = svadd_m(pg, vti, wt1r);
   save_vec(pg, &v[VLEN*(ic2 +   ID3)], vtr);
   save_vec(pg, &v[VLEN*(ic2 + 1+ID3)], vti);

   load_vec(pg, vtr, &v[VLEN*(ic2 +   ID4)]);
   load_vec(pg, vti, &v[VLEN*(ic2 + 1+ID4)]);
   vtr = svadd_m(pg, vtr, wt2i);
   vti = svsub_m(pg, vti, wt2r);
   save_vec(pg, &v[VLEN*(ic2 +   ID4)], vtr);
   save_vec(pg, &v[VLEN*(ic2 + 1+ID4)], vti);
}

inline void set_sp4_zm(Vsimd_t *v, Vsimd_t *wt1, Vsimd_t *wt2)
{
  for(int k = 0; k < VLEN; ++k){
    v[  ID1].v[k] +=  wt1[0].v[k];
    v[1+ID1].v[k] +=  wt1[1].v[k];
    v[  ID2].v[k] +=  wt2[0].v[k];
    v[1+ID2].v[k] +=  wt2[1].v[k];
    v[  ID3].v[k] +=  wt1[1].v[k];
    v[1+ID3].v[k] += -wt1[0].v[k];
    v[  ID4].v[k] += -wt2[1].v[k];
    v[1+ID4].v[k] +=  wt2[0].v[k];
  }
}

inline void set_sp4_zm(svbool_t pg, Vsimd_t *v,
			svreal_t wt1r, svreal_t wt1i,
			svreal_t wt2r, svreal_t wt2i, int ic)
{
   int ic2 = ND * 2 * ic;
   svreal_t vtr, vti;
   load_vec(pg, vtr, &v[ic2 +   ID1].v[0]);
   load_vec(pg, vti, &v[ic2 + 1+ID1].v[0]);
   vtr = svadd_m(pg, vtr, wt1r);
   vti = svadd_m(pg, vti, wt1i);
   save_vec(pg, &v[ic2 +   ID1].v[0], vtr);
   save_vec(pg, &v[ic2 + 1+ID1].v[0], vti);

   load_vec(pg, vtr, &v[ic2 +   ID2].v[0]);
   load_vec(pg, vti, &v[ic2 + 1+ID2].v[0]);
   vtr = svadd_m(pg, vtr, wt2r);
   vti = svadd_m(pg, vti, wt2i);
   save_vec(pg, &v[ic2 +   ID2].v[0], vtr);
   save_vec(pg, &v[ic2 + 1+ID2].v[0], vti);

   load_vec(pg, vtr, &v[ic2 +   ID3].v[0]);
   load_vec(pg, vti, &v[ic2 + 1+ID3].v[0]);
   vtr = svadd_m(pg, vtr, wt1i);
   vti = svsub_m(pg, vti, wt1r);
   save_vec(pg, &v[ic2 +   ID3].v[0], vtr);
   save_vec(pg, &v[ic2 + 1+ID3].v[0], vti);

   load_vec(pg, vtr, &v[ic2 +   ID4].v[0]);
   load_vec(pg, vti, &v[ic2 + 1+ID4].v[0]);
   vtr = svsub_m(pg, vtr, wt2i);
   vti = svadd_m(pg, vti, wt2r);
   save_vec(pg, &v[ic2 +   ID4].v[0], vtr);
   save_vec(pg, &v[ic2 + 1+ID4].v[0], vti);
}

inline void set_sp4_zm(svbool_t pg, real_t *v,
			svreal_t wt1r, svreal_t wt1i,
			svreal_t wt2r, svreal_t wt2i, int ic)
{
   int ic2 = ND * 2 * ic;
   svreal_t vtr, vti;
   load_vec(pg, vtr, &v[VLEN*(ic2 +   ID1)]);
   load_vec(pg, vti, &v[VLEN*(ic2 + 1+ID1)]);
   vtr = svadd_m(pg, vtr, wt1r);
   vti = svadd_m(pg, vti, wt1i);
   save_vec(pg, &v[VLEN*(ic2 +   ID1)], vtr);
   save_vec(pg, &v[VLEN*(ic2 + 1+ID1)], vti);

   load_vec(pg, vtr, &v[VLEN*(ic2 +   ID2)]);
   load_vec(pg, vti, &v[VLEN*(ic2 + 1+ID2)]);
   vtr = svadd_m(pg, vtr, wt2r);
   vti = svadd_m(pg, vti, wt2i);
   save_vec(pg, &v[VLEN*(ic2 +   ID2)], vtr);
   save_vec(pg, &v[VLEN*(ic2 + 1+ID2)], vti);

   load_vec(pg, vtr, &v[VLEN*(ic2 +   ID3)]);
   load_vec(pg, vti, &v[VLEN*(ic2 + 1+ID3)]);
   vtr = svadd_m(pg, vtr, wt1i);
   vti = svsub_m(pg, vti, wt1r);
   save_vec(pg, &v[VLEN*(ic2 +   ID3)], vtr);
   save_vec(pg, &v[VLEN*(ic2 + 1+ID3)], vti);

   load_vec(pg, vtr, &v[VLEN*(ic2 +   ID4)]);
   load_vec(pg, vti, &v[VLEN*(ic2 + 1+ID4)]);
   vtr = svsub_m(pg, vtr, wt2i);
   vti = svadd_m(pg, vti, wt2r);
   save_vec(pg, &v[VLEN*(ic2 +   ID4)], vtr);
   save_vec(pg, &v[VLEN*(ic2 + 1+ID4)], vti);
}

inline void set_sp4_tp_dirac(Vsimd_t *v, Vsimd_t *wt1, Vsimd_t *wt2)
{
  for(int k = 0; k < VLEN; ++k){
    v[  ID3].v[k] +=  wt1[0].v[k];
    v[1+ID3].v[k] +=  wt1[1].v[k];
    v[  ID4].v[k] +=  wt2[0].v[k];
    v[1+ID4].v[k] +=  wt2[1].v[k];
  }
}

inline void set_sp4_tp_dirac(svbool_t pg, real_t *v,
			      svreal_t wt1r, svreal_t wt1i,
			      svreal_t wt2r, svreal_t wt2i, int ic)
{
   int ic2 = ND * 2 * ic;
   svreal_t vtr, vti;
   load_vec(pg, vtr, &v[VLEN*(ic2 +   ID3)]);
   load_vec(pg, vti, &v[VLEN*(ic2 + 1+ID3)]);
   vtr = svadd_m(pg, vtr, wt1r);
   vti = svadd_m(pg, vti, wt1i);
   save_vec(pg, &v[VLEN*(ic2 +   ID3)], vtr);
   save_vec(pg, &v[VLEN*(ic2 + 1+ID3)], vti);

   load_vec(pg, vtr, &v[VLEN*(ic2 +   ID4)]);
   load_vec(pg, vti, &v[VLEN*(ic2 + 1+ID4)]);
   vtr = svadd_m(pg, vtr, wt2r);
   vti = svadd_m(pg, vti, wt2i);
   save_vec(pg, &v[VLEN*(ic2 +   ID4)], vtr);
   save_vec(pg, &v[VLEN*(ic2 + 1+ID4)], vti);
}

inline void set_sp4_tp_dirac(svbool_t pg, Vsimd_t *v,
			      svreal_t wt1r, svreal_t wt1i,
			      svreal_t wt2r, svreal_t wt2i, int ic)
{
   int ic2 = ND * 2 * ic;
   svreal_t vtr, vti;
   load_vec(pg, vtr, &v[ic2 +   ID3].v[0]);
   load_vec(pg, vti, &v[ic2 + 1+ID3].v[0]);
   vtr = svadd_m(pg, vtr, wt1r);
   vti = svadd_m(pg, vti, wt1i);
   save_vec(pg, &v[ic2 +   ID3].v[0], vtr);
   save_vec(pg, &v[ic2 + 1+ID3].v[0], vti);

   load_vec(pg, vtr, &v[ic2 +   ID4].v[0]);
   load_vec(pg, vti, &v[ic2 + 1+ID4].v[0]);
   vtr = svadd_m(pg, vtr, wt2r);
   vti = svadd_m(pg, vti, wt2i);
   save_vec(pg, &v[ic2 +   ID4].v[0], vtr);
   save_vec(pg, &v[ic2 + 1+ID4].v[0], vti);
}

inline void set_sp4_tm_dirac(Vsimd_t *v, Vsimd_t *wt1, Vsimd_t *wt2)
{
  for(int k = 0; k < VLEN; ++k){
    v[  ID1].v[k] +=  wt1[0].v[k];
    v[1+ID1].v[k] +=  wt1[1].v[k];
    v[  ID2].v[k] +=  wt2[0].v[k];
    v[1+ID2].v[k] +=  wt2[1].v[k];
  }
}

inline void set_sp4_tm_dirac(svbool_t pg, Vsimd_t *v,
			      svreal_t wt1r, svreal_t wt1i,
			      svreal_t wt2r, svreal_t wt2i, int ic)
{
   int ic2 = ND * 2 * ic;
   svreal_t vtr, vti;
   load_vec(pg, vtr, &v[ic2 +   ID1].v[0]);
   load_vec(pg, vti, &v[ic2 + 1+ID1].v[0]);
   vtr = svadd_m(pg, vtr, wt1r);
   vti = svadd_m(pg, vti, wt1i);
   save_vec(pg, &v[ic2 +   ID1].v[0], vtr);
   save_vec(pg, &v[ic2 + 1+ID1].v[0], vti);

   load_vec(pg, vtr, &v[ic2 +   ID2].v[0]);
   load_vec(pg, vti, &v[ic2 + 1+ID2].v[0]);
   vtr = svadd_m(pg, vtr, wt2r);
   vti = svadd_m(pg, vti, wt2i);
   save_vec(pg, &v[ic2 +   ID2].v[0], vtr);
   save_vec(pg, &v[ic2 + 1+ID2].v[0], vti);
}

inline void set_sp4_tm_dirac(svbool_t pg, real_t *v,
			      svreal_t wt1r, svreal_t wt1i,
			      svreal_t wt2r, svreal_t wt2i, int ic)
{
   int ic2 = ND * 2 * ic;
   svreal_t vtr, vti;
   load_vec(pg, vtr, &v[VLEN*(ic2 +   ID1)]);
   load_vec(pg, vti, &v[VLEN*(ic2 + 1+ID1)]);
   vtr = svadd_m(pg, vtr, wt1r);
   vti = svadd_m(pg, vti, wt1i);
   save_vec(pg, &v[VLEN*(ic2 +   ID1)], vtr);
   save_vec(pg, &v[VLEN*(ic2 + 1+ID1)], vti);

   load_vec(pg, vtr, &v[VLEN*(ic2 +   ID2)]);
   load_vec(pg, vti, &v[VLEN*(ic2 + 1+ID2)]);
   vtr = svadd_m(pg, vtr, wt2r);
   vti = svadd_m(pg, vti, wt2i);
   save_vec(pg, &v[VLEN*(ic2 +   ID2)], vtr);
   save_vec(pg, &v[VLEN*(ic2 + 1+ID2)], vti);
}


inline void mult_gm5_dirac_vec(Vsimd_t *v, Vsimd_t *w, int Nc)
{
  for(int ivc = 0; ivc < NVC; ++ivc){
    int ivc2 = (ivc % 2) + 2 * ND * (ivc/2);
    for(int k = 0; k < VLEN; ++k){
      v[ID1+ivc2].v[k] = - w[ID3+ivc2].v[k];
      v[ID2+ivc2].v[k] = - w[ID4+ivc2].v[k];
      v[ID3+ivc2].v[k] = - w[ID1+ivc2].v[k];
      v[ID4+ivc2].v[k] = - w[ID2+ivc2].v[k];
    }
  }
}

template <typename REALTYPE>
inline void mult_gm5_dirac_vec(svbool_t pg,
                               REALTYPE* __restrict v,
                               REALTYPE* __restrict w)
{
  //  svbool_t pg = set_predicate();
  svreal_t vt1r, vt1i, vt2r, vt2i;
  svreal_t vt3r, vt3i, vt4r, vt4i;

  load_vec( pg, vt3r, &w[VLEN*(ID1    )]);
  load_vec( pg, vt3i, &w[VLEN*(ID1 + 1)]);
  flip_sign(pg, vt3r);
  flip_sign(pg, vt3i);

  load_vec( pg, vt4r, &w[VLEN*(ID2    )]);
  load_vec( pg, vt4i, &w[VLEN*(ID2 + 1)]);
  flip_sign(pg, vt4r);
  flip_sign(pg, vt4i);

  load_vec( pg, vt1r, &w[VLEN*(ID3    )]);
  load_vec( pg, vt1i, &w[VLEN*(ID3 + 1)]);
  flip_sign(pg, vt1r);
  flip_sign(pg, vt1i);

  load_vec( pg, vt2r, &w[VLEN*(ID4    )]);
  load_vec( pg, vt2i, &w[VLEN*(ID4 + 1)]);
  flip_sign(pg, vt2r);
  flip_sign(pg, vt2i);

  save_vec(pg, &v[VLEN*(ID1    )], vt1r);
  save_vec(pg, &v[VLEN*(ID1 + 1)], vt1i);
  save_vec(pg, &v[VLEN*(ID2    )], vt2r);
  save_vec(pg, &v[VLEN*(ID2 + 1)], vt2i);
  save_vec(pg, &v[VLEN*(ID3    )], vt3r);
  save_vec(pg, &v[VLEN*(ID3 + 1)], vt3i);
  save_vec(pg, &v[VLEN*(ID4    )], vt4r);
  save_vec(pg, &v[VLEN*(ID4 + 1)], vt4i);

}


template <typename REALTYPE>
inline void load_mult_gm5_dirac_vec(Vsimd_t *v, REALTYPE *w, int Nc)
{
  for(int ivc = 0; ivc < NVC; ++ivc){
    int ivc2 = (ivc % 2) + 2 * ND * (ivc/2);
    for(int k = 0; k < VLEN; ++k){
      v[ID1+ivc2].v[k] = - w[k + VLEN*(ID3 + ivc2)];
      v[ID2+ivc2].v[k] = - w[k + VLEN*(ID4 + ivc2)];
      v[ID3+ivc2].v[k] = - w[k + VLEN*(ID1 + ivc2)];
      v[ID4+ivc2].v[k] = - w[k + VLEN*(ID2 + ivc2)];
    }
  }
}

template <typename REALTYPE>
inline void load_mult_gm5_dirac_vec(svbool_t pg,
                               svreal_t& vt1r, svreal_t& vt1i,
                               svreal_t& vt2r, svreal_t& vt2i,
                               svreal_t& vt3r, svreal_t& vt3i,
                               svreal_t& vt4r, svreal_t& vt4i,
                               REALTYPE *w)
{
  //  svbool_t pg = set_predicate();

  load_vec( pg, vt3r, &w[VLEN*(ID1    )]);
  load_vec( pg, vt3i, &w[VLEN*(ID1 + 1)]);
  flip_sign(pg, vt3r);
  flip_sign(pg, vt3i);

  load_vec( pg, vt4r, &w[VLEN*(ID2    )]);
  load_vec( pg, vt4i, &w[VLEN*(ID2 + 1)]);
  flip_sign(pg, vt4r);
  flip_sign(pg, vt4i);

  load_vec( pg, vt1r, &w[VLEN*(ID3    )]);
  load_vec( pg, vt1i, &w[VLEN*(ID3 + 1)]);
  flip_sign(pg, vt1r);
  flip_sign(pg, vt1i);

  load_vec( pg, vt2r, &w[VLEN*(ID4    )]);
  load_vec( pg, vt2i, &w[VLEN*(ID4 + 1)]);
  flip_sign(pg, vt2r);
  flip_sign(pg, vt2i);

}

template <typename REALTYPE>
inline void load_mult_gm5_dirac_vec(Vsimd_t *v,
                                    REALTYPE a, REALTYPE *w, int Nc)
{
  for(int ivc = 0; ivc < NVC; ++ivc){
    int ivc2 = (ivc % 2) + 2 * ND * (ivc/2);
    for(int k = 0; k < VLEN; ++k){
      v[ID1+ivc2].v[k] = - a * w[k + VLEN*(ID3 + ivc2)];
      v[ID2+ivc2].v[k] = - a * w[k + VLEN*(ID4 + ivc2)];
      v[ID3+ivc2].v[k] = - a * w[k + VLEN*(ID1 + ivc2)];
      v[ID4+ivc2].v[k] = - a * w[k + VLEN*(ID2 + ivc2)];
    }
  }
}

} // end of nameless namespace

#endif
//============================================================END=====
