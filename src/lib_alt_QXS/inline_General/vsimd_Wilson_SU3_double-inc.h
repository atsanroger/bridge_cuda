/*!
      @file    vsimd_Wilson_SU3_double-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2422 $
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


  // to be discarded
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

inline void mult_uv(svbool_t pg,
                    svreal_t& yr, svreal_t& yi,
                    svreal_t u0, svreal_t u1,
                    svreal_t u2, svreal_t u3,
                    svreal_t u4, svreal_t u5,
                    svreal_t v0, svreal_t v1,
                    svreal_t v2, svreal_t v3,
                    svreal_t v4, svreal_t v5)
{
  for(int k = 0; k < VLEN; ++k){
    yr.v[k] =
        u0.v[k] * v0.v[k] - u1.v[k] * v1.v[k]
      + u2.v[k] * v2.v[k] - u3.v[k] * v3.v[k]
      + u4.v[k] * v4.v[k] - u5.v[k] * v5.v[k];
    yi.v[k] =
        u0.v[k] * v1.v[k] + u1.v[k] * v0.v[k]
      + u2.v[k] * v3.v[k] + u3.v[k] * v2.v[k]
      + u4.v[k] * v5.v[k] + u5.v[k] * v4.v[k];
  }
}

  // to be discarded
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

inline void mult_udv(svbool_t pg,
                     svreal_t& yr, svreal_t& yi,
                     svreal_t u0, svreal_t u1,
                     svreal_t u2, svreal_t u3,
                     svreal_t u4, svreal_t u5,
                     svreal_t v0, svreal_t v1,
                     svreal_t v2, svreal_t v3,
                     svreal_t v4, svreal_t v5)
{
  for(int k = 0; k < VLEN; ++k){
    yr.v[k] =
        u0.v[k] * v0.v[k] + u1.v[k] * v1.v[k]
      + u2.v[k] * v2.v[k] + u3.v[k] * v3.v[k]
      + u4.v[k] * v4.v[k] + u5.v[k] * v5.v[k];
    yi.v[k] =
        u0.v[k] * v1.v[k] - u1.v[k] * v0.v[k]
      + u2.v[k] * v3.v[k] - u3.v[k] * v2.v[k]
      + u4.v[k] * v5.v[k] - u5.v[k] * v4.v[k];
  }
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

 template <typename REALTYPE>
 inline void set_sp2_xp(svbool_t pg,
                        svreal_t& vt1r, svreal_t& vt1i,
                        svreal_t& vt2r, svreal_t& vt2i,
                        REALTYPE* w, int ic)
 {
   int icr = ND * 2 * ic;
   int ici = ND * 2 * ic + 1;
    for(int k = 0; k < VLEN; ++k){
      vt1r.v[k] = w[k + VLEN*(icr + ID1)] + w[k + VLEN*(ici + ID4)];
      vt1i.v[k] = w[k + VLEN*(ici + ID1)] - w[k + VLEN*(icr + ID4)];
      vt2r.v[k] = w[k + VLEN*(icr + ID2)] + w[k + VLEN*(ici + ID3)];
      vt2i.v[k] = w[k + VLEN*(ici + ID2)] - w[k + VLEN*(icr + ID3)];
    }
 }


  
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
 inline void set_sp2_xm(svbool_t pg,
                        svreal_t& vt1r, svreal_t& vt1i,
                        svreal_t& vt2r, svreal_t& vt2i,
                        REALTYPE* w, int ic)
 {
   int icr = ND * 2 * ic;
   int ici = ND * 2 * ic + 1;
   for(int k = 0; k < VLEN; ++k){
      vt1r.v[k] = w[k + VLEN*(icr + ID1)] - w[k + VLEN*(ici + ID4)];
      vt1i.v[k] = w[k + VLEN*(ici + ID1)] + w[k + VLEN*(icr + ID4)];
      vt2r.v[k] = w[k + VLEN*(icr + ID2)] - w[k + VLEN*(ici + ID3)];
      vt2i.v[k] = w[k + VLEN*(ici + ID2)] + w[k + VLEN*(icr + ID3)];
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

 template <typename REALTYPE>
 inline void set_sp2_yp(svbool_t pg,
                        svreal_t& vt1r, svreal_t& vt1i,
                        svreal_t& vt2r, svreal_t& vt2i,
                        REALTYPE* w, int ic)
 {
    int icr = ND * 2 * ic;
    int ici = ND * 2 * ic + 1;
    for(int k = 0; k < VLEN; ++k){
      vt1r.v[k] = w[k + VLEN*(icr + ID1)] - w[k + VLEN*(icr + ID4)];
      vt1i.v[k] = w[k + VLEN*(ici + ID1)] - w[k + VLEN*(ici + ID4)];
      vt2r.v[k] = w[k + VLEN*(icr + ID2)] + w[k + VLEN*(icr + ID3)];
      vt2i.v[k] = w[k + VLEN*(ici + ID2)] + w[k + VLEN*(ici + ID3)];
    }
 }


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
 inline void set_sp2_ym(svbool_t pg,
                        svreal_t& vt1r, svreal_t& vt1i,
                        svreal_t& vt2r, svreal_t& vt2i,
                        REALTYPE* w, int ic)
 {
    int icr = ND * 2 * ic;
    int ici = ND * 2 * ic + 1;
    for(int k = 0; k < VLEN; ++k){
      vt1r.v[k] = w[k + VLEN*(icr + ID1)] + w[k + VLEN*(icr + ID4)];
      vt1i.v[k] = w[k + VLEN*(ici + ID1)] + w[k + VLEN*(ici + ID4)];
      vt2r.v[k] = w[k + VLEN*(icr + ID2)] - w[k + VLEN*(icr + ID3)];
      vt2i.v[k] = w[k + VLEN*(ici + ID2)] - w[k + VLEN*(ici + ID3)];
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

 template <typename REALTYPE>
 inline void set_sp2_zp(svbool_t pg,
                        svreal_t& vt1r, svreal_t& vt1i,
                        svreal_t& vt2r, svreal_t& vt2i,
                        REALTYPE* w, int ic)
 {
    int icr = ND * 2 * ic;
    int ici = ND * 2 * ic + 1;
    for(int k = 0; k < VLEN; ++k){
      vt1r.v[k] = w[k + VLEN*(icr + ID1)] + w[k + VLEN*(ici + ID3)];
      vt1i.v[k] = w[k + VLEN*(ici + ID1)] - w[k + VLEN*(icr + ID3)];
      vt2r.v[k] = w[k + VLEN*(icr + ID2)] - w[k + VLEN*(ici + ID4)];
      vt2i.v[k] = w[k + VLEN*(ici + ID2)] + w[k + VLEN*(icr + ID4)];
    }
 }

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

template <typename REALTYPE>
 inline void set_sp2_zm(svbool_t pg,
                        svreal_t& vt1r, svreal_t& vt1i,
                        svreal_t& vt2r, svreal_t& vt2i,
                        REALTYPE* w, int ic)
 {
    int icr = ND * 2 * ic;
    int ici = ND * 2 * ic + 1;
    for(int k = 0; k < VLEN; ++k){
      vt1r.v[k] = w[k + VLEN*(icr + ID1)] - w[k + VLEN*(ici + ID3)];
      vt1i.v[k] = w[k + VLEN*(ici + ID1)] + w[k + VLEN*(icr + ID3)];
      vt2r.v[k] = w[k + VLEN*(icr + ID2)] + w[k + VLEN*(ici + ID4)];
      vt2i.v[k] = w[k + VLEN*(ici + ID2)] - w[k + VLEN*(icr + ID4)];
    }
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

template <typename REALTYPE>
 inline void set_sp2_tp_dirac(svbool_t pg,
                        svreal_t& vt1r, svreal_t& vt1i,
                        svreal_t& vt2r, svreal_t& vt2i,
                        REALTYPE* w, int ic)
 {
    int icr = ND * 2 * ic;
    int ici = ND * 2 * ic + 1;
    for(int k = 0; k < VLEN; ++k){
      vt1r.v[k] = 2.0 * w[k + VLEN*(icr + ID3)];
      vt1i.v[k] = 2.0 * w[k + VLEN*(ici + ID3)];
      vt2r.v[k] = 2.0 * w[k + VLEN*(icr + ID4)];
      vt2i.v[k] = 2.0 * w[k + VLEN*(ici + ID4)];
    }
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

template <typename REALTYPE>
 inline void set_sp2_tm_dirac(svbool_t pg,
                        svreal_t& vt1r, svreal_t& vt1i,
                        svreal_t& vt2r, svreal_t& vt2i,
                        REALTYPE* w, int ic)
 {
    int icr = ND * 2 * ic;
    int ici = ND * 2 * ic + 1;
    for(int k = 0; k < VLEN; ++k){
      vt1r.v[k] = 2.0 * w[k + VLEN*(icr + ID1)];
      vt1i.v[k] = 2.0 * w[k + VLEN*(ici + ID1)];
      vt2r.v[k] = 2.0 * w[k + VLEN*(icr + ID2)];
      vt2i.v[k] = 2.0 * w[k + VLEN*(ici + ID2)];
    }
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

 inline void set_sp4_xp(svbool_t pg, Vsimd_t *x,
                        svreal_t wt1r, svreal_t wt1i,
			svreal_t wt2r, svreal_t wt2i, int ic)
 {
    int ic2 = ND * 2 * ic;
    for(int k = 0; k < VLEN; ++k){
      x[ic2 +   ID1].v[k] +=  wt1r.v[k];
      x[ic2 + 1+ID1].v[k] +=  wt1i.v[k];
      x[ic2 +   ID2].v[k] +=  wt2r.v[k];
      x[ic2 + 1+ID2].v[k] +=  wt2i.v[k];
      x[ic2 +   ID3].v[k] += -wt2i.v[k];
      x[ic2 + 1+ID3].v[k] +=  wt2r.v[k];
      x[ic2 +   ID4].v[k] += -wt1i.v[k];
      x[ic2 + 1+ID4].v[k] +=  wt1r.v[k];
    }
 }

 template<typename REALTYPE>
 inline void set_sp4_xp(svbool_t pg, REALTYPE* x,
                        svreal_t wt1r, svreal_t wt1i,
			svreal_t wt2r, svreal_t wt2i, int ic)
 {
    int ic2 = ND * 2 * ic;
    for(int k = 0; k < VLEN; ++k){
      x[k + VLEN*(ic2 +   ID1)] +=  wt1r.v[k];
      x[k + VLEN*(ic2 + 1+ID1)] +=  wt1i.v[k];
      x[k + VLEN*(ic2 +   ID2)] +=  wt2r.v[k];
      x[k + VLEN*(ic2 + 1+ID2)] +=  wt2i.v[k];
      x[k + VLEN*(ic2 +   ID3)] += -wt2i.v[k];
      x[k + VLEN*(ic2 + 1+ID3)] +=  wt2r.v[k];
      x[k + VLEN*(ic2 +   ID4)] += -wt1i.v[k];
      x[k + VLEN*(ic2 + 1+ID4)] +=  wt1r.v[k];
    }
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

 inline void set_sp4_xm(svbool_t pg, Vsimd_t *x,
                        svreal_t wt1r, svreal_t wt1i,
			svreal_t wt2r, svreal_t wt2i, int ic)
 {
    int ic2 = ND * 2 * ic;
    for(int k = 0; k < VLEN; ++k){
      x[ic2 +   ID1].v[k] +=  wt1r.v[k];
      x[ic2 + 1+ID1].v[k] +=  wt1i.v[k];
      x[ic2 +   ID2].v[k] +=  wt2r.v[k];
      x[ic2 + 1+ID2].v[k] +=  wt2i.v[k];
      x[ic2 +   ID3].v[k] +=  wt2i.v[k];
      x[ic2 + 1+ID3].v[k] += -wt2r.v[k];
      x[ic2 +   ID4].v[k] +=  wt1i.v[k];
      x[ic2 + 1+ID4].v[k] += -wt1r.v[k];
    }
 }

 template <typename REALTYPE>
 inline void set_sp4_xm(svbool_t pg, REALTYPE* x,
                        svreal_t wt1r, svreal_t wt1i,
			svreal_t wt2r, svreal_t wt2i, int ic)
 {
    int ic2 = ND * 2 * ic;
    for(int k = 0; k < VLEN; ++k){
      x[k + VLEN*(ic2 +   ID1)] +=  wt1r.v[k];
      x[k + VLEN*(ic2 + 1+ID1)] +=  wt1i.v[k];
      x[k + VLEN*(ic2 +   ID2)] +=  wt2r.v[k];
      x[k + VLEN*(ic2 + 1+ID2)] +=  wt2i.v[k];
      x[k + VLEN*(ic2 +   ID3)] +=  wt2i.v[k];
      x[k + VLEN*(ic2 + 1+ID3)] += -wt2r.v[k];
      x[k + VLEN*(ic2 +   ID4)] +=  wt1i.v[k];
      x[k + VLEN*(ic2 + 1+ID4)] += -wt1r.v[k];
    }
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

 inline void set_sp4_yp(svbool_t pg, Vsimd_t *v,
                        svreal_t wt1r, svreal_t wt1i,
                        svreal_t wt2r, svreal_t wt2i, int ic)
 {
   int ic2 = ND * 2 * ic;
   for(int k = 0; k < VLEN; ++k){
    v[ic2 +   ID1].v[k] +=  wt1r.v[k];
    v[ic2 + 1+ID1].v[k] +=  wt1i.v[k];
    v[ic2 +   ID2].v[k] +=  wt2r.v[k];
    v[ic2 + 1+ID2].v[k] +=  wt2i.v[k];
    v[ic2 +   ID3].v[k] +=  wt2r.v[k];
    v[ic2 + 1+ID3].v[k] +=  wt2i.v[k];
    v[ic2 +   ID4].v[k] += -wt1r.v[k];
    v[ic2 + 1+ID4].v[k] += -wt1i.v[k];
  }
}

 inline void set_sp4_ym(svbool_t pg, Vsimd_t *v,
                        svreal_t wt1r, svreal_t wt1i,
                        svreal_t wt2r, svreal_t wt2i, int ic)
 {
   int ic2 = ND * 2 * ic;
   for(int k = 0; k < VLEN; ++k){
    v[ic2 +   ID1].v[k] +=  wt1r.v[k];
    v[ic2 + 1+ID1].v[k] +=  wt1i.v[k];
    v[ic2 +   ID2].v[k] +=  wt2r.v[k];
    v[ic2 + 1+ID2].v[k] +=  wt2i.v[k];
    v[ic2 +   ID3].v[k] += -wt2r.v[k];
    v[ic2 + 1+ID3].v[k] += -wt2i.v[k];
    v[ic2 +   ID4].v[k] +=  wt1r.v[k];
    v[ic2 + 1+ID4].v[k] +=  wt1i.v[k];
  }
}

 template <typename REALTYPE>
 inline void set_sp4_yp(svbool_t pg, REALTYPE* v,
                        svreal_t wt1r, svreal_t wt1i,
                        svreal_t wt2r, svreal_t wt2i, int ic)
 {
   int ic2 = ND * 2 * ic;
   for(int k = 0; k < VLEN; ++k){
     v[k + VLEN*(ic2 +   ID1)] +=  wt1r.v[k];
     v[k + VLEN*(ic2 + 1+ID1)] +=  wt1i.v[k];
     v[k + VLEN*(ic2 +   ID2)] +=  wt2r.v[k];
     v[k + VLEN*(ic2 + 1+ID2)] +=  wt2i.v[k];
     v[k + VLEN*(ic2 +   ID3)] +=  wt2r.v[k];
     v[k + VLEN*(ic2 + 1+ID3)] +=  wt2i.v[k];
     v[k + VLEN*(ic2 +   ID4)] += -wt1r.v[k];
     v[k + VLEN*(ic2 + 1+ID4)] += -wt1i.v[k];
  }
}

 template <typename REALTYPE>
 inline void set_sp4_ym(svbool_t pg, REALTYPE* v,
                        svreal_t wt1r, svreal_t wt1i,
                        svreal_t wt2r, svreal_t wt2i, int ic)
 {
   int ic2 = ND * 2 * ic;
   for(int k = 0; k < VLEN; ++k){
     v[k + VLEN*(ic2 +   ID1)] +=  wt1r.v[k];
     v[k + VLEN*(ic2 + 1+ID1)] +=  wt1i.v[k];
     v[k + VLEN*(ic2 +   ID2)] +=  wt2r.v[k];
     v[k + VLEN*(ic2 + 1+ID2)] +=  wt2i.v[k];
     v[k + VLEN*(ic2 +   ID3)] += -wt2r.v[k];
     v[k + VLEN*(ic2 + 1+ID3)] += -wt2i.v[k];
     v[k + VLEN*(ic2 +   ID4)] +=  wt1r.v[k];
     v[k + VLEN*(ic2 + 1+ID4)] +=  wt1i.v[k];
  }
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



 inline void set_sp4_zp(svbool_t pg, Vsimd_t *v,
                        svreal_t wt1r, svreal_t wt1i,
                        svreal_t wt2r, svreal_t wt2i, int ic)
 {
  int ic2 = ND * 2 * ic;
  for(int k = 0; k < VLEN; ++k){
    v[ic2 +   ID1].v[k] +=  wt1r.v[k];
    v[ic2 + 1+ID1].v[k] +=  wt1i.v[k];
    v[ic2 +   ID2].v[k] +=  wt2r.v[k];
    v[ic2 + 1+ID2].v[k] +=  wt2i.v[k];
    v[ic2 +   ID3].v[k] += -wt1i.v[k];
    v[ic2 + 1+ID3].v[k] +=  wt1r.v[k];
    v[ic2 +   ID4].v[k] +=  wt2i.v[k];
    v[ic2 + 1+ID4].v[k] += -wt2r.v[k];
  }
}

 inline void set_sp4_zm(svbool_t pg, Vsimd_t *v,
                        svreal_t wt1r, svreal_t wt1i,
                        svreal_t wt2r, svreal_t wt2i, int ic)
 {
  int ic2 = ND * 2 * ic;
  for(int k = 0; k < VLEN; ++k){
    v[ic2 +   ID1].v[k] +=  wt1r.v[k];
    v[ic2 + 1+ID1].v[k] +=  wt1i.v[k];
    v[ic2 +   ID2].v[k] +=  wt2r.v[k];
    v[ic2 + 1+ID2].v[k] +=  wt2i.v[k];
    v[ic2 +   ID3].v[k] +=  wt1i.v[k];
    v[ic2 + 1+ID3].v[k] += -wt1r.v[k];
    v[ic2 +   ID4].v[k] += -wt2i.v[k];
    v[ic2 + 1+ID4].v[k] +=  wt2r.v[k];
  }
}

 template <typename REALTYPE>
 inline void set_sp4_zp(svbool_t pg, REALTYPE* v,
                        svreal_t wt1r, svreal_t wt1i,
                        svreal_t wt2r, svreal_t wt2i, int ic)
 {
  int ic2 = ND * 2 * ic;
  for(int k = 0; k < VLEN; ++k){
    v[k + VLEN * (ic2 +   ID1)] +=  wt1r.v[k];
    v[k + VLEN * (ic2 + 1+ID1)] +=  wt1i.v[k];
    v[k + VLEN * (ic2 +   ID2)] +=  wt2r.v[k];
    v[k + VLEN * (ic2 + 1+ID2)] +=  wt2i.v[k];
    v[k + VLEN * (ic2 +   ID3)] += -wt1i.v[k];
    v[k + VLEN * (ic2 + 1+ID3)] +=  wt1r.v[k];
    v[k + VLEN * (ic2 +   ID4)] +=  wt2i.v[k];
    v[k + VLEN * (ic2 + 1+ID4)] += -wt2r.v[k];
  }
}

 template <typename REALTYPE>
 inline void set_sp4_zm(svbool_t pg, REALTYPE* v,
                        svreal_t wt1r, svreal_t wt1i,
                        svreal_t wt2r, svreal_t wt2i, int ic)
 {
  int ic2 = ND * 2 * ic;
  for(int k = 0; k < VLEN; ++k){
    v[k + VLEN * (ic2 +   ID1)] +=  wt1r.v[k];
    v[k + VLEN * (ic2 + 1+ID1)] +=  wt1i.v[k];
    v[k + VLEN * (ic2 +   ID2)] +=  wt2r.v[k];
    v[k + VLEN * (ic2 + 1+ID2)] +=  wt2i.v[k];
    v[k + VLEN * (ic2 +   ID3)] +=  wt1i.v[k];
    v[k + VLEN * (ic2 + 1+ID3)] += -wt1r.v[k];
    v[k + VLEN * (ic2 +   ID4)] += -wt2i.v[k];
    v[k + VLEN * (ic2 + 1+ID4)] +=  wt2r.v[k];
  }
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

inline void set_sp4_tm_dirac(Vsimd_t *v, Vsimd_t *wt1, Vsimd_t *wt2)
{
  for(int k = 0; k < VLEN; ++k){
    v[  ID1].v[k] +=  wt1[0].v[k];
    v[1+ID1].v[k] +=  wt1[1].v[k];
    v[  ID2].v[k] +=  wt2[0].v[k];
    v[1+ID2].v[k] +=  wt2[1].v[k];
  }
}

 inline void set_sp4_tp_dirac(svbool_t pg, Vsimd_t *v,
                              svreal_t wt1r, svreal_t wt1i,
                              svreal_t wt2r, svreal_t wt2i, int ic)
 {
  int ic2 = ND * 2 * ic;
  for(int k = 0; k < VLEN; ++k){
    v[ic2 +   ID3].v[k] +=  wt1r.v[k];
    v[ic2 + 1+ID3].v[k] +=  wt1i.v[k];
    v[ic2 +   ID4].v[k] +=  wt2r.v[k];
    v[ic2 + 1+ID4].v[k] +=  wt2i.v[k];
  }
 }

 inline void set_sp4_tm_dirac(svbool_t pg, Vsimd_t *v,
                              svreal_t wt1r, svreal_t wt1i,
                              svreal_t wt2r, svreal_t wt2i, int ic)
 {
  int ic2 = ND * 2 * ic;
  for(int k = 0; k < VLEN; ++k){
    v[ic2 +   ID1].v[k] +=  wt1r.v[k];
    v[ic2 + 1+ID1].v[k] +=  wt1i.v[k];
    v[ic2 +   ID2].v[k] +=  wt2r.v[k];
    v[ic2 + 1+ID2].v[k] +=  wt2i.v[k];
  }
 }

 template <typename REALTYPE>
 inline void set_sp4_tp_dirac(svbool_t pg, REALTYPE* v,
                              svreal_t wt1r, svreal_t wt1i,
                              svreal_t wt2r, svreal_t wt2i, int ic)
 {
  int ic2 = ND * 2 * ic;
  for(int k = 0; k < VLEN; ++k){
    v[k + VLEN * (ic2 +   ID3)] +=  wt1r.v[k];
    v[k + VLEN * (ic2 + 1+ID3)] +=  wt1i.v[k];
    v[k + VLEN * (ic2 +   ID4)] +=  wt2r.v[k];
    v[k + VLEN * (ic2 + 1+ID4)] +=  wt2i.v[k];
  }
 }

 template <typename REALTYPE>
 inline void set_sp4_tm_dirac(svbool_t pg, REALTYPE* v,
                              svreal_t wt1r, svreal_t wt1i,
                              svreal_t wt2r, svreal_t wt2i, int ic)
 {
  int ic2 = ND * 2 * ic;
  for(int k = 0; k < VLEN; ++k){
    v[k + VLEN * (ic2 +   ID1)] +=  wt1r.v[k];
    v[k + VLEN * (ic2 + 1+ID1)] +=  wt1i.v[k];
    v[k + VLEN * (ic2 +   ID2)] +=  wt2r.v[k];
    v[k + VLEN * (ic2 + 1+ID2)] +=  wt2i.v[k];
  }
 }

template <typename REALTYPE>
inline void mult_gm5_dirac_vec(svbool_t pg,
                               REALTYPE* v, REALTYPE* w)
{
  for(int k = 0; k < VLEN; ++k){
    v[k + VLEN*(ID3  )] = - w[k + VLEN*(ID1  )];
    v[k + VLEN*(ID3+1)] = - w[k + VLEN*(ID1+1)];
    v[k + VLEN*(ID4  )] = - w[k + VLEN*(ID2  )];
    v[k + VLEN*(ID4+1)] = - w[k + VLEN*(ID2+1)];
    v[k + VLEN*(ID1  )] = - w[k + VLEN*(ID3  )];
    v[k + VLEN*(ID1+1)] = - w[k + VLEN*(ID3+1)];
    v[k + VLEN*(ID2  )] = - w[k + VLEN*(ID4  )];
    v[k + VLEN*(ID2+1)] = - w[k + VLEN*(ID4+1)];
  }    

}



} // end of nameless namespace

#endif
//============================================================END=====
