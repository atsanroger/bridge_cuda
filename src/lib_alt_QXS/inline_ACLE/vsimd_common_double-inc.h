/*!
      @file    vsimd_common_double-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2562 $
*/

#ifndef QXS_VSIMD_COMMON_INCLUDED
#define QXS_VSIMD_COMMON_INCLUDED

namespace {

template <typename REALTYPE>
inline void load_vec(Vsimd_t *vt, REALTYPE *vp, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int k = 0; k < VLEN; ++k){
      vt[in].v[k] = vp[k + VLEN * in];
    }
  }
}

template <typename REALTYPE>
inline void load_vec1_x(REALTYPE *vt, REALTYPE *v, int kx, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int ky = 0; ky < VLENY; ++ky){
      vt[ky + VLENY * in] = v[kx + VLENX * ky + VLEN * in];
    }
  }
}

template <typename REALTYPE>
inline void load_vec1_x(Vsimd_t *vt, int kx, REALTYPE *v, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int ky = 0; ky < VLENY; ++ky){
      vt[in].v[kx + VLENX * ky + VLEN * in] = v[ky + VLENY * in];
    }
  }
}

template <typename REALTYPE>
inline void load_vec1_y(REALTYPE *vt, REALTYPE *v, int ky, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int kx = 0; kx < VLENX; ++kx){
      vt[kx + VLENX * in] = v[kx + VLENX * ky + VLEN * in];
    }
  }
}

template <typename REALTYPE>
inline void save_vec(REALTYPE *x, Vsimd_t *vt, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int k = 0; k < VLEN; ++k){
      x[k + VLEN * in] = vt[in].v[k];
    }
  }
}

/*
template <typename REALTYPE>
inline void save_vec(svbool_t pg, REALTYPE *x, const svreal_t& vt)
{
  svst1(pg, x, vt);
}
*/

inline void save_vec(svbool_t pg, Vsimd_t* x, const svreal_t& vt)
{
  //svst1(pg, &x->v[0], vt);
  save_vec(pg, &x->v[0], vt);
}

template <typename REALTYPE>
inline void save_vec_scatter(svbool_t pg, REALTYPE* vp,
                             svreal_t& vt, svint_t& index)
{
   svst1_scatter_index(pg, vp, index, vt);
 }


template <typename REALTYPE>
inline void save_vec1_x(REALTYPE *x, Vsimd_t *vt, int kx, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int ky = 0; ky < VLENY; ++ky){
      x[ky + VLENY * in] = vt[in].v[kx + VLENX * ky];
    }
  }
}

template <typename REALTYPE>
inline void save_vec1_y(REALTYPE *x, Vsimd_t *vt, int ky, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int kx = 0; kx < VLENX; ++kx){
      x[kx + VLENX * in] = vt[in].v[kx + VLENX * ky];
    }
  }
}

inline void clear_vec(Vsimd_t *vt, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int k = 0; k < VLEN; ++k){
      vt[in].v[k] = 0.0;
    }
  }
}

template <typename REALTYPE>
inline void set_vec(Vsimd_t *vt, REALTYPE a, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int k = 0; k < VLEN; ++k){
      vt[in].v[k] = a;
    }
  }
}

template <typename REALTYPE>
inline void add_vec(REALTYPE *x, Vsimd_t *vt, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int k = 0; k < VLEN; ++k){
      x[k + VLEN * in] += vt[in].v[k];
    }
  }
}

inline void add_vec(Vsimd_t *x, Vsimd_t *y, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int k = 0; k < VLEN; ++k){
      x[in].v[k] += y[in].v[k];
    }
  }
}

inline void add_vec(svbool_t pg, svreal_t& x, svreal_t& y)
{
  x = svadd_m(pg, x, y);
}

inline void add_vec(svbool_t pg, svreal_t& z, svreal_t& x, svreal_t& y)
{
  z = svadd_m(pg, x, y);
}

inline void sub_vec(svbool_t pg, svreal_t& x, svreal_t& y)
{
  x = svsub_m(pg, x, y);
}

inline void sub_vec(svbool_t pg, svreal_t& z, svreal_t& x, svreal_t& y)
{
  z = svsub_m(pg, x, y);
}

inline void copy_vec(Vsimd_t *x, Vsimd_t *y, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int k = 0; k < VLEN; ++k){
      x[in].v[k] = y[in].v[k];
    }
  }
}

template <typename REALTYPE>
inline void set_vec(Vsimd_t *x, REALTYPE a, Vsimd_t *y, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int k = 0; k < VLEN; ++k){
      x[in].v[k] = a * y[in].v[k];
    }
  }
}

inline void set_vec(svbool_t pg, svreal_t& x, real_t a, svreal_t y)
{
  x = svmul_m(pg, y, a);
}

inline void set_vec(svbool_t pg, svreal_t& x, svreal_t y)
{
  x = y;
}

template <typename REALTYPE>
inline void axpy_vec(Vsimd_t *y, REALTYPE a, Vsimd_t *x, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int k = 0; k < VLEN; ++k){
      y[in].v[k] += a * x[in].v[k];
    }
  }
}

inline void axpy_vec(svbool_t pg, svreal_t& y, real_t a, svreal_t x)
{
  y = svmla_m(pg, y, x, a);
}

inline void axpy_vec(svbool_t pg, svreal_t& y, svreal_t a, svreal_t x)
{
  y = svmla_m(pg, y, x, a);
}

inline void ymax_vec(svbool_t pg, svreal_t& y, svreal_t a, svreal_t x)
{
  y = svmls_m(pg, y, x, a);
}

template <typename REALTYPE>
inline void aypx_vec(REALTYPE a, Vsimd_t *x, Vsimd_t *y, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int k = 0; k < VLEN; ++k){
      x[in].v[k] = a * x[in].v[k] + y[in].v[k];
    }
  }
}

inline void aypx_vec(svbool_t pg, real_t a, svreal_t& y, svreal_t x)
{
  y = svmla_m(pg, x, y, a);
}

template <typename REALTYPE>
inline void scal_vec(Vsimd_t *x, REALTYPE a, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int k = 0; k < VLEN; ++k){
      x[in].v[k] *= a;
    }
  }
}

inline void mul_vec(svbool_t pg, svreal_t& x, svreal_t y, svreal_t w)
{
  x = svmul_m(pg, y, w);
}

inline void scal_vec(svbool_t pg, svreal_t& x, svreal_t a)
{
  x = svmul_m(pg, x, a);
}

inline void scal_vec(svbool_t pg, svreal_t& x, real_t a)
{
  x = svmul_m(pg, x, a);
}

inline void sqrt_vec(svbool_t pg, svreal_t& x, svreal_t& y)
{
  x = svsqrt_z(pg, y);
}

inline void div_vec(svbool_t pg, svreal_t& x, svreal_t& y, svreal_t& z)
{
  x = svdiv_m(pg, y, z);
}

template <typename REALTYPE>
  inline void dot_vec(REALTYPE& a, Vsimd_t *x, Vsimd_t *y, int Nin)
  {
    a = REALTYPE(0.0);
    for(int in = 0; in < Nin; ++in){
      for(int k = 0; k < VLEN; ++k){
	a += x[in].v[k] * y[in].v[k];
      }
    }
  }

template <typename REALTYPE>
  inline void norm2_vec(REALTYPE& a, Vsimd_t *x, int Nin)
  {
    a = REALTYPE(0.0);
    for(int in = 0; in < Nin; ++in){
      for(int k = 0; k < VLEN; ++k){
	a += x[in].v[k] * x[in].v[k];
      }
    }
  }

template <typename REALTYPE>
  inline void reduce_vec(REALTYPE& a, Vsimd_t *x, int Nin)
  {
    a = REALTYPE(0.0);
    for(int in = 0; in < Nin; ++in){
      for(int k = 0; k < VLEN; ++k){
	a += x[in].v[k];
      }
    }
  }

template <typename REALTYPE>
  inline void reduce_vec(svbool_t pg, REALTYPE& a, svreal_t& x)
  {
    a = svaddv(pg, x);
  }

inline void add_norm2_vec(svbool_t pg, svreal_t& y, svreal_t& x)
{
    y = svmla_m(pg, y, x, x);
}

inline void add_norm2_vec(Vsimd_t *y, Vsimd_t *x, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 0; k < VLEN; ++k){
        y[in].v[k] += x[in].v[k] * x[in].v[k];
      }
    }
  }

inline void add_dot_vec(svbool_t pg, svreal_t& y, svreal_t& x, svreal_t& w)
{
  y = svmla_m(pg, y, x, w);
}

inline void add_dot_vec(Vsimd_t *y, Vsimd_t *x, Vsimd_t *w, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 0; k < VLEN; ++k){
        y[in].v[k] += x[in].v[k] * w[in].v[k];
      }
    }
  }

inline void sub_dot_vec(svbool_t pg, svreal_t& y, svreal_t& x, svreal_t& w)
  {
    y = svmls_m(pg, y, x, w);
  }

inline void sub_dot_vec(Vsimd_t *y, Vsimd_t *x, Vsimd_t *w, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 0; k < VLEN; ++k){
        y[in].v[k] -= x[in].v[k] * w[in].v[k];
      }
    }
  }

 // setup of index for gather load / scatter store
 inline void set_index_xp(svint_t& svindex_xp)
 {
   int_t index[VLEN];
   for(int iy = 0; iy < VLENY; ++iy){
     index[VLENX-1 + VLENX*iy] = iy;
     for(int ix = 0; ix < VLENX-1; ++ix){
       index[ix + VLENX*iy] = 0;
     }
   }
   svbool_t pg = set_predicate();
   load_svint(pg, svindex_xp, index);
 }

 inline void set_index_xm(svint_t& svindex_xm)
 {
   int_t index[VLEN];
   for(int iy = 0; iy < VLENY; ++iy){
     index[VLENX*iy] = iy;
     for(int ix = 1; ix < VLENX; ++ix){
       index[ix + VLENX*iy] = 0;
     }
   }
   svbool_t pg = set_predicate();
   load_svint(pg, svindex_xm, index);
 }

 inline void set_index_xp_eo(svint_t& svindex_xp)
 {
   int_t index[VLEN];
   for(int iy = 0; iy < VLENY; ++iy){
     index[VLENX-1 + VLENX*iy] = iy/2;
     for(int ix = 0; ix < VLENX-1; ++ix){
       index[ix + VLENX*iy] = 0;
     }
   }
   svbool_t pg = set_predicate();
   load_svint(pg, svindex_xp, index);
 }

 inline void set_index_xm_eo(svint_t& svindex_xm)
 {
   int_t index[VLEN];
   for(int iy = 0; iy < VLENY; ++iy){
     index[VLENX*iy] = iy/2;
     for(int ix = 1; ix < VLENX; ++ix){
       index[ix + VLENX*iy] = 0;
     }
   }
   svbool_t pg = set_predicate();
   load_svint(pg, svindex_xm, index);
 }

 inline void set_index_xp_eo(svuint_t& svindex_xp)
 {
   uint_t index[VLEN];
   for(int iy = 0; iy < VLENY; ++iy){
     index[VLENX-1 + VLENX*iy] = iy/2;
     for(int ix = 0; ix < VLENX-1; ++ix){
       index[ix + VLENX*iy] = 0;
     }
   }
   svbool_t pg = set_predicate();
   load_svuint(pg, svindex_xp, index);
 }

 inline void set_index_xm_eo(svuint_t& svindex_xm)
 {
   uint_t index[VLEN];
   for(int iy = 0; iy < VLENY; ++iy){
     index[VLENX*iy] = iy/2;
     for(int ix = 1; ix < VLENX; ++ix){
       index[ix + VLENX*iy] = 0;
     }
   }
   svbool_t pg = set_predicate();
   load_svuint(pg, svindex_xm, index);
 }

 inline void set_index_yp(svint_t& svindex_yp)
 {
   int_t index[VLEN];
   for(int ix = 0; ix < VLENX; ++ix){
     for(int iy = 0; iy < VLENY-1; ++iy){
       index[ix + VLENX*iy] = 0;
     }
     index[ix + VLENX*(VLENY-1)] = ix;
   }
   svbool_t pg = set_predicate();
   load_svint(pg, svindex_yp, index);
 }

 inline void set_index_ym(svint_t& svindex_ym)
 {
   int_t index[VLEN];
   for(int ix = 0; ix < VLENX; ++ix){
     index[ix] = ix;
     for(int iy = 1; iy < VLENY; ++iy){
       index[ix + VLENX*iy] = 0;
     }
   }
   svbool_t pg = set_predicate();
   load_svint(pg, svindex_ym, index);
 }

template <typename REALTYPE>
  inline void shift_vec(svbool_t pg, svuint_t idx, 
                        svreal_t& v, 
                        const REALTYPE* __restrict xc,
                        const REALTYPE* __restrict xn)
  {
    svbool_t pg0 = set_predicate();
    svreal_t vc, vn;
    load_vec(pg0, vc, xc);
    load_vec(pg0, vn, xn);
    svreal_t vv = svsel(pg, vn, vc);
    v = svtbl(vv, idx);
  }

template <typename REALTYPE>
inline void shift_vec_xbw(svbool_t& pg1, svbool_t& pg2,
                          svreal_t& v, REALTYPE *wx, REALTYPE *wn)
{
  load_vec(pg1, v, &wx[1]);
  load_add(pg2, v, &wn[-VLENX+1]);
}

template <typename REALTYPE>
inline void shift_vec_xbw(svbool_t& pg1, svbool_t& pg2, svbool_t& pg3,
                          svreal_t& v, REALTYPE *wx, REALTYPE *wn)
{
  load_vec(pg3, v, &wx[0]);
  load_add(pg1, v, &wx[1]);
  load_add(pg2, v, &wn[-VLENX+1]);
}

template <typename REALTYPE>
inline void shift_vec_xfw(svbool_t& pg1, svbool_t& pg2,
                          svreal_t& v, REALTYPE *wx, REALTYPE *wn)
{
  load_vec(pg1, v, &wx[-1]);
  load_add(pg2, v, &wn[VLENX-1]);
}

template <typename REALTYPE>
inline void shift_vec_xfw(svbool_t& pg1, svbool_t& pg2, svbool_t& pg3,
                          svreal_t& v, REALTYPE *wx, REALTYPE *wn)
{
  load_vec(pg3, v, &wx[0]);
  load_add(pg1, v, &wx[-1]);
  load_add(pg2, v, &wn[VLENX-1]);
}

template <typename REALTYPE>
inline void shift_vec_xfw(svbool_t& pg1, svbool_t& pg2,
                          Vsimd_t *x, REALTYPE *wx, REALTYPE *wn,
                          int Nin)
{
  svbool_t pg = set_predicate(); 
  for(int in = 0; in < Nin; ++in){
    svreal_t vt;
    load_vec(pg1, vt, &wx[VLEN*in -1]);
    load_add(pg2, vt, &wn[VLEN*in + VLENX-1]);
    //svst1(pg, &x[in].v[0], vt);
    save_vec(pg, &x[in].v[0], vt);
  }
}

template <typename REALTYPE>
inline void shift_vec_ybw(svbool_t& pg1, svbool_t& pg2,
                          svreal_t& v, REALTYPE *wx, REALTYPE *wn)
{
  load_vec(pg1, v, &wx[VLENX]);
  load_add(pg2, v, &wn[-VLENX*(VLENY-1)]);
}

template <typename REALTYPE>
inline void shift_vec_yfw(svbool_t& pg1, svbool_t& pg2,
                          svreal_t& v, REALTYPE *wx, REALTYPE *wn)
{
  load_vec(pg1, v, &wx[-VLENX]);
  load_add(pg2, v, &wn[VLENX*(VLENY-1)]);
}

template <typename REALTYPE>
inline void shift_vec_ybw(svreal_t& v, REALTYPE *wx, REALTYPE *wn)
{
  svbool_t pg = set_predicate();
  svreal_t v1, v2;
  load_vec(pg, v1, &wx[0]);
  load_vec(pg, v2, &wn[0]);
  v = svext(v1, v2, VLENX);
}

template <typename REALTYPE>
inline void shift_vec_yfw(svreal_t& v, REALTYPE *wx, REALTYPE *wn)
{
  svbool_t pg = set_predicate();
  svreal_t v1, v2;
  load_vec(pg, v1, &wx[0]);
  load_vec(pg, v2, &wn[0]);
  v = svext(v2, v1, VLENX*(VLENY-1));
}

template <typename REALTYPE>
inline void shift_vec_yfw(svbool_t& pg1, svbool_t& pg2,
                          Vsimd_t *x, REALTYPE *wx, REALTYPE *wn,
                          int Nin)
{
  svbool_t pg = set_predicate(); 
  for(int in = 0; in < Nin; ++in){
    svreal_t vt;
    load_vec(pg1, vt, &wx[VLEN*in - VLENX]);
    load_add(pg2, vt, &wn[VLEN*in + VLENX*(VLENY-1)]);
    //svst1(pg, &x[in].v[0], vt);
    save_vec(pg, &x[in].v[0], vt);
  }
}

template <typename REALTYPE>
inline void shift_vec0_xbw(REALTYPE *v, REALTYPE *w, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int kx = 0; kx < VLENX-1; ++kx){
      for(int ky = 0; ky < VLENY; ++ky){
        v[kx + VLENX * ky + VLEN * in] = w[kx+1 + VLENX * ky + VLEN * in];
      }
    }
    int kx = VLENX - 1;
    for(int ky = 0; ky < VLENY; ++ky){
      v[kx + VLENX * ky + VLEN * in] = 0.0;
    }
  }
}

template <typename REALTYPE>
inline void shift_vec0_xfw(REALTYPE *v, REALTYPE *w, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int kx = 1; kx < VLENX; ++kx){
      for(int ky = 0; ky < VLENY; ++ky){
        v[kx + VLENX * ky + VLEN * in] = w[kx-1 + VLENX * ky + VLEN * in];
      }
    }
    for(int ky = 0; ky < VLENY; ++ky){
      v[0 + VLENX * ky + VLEN * in] = 0.0;
    }
  }
}

template <typename REALTYPE>
inline void shift_vec0_ybw(REALTYPE *v, REALTYPE *w, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int kx = 0; kx < VLENX; ++kx){
      for(int ky = 0; ky < VLENY-1; ++ky){
        v[kx + VLENX * ky + VLEN * in] = w[kx + VLENX * (ky+1) + VLEN * in];
      }
    }
    int ky = VLENY - 1;
    for(int kx = 0; kx < VLENX; ++kx){
      v[kx + VLENX * ky + VLEN * in] = 0.0;
    }
  }
}

template <typename REALTYPE>
inline void shift_vec0_yfw(REALTYPE *v, REALTYPE *w, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int kx = 0; kx < VLENX; ++kx){
      for(int ky = 1; ky < VLENY; ++ky){
        v[kx + VLENX * ky + VLEN * in] = w[kx + VLENX * (ky-1) + VLEN * in];
      }
    }
    int ky = 0;
    for(int kx = 0; kx < VLENX; ++kx){
      v[kx + VLENX * ky + VLEN * in] = 0.0;
    }
  }
}


template <typename REALTYPE>
inline void shift_vec1_xbw(Vsimd_t *x, REALTYPE *buf, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int kx = 0; kx < VLENX-1; ++kx){
      for(int ky = 0; ky < VLENY; ++ky){
        x[in].v[kx + VLENX * ky] = 0.0;
      }
    }
    int kx = VLENX - 1;
    for(int ky = 0; ky < VLENY; ++ky){
      x[in].v[kx + VLENX * ky] = buf[ky + VLENY * in];
    }
  }
}

template <typename REALTYPE>
inline void shift_vec1_xfw(Vsimd_t *x, REALTYPE *buf, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int kx = 1; kx < VLENX; ++kx){
      for(int ky = 0; ky < VLENY; ++ky){
        x[in].v[kx + VLENX * ky] = 0.0;
      }
    }
    for(int ky = 0; ky < VLENY; ++ky){
      x[in].v[0 + VLENX * ky] = buf[ky + VLENY * in];
    }
  }
}


template <typename REALTYPE>
inline void shift_vec1_ybw(Vsimd_t *v, REALTYPE *buf, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int kx = 0; kx < VLENX; ++kx){
      for(int ky = 0; ky < VLENY-1; ++ky){
        v[in].v[kx + VLENX * ky] = 0.0;
      }
    }
    int ky = VLENY - 1;
    for(int kx = 0; kx < VLENX; ++kx){
      v[in].v[kx + VLENX * ky] = buf[kx + VLENX * in];
    }
  }
}

template <typename REALTYPE>
inline void shift_vec1_yfw(Vsimd_t *v, REALTYPE *buf, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int kx = 0; kx < VLENX; ++kx){
      for(int ky = 1; ky < VLENY; ++ky){
        v[in].v[kx + VLENX * ky] = 0.0;
      }
    }
    int ky = 0;
    for(int kx = 0; kx < VLENX; ++kx){
      v[in].v[kx + VLENX * ky] = buf[kx + VLENX * in];
    }
  }
}

template <typename REALTYPE>
inline void shift_vec2_xbw(REALTYPE *v, REALTYPE *w, REALTYPE *y, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int kx = 0; kx < VLENX-1; ++kx){
      for(int ky = 0; ky < VLENY; ++ky){
        v[kx + VLENX * ky + VLEN * in] = w[kx+1 + VLENX * ky + VLEN * in];
      }
    }
    int kx = VLENX - 1;
    for(int ky = 0; ky < VLENY; ++ky){
      v[kx + VLENX * ky + VLEN * in] = y[0 + VLENX * ky + VLEN * in];
    }
  }
}

template <typename REALTYPE>
inline void shift_vec2_xfw(REALTYPE *v, REALTYPE *w, REALTYPE *y, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int kx = 1; kx < VLENX; ++kx){
      for(int ky = 0; ky < VLENY; ++ky){
        v[kx + VLENX * ky + VLEN * in] = w[kx-1 + VLENX * ky + VLEN * in];
      }
    }
    for(int ky = 0; ky < VLENY; ++ky){
      v[0 + VLENX * ky + VLEN * in] = y[VLENX-1 + VLENX * ky + VLEN * in];
    }
  }
}

template <typename REALTYPE>
inline void shift_vec2_xbw(Vsimd_t *v, REALTYPE *w, REALTYPE *y, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int kx = 0; kx < VLENX-1; ++kx){
      for(int ky = 0; ky < VLENY; ++ky){
        v[in].v[kx + VLENX * ky] = w[kx+1 + VLENX * ky + VLEN * in];
      }
    }
    int kx = VLENX - 1;
    for(int ky = 0; ky < VLENY; ++ky){
      v[in].v[kx + VLENX * ky] = y[0 + VLENX * ky + VLEN * in];
    }
  }
}

template <typename REALTYPE>
inline void shift_vec2_xbw_eo(Vsimd_t *v, REALTYPE *w, REALTYPE *y,
                              int ieo, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int ky = 0; ky < VLENY; ++ky){
      if((ky % 2) == ieo){
        for(int kx = 0; kx < VLENX; ++kx){
          v[in].v[kx + VLENX * ky] = w[kx + VLENX * ky + VLEN * in];
        }
      }else{
        for(int kx = 0; kx < VLENX-1; ++kx){
          v[in].v[kx + VLENX * ky] = w[kx+1 + VLENX * ky + VLEN * in];
        }
        v[in].v[VLENX-1 + VLENX * ky] = y[0 + VLENX * ky + VLEN * in];
      }
    }
  }
}

template <typename REALTYPE>
inline void shift_vec2_xfw_eo(Vsimd_t *v, REALTYPE *w, REALTYPE *y,
                              int ieo, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int ky = 0; ky < VLENY; ++ky){
      if((ky % 2) == ieo){
        for(int kx = 1; kx < VLENX; ++kx){
          v[in].v[kx + VLENX * ky] = w[kx-1 + VLENX * ky + VLEN * in];
        }
        v[in].v[0 + VLENX * ky] = y[VLENX-1 + VLENX * ky + VLEN * in];
      }else{
        for(int kx = 0; kx < VLENX; ++kx){
          v[in].v[kx + VLENX * ky] = w[kx + VLENX * ky + VLEN * in];
        }
      }
    }
  }
}

template <typename REALTYPE>
inline void shift_vec2_xfw(Vsimd_t *v, REALTYPE *w, REALTYPE *y, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int kx = 1; kx < VLENX; ++kx){
      for(int ky = 0; ky < VLENY; ++ky){
        v[in].v[kx + VLENX * ky] = w[kx-1 + VLENX * ky + VLEN * in];
      }
    }
    for(int ky = 0; ky < VLENY; ++ky){
      v[in].v[0 + VLENX * ky] = y[VLENX-1 + VLENX * ky + VLEN * in];
    }
  }
}

template <typename REALTYPE>
inline void shift_vec2_ybw(REALTYPE *v, REALTYPE *w, REALTYPE *y, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int kx = 0; kx < VLENX; ++kx){
      for(int ky = 0; ky < VLENY-1; ++ky){
        v[kx + VLENX * ky + VLEN * in] = w[kx + VLENX * (ky+1) + VLEN * in];
      }
    }
    int ky = VLENY - 1;
    for(int kx = 0; kx < VLENX; ++kx){
      v[kx + VLENX * ky + VLEN * in] = y[kx + VLENX * 0 + VLEN * in];
    }
  }
}

template <typename REALTYPE>
inline void shift_vec2_yfw(REALTYPE *v, REALTYPE *w, REALTYPE *y, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int kx = 0; kx < VLENX; ++kx){
      for(int ky = 1; ky < VLENY; ++ky){
        v[kx + VLENX * ky + VLEN * in] = w[kx + VLENX * (ky-1) + VLEN * in];
      }
    }
    int ky = 0;
    for(int kx = 0; kx < VLENX; ++kx){
      v[kx + VLENX * ky + VLEN * in] = y[kx + VLENX * (VLENY-1) + VLEN * in];
    }
  }
}

template <typename REALTYPE>
inline void shift_vec2_ybw(Vsimd_t *v, REALTYPE *w, REALTYPE *y, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int kx = 0; kx < VLENX; ++kx){
      for(int ky = 0; ky < VLENY-1; ++ky){
        v[in].v[kx + VLENX * ky] = w[kx + VLENX * (ky+1) + VLEN * in];
      }
    }
    int ky = VLENY - 1;
    for(int kx = 0; kx < VLENX; ++kx){
      v[in].v[kx + VLENX * ky] = y[kx + VLENX * 0 + VLEN * in];
    }
  }
}

template <typename REALTYPE>
inline void shift_vec2_yfw(Vsimd_t *v, REALTYPE *w, REALTYPE *y, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int kx = 0; kx < VLENX; ++kx){
      for(int ky = 1; ky < VLENY; ++ky){
        v[in].v[kx + VLENX * ky] = w[kx + VLENX * (ky-1) + VLEN * in];
      }
    }
    int ky = 0;
    for(int kx = 0; kx < VLENX; ++kx){
      v[in].v[kx + VLENX * ky] = y[kx + VLENX * (VLENY-1) + VLEN * in];
    }
  }
}


// the following definitions are to be discarded.

template <typename REALTYPE>
inline void shift_vec0_bw(REALTYPE *v, REALTYPE *w, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int k = 0; k < VLEN-1; ++k){
      v[k + VLEN * in] = w[k+1 + VLEN * in];
    }
    v[VLEN-1 + VLEN*in] = 0.0;
  }
}

template <typename REALTYPE>
inline void shift_vec0_fw(REALTYPE *v, REALTYPE *w, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int k = 1; k < VLEN; ++k){
      v[k + VLEN * in] = w[k-1 + VLEN * in];
    }
    v[0 + VLEN*in] = 0.0;
  }
}

template <typename REALTYPE>
inline void shift_vec1_bw(Vsimd_t *x, REALTYPE *buf, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int k = 0; k < VLEN-1; ++k){
      x[in].v[k] = 0.0;
    }
    x[in].v[VLEN-1] = buf[in];
  }
}

template <typename REALTYPE>
inline void shift_vec1_fw(Vsimd_t *x, REALTYPE *buf, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int k = 1; k < VLEN; ++k){
      x[in].v[k] = 0.0;
    }
    x[in].v[0] = buf[in];
  }
}

template <typename REALTYPE>
inline void shift_vec1_bw(REALTYPE *v, REALTYPE *w, REALTYPE *buf, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int k = 0; k < VLEN-1; ++k){
      v[k + VLEN * in] = w[k+1 + VLEN * in];
      v[k + VLEN * in] = w[k+1 + VLEN * in];
    }
    v[VLEN-1 + VLEN * in] = buf[in];
  }
}

template <typename REALTYPE>
inline void shift_vec1_fw(REALTYPE *v, REALTYPE *w, REALTYPE *buf, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int k = 1; k < VLEN; ++k){
      v[k + VLEN * in] = w[k-1 + VLEN * in];
    }
    v[0 + VLEN * in] = buf[in];
  }
}

template <typename REALTYPE>
inline void shift_vec2_bw(REALTYPE *v, REALTYPE *w, REALTYPE *y, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int k = 0; k < VLEN-1; ++k){
      v[k + VLEN * in] = w[k+1 + VLEN * in];
    }
    v[VLEN-1 + VLEN * in] = y[0 + VLEN * in];
  }
}

template <typename REALTYPE>
inline void shift_vec2_fw(REALTYPE *v, REALTYPE *w, REALTYPE *y, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int k = 1; k < VLEN; ++k){
      v[k + VLEN*in] = w[k-1 + VLEN * in];
    }
    v[0 + VLEN*in] = y[VLEN-1 + VLEN * in];
  }
}

template <typename REALTYPE>
inline void shift_vec2_bw(Vsimd_t *v, REALTYPE *w, REALTYPE *y, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int k = 0; k < VLEN-1; ++k){
      v[in].v[k] = w[k+1 + VLEN * in];
    }
    v[in].v[VLEN-1] = y[0 + VLEN * in];
  }
}

template <typename REALTYPE>
inline void shift_vec2_fw(Vsimd_t *v, REALTYPE *w, REALTYPE *y, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int k = 1; k < VLEN; ++k){
      v[in].v[k] = w[k-1 + VLEN * in];
    }
    v[in].v[0] = y[VLEN-1 + VLEN * in];
  }
}


template <typename REALTYPE>
inline void load_vec1(REALTYPE *vt, REALTYPE *v, int k, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    vt[in] = v[k + VLEN * in];
  }
}

template <typename REALTYPE>
inline void save_vec1(REALTYPE *x, Vsimd_t *vt, int k, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    x[in] = vt[in].v[k];
  }
}

inline svreal_t compact_vec(svbool_t pg, svreal_t& yt)
{
  return svcompact(pg, yt);
}

template <typename REALTYPE>
inline void load_add_gather(svbool_t pg2, svreal_t& vt, REALTYPE* v,
			    svuint_t& index, int skip)
 {
   svbool_t pg1 = set_predicate_whilelt(skip);
   svreal_t v1, v2;
   load_vec(pg1, v1, v);
   v2 = svtbl(v1, index);
   vt = svsel(pg2, v2, vt);
 }

} // end of nameless namespace

#endif
//============================================================END=====
