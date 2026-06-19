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

  inline svbool_t set_predicate(){
    svbool_t pg;
    for(int i = 0; i < VLEN; ++i){
      pg.v[i] = true;
    }
    return pg;
  }

  inline svbool_t set_predicate_false(){
    svbool_t pg;
    for(int i = 0; i < VLEN; ++i){
      pg.v[i] = false;
    }
    return pg;
  }

  inline svbool_t set_predicate_not(const svbool_t pg1){
    svbool_t pg2;
    for(int i = 0; i < VLEN; ++i){
      if(pg1.v[i] == true){
	pg2.v[i] = false;
      }else{
	pg2.v[i] = true;
      }
    }
    return pg2;
  }

  inline svbool_t sveor_z(svbool_t pg, svbool_t pg1, svbool_t pg2){
    svbool_t pg3;
    for(int i = 0; i < VLEN; ++i){
      pg3.v[i] = pg1.v[i] ^ pg2.v[i];
    }
    return pg3;
  }

  inline void set_predicate_xp(svbool_t& pg1, svbool_t& pg2){
    for(int k = 0; k < VLEN; ++k){
      pg1.v[k] = true;
      pg2.v[k] = false;
    }
    for(int iy = 0; iy < VLENY; ++iy){
      pg1.v[VLENX-1 + VLENX * iy] = false;
      pg2.v[VLENX-1 + VLENX * iy] = true;
    }
  }

  inline void set_predicate_xm(svbool_t& pg1, svbool_t& pg2){
    for(int k = 0; k < VLEN; ++k){
      pg1.v[k] = true;
      pg2.v[k] = false;
    }
    for(int iy = 0; iy < VLENY; ++iy){
      pg1.v[VLENX * iy] = false;
      pg2.v[VLENX * iy] = true;
    }
  }

  inline void set_predicate_xp_eo(svbool_t& pg1, svbool_t& pg2,
                                  svbool_t& pg3, int ieo) {
    for(int k = 0; k < VLEN; ++k){
      pg1.v[k] = false;
      pg2.v[k] = false;
      pg3.v[k] = false;
    }
    for(int iy = 0; iy < VLENY; ++iy){
      if(iy % 2 != ieo){
        for(int ix = 0; ix < VLENX-1; ++ix){
          pg1.v[ix + VLENX * iy] = true;
        }
        pg2.v[VLENX-1 + VLENX * iy] = true;
      }else{
        for(int ix = 0; ix < VLENX; ++ix){
          pg3.v[ix + VLENX * iy] = true;
        }
      }
    }
  }

  inline void set_predicate_xm_eo(svbool_t& pg1, svbool_t& pg2,
                                  svbool_t& pg3, int ieo) {
    for(int k = 0; k < VLEN; ++k){
      pg1.v[k] = false;
      pg2.v[k] = false;
      pg3.v[k] = false;
    }
    for(int iy = 0; iy < VLENY; ++iy){
      if(iy % 2 != ieo){
        for(int ix = 0; ix < VLENX; ++ix){
          pg3.v[ix + VLENX * iy] = true;
        }
      }else{
        for(int ix = 1; ix < VLENX; ++ix){
          pg1.v[ix + VLENX * iy] = true;
        }
        pg2.v[0 + VLENX * iy] = true;
      }
    }
  }

  inline void set_predicate_yp(svbool_t& pg1, svbool_t& pg2){
    for(int k = 0; k < VLEN; ++k){
      pg1.v[k] = true;
      pg2.v[k] = false;
    }
    for(int ix = 0; ix < VLENX; ++ix){
      pg1.v[ix + VLENX * (VLENY-1)] = false;
      pg2.v[ix + VLENX * (VLENY-1)] = true;
    }
  }

  inline void set_predicate_ym(svbool_t& pg1, svbool_t& pg2){
    for(int k = 0; k < VLEN; ++k){
      pg1.v[k] = true;
      pg2.v[k] = false;
    }
    for(int ix = 0; ix < VLENX; ++ix){
      pg1.v[ix] = false;
      pg2.v[ix] = true;
    }
  }

inline void set_index_xp(svint_t& svindex_xp)
{
  for(int iy = 0; iy < VLENY; ++iy){
    svindex_xp.v[VLENX-1 + VLENX*iy] = iy;
    for(int ix = 0; ix < VLENX-1; ++ix){
      svindex_xp.v[ix + VLENX*iy] = 0;
    }
  }
}

inline void set_index_xm(svint_t& svindex_xm)
{
  for(int iy = 0; iy < VLENY; ++iy){
    svindex_xm.v[VLENX*iy] = iy;
    for(int ix = 1; ix < VLENX; ++ix){
      svindex_xm.v[ix + VLENX*iy] = 0;
    }
  }
 }

 inline void set_index_xp_eo(svint_t& svindex_xp)
 {
   for(int iy = 0; iy < VLENY; ++iy){
     svindex_xp.v[VLENX-1 + VLENX*iy] = iy/2;
     for(int ix = 0; ix < VLENX-1; ++ix){
       svindex_xp.v[ix + VLENX*iy] = 0;
     }
   }
 }

 inline void set_index_xm_eo(svint_t& svindex_xm)
 {
   for(int iy = 0; iy < VLENY; ++iy){
     svindex_xm.v[VLENX*iy] = iy/2;
     for(int ix = 1; ix < VLENX; ++ix){
       svindex_xm.v[ix + VLENX*iy] = 0;
     }
   }
 }

 inline void set_index_xp_eo(svuint_t& svindex_xp)
 {
   for(int iy = 0; iy < VLENY; ++iy){
     svindex_xp.v[VLENX-1 + VLENX*iy] = iy/2;
     for(int ix = 0; ix < VLENX-1; ++ix){
       svindex_xp.v[ix + VLENX*iy] = 0;
     }
   }
 }

 inline void set_index_xm_eo(svuint_t& svindex_xm)
 {
   for(int iy = 0; iy < VLENY; ++iy){
     svindex_xm.v[VLENX*iy] = iy/2;
     for(int ix = 1; ix < VLENX; ++ix){
       svindex_xm.v[ix + VLENX*iy] = 0;
     }
   }
 }

 inline void set_index_yp(svint_t& svindex_yp)
 {
   for(int ix = 0; ix < VLENX; ++ix){
     for(int iy = 0; iy < VLENY-1; ++iy){
       svindex_yp.v[ix + VLENX*iy] = 0;
     }
     svindex_yp.v[ix + VLENX*(VLENY-1)] = ix;
   }
 }

 inline void set_index_ym(svint_t& svindex_ym)
 {
   for(int ix = 0; ix < VLENX; ++ix){
     svindex_ym.v[ix] = ix;
     for(int iy = 1; iy < VLENY; ++iy){
       svindex_ym.v[ix + VLENX*iy] = 0;
     }
   }
 }

// to be discarded
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
inline void load_vec(svbool_t pg, Vsimd_t& vt, REALTYPE *vp)
{
  for(int k = 0; k < VLEN; ++k){
    if(pg.v[k]){
      vt.v[k] = vp[k];
    }else{
      vt.v[k] = 0.0;
    }
  }
}

template <typename REALTYPE>
inline void load_vec_gather(svbool_t pg, Vsimd_t& v, REALTYPE* vp,
                              svint_t index)
{
  for(int k = 0; k < VLEN; ++k){
    if(pg.v[k]) v.v[k] = vp[index.v[k]];
  }
}

template <typename REALTYPE>
inline void save_vec_scatter(svbool_t pg, REALTYPE* vp, Vsimd_t& vt,
                             svint_t index)
{
  for(int k = 0; k < VLEN; ++k){
    if(pg.v[k]) vp[index.v[k]] = vt.v[k];
  }
}

template <typename REALTYPE>
inline void load_add(svbool_t pg, Vsimd_t& vt, REALTYPE *vp)
{
  for(int k = 0; k < VLEN; ++k){
    if(pg.v[k]) vt.v[k] = vp[k];
  }
}

template <typename REALTYPE>
inline void load_add_gather(svbool_t pg, Vsimd_t& v, REALTYPE* vp,
                            Isimd_t& index)
{
  for(int k = 0; k < VLEN; ++k){
    if(pg.v[k]) v.v[k] = vp[index.v[k]];
  }
}

template <typename REALTYPE>
inline void load_add_gather(svbool_t pg, Vsimd_t& vt, REALTYPE* vp,
                            Usimd_t& index, int skip)
{
  for(int k = 0; k < VLEN; ++k){
    if(pg.v[k]) vt.v[k] = vp[index.v[k]];
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
      vt.v[kx + VLENX * ky + VLEN * in] = v[ky + VLENY * in];
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

template <typename REALTYPE>
inline void save_vec(svbool_t pg, REALTYPE *x, Vsimd_t& vt)
{
  for(int k = 0; k < VLEN; ++k){
    if(pg.v[k]) x[k] = vt.v[k];
  }
}

template <typename REALTYPE>
inline void svst1_scatter_index(svbool_t pg, REALTYPE *x,
				Isimd_t index, Vsimd_t& vt)
{
  for(int k = 0; k < VLEN; ++k){
    if(pg.v[k]) x[index.v[k]] = vt.v[k];
  }
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

template <typename REALTYPE>
  inline void shift_vec(svbool_t pg, svuint_t idx,
                        svreal_t& v, REALTYPE* xc, REALTYPE* xn)
{
  svreal_t vt;
  for(int k = 0; k < VLEN; ++k){
    if(pg.v[k]){
      vt.v[k] = xn[k];
    }else{
      vt.v[k] = xc[k];
    }
  }
  for(int k = 0; k < VLEN; ++k){
    v.v[k] = vt.v[idx.v[k]];
  }
}

template <typename REALTYPE>
inline void shift_vec_xbw(svbool_t& pg1, svbool_t& pg2,
                          svreal_t& v, REALTYPE* wx, REALTYPE* wn)
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
                          svreal_t& v, REALTYPE* wx, REALTYPE* wn)
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
inline void shift_vec_ybw(svbool_t& pg1, svbool_t& pg2,
                          svreal_t& v, REALTYPE* wx, REALTYPE* wn)
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
inline void shift_vec_ybw(svreal_t& v, REALTYPE* wx, REALTYPE* wn)
{
  for(int k = 0; k < VLEN - VLENX; ++k){
    v.v[k] = wx[k + VLENX];
  }
  for(int k = 0; k < VLENX; ++k){
    v.v[k + VLEN - VLENX] = wn[k];
  }
}

template <typename REALTYPE>
inline void shift_vec_yfw(svreal_t& v, REALTYPE *wx, REALTYPE *wn)
{
  for(int k = 0; k < VLENX; ++k){
    v.v[k] = wn[VLEN - VLENX + k];
  }
  for(int k = VLENX; k < VLEN; ++k){
    v.v[k] = wx[k - VLENX];
  }
}

inline void clear_vec(svbool_t pg, Vsimd_t& vt)
{
  for(int k = 0; k < VLEN; ++k){
    if(pg.v[k]) vt.v[k] = 0.0;
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

// to be discarded
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

inline void add_vec(svbool_t pg, Vsimd_t& x, Vsimd_t& y)
 {
   for(int k = 0; k < VLEN; ++k){
     if(pg.v[k]) x.v[k] += y.v[k];
   }
 }

inline void add_vec(svbool_t pg, Vsimd_t& z, Vsimd_t& x, Vsimd_t& y)
 {
   for(int k = 0; k < VLEN; ++k){
     if(pg.v[k]) z.v[k] = x.v[k] + y.v[k];
   }
 }

inline void sub_vec(svbool_t pg, Vsimd_t& x, Vsimd_t& y)
 {
   for(int k = 0; k < VLEN; ++k){
     if(pg.v[k]) x.v[k] -= y.v[k];
   }
 }

inline void sub_vec(svbool_t pg, Vsimd_t& z, Vsimd_t& x, Vsimd_t& y)
 {
   for(int k = 0; k < VLEN; ++k){
     if(pg.v[k]) z.v[k] = x.v[k] - y.v[k];
   }
 }

inline void mul_vec(svbool_t pg, Vsimd_t& w, Vsimd_t& x, Vsimd_t& y)
 {
   for(int k = 0; k < VLEN; ++k){
     if(pg.v[k]) w.v[k] = x.v[k] * y.v[k];
   }
 }

inline void copy_vec(Vsimd_t *x, Vsimd_t *y, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int k = 0; k < VLEN; ++k){
      x[in].v[k] = y[in].v[k];
    }
  }
}

// to be discarded
template <typename REALTYPE>
inline void set_vec(Vsimd_t *x, REALTYPE a, Vsimd_t *y, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int k = 0; k < VLEN; ++k){
      x[in].v[k] = a * y[in].v[k];
    }
  }
}

template <typename REALTYPE>
inline void set_vec(svbool_t pg, Vsimd_t& x, REALTYPE a)
{
  for(int k = 0; k < VLEN; ++k){
    if(pg.v[k]) x.v[k] = a;
  }
}

// to be discarded
template <typename REALTYPE>
inline void axpy_vec(Vsimd_t *y, REALTYPE a, Vsimd_t *x, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int k = 0; k < VLEN; ++k){
      y[in].v[k] += a * x[in].v[k];
    }
  }
}

template <typename REALTYPE>
inline void axpy_vec(svbool_t pg, Vsimd_t& y, REALTYPE a, Vsimd_t& x)
{
  for(int k = 0; k < VLEN; ++k){
    if(pg.v[k]) y.v[k] += a * x.v[k];
  }
}

inline void axpy_vec(svbool_t pg, Vsimd_t& y, Vsimd_t& a, Vsimd_t& x)
{
  for(int k = 0; k < VLEN; ++k){
    if(pg.v[k]) y.v[k] += a.v[k] * x.v[k];
  }
}

inline void ymax_vec(svbool_t pg, Vsimd_t& y, Vsimd_t& a, Vsimd_t& x)
{
  for(int k = 0; k < VLEN; ++k){
    if(pg.v[k]) y.v[k] -= a.v[k] * x.v[k];
  }
}

// to be discarded
template <typename REALTYPE>
inline void aypx_vec(REALTYPE a, Vsimd_t *x, Vsimd_t *y, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int k = 0; k < VLEN; ++k){
      x[in].v[k] = a * x[in].v[k] + y[in].v[k];
    }
  }
}

template <typename REALTYPE>
inline void aypx_vec(svbool_t pg, REALTYPE a, Vsimd_t& x, Vsimd_t& y)
{
  for(int k = 0; k < VLEN; ++k){
    if(pg.v[k]) x.v[k] = a * x.v[k] + y.v[k];
  }
}

template <typename REALTYPE>
inline void scal_vec(svbool_t pg, Vsimd_t& x, REALTYPE a)
{
  for(int k = 0; k < VLEN; ++k){
    if(pg.v[k]) x.v[k] *= a;
  }
}

inline void scal_vec(svbool_t pg, Vsimd_t& x, Vsimd_t& a)
{
  for(int k = 0; k < VLEN; ++k){
    if(pg.v[k]) x.v[k] *= a.v[k];
  }
}

// to be discarded
template <typename REALTYPE>
inline void scal_vec(Vsimd_t *x, REALTYPE a, int Nin)
{
  for(int in = 0; in < Nin; ++in){
    for(int k = 0; k < VLEN; ++k){
      x[in].v[k] *= a;
    }
  }
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

// to be discarded
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
inline void reduce_vec(svbool_t pg, REALTYPE& a, Vsimd_t& x)
  {
    a = REALTYPE(0.0);
    for(int k = 0; k < VLEN; ++k){
      if(pg.v[k]) a += x.v[k];
    }
  }

// to be discarded
inline void add_norm2_vec(Vsimd_t *y, Vsimd_t *x, int Nin)
 {
   for(int in = 0; in < Nin; ++in){
     for(int k = 0; k < VLEN; ++k){
       y[in].v[k] += x[in].v[k] * x[in].v[k];
     }
   }
 }

inline void add_norm2_vec(svbool_t pg, Vsimd_t& y, Vsimd_t& x)
 {
   for(int k = 0; k < VLEN; ++k){
     if(pg.v[k]) y.v[k] += x.v[k] * x.v[k];
   }
 }

// to be discaeded
inline void add_dot_vec(Vsimd_t *y, Vsimd_t *x, Vsimd_t *w, int Nin)
 {
   for(int in = 0; in < Nin; ++in){
     for(int k = 0; k < VLEN; ++k){
       y[in].v[k] += x[in].v[k] * w[in].v[k];
     }
   }
 }

inline void add_dot_vec(svbool_t pg, Vsimd_t& y, Vsimd_t& x, Vsimd_t& w)
 {
   for(int k = 0; k < VLEN; ++k){
     if(pg.v[k]) y.v[k] += x.v[k] * w.v[k];
   }
 }

// to be discarded
inline void sub_dot_vec(Vsimd_t *y, Vsimd_t *x, Vsimd_t *w, int Nin)
 {
   for(int in = 0; in < Nin; ++in){
     for(int k = 0; k < VLEN; ++k){
       y[in].v[k] -= x[in].v[k] * w[in].v[k];
     }
   }
 }

inline void sub_dot_vec(svbool_t pg, Vsimd_t& y, Vsimd_t& x, Vsimd_t& w)
 {
   for(int k = 0; k < VLEN; ++k){
     if(pg.v[k]) y.v[k] -= x.v[k] * w.v[k];
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

inline void set_vec(svbool_t pg, svreal_t& x, real_t a, svreal_t& y)
{
  for(int k = 0; k < VLEN; ++k){
    if(pg.v[k]) x.v[k] = a * y.v[k];
  }
}
  
inline void set_vec(svbool_t pg, svreal_t& x, svreal_t& y)
{
  for(int k = 0; k < VLEN; ++k){
    if(pg.v[k]) x.v[k] = y.v[k];
  }
}
  
inline void sqrt_vec(svbool_t pg, svreal_t& x, svreal_t& y)
{
  for(int k = 0; k < VLEN; ++k){
    if(pg.v[k]) x.v[k] = sqrt(y.v[k]);
  }
}
  
inline void div_vec(svbool_t pg, svreal_t& x, svreal_t& y, svreal_t& z)
{
  for(int k = 0; k < VLEN; ++k){
    if(pg.v[k]) x.v[k] = y.v[k] / z.v[k];
  }
}
  
inline svreal_t svadd_z(svbool_t pg, svreal_t wt1, svreal_t wt2)
{
  svreal_t vt;
  for(int k = 0; k < VLEN; ++k){
    if(pg.v[k]){
      vt.v[k] = wt1.v[k] + wt2.v[k];
    }else{
      vt.v[k] = 0.0;
    }
  }
  return vt;
}

inline svreal_t svsub_z(svbool_t pg, svreal_t wt1, svreal_t wt2)
{
  svreal_t vt;
  for(int k = 0; k < VLEN; ++k){
    if(pg.v[k]){
      vt.v[k] = wt1.v[k] - wt2.v[k];
    }else{
      vt.v[k] = 0.0;
    }
  }
  return vt;
}

inline svreal_t svadd_m(svbool_t pg, svreal_t wt1, svreal_t wt2)
{
  svreal_t vt;
  for(int k = 0; k < VLEN; ++k){
    if(pg.v[k]) vt.v[k] = wt1.v[k] + wt2.v[k];
  }
  return vt;
}

inline svreal_t svsub_m(svbool_t pg, svreal_t wt1, svreal_t wt2)
{
  svreal_t vt;
  for(int k = 0; k < VLEN; ++k){
    if(pg.v[k]) vt.v[k] = wt1.v[k] - wt2.v[k];
  }
  return vt;
}

inline void flip_sign(svbool_t pg, svreal_t& vt)
{
  for(int k = 0; k < VLEN; ++k){
    if(pg.v[k]) vt.v[k] = - vt.v[k];
  }
}

 template<typename REALTYPE>
 inline svreal_t svmla_m(svbool_t pg, svreal_t wt1, svreal_t wt2,
                          REALTYPE a)
 {
  svreal_t vt;
  for(int k = 0; k < VLEN; ++k){
    if(pg.v[k]) vt.v[k] = wt1.v[k] + a * wt2.v[k];
  }
  return vt;
 }

inline svreal_t compact_vec(svbool_t pg, svreal_t& yt)
{
  svreal_t vt;
  int j = 0;
  for(int k = 0; k < VLEN; ++k){
    if(pg.v[k]){
      vt.v[j] = yt.v[k];
      ++j;
    }
  }
  for(int k = j; k < VLEN; ++k){
    vt.v[j] = 0.0;
  }
  return vt;
}

inline svbool_t set_predicate_whilelt(int n)
{
  svbool_t pg;
  for(int k = 0; k < n; ++k){
    pg.v[k] = true;
  }
  for(int k = n; k < VLEN; ++k){
    pg.v[k] = false;
  }
  return pg;
}

inline void set1_at(const int i, svbool_t &pg){
  if(pg.v[i] == true){
    pg.v[i] = false;
  }else{
    pg.v[i] = true;
  }
}

inline void rot1_R(uint_t *u)
{
  uint_t tmp = u[VLENX-1]; // tail
  for(int i = VLENX-1; i >= 1; --i){
    u[i] = u[i-1];
  }
  u[0] = tmp;
}

inline void rot1_L(uint_t *u)
{
  uint_t tmp = u[0]; // head                                                
  for(int i = 0; i < VLENX-1; ++i){
    u[i] = u[i+1];
  }
  u[VLENX-1] = tmp;
}

inline void set_idx_predicate_xp_eo(svbool_t& pg, svuint_t& idx,
                                    const int ieo)
{
  uint_t u[VLEN];
  for(int i = 0; i < VLEN; ++i){
    u[i] = i;
  }
  pg = set_predicate_false();
  if(ieo == 0){
    // L-shift odd rows                                                       
    for(int i = VLENX; i < VLEN; i += 2*VLENX){
      set1_at(i, pg);
      rot1_L(u+i);
    }
  }
  if(ieo == 1){
    // L-shift env rows                                                       
    for(int i = 0; i < VLEN; i += 2*VLENX){
      set1_at(i, pg);
      rot1_L(u+i);
    }
  }
  for(int k = 0; k < VLEN; ++k){
    idx.v[k] = u[k];
  }
}

inline void set_idx_predicate_xm_eo(svbool_t& pg, svuint_t& idx,
                                    const int ieo)
{
  uint_t u[VLEN];
  for(int i = 0; i<VLEN; i++){
    u[i] = i;
  }
  pg = set_predicate_false();
  if(ieo == 0){
    // R-shift evn rows                                                       
    for(int i = 0; i < VLEN; i += 2*VLENX){
      set1_at(i + VLENX-1, pg); // 3, 11
      rot1_R(u+i);
    }
  }
  if(ieo == 1){
    // R-shift odd rows                                                       
    for(int i = VLENX; i < VLEN; i += 2*VLENX){
      set1_at(i + VLENX-1, pg); // 7, 15                                     
      rot1_R(u+i);
    }
  }
  for(int k = 0; k < VLEN; ++k){
    idx.v[k] = u[k];
  }
}

inline void set_idx_predicate_yp(svbool_t& pg1, svuint_t& idx)
{
  pg1 = set_predicate_whilelt(VLENX);
  uint_t u[VLEN];
  for(int i = 0; i < VLENX * (VLENY-1); ++i){
    u[i] = i + VLENX;
  }
  for(int i = 0; i < VLENX; ++i){
    u[i + VLENX * (VLENY-1)] = i;
  }
  for(int k = 0; k < VLEN; ++k){
    idx.v[k] = u[k];
  }
}

inline void set_idx_predicate_ym(svbool_t& pg1, svuint_t& idx)
{
  svbool_t pg2 = set_predicate_whilelt(VLENX * (VLENY-1));
  pg1 = set_predicate_not(pg2);
  uint_t u[VLEN];
  for(int i = 0; i < VLENX; ++i){
    u[i] = i + VLENX * (VLENY-1);
  }
  for(int i = VLENX; i < VLEN; ++i){
    u[i] = i - VLENX;
  }
  for(int k = 0; k < VLEN; ++k){
    idx.v[k] = u[k];
  }
}




} // end of nameless namespace

#endif
//============================================================END=====
