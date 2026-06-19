/*!
      @file    vsimd_Domainwall_SU3_double-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2581 $
*/

#ifndef QXS_VSIMD_DOMAINWALL_SU3_DOUBLE_INC_INCLUDED
#define QXS_VSIMD_DOMAINWALL_SU3_DOUBLE_INC_INCLUDED

/*!
  CAUTION: this file is not modified from that for QXS impl.


  These inline functions are used in the 5din-type implementation
  of the domain-wall fermion operators for QXS implementation.
  While the gamma_5 matrix of QXS is the same as Bridsge++'s,
  the spinor fields are multipled by gamma_4 before and after
  applying QXS code to compensate the negative sign of the
  spatial gamma matrices. This requires the change of sign of
  the gamma_5 in the domain-wall ferimion implemtation so as to
  give the same propagator as the Bridge++'s. This change is
  incorporated in the following inline functions.
                                         [H.Matsufuru, 5 Dec 2019]
 */

namespace {


template<typename REALTYPE>
inline void set_aPp5_dirac_vec(Vsimd_t *v,
                               REALTYPE a, Vsimd_t *w, int Nc)
{
  for(int ivc = 0; ivc < NVC; ++ivc){
    int ivc2 = (ivc % 2) + 2 * ND * (ivc/2);
    for(int k = 0; k < VLEN; ++k){
      v[ID1+ivc2].v[k] = a * (w[ID1+ivc2].v[k] - w[ID3+ivc2].v[k]);
      v[ID2+ivc2].v[k] = a * (w[ID2+ivc2].v[k] - w[ID4+ivc2].v[k]);
      v[ID3+ivc2].v[k] = a * (w[ID3+ivc2].v[k] - w[ID1+ivc2].v[k]);
      v[ID4+ivc2].v[k] = a * (w[ID4+ivc2].v[k] - w[ID2+ivc2].v[k]);
    }
  }
}

template<typename REALTYPE>
inline void set_aPm5_dirac_vec(Vsimd_t *v,
                               REALTYPE a, Vsimd_t *w, int Nc)
{
  for(int ivc = 0; ivc < NVC; ++ivc){
    int ivc2 = (ivc % 2) + 2 * ND * (ivc/2);
    for(int k = 0; k < VLEN; ++k){
      v[ID1+ivc2].v[k] = a * (w[ID1+ivc2].v[k] + w[ID3+ivc2].v[k]);
      v[ID2+ivc2].v[k] = a * (w[ID2+ivc2].v[k] + w[ID4+ivc2].v[k]);
      v[ID3+ivc2].v[k] = a * (w[ID3+ivc2].v[k] + w[ID1+ivc2].v[k]);
      v[ID4+ivc2].v[k] = a * (w[ID4+ivc2].v[k] + w[ID2+ivc2].v[k]);
    }
  }
}

template<typename REALTYPE>
inline void add_aPp5_dirac_vec(Vsimd_t *v,
                               REALTYPE a, Vsimd_t *w, int Nc)
{
  for(int ivc = 0; ivc < NVC; ++ivc){
    int ivc2 = (ivc % 2) + 2 * ND * (ivc/2);
    for(int k = 0; k < VLEN; ++k){
      v[ID1+ivc2].v[k] += a * (w[ID1+ivc2].v[k] - w[ID3+ivc2].v[k]);
      v[ID2+ivc2].v[k] += a * (w[ID2+ivc2].v[k] - w[ID4+ivc2].v[k]);
      v[ID3+ivc2].v[k] += a * (w[ID3+ivc2].v[k] - w[ID1+ivc2].v[k]);
      v[ID4+ivc2].v[k] += a * (w[ID4+ivc2].v[k] - w[ID2+ivc2].v[k]);
    }
  }
}

template<typename REALTYPE>
inline void add_aPm5_dirac_vec(Vsimd_t *v,
                               REALTYPE a, Vsimd_t *w, int Nc)
{
  for(int ivc = 0; ivc < NVC; ++ivc){
    int ivc2 = (ivc % 2) + 2 * ND * (ivc/2);
    for(int k = 0; k < VLEN; ++k){
      v[ID1+ivc2].v[k] += a * (w[ID1+ivc2].v[k] + w[ID3+ivc2].v[k]);
      v[ID2+ivc2].v[k] += a * (w[ID2+ivc2].v[k] + w[ID4+ivc2].v[k]);
      v[ID3+ivc2].v[k] += a * (w[ID3+ivc2].v[k] + w[ID1+ivc2].v[k]);
      v[ID4+ivc2].v[k] += a * (w[ID4+ivc2].v[k] + w[ID2+ivc2].v[k]);
    }
  }
}



template<typename REALTYPE>
inline void set_aPm5_dirac_vec(svreal_t& vt1r, svreal_t& vt1i,
                               svreal_t& vt2r, svreal_t& vt2i,
                               svreal_t& vt3r, svreal_t& vt3i,
                               svreal_t& vt4r, svreal_t& vt4i,
                               REALTYPE a, REALTYPE *w, int is, int ic)
{
  svbool_t pg = set_predicate();
  int off_up = 2*ND * ic + NVCD * is;
  svreal_t w3r, w3i, w4r, w4i;

  load_vec(pg, vt1r, &w[VLEN*(ID1   + off_up)]);
  load_vec(pg, vt1i, &w[VLEN*(ID1+1 + off_up)]);
  load_vec(pg, vt2r, &w[VLEN*(ID2   + off_up)]);
  load_vec(pg, vt2i, &w[VLEN*(ID2+1 + off_up)]);

  load_vec(pg, w3r, &w[VLEN*(ID3   + off_up)]);
  load_vec(pg, w3i, &w[VLEN*(ID3+1 + off_up)]);
  load_vec(pg, w4r, &w[VLEN*(ID4   + off_up)]);
  load_vec(pg, w4i, &w[VLEN*(ID4+1 + off_up)]);

  add_vec( pg, vt1r, w3r);
  add_vec( pg, vt1i, w3i);
  scal_vec(pg, vt1r, a);
  scal_vec(pg, vt1i, a);

  vt3r = vt1r;
  vt3i = vt1i;

  add_vec( pg, vt2r, w4r);
  add_vec( pg, vt2i, w4i);
  scal_vec(pg, vt2r, a);
  scal_vec(pg, vt2i, a);

  vt4r = vt2r;
  vt4i = vt2i;

}

template<typename REALTYPE>
inline void add_aPp5_dirac_vec(svreal_t& vt1r, svreal_t& vt1i,
                               svreal_t& vt2r, svreal_t& vt2i,
                               svreal_t& vt3r, svreal_t& vt3i,
                               svreal_t& vt4r, svreal_t& vt4i,
                               REALTYPE a, REALTYPE *w, int is, int ic)
{
  svbool_t pg = set_predicate();
  int off_up = 2*ND * ic + NVCD * is;
  svreal_t w1r, w1i, w2r, w2i, w3r, w3i, w4r, w4i;

  load_vec(pg, w1r, &w[VLEN*(ID1   + off_up)]);
  load_vec(pg, w1i, &w[VLEN*(ID1+1 + off_up)]);
  load_vec(pg, w2r, &w[VLEN*(ID2   + off_up)]);
  load_vec(pg, w2i, &w[VLEN*(ID2+1 + off_up)]);
  load_vec(pg, w3r, &w[VLEN*(ID3   + off_up)]);
  load_vec(pg, w3i, &w[VLEN*(ID3+1 + off_up)]);
  load_vec(pg, w4r, &w[VLEN*(ID4   + off_up)]);
  load_vec(pg, w4i, &w[VLEN*(ID4+1 + off_up)]);

  sub_vec( pg, w1r, w3r);
  sub_vec( pg, w1i, w3i);
  sub_vec( pg, w2r, w4r);
  sub_vec( pg, w2i, w4i);

  axpy_vec(pg, vt1r, a, w1r);
  axpy_vec(pg, vt1i, a, w1i);
  axpy_vec(pg, vt3r, -a, w1r);
  axpy_vec(pg, vt3i, -a, w1i);

  axpy_vec(pg, vt2r, a, w2r);
  axpy_vec(pg, vt2i, a, w2i);
  axpy_vec(pg, vt4r, -a, w2r);
  axpy_vec(pg, vt4i, -a, w2i);

}

template<typename REALTYPE>
inline void add_aPm5_dirac_vec(svreal_t& vt1r, svreal_t& vt1i,
                               svreal_t& vt2r, svreal_t& vt2i,
                               svreal_t& vt3r, svreal_t& vt3i,
                               svreal_t& vt4r, svreal_t& vt4i,
                               REALTYPE a, REALTYPE *w, int is, int ic)
{
  svbool_t pg = set_predicate();
  int off_up = 2*ND * ic + NVCD * is;
  svreal_t w1r, w1i, w2r, w2i, w3r, w3i, w4r, w4i;

  load_vec(pg, w1r, &w[VLEN*(ID1   + off_up)]);
  load_vec(pg, w1i, &w[VLEN*(ID1+1 + off_up)]);
  load_vec(pg, w2r, &w[VLEN*(ID2   + off_up)]);
  load_vec(pg, w2i, &w[VLEN*(ID2+1 + off_up)]);
  load_vec(pg, w3r, &w[VLEN*(ID3   + off_up)]);
  load_vec(pg, w3i, &w[VLEN*(ID3+1 + off_up)]);
  load_vec(pg, w4r, &w[VLEN*(ID4   + off_up)]);
  load_vec(pg, w4i, &w[VLEN*(ID4+1 + off_up)]);

  add_vec( pg, w1r, w3r);
  add_vec( pg, w1i, w3i);
  add_vec( pg, w2r, w4r);
  add_vec( pg, w2i, w4i);

  axpy_vec(pg, vt1r, a, w1r);
  axpy_vec(pg, vt1i, a, w1i);
  axpy_vec(pg, vt3r, a, w1r);
  axpy_vec(pg, vt3i, a, w1i);

  axpy_vec(pg, vt2r, a, w2r);
  axpy_vec(pg, vt2i, a, w2i);
  axpy_vec(pg, vt4r, a, w2r);
  axpy_vec(pg, vt4i, a, w2i);
}

template<typename REALTYPE>
inline void add_aPp5_dirac_vec(svbool_t pg,
                               svreal_t& vt1r, svreal_t& vt1i,
                               svreal_t& vt2r, svreal_t& vt2i,
                               svreal_t& vt3r, svreal_t& vt3i,
                               svreal_t& vt4r, svreal_t& vt4i,
                               REALTYPE a,
                               svreal_t& xt1r, svreal_t& xt1i,
                               svreal_t& xt2r, svreal_t& xt2i,
                               svreal_t& xt3r, svreal_t& xt3i,
                               svreal_t& xt4r, svreal_t& xt4i)
{
  svreal_t yt1r, yt1i, yt2r, yt2i;
  yt1r = svsub_m(pg, xt1r, xt3r);
  yt1i = svsub_m(pg, xt1i, xt3i);
  yt2r = svsub_m(pg, xt2r, xt4r);
  yt2i = svsub_m(pg, xt2i, xt4i);
  axpy_vec(pg, vt1r,  a, yt1r);
  axpy_vec(pg, vt1i,  a, yt1i);
  axpy_vec(pg, vt2r,  a, yt2r);
  axpy_vec(pg, vt2i,  a, yt2i);
  axpy_vec(pg, vt3r, -a, yt1r);
  axpy_vec(pg, vt3i, -a, yt1i);
  axpy_vec(pg, vt4r, -a, yt2r);
  axpy_vec(pg, vt4i, -a, yt2i);
}

template<typename REALTYPE>
inline void add_aPm5_dirac_vec(svbool_t pg,
                               svreal_t& vt1r, svreal_t& vt1i,
                               svreal_t& vt2r, svreal_t& vt2i,
                               svreal_t& vt3r, svreal_t& vt3i,
                               svreal_t& vt4r, svreal_t& vt4i,
                               REALTYPE a,
                               svreal_t& xt1r, svreal_t& xt1i,
                               svreal_t& xt2r, svreal_t& xt2i,
                               svreal_t& xt3r, svreal_t& xt3i,
                               svreal_t& xt4r, svreal_t& xt4i)
{
  svreal_t yt1r, yt1i, yt2r, yt2i;
  yt1r = svadd_m(pg, xt1r, xt3r);
  yt1i = svadd_m(pg, xt1i, xt3i);
  yt2r = svadd_m(pg, xt2r, xt4r);
  yt2i = svadd_m(pg, xt2i, xt4i);
  axpy_vec(pg, vt1r, a, yt1r);
  axpy_vec(pg, vt1i, a, yt1i);
  axpy_vec(pg, vt2r, a, yt2r);
  axpy_vec(pg, vt2i, a, yt2i);
  axpy_vec(pg, vt3r, a, yt1r);
  axpy_vec(pg, vt3i, a, yt1i);
  axpy_vec(pg, vt4r, a, yt2r);
  axpy_vec(pg, vt4i, a, yt2i);
}

template<typename REALTYPE>
inline void set_aPp5_dirac_vec(svbool_t pg,
                               svreal_t& vt1r, svreal_t& vt1i,
                               svreal_t& vt2r, svreal_t& vt2i,
                               svreal_t& vt3r, svreal_t& vt3i,
                               svreal_t& vt4r, svreal_t& vt4i,
                               REALTYPE a,
                               svreal_t& xt1r, svreal_t& xt1i,
                               svreal_t& xt2r, svreal_t& xt2i,
                               svreal_t& xt3r, svreal_t& xt3i,
                               svreal_t& xt4r, svreal_t& xt4i)
{
  vt1r = svsub_m(pg, xt1r, xt3r);
  vt1i = svsub_m(pg, xt1i, xt3i);
  vt2r = svsub_m(pg, xt2r, xt4r);
  vt2i = svsub_m(pg, xt2i, xt4i);
  vt3r = vt1r;
  vt3i = vt1i;
  vt4r = vt2r;
  vt4i = vt2i;

  scal_vec(pg, vt1r,  a);
  scal_vec(pg, vt1i,  a);
  scal_vec(pg, vt2r,  a);
  scal_vec(pg, vt2i,  a);
  scal_vec(pg, vt3r, -a);
  scal_vec(pg, vt3i, -a);
  scal_vec(pg, vt4r, -a);
  scal_vec(pg, vt4i, -a);
}

template<typename REALTYPE>
inline void set_aPm5_dirac_vec(svbool_t pg,
                               svreal_t& vt1r, svreal_t& vt1i,
                               svreal_t& vt2r, svreal_t& vt2i,
                               svreal_t& vt3r, svreal_t& vt3i,
                               svreal_t& vt4r, svreal_t& vt4i,
                               REALTYPE a,
                               svreal_t& xt1r, svreal_t& xt1i,
                               svreal_t& xt2r, svreal_t& xt2i,
                               svreal_t& xt3r, svreal_t& xt3i,
                               svreal_t& xt4r, svreal_t& xt4i)
{
  vt1r = svadd_m(pg, xt1r, xt3r);
  vt1i = svadd_m(pg, xt1i, xt3i);
  vt2r = svadd_m(pg, xt2r, xt4r);
  vt2i = svadd_m(pg, xt2i, xt4i);
  scal_vec(pg, vt1r, a);
  scal_vec(pg, vt1i, a);
  scal_vec(pg, vt2r, a);
  scal_vec(pg, vt2i, a);
  vt3r = vt1r;
  vt3i = vt1i;
  vt4r = vt2r;
  vt4i = vt2i;
}



template <typename REALTYPE>
inline void load_mult_gm5_dirac_vec(svbool_t pg,
                               svreal_t& vt1r, svreal_t& vt1i,
                               svreal_t& vt2r, svreal_t& vt2i,
                               svreal_t& vt3r, svreal_t& vt3i,
                               svreal_t& vt4r, svreal_t& vt4i,
                               REALTYPE *w)
{
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

template<typename REALTYPE>
inline void dw_5dir_axpy(svbool_t pg, REALTYPE *v,
                         REALTYPE *y, REALTYPE *w,
                         REALTYPE a1, REALTYPE a2,
                         REALTYPE b1, REALTYPE b2,
                         svreal_t zt, int index)
{
  //v[i]=a1*w[i]+a2*zt
  //y[i]=-0.5*(b1*w[i]+b2*zt)
  svreal_t vt, wt, yt;
  load_vec(pg, wt, &w[VLEN*index]);
  set_vec( pg, vt, a1, wt);
  axpy_vec(pg, vt, a2, zt);
  save_vec(pg, &v[VLEN*index], vt);

  set_vec( pg, yt, -0.5*b1, wt);
  axpy_vec(pg, yt, -0.5*b2, zt);
  save_vec(pg, &y[VLEN*index], yt);
}

template<typename REALTYPE>
inline void accum_dw_5dir_axpy(svbool_t pg, REALTYPE *v,
                                 REALTYPE *y, REALTYPE *w,
                                 REALTYPE a1, REALTYPE a2,
                                 REALTYPE b1, REALTYPE b2,
                                 svreal_t zt, int index)
{
    //v[i]+=a1*w[i]+a2*zt
    //y[i]=-0.5*(b1*w[i]+b2*zt)
    svreal_t vt, wt, yt;
    load_vec(pg, vt, &v[VLEN * index]);
    load_vec(pg, wt, &w[VLEN * index]);
    axpy_vec(pg, vt, a2, zt);
    axpy_vec(pg, vt, a1, wt);
    save_vec(pg, &v[VLEN * index], vt);

    set_vec(pg, yt, -0.5 * b1, wt);
    axpy_vec(pg, yt, -0.5 * b2, zt);
    save_vec(pg, &y[VLEN * index], yt);
}

template<typename REALTYPE>
inline void dw_5dir_dag(svbool_t pg,
                        svreal_t& vt1r, svreal_t& vt1i,
                        svreal_t& vt2r, svreal_t& vt2i,
                        svreal_t& vt3r, svreal_t& vt3i,
                        svreal_t& vt4r, svreal_t& vt4i,
                        REALTYPE *w, REALTYPE *y,
                        REALTYPE a1, REALTYPE a2, int index)
{
  load_vec(pg, vt1r, &w[VLEN*(ID3   + index)]);
  load_vec(pg, vt1i, &w[VLEN*(ID3+1 + index)]);
  load_vec(pg, vt2r, &w[VLEN*(ID4   + index)]);
  load_vec(pg, vt2i, &w[VLEN*(ID4+1 + index)]);
  load_vec(pg, vt3r, &w[VLEN*(ID1   + index)]);
  load_vec(pg, vt3i, &w[VLEN*(ID1+1 + index)]);
  load_vec(pg, vt4r, &w[VLEN*(ID2   + index)]);
  load_vec(pg, vt4i, &w[VLEN*(ID2+1 + index)]);

  svreal_t yt1r, yt1i;
  load_vec(pg, yt1r, &y[VLEN*(ID1   + index)]);
  load_vec(pg, yt1i, &y[VLEN*(ID1+1 + index)]);
  scal_vec(pg, vt1r, -a1);
  scal_vec(pg, vt1i, -a1);
  axpy_vec(pg, vt1r, a2, yt1r);
  axpy_vec(pg, vt1i, a2, yt1i);

  svreal_t yt2r, yt2i;
  load_vec(pg, yt2r, &y[VLEN*(ID2   + index)]);
  load_vec(pg, yt2i, &y[VLEN*(ID2+1 + index)]);
  scal_vec(pg, vt2r, -a1);
  scal_vec(pg, vt2i, -a1);
  axpy_vec(pg, vt2r, a2, yt2r);
  axpy_vec(pg, vt2i, a2, yt2i);

    svreal_t yt3r, yt3i;
  load_vec(pg, yt3r, &y[VLEN*(ID3   + index)]);
  load_vec(pg, yt3i, &y[VLEN*(ID3+1 + index)]);
  scal_vec(pg, vt3r, -a1);
  scal_vec(pg, vt3i, -a1);
  axpy_vec(pg, vt3r, a2, yt3r);
  axpy_vec(pg, vt3i, a2, yt3i);

  svreal_t yt4r, yt4i;
  load_vec(pg, yt4r, &y[VLEN*(ID4   + index)]);
  load_vec(pg, yt4i, &y[VLEN*(ID4+1 + index)]);
  scal_vec(pg, vt4r, -a1);
  scal_vec(pg, vt4i, -a1);
  axpy_vec(pg, vt4r, a2, yt4r);
  axpy_vec(pg, vt4i, a2, yt4i);
}





} // nameless namespace end

#endif
//============================================================END=====
