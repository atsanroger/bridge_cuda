/*!
        @file    afield_Spinor-inc.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-10-17 13:51:32 #$
        @version $LastChangedRevision: 2550 $
*/

#ifndef QXS_AFIELD_SPINOR_INC_INCLUDED
#define QXS_AFIELD_SPINOR_INC_INCLUDED

#include <cstdlib> 

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/inline/afield_th-inc.h"

#include "lib_alt_QXS/Field/afield_Gauge_index-inc.h"

namespace {
  inline int idxSr(int ic, int id){ return     2*(id + ND * ic); }
  inline int idxSi(int ic, int id){ return 1 + 2*(id + ND * ic); }
}

namespace QXS_Spinor{

//====================================================================
template <typename REALTYPE>
void mult_Gnv(AField<REALTYPE,QXS>& v,       const int exv,
              const AField<REALTYPE,QXS>& u, const int exu,
              const AField<REALTYPE,QXS>& w, const int exw)
{
#pragma omp barrier

  int Nst  = w.nvol();
  int Nstv = Nst/VLEN;

  typedef AField<REALTYPE,QXS> AFIELD;

  REALTYPE* vp = v.ptr(NVCD * Nst * exv);
  REALTYPE* up = const_cast<AFIELD*>(&u)->ptr(NDF * Nst * exu);
  REALTYPE* wp = const_cast<AFIELD*>(&w)->ptr(NVCD * Nst * exw);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstv);

  svbool_t pg = set_predicate();

  for(int site = is; site < ns; ++site){
    REALTYPE* vpt = &vp[VLEN * NVCD * site];
    REALTYPE* upt = &up[VLEN * NDF  * site];
    REALTYPE* wpt = &wp[VLEN * NVCD * site];

    svreal_t u00r, u00i, u01r, u01i, u02r, u02i;
    load_vec(pg, u00r, &upt[VLEN * idxGr(0,0)]);
    load_vec(pg, u00i, &upt[VLEN * idxGi(0,0)]);
    load_vec(pg, u01r, &upt[VLEN * idxGr(0,1)]);
    load_vec(pg, u01i, &upt[VLEN * idxGi(0,1)]);
    load_vec(pg, u02r, &upt[VLEN * idxGr(0,2)]);
    load_vec(pg, u02i, &upt[VLEN * idxGi(0,2)]);

    svreal_t u10r, u10i, u11r, u11i, u12r, u12i;
    load_vec(pg, u10r, &upt[VLEN * idxGr(1,0)]);
    load_vec(pg, u10i, &upt[VLEN * idxGi(1,0)]);
    load_vec(pg, u11r, &upt[VLEN * idxGr(1,1)]);
    load_vec(pg, u11i, &upt[VLEN * idxGi(1,1)]);
    load_vec(pg, u12r, &upt[VLEN * idxGr(1,2)]);
    load_vec(pg, u12i, &upt[VLEN * idxGi(1,2)]);

    svreal_t u20r, u20i, u21r, u21i, u22r, u22i;
    load_vec(pg, u20r, &upt[VLEN * idxGr(2,0)]);
    load_vec(pg, u20i, &upt[VLEN * idxGi(2,0)]);
    load_vec(pg, u21r, &upt[VLEN * idxGr(2,1)]);
    load_vec(pg, u21i, &upt[VLEN * idxGi(2,1)]);
    load_vec(pg, u22r, &upt[VLEN * idxGr(2,2)]);
    load_vec(pg, u22i, &upt[VLEN * idxGi(2,2)]);

    for(int id = 0; id < ND; ++id){

      svreal_t w0r, w0i, w1r, w1i, w2r, w2i;
      load_vec(pg, w0r, &wpt[VLEN * idxSr(0,id)]);
      load_vec(pg, w0i, &wpt[VLEN * idxSi(0,id)]);
      load_vec(pg, w1r, &wpt[VLEN * idxSr(1,id)]);
      load_vec(pg, w1i, &wpt[VLEN * idxSi(1,id)]);
      load_vec(pg, w2r, &wpt[VLEN * idxSr(2,id)]);
      load_vec(pg, w2i, &wpt[VLEN * idxSi(2,id)]);

      svreal_t vt0r, vt0i, vt1r, vt1i, vt2r, vt2i;
      mul_vec( pg, vt0r, u00r, w0r);
      mul_vec( pg, vt0i, u00r, w0i);
      mul_vec( pg, vt1r, u10r, w0r);
      mul_vec( pg, vt1i, u10r, w0i);
      mul_vec( pg, vt2r, u20r, w0r);
      mul_vec( pg, vt2i, u20r, w0i);

      ymax_vec(pg, vt0r, u00i, w0i);
      axpy_vec(pg, vt0i, u00i, w0r);
      ymax_vec(pg, vt1r, u10i, w0i);
      axpy_vec(pg, vt1i, u10i, w0r);
      ymax_vec(pg, vt2r, u20i, w0i);
      axpy_vec(pg, vt2i, u20i, w0r);

      axpy_vec(pg, vt0r, u01r, w1r);
      axpy_vec(pg, vt0i, u01r, w1i);
      axpy_vec(pg, vt1r, u11r, w1r);
      axpy_vec(pg, vt1i, u11r, w1i);
      axpy_vec(pg, vt2r, u21r, w1r);
      axpy_vec(pg, vt2i, u21r, w1i);

      ymax_vec(pg, vt0r, u01i, w1i);
      axpy_vec(pg, vt0i, u01i, w1r);
      ymax_vec(pg, vt1r, u11i, w1i);
      axpy_vec(pg, vt1i, u11i, w1r);
      ymax_vec(pg, vt2r, u21i, w1i);
      axpy_vec(pg, vt2i, u21i, w1r);

      axpy_vec(pg, vt0r, u02r, w2r);
      axpy_vec(pg, vt0i, u02r, w2i);
      axpy_vec(pg, vt1r, u12r, w2r);
      axpy_vec(pg, vt1i, u12r, w2i);
      axpy_vec(pg, vt2r, u22r, w2r);
      axpy_vec(pg, vt2i, u22r, w2i);

      ymax_vec(pg, vt0r, u02i, w2i);
      axpy_vec(pg, vt0i, u02i, w2r);
      ymax_vec(pg, vt1r, u12i, w2i);
      axpy_vec(pg, vt1i, u12i, w2r);
      ymax_vec(pg, vt2r, u22i, w2i);
      axpy_vec(pg, vt2i, u22i, w2r);

      save_vec(pg, &vpt[VLEN * idxSr(0,id)], vt0r);
      save_vec(pg, &vpt[VLEN * idxSi(0,id)], vt0i);
      save_vec(pg, &vpt[VLEN * idxSr(1,id)], vt1r);
      save_vec(pg, &vpt[VLEN * idxSi(1,id)], vt1i);
      save_vec(pg, &vpt[VLEN * idxSr(2,id)], vt2r);
      save_vec(pg, &vpt[VLEN * idxSi(2,id)], vt2i);
    }

  }

#pragma omp barrier
}

//====================================================================
template <typename REALTYPE>
void mult_Gdv(AField<REALTYPE,QXS>& v,       const int exv,
              const AField<REALTYPE,QXS>& u, const int exu,
              const AField<REALTYPE,QXS>& w, const int exw)
{
#pragma omp barrier

  int Nst  = w.nvol();
  int Nstv = Nst/VLEN;

  typedef AField<REALTYPE,QXS> AFIELD;

  REALTYPE* vp = v.ptr(NVCD * Nst * exv);
  REALTYPE* up = const_cast<AFIELD*>(&u)->ptr(NDF * Nst * exu);
  REALTYPE* wp = const_cast<AFIELD*>(&w)->ptr(NVCD * Nst * exw);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstv);

  svbool_t pg = set_predicate();

  for(int site = is; site < ns; ++site){
    REALTYPE* vpt = &vp[VLEN * NVCD * site];
    REALTYPE* upt = &up[VLEN * NDF  * site];
    REALTYPE* wpt = &wp[VLEN * NVCD * site];

    svreal_t u00r, u00i, u01r, u01i, u02r, u02i;
    load_vec(pg, u00r, &upt[VLEN * idxGr(0,0)]);
    load_vec(pg, u00i, &upt[VLEN * idxGi(0,0)]);
    load_vec(pg, u01r, &upt[VLEN * idxGr(0,1)]);
    load_vec(pg, u01i, &upt[VLEN * idxGi(0,1)]);
    load_vec(pg, u02r, &upt[VLEN * idxGr(0,2)]);
    load_vec(pg, u02i, &upt[VLEN * idxGi(0,2)]);

    svreal_t u10r, u10i, u11r, u11i, u12r, u12i;
    load_vec(pg, u10r, &upt[VLEN * idxGr(1,0)]);
    load_vec(pg, u10i, &upt[VLEN * idxGi(1,0)]);
    load_vec(pg, u11r, &upt[VLEN * idxGr(1,1)]);
    load_vec(pg, u11i, &upt[VLEN * idxGi(1,1)]);
    load_vec(pg, u12r, &upt[VLEN * idxGr(1,2)]);
    load_vec(pg, u12i, &upt[VLEN * idxGi(1,2)]);

    svreal_t u20r, u20i, u21r, u21i, u22r, u22i;
    load_vec(pg, u20r, &upt[VLEN * idxGr(2,0)]);
    load_vec(pg, u20i, &upt[VLEN * idxGi(2,0)]);
    load_vec(pg, u21r, &upt[VLEN * idxGr(2,1)]);
    load_vec(pg, u21i, &upt[VLEN * idxGi(2,1)]);
    load_vec(pg, u22r, &upt[VLEN * idxGr(2,2)]);
    load_vec(pg, u22i, &upt[VLEN * idxGi(2,2)]);

    for(int id = 0; id < ND; ++id){

      svreal_t w0r, w0i, w1r, w1i, w2r, w2i;
      load_vec(pg, w0r, &wpt[VLEN * idxSr(0,id)]);
      load_vec(pg, w0i, &wpt[VLEN * idxSi(0,id)]);
      load_vec(pg, w1r, &wpt[VLEN * idxSr(1,id)]);
      load_vec(pg, w1i, &wpt[VLEN * idxSi(1,id)]);
      load_vec(pg, w2r, &wpt[VLEN * idxSr(2,id)]);
      load_vec(pg, w2i, &wpt[VLEN * idxSi(2,id)]);

      svreal_t vt0r, vt0i, vt1r, vt1i, vt2r, vt2i;
      mul_vec( pg, vt0r, u00r, w0r);
      mul_vec( pg, vt0i, u00r, w0i);
      mul_vec( pg, vt1r, u01r, w0r);
      mul_vec( pg, vt1i, u01r, w0i);
      mul_vec( pg, vt2r, u02r, w0r);
      mul_vec( pg, vt2i, u02r, w0i);

      axpy_vec(pg, vt0r, u00i, w0i);
      ymax_vec(pg, vt0i, u00i, w0r);
      axpy_vec(pg, vt1r, u01i, w0i);
      ymax_vec(pg, vt1i, u01i, w0r);
      axpy_vec(pg, vt2r, u02i, w0i);
      ymax_vec(pg, vt2i, u02i, w0r);

      axpy_vec(pg, vt0r, u10r, w1r);
      axpy_vec(pg, vt0i, u10r, w1i);
      axpy_vec(pg, vt1r, u11r, w1r);
      axpy_vec(pg, vt1i, u11r, w1i);
      axpy_vec(pg, vt2r, u12r, w1r);
      axpy_vec(pg, vt2i, u12r, w1i);

      axpy_vec(pg, vt0r, u10i, w1i);
      ymax_vec(pg, vt0i, u10i, w1r);
      axpy_vec(pg, vt1r, u11i, w1i);
      ymax_vec(pg, vt1i, u11i, w1r);
      axpy_vec(pg, vt2r, u12i, w1i);
      ymax_vec(pg, vt2i, u12i, w1r);

      axpy_vec(pg, vt0r, u20r, w2r);
      axpy_vec(pg, vt0i, u20r, w2i);
      axpy_vec(pg, vt1r, u21r, w2r);
      axpy_vec(pg, vt1i, u21r, w2i);
      axpy_vec(pg, vt2r, u22r, w2r);
      axpy_vec(pg, vt2i, u22r, w2i);

      axpy_vec(pg, vt0r, u20i, w2i);
      ymax_vec(pg, vt0i, u20i, w2r);
      axpy_vec(pg, vt1r, u21i, w2i);
      ymax_vec(pg, vt1i, u21i, w2r);
      axpy_vec(pg, vt2r, u22i, w2i);
      ymax_vec(pg, vt2i, u22i, w2r);

      save_vec(pg, &vpt[VLEN * idxSr(0,id)], vt0r);
      save_vec(pg, &vpt[VLEN * idxSi(0,id)], vt0i);
      save_vec(pg, &vpt[VLEN * idxSr(1,id)], vt1r);
      save_vec(pg, &vpt[VLEN * idxSi(1,id)], vt1i);
      save_vec(pg, &vpt[VLEN * idxSr(2,id)], vt2r);
      save_vec(pg, &vpt[VLEN * idxSi(2,id)], vt2i);
    }

  }

#pragma omp barrier
}

//====================================================================
template <typename REALTYPE>
void tensorProd(AField<REALTYPE,QXS>& u, const int mu,
                const AField<REALTYPE,QXS>& v,
                const AField<REALTYPE,QXS>& w)
{
#pragma omp barrier

  int Nst  = w.nvol();
  int Nstv = Nst/VLEN;

  typedef AField<REALTYPE,QXS> AFIELD;

  REALTYPE* up = u.ptr(NDF * Nst * mu);
  REALTYPE* vp = const_cast<AFIELD*>(&v)->ptr(0);
  REALTYPE* wp = const_cast<AFIELD*>(&w)->ptr(0);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstv);

  svbool_t pg = set_predicate();

  for(int site = is; site < ns; ++site){
    REALTYPE* upt = &up[VLEN * NDF  * site];
    REALTYPE* vpt = &vp[VLEN * NVCD * site];
    REALTYPE* wpt = &wp[VLEN * NVCD * site];

    for(int ic1 = 0; ic1 < NC; ++ic1){
      for(int ic2 = 0; ic2 < NC; ++ic2){
        svreal_t ur, ui;
        clear_vec(pg, ur);
        clear_vec(pg, ui);
        for(int id = 0; id < ND; ++id){
          svreal_t wr, wi;
          load_vec(pg, wr, &wpt[VLEN * idxSr(ic1,id)]);
          load_vec(pg, wi, &wpt[VLEN * idxSi(ic1,id)]);
          svreal_t vr, vi;
          load_vec(pg, vr, &vpt[VLEN * idxSr(ic2,id)]);
          load_vec(pg, vi, &vpt[VLEN * idxSi(ic2,id)]);
          axpy_vec(pg, ur, vr, wr);
          axpy_vec(pg, ui, vr, wi);
          axpy_vec(pg, ur, vi, wi);
          ymax_vec(pg, ui, vi, wr);
	}
        save_vec(pg, &upt[VLEN * idxGr(ic1,ic2)], ur);
        save_vec(pg, &upt[VLEN * idxGi(ic1,ic2)], ui);
      }
    }

  }

#pragma omp barrier
}


} // namespace QXS_Gauge

//============================================================END=====
#endif
