/*!
        @file    afield-inc.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-10-02 13:50:31 #$
        @version $LastChangedRevision: 2543 $
*/

#ifndef AFIELD_VECTOR_INC_INCLUDED
#define AFIELD_VECTOR_INC_INCLUDED

#include <cstdlib> 

#include "complexTraits.h"
#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_Vector/inline/afield_th-inc.h" 
#include "lib_alt_Vector/Field/aindex_lex.h"
#include "lib_alt_Vector/BridgeVec/bridgeVec_afield.h"

//====================================================================
template <typename REALTYPE>
void copy(AField<REALTYPE, VECTOR>& v, const AField<REALTYPE, VECTOR> &w)
{
  v.copy(w);
}

//====================================================================
template <typename REALTYPE>
void copy(AField<REALTYPE, VECTOR>& v, const int ex,
          const AField<REALTYPE, VECTOR> &w, const int ex_w)
{
  v.copy(ex, w, ex_w);
}

//====================================================================
template <typename REALTYPE>
void axpy(AField<REALTYPE, VECTOR>& v, const REALTYPE a, const AField<REALTYPE, VECTOR> &w)
{
  v.axpy(a, w);
}

//====================================================================
template <typename REALTYPE>
void axpy(AField<REALTYPE, VECTOR>& v,
          const typename ComplexTraits<REALTYPE>::complex_t a,
          const AField<REALTYPE, VECTOR> &w)
{
  v.axpy(a, w);
}

//====================================================================
template <typename REALTYPE>
void axpy(AField<REALTYPE, VECTOR>& v, const int ex,
          const REALTYPE a, 
          const AField<REALTYPE, VECTOR> &w, const int ex_w)
{
  v.axpy(ex, a, w, ex_w);
}

//====================================================================
template <typename REALTYPE>
void aypx(const REALTYPE a, AField<REALTYPE, VECTOR>& v,
          const AField<REALTYPE, VECTOR>& w)
{
  v.aypx(a, w);
}

//====================================================================
template <typename REALTYPE>
void aypx(const typename ComplexTraits<REALTYPE>::complex_t a,
          AField<REALTYPE, VECTOR>& v,
          const AField<REALTYPE, VECTOR>& w)
{
  v.aypx(a, w);
}

//====================================================================
template <typename REALTYPE>
void scal(AField<REALTYPE, VECTOR>& v, REALTYPE a)
{
  v.scal(a);
}

//====================================================================
template <typename REALTYPE>
REALTYPE dot(AField<REALTYPE, VECTOR>& v, AField<REALTYPE, VECTOR>& w)
{
  return v.dot(w);
}

//====================================================================
template <typename REALTYPE>
typename ComplexTraits<REALTYPE>::complex_t dotc(
                                  const AField<REALTYPE, VECTOR>& v,
                                  const AField<REALTYPE, VECTOR>& w)
{
  REALTYPE vw_r, vw_i;
  v.dotc(vw_r, vw_i, w);
  return cmplx(vw_r, vw_i);
}

//====================================================================
template <typename REALTYPE>
REALTYPE norm2(const AField<REALTYPE, VECTOR>& v)
{
  return v.norm2();
}

//====================================================================
template <typename REALTYPE>
void xI(AField<REALTYPE, VECTOR>& v)
{
  v.xI();
}

//====================================================================
template <typename REALTYPE>
void conjg(AField<REALTYPE, VECTOR>& v)
{
  v.conjg();
}

//====================================================================
template <class INDEX, class AFIELD>
void convert(INDEX& index, AFIELD& v, const Field& w) 
{
  int Nin  = w.nin();
  int Nvol = w.nvol();
  int Nex  = w.nex();
  assert(v.check_size(Nin, Nvol, Nex));

  //  v.clear();

  typename AFIELD::real_t *vp = v.ptr(0);

  size_t Nvol_pad = v.nvol_pad();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol_pad);

  for(int ex = 0; ex < Nex; ++ex){   
   for(int site = is; site < ns; ++site){   

     if(site < Nvol){
       for(int in = 0; in < Nin; ++in){   
         int iw = in + Nin * (site + Nvol * ex);
         int iv = IDXV_BARE(Nin, Nvol_pad, in, site, ex);
         vp[iv] = w.cmp(iw);
       }
     }else{
       for(int in = 0; in < Nin; ++in){   
         int iv = IDXV_BARE(Nin, Nvol_pad, in, site, ex);
         vp[iv] = 0.0;
       }

     }

   }
  }

  vout.detailed("convert finished.\n");

#pragma omp barrier

}  

//====================================================================
template <class INDEX, class AFIELD>
void convert_spinor(INDEX& index, AFIELD& v, const Field& w) 
{
  int Nc = CommonParameters::Nc();
  int Nd = CommonParameters::Nd();
  assert(w.nin() == 2*Nc*Nd);
  assert(v.nin() == 2*Nc*Nd);

  convert(index, v, w);

}  

//====================================================================
template <class INDEX, class AFIELD>
void convert_gauge(INDEX& index, AFIELD& v, const Field& w) 
{
  int Nc = CommonParameters::Nc();
  assert(w.nin() == 2*Nc*Nc);
  assert(v.nin() == 2*Nc*Nc);

  convert(index, v, w);

}  

//====================================================================
template<class INDEX2, class AFIELD2, class INDEX1, class AFIELD1>
void convert(INDEX2& index2, AFIELD2& v2,
             const INDEX1& index1, const AFIELD1& v1) 
{
  int Nin  = v1.nin();
  int Nvol = v1.nvol();
  int Nex  = v1.nex();
  assert(v2.check_size(Nin, Nvol, Nex));

  typename AFIELD1::real_t *v1p = const_cast<AFIELD1*>(&v1)->ptr(0);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

#pragma omp barrier

  typename AFIELD2::real_t *v2p = v2.ptr(0);

  for(int ex = 0; ex < Nex; ++ex){   
   for(int site = is; site < ns; ++site){   
    for(int in = 0; in < Nin; ++in){   
      int iv1 = index1.idx(in, Nin, site, ex);
      int iv2 = index2.idx(in, Nin, site, ex);
      v2p[iv2] = v1.cmp(iv1);
    }
   }
  }

#pragma omp barrier

}  

//====================================================================
template <class INDEX2, class AFIELD2, class INDEX1, class AFIELD1>
void convert_h(INDEX2& index2, AFIELD2& v2,
               const INDEX1& index1, const AFIELD1& v1) 
{
  int Nin  = v1.nin();
  int Nvol = v1.nvol();
  int Nex  = v1.nex();
  assert(v2.check_size(Nin, Nvol, Nex));

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

#pragma omp barrier

  typename AFIELD1::real_t *v1p = const_cast<AFIELD1*>(&v1)->ptr(0);
  typename AFIELD2::real_t *v2p = v2.ptr(0);

  for(int ex = 0; ex < Nex; ++ex){   
   for(int site = is; site < ns; ++site){   
    for(int in = 0; in < Nin; ++in){   
      int iv1 = index1.idxh(in, Nin, site, ex);
      int iv2 = index2.idxh(in, Nin, site, ex);
      v2p[iv2] = v1.cmp(iv1);
    }
   }
  }

#pragma omp barrier

}  

//====================================================================
template <class INDEX, class AFIELD>
void reverse(INDEX& index, Field& v, const AFIELD& w) 
{
  int Nin  = v.nin();
  int Nvol = v.nvol();
  int Nex  = v.nex();
  w.check_size(Nin, Nvol, Nex);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  typedef typename AFIELD::real_t real_t;
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);

#pragma omp barrier

  for(int ex = 0; ex < Nex; ++ex){   
    for(int site = is; site < ns; ++site){   
      for(int in = 0; in < Nin; ++in){   
        int iv = in + Nin * (site + Nvol * ex);
        int iw = index.idx(in, Nin, site, ex);
        v.set(iv, double(w.cmp(iw)));
      }
    }
  }

  vout.detailed("reverse finished.\n");

#pragma omp barrier

}  

//====================================================================
/*
template <class AFIELD>
void reverse(AIndex_lex<typename AFIELD::real_t,VECTOR>& index,
             Field& v, const AFIELD& w) 
{
  int Nin  = v.nin();
  int Nvol = v.nvol();
  int Nex  = v.nex();
  assert(w.check_size(Nin, Nvol, Nex));

  typedef typename AFIELD::real_t real_t;
  double *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);

  //  int Nvol_pad = w.nvol_pad();

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){
    for(int ex = 0; ex < Nex; ++ex){
      double *vp2 = &vp[Nin * Nvol * ex];
      real_t *wp2 = &wp[Nin * Nvol_pad * ex];
      reverse_dev(vp2, wp2, Nin, Nvol, Nvol_pad);
    }
  }

  vout.detailed("reverse from Vector AField finished.\n");

#pragma omp barrier

}  
*/
//====================================================================
template <class INDEX, class AFIELD>
void reverse_spinor(INDEX& index, Field& v, const AFIELD& w) 
{
  int Nc = CommonParameters::Nc();
  int Nd = CommonParameters::Nd();
  assert(w.nin() == 2*Nc*Nd);
  assert(v.nin() == 2*Nc*Nd);

  reverse(index, v, w);

}  

//====================================================================
template <class INDEX, class AFIELD>
void reverse_gauge(INDEX& index, Field& v, const AFIELD& w) 
{
  int Nc = CommonParameters::Nc();
  assert(w.nin() == 2*Nc*Nc);
  assert(v.nin() == 2*Nc*Nc);

  reverse(index, v, w);

}  

//============================================================END=====
#endif
