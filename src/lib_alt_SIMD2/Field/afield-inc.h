/*!
        @file    afield-inc.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
        @version $LastChangedRevision: 2160 $
*/

#ifndef AFIELD_SIMD2_INC_INCLUDED
#define AFIELD_SIMD2_INC_INCLUDED

#include <cstdlib> 

#include "lib/ResourceManager/threadManager_OpenMP.h"

#include "lib_alt_SIMD2/inline/afield_th-inc.h" 

//====================================================================
template <typename REALTYPE>
void copy(AField<REALTYPE,SIMD2>& v, const AField<REALTYPE,SIMD2> &w)
{
  v.copy(w);
}

//====================================================================
template <typename REALTYPE>
void copy(AField<REALTYPE,SIMD2>& v, const int ex,
          const AField<REALTYPE,SIMD2> &w, const int ex_w)
{
  v.copy(ex, w, ex_w);
}

//====================================================================
template <typename REALTYPE>
void axpy(AField<REALTYPE,SIMD2>& v, const REALTYPE a, const AField<REALTYPE,SIMD2> &w)
{
  v.axpy(a, w);
}

//====================================================================
template <typename REALTYPE>
void axpy(AField<REALTYPE,SIMD2>& v,
          const typename ComplexTraits<REALTYPE>::complex_t a,
          const AField<REALTYPE,SIMD2> &w)
{
  v.axpy(a, w);
}

//====================================================================
template <typename REALTYPE>
void axpy(AField<REALTYPE,SIMD2>& v, const int ex,
          const REALTYPE a, 
          const AField<REALTYPE,SIMD2> &w, const int ex_w)
{
  v.axpy(ex, a, w, ex_w);
}

//====================================================================
template <typename REALTYPE>
void aypx(REALTYPE a, AField<REALTYPE,SIMD2>& v, AField<REALTYPE,SIMD2> &w)
{
  v.aypx(a, w);
}

//====================================================================
template <typename REALTYPE>
void aypx(const typename ComplexTraits<REALTYPE>::complex_t a,
          AField<REALTYPE,SIMD2>& v,
          const AField<REALTYPE,SIMD2> &w)
{
  v.aypx(a, w);
}

//====================================================================
template <typename REALTYPE>
void scal(AField<REALTYPE,SIMD2>& v, const REALTYPE a)
{
  v.scal(a);
}

//====================================================================
template <typename REALTYPE>
void scal(AField<REALTYPE,SIMD2>& v,
          const typename ComplexTraits<REALTYPE>::complex_t a)
{
  v.scal(a);
}

//====================================================================
template <typename REALTYPE>
REALTYPE dot(AField<REALTYPE,SIMD2>& v, AField<REALTYPE,SIMD2>& w)
{
  return v.dot(w);
}

//====================================================================
template <typename REALTYPE>
typename ComplexTraits<REALTYPE>::complex_t dotc(
                                   const AField<REALTYPE, SIMD2>& v,
                                   const AField<REALTYPE, SIMD2>& w)
{
  REALTYPE vw_r, vw_i;
  v.dotc(vw_r, vw_i, w);
  return cmplx(vw_r, vw_i);
}

//====================================================================
template <typename REALTYPE>
void xI(AField<REALTYPE,SIMD2>& v)
{
  v.xI();
}

//====================================================================
template <typename REALTYPE>
void conjg(AField<REALTYPE,SIMD2>& v)
{
  v.conjg();
}

//====================================================================
template <class INDEX, class FIELD>
void convert(INDEX& index, FIELD& v, const Field& w) 
{
  int Nin  = w.nin();
  int Nvol = w.nvol();
  int Nex  = w.nex();
  assert(v.check_size(Nin, Nvol, Nex));

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  for(int ex = 0; ex < Nex; ++ex){   
   for(int site = is; site < ns; ++site){   
    for(int in = 0; in < Nin; ++in){   
      int iw = in + Nin * (site + Nvol * ex);
      int iv = index.idx(in, Nin, site, ex);
      v.e(iv) = w.cmp(iw);
    }
   }
  }

}  

//====================================================================
template <class INDEX, class FIELD>
void convert_spinor(INDEX& index, FIELD& v, const Field& w) 
{
  int Nin  = w.nin();
  int Nvol = w.nvol();
  int Nex  = w.nex();
  assert(v.check_size(Nin, Nvol, Nex));

  int Nc = CommonParameters::Nc();
  int Nd = CommonParameters::Nd();
  assert(Nin == 2*Nc*Nd);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  for(int ex = 0; ex < Nex; ++ex){   
   for(int site = is; site < ns; ++site){   
    for(int id = 0; id < Nd; ++id){   
     for(int ic = 0; ic < Nc; ++ic){   
      int iwr =    2*(ic + Nc*id) + Nin * (site + Nvol * ex);
      int iwi = 1+ 2*(ic + Nc*id) + Nin * (site + Nvol * ex);
      int ivr = index.idx_SPr(ic, id, site, ex);
      int ivi = index.idx_SPi(ic, id, site, ex);
      v.e(ivr) = w.cmp(iwr);
      v.e(ivi) = w.cmp(iwi);
     }
    }
   }
  }

}  

//====================================================================
template <class INDEX, class FIELD>
void convert_gauge(INDEX& index, FIELD& v, const Field& w) 
{
  int Nin  = w.nin();
  int Nvol = w.nvol();
  int Nex  = w.nex();
  assert(v.check_size(Nin, Nvol, Nex));

  int Nc = CommonParameters::Nc();
  assert(Nin == 2*Nc*Nc);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  for(int ex = 0; ex < Nex; ++ex){   
   for(int site = is; site < ns; ++site){   
    for(int ic2 = 0; ic2 < Nc; ++ic2){   
     for(int ic1 = 0; ic1 < Nc; ++ic1){   
      int iwr =    2*(ic1 + Nc*ic2) + Nin * (site + Nvol * ex);
      int iwi = 1+ 2*(ic1 + Nc*ic2) + Nin * (site + Nvol * ex);
      int ivr = index.idx_Gr(ic1, ic2, site, ex);
      int ivi = index.idx_Gi(ic1, ic2, site, ex);
      v.e(ivr) = w.cmp(iwr);
      v.e(ivi) = w.cmp(iwi);
     }
    }
   }
  }

}  

//====================================================================
template <class INDEX2, class FIELD2, class INDEX1, class FIELD1>
void convert(INDEX2& index2, FIELD2& v2,
             const INDEX1& index1, const FIELD1& v1) 
{
  int Nin  = v1.nin();
  int Nvol = v1.nvol();
  int Nex  = v1.nex();
  assert(v2.check_size(Nin, Nvol, Nex));

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  for(int ex = 0; ex < Nex; ++ex){   
   for(int site = is; site < ns; ++site){   
    for(int in = 0; in < Nin; ++in){   
      int iv1 = index1.idx(in, Nin, site, ex);
      int iv2 = index2.idx(in, Nin, site, ex);
      v2.e(iv2) = v1.cmp(iv1);
    }
   }
  }

}  

//====================================================================
template <class INDEX2, class FIELD2, class INDEX1, class FIELD1>
void convert_h(INDEX2& index2, FIELD2& v2,
               const INDEX1& index1, const FIELD1& v1) 
{
  int Nin  = v1.nin();
  int Nvol = v1.nvol();
  int Nex  = v1.nex();
  assert(v2.check_size(Nin, Nvol, Nex));

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  for(int ex = 0; ex < Nex; ++ex){   
   for(int site = is; site < ns; ++site){   
    for(int in = 0; in < Nin; ++in){   
      int iv1 = index1.idxh(in, Nin, site, ex);
      int iv2 = index2.idxh(in, Nin, site, ex);
      v2.e(iv2) = v1.cmp(iv1);
    }
   }
  }

}  

//====================================================================
template <class INDEX, class FIELD>
void reverse(INDEX& index, Field& v, FIELD& w) 
{
  int Nin  = v.nin();
  int Nvol = v.nvol();
  int Nex  = v.nex();
  w.check_size(Nin, Nvol, Nex);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  for(int ex = 0; ex < Nex; ++ex){   
   for(int site = is; site < ns; ++site){   
    for(int in = 0; in < Nin; ++in){   
      int iv = in + Nin * (site + Nvol * ex);
      int iw = index.idx(in, Nin, site, ex);
      v.set(iv, double(w.cmp(iw)));
    }
   }
  }

}  

//====================================================================
template <class INDEX, class FIELD>
void reverse_spinor(INDEX& index, Field& v, FIELD& w) 
{
  int Nin  = v.nin();
  int Nvol = v.nvol();
  int Nex  = v.nex();
  w.check_size(Nin, Nvol, Nex);

  int Nc = CommonParameters::Nc();
  int Nd = CommonParameters::Nd();
  assert(Nin == 2*Nc*Nd);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  for(int ex = 0; ex < Nex; ++ex){   
   for(int site = is; site < ns; ++site){   
    for(int id = 0; id < Nd; ++id){   
     for(int ic = 0; ic < Nc; ++ic){   
      int ivr =    2*(ic + Nc*id) + Nin * (site + Nvol * ex);
      int ivi = 1+ 2*(ic + Nc*id) + Nin * (site + Nvol * ex);
      int iwr = index.idx_SPr(ic, id, site, ex);
      int iwi = index.idx_SPi(ic, id, site, ex);
      v.set(ivr, double(w.cmp(iwr)));
      v.set(ivi, double(w.cmp(iwi)));
     }
    }
   }
  }

}  

//====================================================================
template <class INDEX, class FIELD>
void reverse_gauge(INDEX& index, Field& v, FIELD& w) 
{
  int Nin  = v.nin();
  int Nvol = v.nvol();
  int Nex  = v.nex();
  w.check_size(Nin, Nvol, Nex);

  int Nc = CommonParameters::Nc();
  assert(Nin == 2*Nc*Nc);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  for(int ex = 0; ex < Nex; ++ex){   
   for(int site = is; site < ns; ++site){   
    for(int ic2 = 0; ic2 < Nc; ++ic2){   
     for(int ic1 = 0; ic1 < Nc; ++ic1){   
      int ivr =    2*(ic1 + Nc*ic2) + Nin * (site + Nvol * ex);
      int ivi = 1+ 2*(ic1 + Nc*ic2) + Nin * (site + Nvol * ex);
      int iwr = index.idx_Gr(ic1, ic2, site, ex);
      int iwi = index.idx_Gi(ic1, ic2, site, ex);
      v.set(ivr, double(w.cmp(iwr)));
      v.set(ivi, double(w.cmp(iwi)));
     }
    }
   }
  }

}  

//====================================================================
template <typename REALTYPE>
REALTYPE norm2(const AField<REALTYPE,SIMD2>& v)
{
  AField<REALTYPE,SIMD2> *vp = const_cast<AField<REALTYPE,SIMD2>*>(&v);
  return vp->norm2();
}

//============================================================END=====
#endif
