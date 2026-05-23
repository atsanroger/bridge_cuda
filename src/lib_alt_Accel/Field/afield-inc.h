/*!
        @file    afield-inc.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2024-11-05 17:41:16 #$
        @version $LastChangedRevision: 2619 $
*/

#ifndef ACCEL_AFIELD_INC_INCLUDED
#define ACCEL_AFIELD_INC_INCLUDED

#include <cstdlib> 

#include "complexTraits.h"
#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_Accel/inline/afield_th-inc.h" 
#include "lib_alt_Accel/Field/aindex_lex.h"
#include "lib_alt_Accel/Field/aindex_eo.h"

#include "lib_alt_Accel/BridgeACC/bridgeACC_AField.h"

//====================================================================
template <typename REALTYPE>
void update_host(AField<REALTYPE, ACCEL>& v)
{
  v.update_host();
}

//====================================================================
template <typename REALTYPE>
void update_device(AField<REALTYPE, ACCEL>& v)
{
  v.update_device();
}

//====================================================================
template <typename REALTYPE>
void copy(AField<REALTYPE, ACCEL>& v, const AField<REALTYPE, ACCEL> &w)
{
  v.copy(w);
}

//====================================================================
template <typename REALTYPE>
void copy(AField<REALTYPE, ACCEL>& v, const int ex,
          const AField<REALTYPE, ACCEL> &w, const int ex_w)
{
  v.copy(ex, w, ex_w);
}

//====================================================================
template <typename REALTYPE>
void axpy(AField<REALTYPE, ACCEL>& v, const REALTYPE a, const AField<REALTYPE, ACCEL> &w)
{
  v.axpy(a, w);
}

//====================================================================
template <typename REALTYPE>
void axpy(AField<REALTYPE, ACCEL>& v,
          const typename ComplexTraits<REALTYPE>::complex_t a,
          const AField<REALTYPE, ACCEL> &w)
{
  v.axpy(a, w);
}

//====================================================================
template <typename REALTYPE>
void axpy(AField<REALTYPE, ACCEL>& v, const int ex,
          const REALTYPE a, 
          const AField<REALTYPE, ACCEL> &w, const int ex_w)
{
  v.axpy(ex, a, w, ex_w);
}

//====================================================================
template <typename REALTYPE>
void aypx(const REALTYPE a, AField<REALTYPE, ACCEL>& v,
          const AField<REALTYPE, ACCEL>& w, const int mode = 0)
{
  v.aypx(a, w, mode);
}

//====================================================================
template <typename REALTYPE>
void aypx(const typename ComplexTraits<REALTYPE>::complex_t a,
          AField<REALTYPE, ACCEL>& v,
          const AField<REALTYPE, ACCEL>& w, const int mode = 0)
{
  v.aypx(a, w, mode);
}

//====================================================================
template <typename REALTYPE>
void scal(AField<REALTYPE, ACCEL>& v, REALTYPE a)
{
  v.scal(a);
}

//====================================================================
template <typename REALTYPE>
void scal(AField<REALTYPE, ACCEL>& v,
          const typename ComplexTraits<REALTYPE>::complex_t a)
{
  v.scal(a);
}

//====================================================================
template <typename REALTYPE>
REALTYPE dot(AField<REALTYPE, ACCEL>& v, AField<REALTYPE, ACCEL>& w)
{
  return v.dot(w);
}

//====================================================================
template <typename REALTYPE>
typename ComplexTraits<REALTYPE>::complex_t dotc(
                                  const AField<REALTYPE, ACCEL>& v,
                                  const AField<REALTYPE, ACCEL>& w)
{
  REALTYPE vw_r, vw_i;
  v.dotc(vw_r, vw_i, w);
  return cmplx(vw_r, vw_i);
}

//====================================================================
template <typename REALTYPE>
REALTYPE norm2(const AField<REALTYPE, ACCEL>& v)
{
  return v.norm2();
}

//====================================================================
template <typename REALTYPE>
void xI(AField<REALTYPE, ACCEL>& v)
{
  v.xI();
}

//====================================================================
template <typename REALTYPE>
void conjg(AField<REALTYPE, ACCEL>& v)
{
  v.conjg();
}

//====================================================================
template <class INDEX, class AFIELD>
void convert(INDEX& index, AFIELD& v, const Field& w) 
{
#pragma omp barrier

  int Nin  = w.nin();
  int Nvol = w.nvol();
  int Nex  = w.nex();
  assert(v.check_size(Nin, Nvol, Nex));

  vout.paranoiac("convert to ACCEL from Field start.\n");
  vout.paranoiac("  AFIELD = %s\n", AFIELD::class_name.c_str());
  vout.paranoiac("  Nin = %d Nvol = %d Nex = %d\n", Nin, Nvol, Nex);

  typedef typename AFIELD::real_t real_t;
  real_t* vp = v.ptr(0);

  int Nvol_pad = v.nvol_pad();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol_pad);

  for(int ex = 0; ex < Nex; ++ex){   
    for(int site = is; site < ns; ++site){   

     if(site < Nvol){
       for(int in = 0; in < Nin; ++in){   
         int iw = in + Nin * (site + Nvol * ex);
         int iv = index.idx(in, Nin, site, ex);
         vp[iv] = w.cmp(iw);
       }
     }else{
       for(int in = 0; in < Nin; ++in){   
         int iv = index.idx(in, Nin, site, ex);
         vp[iv] = real_t(0.0);
       }
     }

   }
  }

#pragma omp barrier

  if(ith == 0){
    size_t nv = Nin * Nvol_pad * Nex;
    BridgeACC::copy_to_device(vp, nv);
  }

  vout.paranoiac("convert to ACCEL from Field finished.\n");

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

  if(ith == 0){
    int nv = Nin * Nvol * Nex;
    BridgeACC::copy_from_device(v1p, nv);
  }

#pragma omp barrier

  typename AFIELD2::real_t *v2p = v2.ptr(0);

  for(int ex = 0; ex < Nex; ++ex){   
   for(int site = is; site < ns; ++site){   
    for(int in = 0; in < Nin; ++in){   
      int iv1 = index1.idx(in, Nin, site, ex);
      int iv2 = index2.idx(in, Nin, site, ex);
      v2p[iv2] = v1p[iv1];
      //v2p[iv2] = v1.cmp(iv1);
    }
   }
  }

#pragma omp barrier

  if(ith == 0){
    int nv = Nin * Nvol * Nex;
    BridgeACC::copy_to_device(v2p, nv);
  }

#pragma omp barrier

}

//====================================================================
template<>
void convert(AIndex_lex<double,ACCEL>& index2,
             AField<double,ACCEL>& v2,
             const AIndex_lex<float,ACCEL>& index1,
             const AField<float,ACCEL>& v1)
{
  int Nin  = v1.nin();
  int Nvol = v1.nvol();
  int Nex  = v1.nex();
  assert(v2.check_size(Nin, Nvol, Nex));

  double *v2p = v2.ptr(0);
  float  *v1p = const_cast<AField<float,ACCEL>* >(&v1)->ptr(0);

  int ith, nth;
  set_thread(ith, nth);

#pragma omp barrier

  if(ith == 0){
    int nvx = Nvol * Nex;
    BridgeACC::copy(v2p, 0, v1p, 0, Nin, nvx);
  }

#pragma omp barrier

}  

//====================================================================
template<>
void convert(AIndex_lex<float,ACCEL>& index2,
             AField<float,ACCEL>& v2,
             const AIndex_lex<double,ACCEL>& index1,
             const AField<double,ACCEL>& v1)
{
  int Nin  = v1.nin();
  int Nvol = v1.nvol();
  int Nex  = v1.nex();
  assert(v2.check_size(Nin, Nvol, Nex));

  float  *v2p = v2.ptr(0);
  double *v1p = const_cast<AField<double,ACCEL>* >(&v1)->ptr(0);

  int ith, nth;
  set_thread(ith, nth);

#pragma omp barrier

  if(ith == 0){
    int nvx = Nvol * Nex;
    BridgeACC::copy(v2p, 0, v1p, 0, Nin, nvx);
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

  if(ith == 0){
    int nv = Nin * Nvol * Nex;
    BridgeACC::copy_from_device(v1p, nv);
  }

#pragma omp barrier

  typename AFIELD2::real_t *v2p = v2.ptr(0);

  for(int ex = 0; ex < Nex; ++ex){
   for(int site = is; site < ns; ++site){
    for(int in = 0; in < Nin; ++in){
      int iv1 = index1.idxh(in, Nin, site, ex);
      int iv2 = index2.idxh(in, Nin, site, ex);
      v2p[iv2] = v1p[iv1];
      //v2p[iv2] = v1.cmp(iv1);
    }
   }
  }

#pragma omp barrier

  if(ith == 0){
    int nv = Nin * Nvol * Nex;
    BridgeACC::copy_to_device(v2p, nv);
  }

#pragma omp barrier

}


//====================================================================
template<>
void convert_h(AIndex_eo<double,ACCEL>& index2,
               AField<double,ACCEL>& v2,
               const AIndex_eo<float,ACCEL>& index1,
               const AField<float,ACCEL>& v1)
{
  int Nin  = v1.nin();
  int Nvol = v1.nvol();
  int Nex  = v1.nex();
  assert(v2.check_size(Nin, Nvol, Nex));

  double *v2p = v2.ptr(0);
  float  *v1p = const_cast<AField<float,ACCEL>* >(&v1)->ptr(0);

  int ith, nth;
  set_thread(ith, nth);

#pragma omp barrier

  if(ith == 0){
    int nvx = Nvol * Nex;
    BridgeACC::copy(v2p, 0, v1p, 0, Nin, nvx);
  }

#pragma omp barrier

}

//====================================================================
template<>
void convert_h(AIndex_eo<float,ACCEL>& index2,
               AField<float,ACCEL>& v2,
               const AIndex_eo<double,ACCEL>& index1,
               const AField<double,ACCEL>& v1)
{
  int Nin  = v1.nin();
  int Nvol = v1.nvol();
  int Nex  = v1.nex();
  assert(v2.check_size(Nin, Nvol, Nex));

  float  *v2p = v2.ptr(0);
  double *v1p = const_cast<AField<double,ACCEL>* >(&v1)->ptr(0);

  int ith, nth;
  set_thread(ith, nth);

#pragma omp barrier

  if(ith == 0){
    int nvx = Nvol * Nex;
    BridgeACC::copy(v2p, 0, v1p, 0, Nin, nvx);
  }

#pragma omp barrier

}

//====================================================================
template <class INDEX, class AFIELD>
void reverse(INDEX& index, Field& v, const AFIELD& w) 
{
#pragma omp barrier

  int Nin  = v.nin();
  int Nvol = v.nvol();
  int Nex  = v.nex();
  w.check_size(Nin, Nvol, Nex);

  vout.paranoiac("reverse to Field from ACCEL start.\n");

  int Nvol_pad = w.nvol_pad();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  typedef typename AFIELD::real_t real_t;
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);

  if(ith == 0){
    int nv = Nin * Nvol_pad * Nex;
    BridgeACC::copy_from_device(wp, nv);
  }

#pragma omp barrier

  for(int ex = 0; ex < Nex; ++ex){   
    for(int site = is; site < ns; ++site){   
      for(int in = 0; in < Nin; ++in){   
        int iv = in + Nin * (site + Nvol * ex);
        int iw = index.idx(in, Nin, site, ex);
        v.set(iv, double(wp[iw]));
        //v.set(iv, double(w.cmp(iw)));
      }
    }
  }

  vout.paranoiac("reverse to Field from ACCEL finished.\n");

#pragma omp barrier

}  

//====================================================================
template <class AFIELD>
void reverse(AIndex_lex<typename AFIELD::real_t,ACCEL>& index,
             Field& v, const AFIELD& w) 
{
#pragma omp barrier

  int Nin  = v.nin();
  int Nvol = v.nvol();
  int Nex  = v.nex();
  assert(w.check_size(Nin, Nvol, Nex));

  vout.paranoiac("reverse to Field from ACCEL start.\n");

  typedef typename AFIELD::real_t real_t;
  double *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);

  int Nvol_pad = w.nvol_pad();

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){
    for(int ex = 0; ex < Nex; ++ex){
      double *vp2 = &vp[Nin * Nvol * ex];
      real_t *wp2 = &wp[Nin * Nvol_pad * ex];
      BridgeACC::reverse(vp2, wp2, Nin, Nvol, Nvol_pad);
    }
  }

  vout.paranoiac("reverse to Field from ACCEL finished.\n");

#pragma omp barrier

}  

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
