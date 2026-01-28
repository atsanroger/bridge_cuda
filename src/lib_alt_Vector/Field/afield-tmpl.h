/*!
      @file    afield-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-10-02 13:50:31 #$
      @version $LastChangedRevision: 2543 $
*/

#include "lib_alt_Vector/inline/define_params.h"
#include "lib_alt_Vector/inline/define_index.h"

template<typename REALTYPE>
const std::string AField<REALTYPE,VECTOR>::class_name
                                      = "AField<REALTYPE,VECTOR>";

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,VECTOR>::init(const int nin,
                                    const int nvol,
                                    const int nex,
                                    const element_type cmpl)
{
  int ith, nth;
  set_thread(ith, nth);
  if(nth > 1){
    vout.crucial("%s: init must be called outisde parallel region\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  m_nin  = nin;
  m_nvol = nvol;
  m_nex  = nex;
  m_element_type = cmpl;                                                    
 
  m_nsize = m_nin * m_nvol * m_nex;

  //  m_nvol_pad = ceil_nwp(m_nvol);
  m_nvol_pad = m_nvol + VPAD;
  m_nsize_pad = m_nin * m_nvol_pad * m_nex;

  if(m_nsize > 0){
    int nodeid = Communicator::nodeid();
    size_t size_base = m_nsize_pad + VSHIFT * nodeid;
    m_field_base = (real_t*)malloc(size_base * sizeof(real_t));
    m_field = &m_field_base[VSHIFT * nodeid];
    vout.paranoiac("%s: data memory allocated\n", class_name.c_str());
  }else{
    m_field = NULL;
    vout.paranoiac("%s: null data size\n", class_name.c_str());
  }

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,VECTOR>::tidyup()
{
  int ith, nth;
  set_thread(ith, nth);
  if(nth > 1){
    vout.crucial("%s: tidyup must be called outisde parallel region\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  if(m_field != NULL) free(m_field_base);

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,VECTOR>::reset(const int nin,
                                     const int nvol,
                                     const int nex,
                                     const element_type cmpl)
{
  int ith, nth;
  set_thread(ith, nth);
  if(nth > 1){
    vout.crucial("%s: reset must be called outisde parallel region\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  if(check_size(nin,nvol,nex) && m_element_type == cmpl) return;

  std::size_t nsize_pad_prev = m_nsize_pad;

  m_nin  = nin;
  m_nvol = nvol;
  m_nex  = nex;
  m_element_type = cmpl;

  m_nsize = m_nin * m_nvol * m_nex;

  m_nvol_pad = m_nvol + VPAD;
  m_nsize_pad = m_nin * m_nvol_pad * m_nex;

  if(m_nsize_pad != nsize_pad_prev){
    if(m_field != NULL){
      free(m_field_base);
    }
    if(m_nsize > 0){
      int nodeid = Communicator::nodeid();
      size_t size_base = m_nsize_pad + VSHIFT * nodeid;
      m_field_base = (real_t*)malloc(size_base * sizeof(real_t));
      m_field = &m_field_base[VSHIFT * nodeid];
    }else{
      m_field = NULL;
    }
    vout.paranoiac("%s: data resized\n", class_name.c_str());
  }

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,VECTOR>::update_host()
{
}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,VECTOR>::update_device()
{
}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,VECTOR>::set(const REALTYPE a)
{

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_nvol_pad);

  for(int ex = 0; ex < m_nex; ++ex){
    for(int ist = is; ist < ns; ++ist){

      if(ist < m_nvol){
        for(int in = 0; in < m_nin; ++in){
          int idx2 = IDXV_BARE(m_nin, m_nvol_pad, in, ist, ex);
          m_field[idx2] = a;
        }
      }else{
        for(int in = 0; in < m_nin; ++in){
          int idx2 = IDXV_BARE(m_nin, m_nvol_pad, in, ist, ex);
          m_field[idx2] = 0.0;
        }
      }

    }
  }
#pragma omp barrier

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,VECTOR>::set_host(const REALTYPE a)
{
  set(a);
}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,VECTOR>::set(const int index, const REALTYPE a)
{
  m_field[index] = a;
}

//====================================================================
template<typename REALTYPE>
REALTYPE AField<REALTYPE,VECTOR>::cmp(const int index) const
{
  return m_field[index];
}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,VECTOR>::clear()
{
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_nsize_pad);

  for(int i = is; i < ns; ++is){
    m_field[i] = 0.0;
  }

#pragma omp barrier

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,VECTOR>::copy(const Field& w)
{
  assert(check_size(w));

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_nvol_pad);

  for(int ex = 0; ex < m_nex; ++ex){
    for(int ist = is; ist < ns; ++ist){

      if(ist < m_nvol){
        for(int in = 0; in < m_nin; ++in){
          int idx2 = IDXV_BARE(m_nin, m_nvol_pad, in, ist, ex);
          m_field[idx2] = REALTYPE(w.cmp(in, ist, ex));
        }
      }else{
        for(int in = 0; in < m_nin; ++in){
          int idx2 = IDXV_BARE(m_nin, m_nvol_pad, in, ist, ex);
          m_field[idx2] = 0.0;
        }
      }

    }
  }
#pragma omp barrier

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,VECTOR>::copy(const AField<REALTYPE,VECTOR>& w)
{
  assert(check_size(w));

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_nvol_pad);

  for(int ex = 0; ex < m_nex; ++ex){
    for(int ist = is; ist < ns; ++ist){

      if(ist < m_nvol){
        for(int in = 0; in < m_nin; ++in){
          int idx2 = IDXV_BARE(m_nin, m_nvol_pad, in, ist, ex);
          m_field[idx2] = w.cmp(idx2);
        }
      }else{
        for(int in = 0; in < m_nin; ++in){
          int idx2 = IDXV_BARE(m_nin, m_nvol_pad, in, ist, ex);
          m_field[idx2] = 0.0;
        }
      }

    }
  }
#pragma omp barrier

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,VECTOR>::copy(const int ex,
                                    const AField<REALTYPE,VECTOR>& w,
                                    const int ex_w)
{
  assert( nin() == w.nin());
  assert(nvol() == w.nvol());
  assert(ex < nex());
  assert(ex_w < w.nex());

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_nvol_pad);

  for(int ist = is; ist < ns; ++ist){

    if(ist < m_nvol){
      for(int in = 0; in < m_nin; ++in){
        int idx2 = IDXV_BARE(m_nin, m_nvol_pad, in, ist, ex);
        m_field[idx2] = w.cmp(IDXV(m_nin, m_nvol, in, ist, ex_w));
      }
    }else{
      for(int in = 0; in < m_nin; ++in){
        int idx2 = IDXV_BARE(m_nin, m_nvol_pad, in, ist, ex);
        m_field[idx2] = 0.0;
      }
    }

  }
#pragma omp barrier

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,VECTOR>::axpy(const REALTYPE a,
                                   const AField<REALTYPE,VECTOR>& w)
{
  assert(check_size(w));

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_nvol);

  for(int ex = 0; ex < m_nex; ++ex){
    for(int in = 0; in < m_nin; ++in){
      for(int ist = is; ist < ns; ++ist){
        int idx2 = IDXV(m_nin, m_nvol, in, ist, ex);
        m_field[idx2] += a * w.m_field[idx2];
      }
    }
  }

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,VECTOR>::axpy(const int ex, const REALTYPE a,
                                    const AField<REALTYPE,VECTOR>& w,
                                    const int ex_w)
{
  assert( nin() == w.nin());
  assert(nvol() == w.nvol());
  assert(ex < nex());
  assert(ex_w < w.nex());

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_nvol);

  for(int in = 0; in < m_nin; ++in){
    for(int ist = is; ist < ns; ++ist){
      int idx1 = IDXV(m_nin, m_nvol, in, ist, ex_w);
      int idx2 = IDXV(m_nin, m_nvol, in, ist, ex);
      m_field[idx2] += a * w.m_field[idx1];
    }
  }

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,VECTOR>::axpy(const REALTYPE ar, const REALTYPE ai,
                                    const AField<REALTYPE,VECTOR>& w)
{
  assert(check_size(w));

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_nvol);

  for(int ex = 0; ex < m_nex; ++ex){
    for(int in2 = 0; in2 < m_nin/2; ++in2){
      int inr = 2*in2;
      int ini = 2*in2+1;
      for(int ist = is; ist < ns; ++ist){
        int idxr = IDXV(m_nin, m_nvol, inr, ist, ex);
        int idxi = IDXV(m_nin, m_nvol, ini, ist, ex);
        REALTYPE wr = w.m_field[IDXV(m_nin, m_nvol, inr, ist, ex)];
        REALTYPE wi = w.m_field[IDXV(m_nin, m_nvol, ini, ist, ex)];
        m_field[idxr] += ar * wr - ai * wi;
        m_field[idxr] += ar * wi + ai * wr;
      }
    }
  }

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,VECTOR>::axpy(const int ex,
                                    const REALTYPE ar, const REALTYPE ai,
                                    const AField<REALTYPE,VECTOR>& w,
                                    const int ex_w)
{
  assert( nin() == w.nin());
  assert(nvol() == w.nvol());
  assert(ex < nex());
  assert(ex_w < w.nex());

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_nvol);

  for(int in2 = 0; in2 < m_nin/2; ++in2){
    int inr = 2*in2;
    int ini = 2*in2+1;
    for(int ist = is; ist < ns; ++ist){
      int idxr = IDXV(m_nin, m_nvol, inr, ist, ex);
      int idxi = IDXV(m_nin, m_nvol, ini, ist, ex);
      REALTYPE wr = w.m_field[IDXV(m_nin, m_nvol, inr, ist, ex_w)];
      REALTYPE wi = w.m_field[IDXV(m_nin, m_nvol, ini, ist, ex_w)];
      m_field[idxr] += ar * wr - ai * wi;
      m_field[idxr] += ar * wi + ai * wr;
    }
  }

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,VECTOR>::axpy(
                 const typename ComplexTraits<REALTYPE>::complex_t a,
                 const AField<REALTYPE,VECTOR>& w)
{
  assert(check_size(w));

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_nvol);

  REALTYPE ar = real(a);
  REALTYPE ai = imag(a);

  for(int ex = 0; ex < m_nex; ++ex){
    for(int in2 = 0; in2 < m_nin/2; ++in2){
      int inr = 2*in2;
      int ini = 2*in2+1;
      for(int ist = is; ist < ns; ++ist){
        int idxr = IDXV(m_nin, m_nvol, inr, ist, ex);
        int idxi = IDXV(m_nin, m_nvol, ini, ist, ex);
        REALTYPE wr = w.m_field[IDXV(m_nin, m_nvol, inr, ist, ex)];
        REALTYPE wi = w.m_field[IDXV(m_nin, m_nvol, ini, ist, ex)];
        m_field[idxr] += ar * wr - ai * wi;
        m_field[idxr] += ar * wi + ai * wr;
      }
    }
  }

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,VECTOR>::axpy(
                 const int ex,
                 const typename ComplexTraits<REALTYPE>::complex_t a,
                 const AField<REALTYPE,VECTOR>& w,
                 const int ex_w)
{
  assert( nin() == w.nin());
  assert(nvol() == w.nvol());
  assert(ex < nex());
  assert(ex_w < w.nex());

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_nvol);

  REALTYPE ar = real(a);
  REALTYPE ai = imag(a);

  for(int in2 = 0; in2 < m_nin/2; ++in2){
    int inr = 2*in2;
    int ini = 2*in2+1;
    for(int ist = is; ist < ns; ++ist){
      int idxr = IDXV(m_nin, m_nvol, inr, ist, ex);
      int idxi = IDXV(m_nin, m_nvol, ini, ist, ex);
      REALTYPE wr = w.m_field[IDXV(m_nin, m_nvol, inr, ist, ex_w)];
      REALTYPE wi = w.m_field[IDXV(m_nin, m_nvol, ini, ist, ex_w)];
      m_field[idxr] += ar * wr - ai * wi;
      m_field[idxr] += ar * wi + ai * wr;
    }
  }

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,VECTOR>::aypx(const REALTYPE a,
                                    const AField<REALTYPE,VECTOR>& w)
{
  assert(check_size(w));

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_nvol);

  for(int ex = 0; ex < m_nex; ++ex){
    for(int in = 0; in < m_nin; ++in){
      for(int ist = is; ist < ns; ++ist){
        int idx2 = IDXV(m_nin, m_nvol, in, ist, ex);
        m_field[idx2] = a * m_field[idx2] + w.m_field[idx2];
      }
    }
  }

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,VECTOR>::aypx(const int ex,
                                    const REALTYPE a,
                                    const AField<REALTYPE,VECTOR>& w,
                                    const int ex_w)
{
  assert( nin() == w.nin());
  assert(nvol() == w.nvol());
  assert(ex < nex());
  assert(ex_w < w.nex());

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_nvol);

  for(int ist = is; ist < ns; ++ist){
    for(int in = 0; in < m_nin; ++in){
      int idx1 = IDXV(m_nin, m_nvol, in, ist, ex_w);
      int idx2 = IDXV(m_nin, m_nvol, in, ist, ex);
      m_field[idx2] = a * m_field[idx2] + w.m_field[idx1];
    }
  }

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,VECTOR>::aypx(const REALTYPE ar, const REALTYPE ai,
                                    const AField<REALTYPE,VECTOR>& w)
{
  assert(check_size(w));

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_nvol);

  for(int ex = 0; ex < m_nex; ++ex){
    for(int in2 = 0; in2 < m_nin/2; ++in2){
      int inr = 2*in2;
      int ini = 2*in2+1;
      for(int ist = is; ist < ns; ++ist){
        int idxr = IDXV(m_nin, m_nvol, inr, ist, ex);
        int idxi = IDXV(m_nin, m_nvol, ini, ist, ex);
        REALTYPE vr = m_field[idxr];
        REALTYPE vi = m_field[idxi];
        REALTYPE wr = w.m_field[idxr];
        REALTYPE wi = w.m_field[idxi];
        m_field[idxr] = ar * vr - ai * vi + wr;
        m_field[idxr] = ar * vi + ai * vr + wi;
      }
    }
  }

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,VECTOR>::aypx(
                 const typename ComplexTraits<REALTYPE>::complex_t a,
                 const AField<REALTYPE,VECTOR>& w)
{
  assert(check_size(w));

  REALTYPE ar = real(a);
  REALTYPE ai = imag(a);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_nvol);

  for(int ex = 0; ex < m_nex; ++ex){
    for(int in2 = 0; in2 < m_nin/2; ++in2){
      int inr = 2*in2;
      int ini = 2*in2+1;
      for(int ist = is; ist < ns; ++ist){
        int idxr = IDXV(m_nin, m_nvol, inr, ist, ex);
        int idxi = IDXV(m_nin, m_nvol, ini, ist, ex);
        REALTYPE vr = m_field[idxr];
        REALTYPE vi = m_field[idxi];
        REALTYPE wr = w.m_field[idxr];
        REALTYPE wi = w.m_field[idxi];
        m_field[idxr] = ar * vr - ai * vi + wr;
        m_field[idxr] = ar * vi + ai * vr + wi;
      }
    }
  }

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,VECTOR>::aypx(
                 const int ex,
                 const typename ComplexTraits<REALTYPE>::complex_t a,
                 const AField<REALTYPE,VECTOR>& w,
                 const int ex_w)
{
  assert( nin() == w.nin());
  assert(nvol() == w.nvol());
  assert(ex < nex());
  assert(ex_w < w.nex());

  REALTYPE ar = real(a);
  REALTYPE ai = imag(a);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_nvol);

  for(int in2 = 0; in2 < m_nin/2; ++in2){
    int inr = 2*in2;
    int ini = 2*in2+1;
    for(int ist = is; ist < ns; ++ist){
      int idxr = IDXV(m_nin, m_nvol, inr, ist, ex);
      int idxi = IDXV(m_nin, m_nvol, ini, ist, ex);
      REALTYPE vr =  m_field[idxr];
      REALTYPE vi =  m_field[idxi];
      REALTYPE wr = w.m_field[IDXV(m_nin, m_nvol, inr, ist, ex_w)];
      REALTYPE wi = w.m_field[IDXV(m_nin, m_nvol, ini, ist, ex_w)];
      m_field[idxr] = ar * vr - ai * vi + wr;
      m_field[idxr] = ar * vi + ai * vr + wi;
    }
  }

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,VECTOR>::scal(const REALTYPE a)
{
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_nvol);

  for(int ex = 0; ex < m_nex; ++ex){
    for(int in = 0; in < m_nin; ++in){
      for(int ist = is; ist < ns; ++ist){
        int idx2 = IDXV(m_nin, m_nvol, in, ist, ex);
        m_field[idx2] *= a;
      }
    }
  }

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,VECTOR>::scal(const REALTYPE a, const int ex)
{
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_nvol);

  for(int in = 0; in < m_nin; ++in){
    for(int ist = is; ist < ns; ++ist){
      int idx2 = IDXV(m_nin, m_nvol, in, ist, ex);
      m_field[idx2] *= a;
    }
  }

}

//====================================================================
template<typename REALTYPE>
REALTYPE AField<REALTYPE,VECTOR>::dot(const AField<REALTYPE,VECTOR>& w)
{
  assert(check_size(w));

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_nvol);

  real_t a = 0.0;

  for(int ex = 0; ex < m_nex; ++ex){
    for(int in = 0; in < m_nin; ++in){
      for(int ist = is; ist < ns; ++ist){
        int idx2 = IDXV(m_nin, m_nvol, in, ist, ex);
        a += m_field[idx2] * w.m_field[idx2];
      }
    }
  }

#pragma omp barrier

  ThreadManager_OpenMP::reduce_sum_global(a, ith, nth);

  return a;

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,VECTOR>::dotc(REALTYPE& ar, REALTYPE& ai,
                                    const AField<REALTYPE,VECTOR>& w) const
{
  assert(check_size(w));

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_nvol);

  real_t ar2 = 0.0;
  real_t ai2 = 0.0;

  for(int ex = 0; ex < m_nex; ++ex){
    for(int in2 = 0; in2 < m_nin/2; ++in2){
      int inr = 2*in2;
      int ini = 2*in2+1;
      for(int ist = is; ist < ns; ++ist){
        int idxr = IDXV(m_nin, m_nvol, inr, ist, ex);
        int idxi = IDXV(m_nin, m_nvol, ini, ist, ex);
        REALTYPE vr = m_field[idxr];
        REALTYPE vi = m_field[idxi];
        REALTYPE wr = w.m_field[idxr];
        REALTYPE wi = w.m_field[idxi];
        ar2 += vr * wr + vi * wi;
        ai2 += vr * wi - vi * wr;
      }
    }
  }

#pragma omp barrier

  ThreadManager_OpenMP::reduce_sum_global(ar2, ith, nth);
  ThreadManager_OpenMP::reduce_sum_global(ai2, ith, nth);

  ar = ar2;
  ai = ai2;

}

//====================================================================
template<typename REALTYPE>
REALTYPE AField<REALTYPE,VECTOR>::norm2() const
{
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_nvol);

  real_t a = 0.0;

  for(int ex = 0; ex < m_nex; ++ex){
    for(int in = 0; in < m_nin; ++in){
      for(int ist = is; ist < ns; ++ist){
        int idx2 = IDXV(m_nin, m_nvol, in, ist, ex);
        REALTYPE v =  m_field[idx2] ;
        a += v * v;
      }
    }
  }

#pragma omp barrier

  ThreadManager_OpenMP::reduce_sum_global(a, ith, nth);

  return a;

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,VECTOR>::xI()
{
  if(field_element_type() != Element_type::COMPLEX){
    vout.general("%s: xI is not relevant opearation\n",
                 class_name.c_str());
    return;
  }

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_nvol);

  for(int ex = 0; ex < m_nex; ++ex){
    for(int in2 = 0; in2 < m_nin/2; ++in2){
      int inr = 2*in2;
      int ini = 2*in2+1;
      for(int ist = is; ist < ns; ++ist){
        int idxr = IDXV(m_nin, m_nvol, inr, ist, ex);
        int idxi = IDXV(m_nin, m_nvol, ini, ist, ex);
        REALTYPE vr = m_field[idxr];
        REALTYPE vi = m_field[idxi];
        m_field[idxr] = - vi;
        m_field[idxi] =   vr;
      }
    }
  }

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,VECTOR>::conjg()
{
  if(field_element_type() == Element_type::REAL) return;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_nvol);

  for(int ex = 0; ex < m_nex; ++ex){
    for(int in2 = 0; in2 < m_nin/2; ++in2){
      int ini = 2*in2+1;
      for(int ist = is; ist < ns; ++ist){
        int idxi = IDXV(m_nin, m_nvol, ini, ist, ex);
        m_field[idxi] = - m_field[idxi];
      }
    }
  }

}

//============================================================END=====
