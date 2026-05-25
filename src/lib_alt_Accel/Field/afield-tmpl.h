/*!
      @file    afield-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2024-03-08 16:00:51 #$
      @version $LastChangedRevision: 2574 $
*/

#include "lib_alt_Accel/inline/define_params.h"
#include "lib_alt_Accel/inline/define_index.h"

#ifdef USE_ACCEL_CUDA
#include <cuda_runtime.h>
#endif

template<typename REALTYPE>
const std::string AField<REALTYPE,ACCEL>::class_name
                                     = "AField<REALTYPE,ACCEL>";

template<typename REALTYPE>
int AField<REALTYPE,ACCEL>::m_num_instance = 0;

template<typename REALTYPE>
REALTYPE* AField<REALTYPE,ACCEL>::m_red1 = 0;

template<typename REALTYPE>
REALTYPE* AField<REALTYPE,ACCEL>::m_red2 = 0;

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,ACCEL>::init(const int nin,
                                    const int nvol,
                                    const int nex,
                                    const element_type cmpl)
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();

  m_nin  = nin;
  m_nvol = nvol;
  m_nex  = nex;
  m_element_type = cmpl;                                                    
 
  m_nsize = m_nin * m_nvol * m_nex;

  m_nvol_pad = ceil_nwp(m_nvol);
  m_nsize_pad = m_nin * m_nvol_pad * m_nex;

  if(m_nsize > 0){
    m_field = (real_t*)malloc(m_nsize_pad*sizeof(real_t));
    BridgeACC::afield_init(m_field, m_nsize_pad);
    vout.paranoiac("%s: data memory allocated\n", class_name.c_str());
  }else{
    m_field = NULL;
    vout.paranoiac("%s: null data size\n", class_name.c_str());
  }

#ifndef USE_ACCEL_OPENACC
  // setup fields for resuction
  if(m_num_instance == 0){
    int Nvol_pad = ceil_nwp(CommonParameters::Nvol());

    m_red1 = (real_t*)malloc(Nvol_pad * sizeof(real_t));
    BridgeACC::afield_init(m_red1, Nvol_pad);
    
    m_red2 = (real_t*)malloc(Nvol_pad * sizeof(real_t));
    BridgeACC::afield_init(m_red2, Nvol_pad);
  }
#endif

  ++m_num_instance;

  vout.paranoiac("%s: construction finished\n", class_name.c_str());

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,ACCEL>::tidyup()
{
  ThreadManager::assert_single_thread(class_name);

  if(m_field != NULL){
    BridgeACC::afield_tidyup(m_field, m_nsize_pad);
    free(m_field);
  }

  --m_num_instance;

#ifndef USE_ACCEL_OPENACC
  // tidyup fields for reduction
  if(m_num_instance == 0){
    int Nvol_pad = ceil_nwp(CommonParameters::Nvol());
    BridgeACC::exit_data_delete(m_red1, Nvol_pad);
    free(m_red1);
    BridgeACC::exit_data_delete(m_red2, Nvol_pad);
    free(m_red2);
  }
#endif
}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,ACCEL>::reset(const int nin,
                                     const int nvol,
                                     const int nex,
                                     const element_type cmpl)
{
  ThreadManager::assert_single_thread(class_name);

  vout.paranoiac("%s: reset called\n", class_name.c_str());

  if(check_size(nin,nvol,nex) && m_element_type == cmpl) return;

  std::size_t nsize_pad_prev = m_nsize_pad;

  m_nin  = nin;
  m_nvol = nvol;
  m_nex  = nex;
  m_element_type = cmpl;

  m_nsize = m_nin * m_nvol * m_nex;

  m_nvol_pad = ceil_nwp(m_nvol);
  m_nsize_pad = m_nin * m_nvol_pad * m_nex;

  if(m_nsize_pad != nsize_pad_prev){
    if(m_field != NULL){
      BridgeACC::afield_tidyup(m_field, nsize_pad_prev);
      free(m_field);
    }
    if(m_nsize > 0){
      m_field = (real_t*)malloc(m_nsize_pad*sizeof(real_t));
      BridgeACC::afield_init(m_field, m_nsize_pad);
    }else{
      m_field = NULL;
    }
    vout.paranoiac("%s: data resized\n", class_name.c_str());
  }

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,ACCEL>::update_host() const
{
#pragma omp barrier

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){
    BridgeACC::copy_from_device(m_field, m_nsize_pad);
  }

#pragma omp barrier
}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,ACCEL>::update_device()
{
#pragma omp barrier

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){
    BridgeACC::copy_to_device(m_field, m_nsize_pad);
  }

#pragma omp barrier
}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,ACCEL>::set(const REALTYPE a)
{
  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){
    int nsize2 = m_nsize_pad/m_nin;
    BridgeACC::afield_set(m_field, a, m_nin, nsize2);
  }
#pragma omp barrier

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,ACCEL>::set(const int index, const REALTYPE a)
{
  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){
    m_field[index] = a;
    //BridgeACC::copy_to_device(&m_field[index], 1);
    BridgeACC::copy_to_device(m_field, index, 1);
  }
#pragma omp barrier

}

//====================================================================
template<typename REALTYPE>
REALTYPE AField<REALTYPE,ACCEL>::cmp(const int index) const
{
  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){
    //BridgeACC::copy_from_device(&m_field[index], 1);
    BridgeACC::copy_from_device(m_field, index, 1);
  }
#pragma omp barrier

  return m_field[index];

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,ACCEL>::set_host(const REALTYPE a)
{
#pragma omp barrier

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_nsize_pad);

  for (int i = is; i < ns; ++i) {
    m_field[i] = a;
  }

#pragma omp barrier

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,ACCEL>::copy(const Field& w)
{
  assert(check_size(w));

  int size2 = m_nin * m_nvol_pad;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_nvol_pad);

  for(int ex = 0; ex < m_nex; ++ex){
   for(int ist = is; ist < ns; ++ist){

     if(ist < m_nvol){
       for(int in = 0; in < m_nin; ++in){
         int idx2 = IDX2(m_nin, in, ist) + size2 * ex;
         m_field[idx2] = REALTYPE(w.cmp(in, ist, ex));
       }
     }else{
       for(int in = 0; in < m_nin; ++in){
         int idx2 = IDX2(m_nin, in, ist) + size2 * ex;
         m_field[idx2] = 0.0;
       }
     }

   }
  }
#pragma omp barrier

  if(ith == 0){
    BridgeACC::copy_to_device(m_field, m_nsize_pad);
  }

#pragma omp barrier

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,ACCEL>::copy(const AField<REALTYPE,ACCEL>& w)
{
  assert(check_size(w));

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){
    if (m_nin == w.nin()) {
      int nvex = m_nsize_pad/m_nin;
      BridgeACC::copy(m_field, w.m_field, m_nin, nvex);
    } else if (m_nin == 2 * w.nin()) {
      BridgeACC::copy_to_qdw(m_field, w.m_field, m_nvol_pad, m_nin / 4);
    } else if (m_nin * 2 == w.nin()) {
      BridgeACC::copy_from_qdw(m_field, w.m_field, m_nvol_pad, w.m_nin / 4);
    } else {
      assert(false);
    }
  }
#pragma omp barrier

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,ACCEL>::copy(const int ex,
                                    const AField<REALTYPE,ACCEL>& w,
                                    const int ex_w)
{
  assert( nin() == w.nin());
  assert(nvol() == w.nvol());
  assert(ex < nex());
  assert(ex_w < w.nex());

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){
    int size2 = m_nin * m_nvol_pad;
    BridgeACC::copy(m_field, size2 * ex, w.m_field, size2 * ex_w,
                    m_nin, m_nvol_pad);
  }
#pragma omp barrier

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,ACCEL>::axpy(const REALTYPE a,
                                    const AField<REALTYPE,ACCEL>& w,
                                    const int mode)
{
  assert(check_size(w));

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){
    int nvex = m_nsize_pad/m_nin;
    BridgeACC::axpy(m_field, 0, a, w.m_field, 0, m_nin, nvex, mode);
  }
#pragma omp barrier

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,ACCEL>::axpy(const int ex, const REALTYPE a,
                                    const AField<REALTYPE,ACCEL>& w,
                                    const int ex_w,
                                    const int mode)
{
  assert( nin() == w.nin());
  assert(nvol() == w.nvol());
  assert(ex < nex());
  assert(ex_w < w.nex());

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){
    int size2 = m_nin * m_nvol_pad;
    BridgeACC::axpy(m_field, size2 * ex, a, w.m_field, size2 * ex_w,
                    m_nin, m_nvol_pad, mode);
  }
#pragma omp barrier

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,ACCEL>::axpy(const REALTYPE ar, const REALTYPE ai,
                                    const AField<REALTYPE,ACCEL>& w)
{
  assert(check_size(w));

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){
    int nvex = m_nsize_pad/m_nin;
    BridgeACC::axpy(m_field, 0, ar, ai, w.m_field, 0, m_nin, nvex);
  }
#pragma omp barrier

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,ACCEL>::axpy(const int ex,
                                    const REALTYPE ar, const REALTYPE ai,
                                    const AField<REALTYPE,ACCEL>& w,
                                    const int ex_w)
{
  assert( nin() == w.nin());
  assert(nvol() == w.nvol());
  assert(ex < nex());
  assert(ex_w < w.nex());

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){
    int size2 = m_nin * m_nvol_pad;
    BridgeACC::axpy(m_field, size2 * ex, ar, ai, w.m_field, size2 * ex_w,
                    m_nin, m_nvol_pad);
  }
#pragma omp barrier

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,ACCEL>::axpy(
                 const typename ComplexTraits<REALTYPE>::complex_t a,
                 const AField<REALTYPE,ACCEL>& w)
{
  assert(check_size(w));

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){
    REALTYPE ar = real(a);
    REALTYPE ai = imag(a);
    int nvex = m_nsize_pad/m_nin;
    BridgeACC::axpy(m_field, 0, ar, ai, w.m_field, 0, m_nin, nvex);
  }
#pragma omp barrier

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,ACCEL>::axpy(
                 const int ex,
                 const typename ComplexTraits<REALTYPE>::complex_t a,
                 const AField<REALTYPE,ACCEL>& w,
                 const int ex_w)
{
  assert( nin() == w.nin());
  assert(nvol() == w.nvol());
  assert(ex < nex());
  assert(ex_w < w.nex());

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){
    int size2 = m_nin * m_nvol_pad;
    REALTYPE ar = real(a);
    REALTYPE ai = imag(a);
    BridgeACC::axpy(m_field, size2 * ex, ar, ai, w.m_field, size2 * ex_w,
                    m_nin, m_nvol_pad);
  }
#pragma omp barrier

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,ACCEL>::aypx(const REALTYPE a,
                                  const AField<REALTYPE,ACCEL>& w, const int mode)
{
  assert(check_size(w));

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){
    int nvex = m_nsize_pad/m_nin;
    BridgeACC::aypx(a, m_field, 0, w.m_field, 0, m_nin, nvex, mode);
  }
#pragma omp barrier

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,ACCEL>::aypx(const int ex,
                                  const REALTYPE a,
                                  const AField<REALTYPE,ACCEL>& w,
                                  const int ex_w, const int mode)
{
  assert( nin() == w.nin());
  assert(nvol() == w.nvol());
  assert(ex < nex());
  assert(ex_w < w.nex());

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){
    int size2 = m_nin * m_nvol_pad;
    BridgeACC::aypx(a, m_field, size2 * ex, w.m_field, size2 * ex_w,
                    m_nin, m_nvol_pad, mode);
  }
#pragma omp barrier

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,ACCEL>::aypx(const REALTYPE ar, const REALTYPE ai,
                                    const AField<REALTYPE,ACCEL>& w, const int mode)
{
  assert(check_size(w));

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){
    int nvex = m_nsize_pad/m_nin;
    BridgeACC::aypx(ar, ai, m_field, 0, w.m_field, 0, m_nin, nvex, mode);
  }
#pragma omp barrier

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,ACCEL>::aypx(
                 const typename ComplexTraits<REALTYPE>::complex_t a,
                 const AField<REALTYPE,ACCEL>& w, const int mode)
{
  assert(check_size(w));

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){
    REALTYPE ar = real(a);
    REALTYPE ai = imag(a);
    int nvex = m_nsize_pad/m_nin;
    BridgeACC::aypx(ar, ai, m_field, 0, w.m_field, 0, m_nin, nvex, mode);
  }
#pragma omp barrier

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,ACCEL>::aypx(
                 const int ex,
                 const typename ComplexTraits<REALTYPE>::complex_t a,
                 const AField<REALTYPE,ACCEL>& w,
                 const int ex_w, const int mode)
{
  assert( nin() == w.nin());
  assert(nvol() == w.nvol());
  assert(ex < nex());
  assert(ex_w < w.nex());

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){
    int size2 = m_nin * m_nvol_pad;
    REALTYPE ar = real(a);
    REALTYPE ai = imag(a);
    BridgeACC::aypx(ar, ai, m_field, size2 * ex, w.m_field, size2 * ex_w,
                    m_nin, m_nvol_pad, mode);
  }

#pragma omp barrier

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,ACCEL>::scal(const REALTYPE a)
{
  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){
    int nvex = m_nsize_pad/m_nin;
    BridgeACC::scal(m_field, 0, a, m_nin, nvex);
  }
#pragma omp barrier

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,ACCEL>::scal(const REALTYPE a, const int ex)
{
  int size2 = m_nin * m_nvol_pad;

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){
    BridgeACC::scal(m_field, size2 * ex, a, m_nin, m_nvol_pad);
  }

#pragma omp barrier

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,ACCEL>::scal(
                const typename ComplexTraits<REALTYPE>::complex_t a)
{
  REALTYPE ar = real(a);
  REALTYPE ai = imag(a);

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){
    int nvex = m_nsize_pad/m_nin;
    BridgeACC::scal(m_field, 0, ar, ai, m_nin, nvex);
  }

#pragma omp barrier

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,ACCEL>::scal(
                const typename ComplexTraits<REALTYPE>::complex_t a,
                const int ex)
{
  int size2 = m_nin * m_nvol_pad;

  REALTYPE ar = real(a);
  REALTYPE ai = imag(a);

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){
    BridgeACC::scal(m_field, size2 * ex, ar, ai, m_nin, m_nvol_pad);
  }
#pragma omp barrier

}

//====================================================================
template<typename REALTYPE>
REALTYPE AField<REALTYPE,ACCEL>::dot(const AField<REALTYPE,ACCEL>& w,
                                       const int mode)
{
  assert(check_size(w));

  int ith, nth;
  set_thread(ith, nth);

  real_t a = 0.0;

  if(ith == 0){
#ifdef USE_ACCEL_OPENACC
    int nvex = m_nsize_pad/m_nin;
    a = BridgeACC::dot(m_field, w.m_field, m_nin, nvex);
#else
    a = BridgeACC::dot(m_field, w.m_field,
                       m_red1, m_nin, m_nvol_pad, m_nex, mode);
#endif
  }
#pragma omp barrier

  ThreadManager::reduce_sum_global(a, ith, nth);

  return a;

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,ACCEL>::dotc(REALTYPE& ar, REALTYPE& ai,
                                  const AField<REALTYPE,ACCEL>& w,
                                  const int mode) const
{
  assert(check_size(w));

  int ith, nth;
  set_thread(ith, nth);

  real_t ar2 = 0.0;
  real_t ai2 = 0.0;

  if(ith == 0){
#ifdef USE_ACCEL_OPENACC
    int nvex = m_nvol_pad * m_nex;
    BridgeACC::dotc(&ar2, &ai2, m_field, w.m_field, m_nin, nvex);
#else
    BridgeACC::dotc(&ar2, &ai2, m_field, w.m_field,
                    m_red1, m_red2, m_nin, m_nvol_pad, m_nex, mode);
#endif
  }

#pragma omp barrier

  real_t sum[2] = { ar2, ai2 };
  ThreadManager::reduce_sum_global(sum, 2, ith, nth);
  ar = sum[0];
  ai = sum[1];

}

//====================================================================
template<typename REALTYPE>
REALTYPE AField<REALTYPE,ACCEL>::norm2(const int mode) const
{
  int ith, nth;
  set_thread(ith, nth);

  real_t a = 0.0;

  if(ith == 0){
#ifdef USE_ACCEL_OPENACC
    int nvex = m_nvol_pad * m_nex;
    a = BridgeACC::norm2(m_field, m_nin, nvex);
#else
    a = BridgeACC::norm2(m_field, m_red1, m_nin, m_nvol_pad, m_nex, mode);
#endif
  }
#pragma omp barrier

  ThreadManager::reduce_sum_global(a, ith, nth);

  return a;

}

//====================================================================
template<typename REALTYPE>
REALTYPE AField<REALTYPE,ACCEL>::norm2_host() const
{
  int ith, nth;
  set_thread(ith, nth);

  real_t a = 0.0;

  if(ith == 0){
    int nvex = m_nvol_pad * m_nex;
    for(int ivex = 0; ivex < nvex; ++ivex){
      real_t at = 0.0;
      for(int in = 0; in < m_nin; ++in){
        real_t ft = m_field[IDX2(m_nin, in, ivex)];
        at += ft * ft;
      }
      a += at;
    }
  }
#pragma omp barrier

  ThreadManager::reduce_sum_global(a, ith, nth);

  return a;

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,ACCEL>::normalize(const int mode)
{
  if (!mode) return;

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){
    int nvol_pad_ex = m_nvol_pad * m_nex;
    BridgeACC::normalize(m_field, m_nin, nvol_pad_ex, mode);
  }
#pragma omp barrier

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,ACCEL>::xI()
{
  if(field_element_type() != Element_type::COMPLEX){
    vout.general("%s: xI is not relevant opearation\n",
                 class_name.c_str());
    return;
  }

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){
    int nvex = m_nvol_pad * m_nex;
    BridgeACC::xI(m_field, m_nin, nvex);
  }
#pragma omp barrier

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,ACCEL>::conjg()
{
  if(field_element_type() == Element_type::REAL) return;

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){
    int nvex = m_nvol_pad * m_nex;
    BridgeACC::conjg(m_field, m_nin, nvex);
  }
#pragma omp barrier

}

//============================================================END=====
