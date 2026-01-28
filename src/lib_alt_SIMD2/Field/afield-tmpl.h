/*!
      @file    afield-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
      @version $LastChangedRevision: 2160 $
*/

template<typename REALTYPE>
const std::string AField<REALTYPE,SIMD2>::class_name = "AField<REALTYPE,SIMD2>";

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD2>::init(const int nin, const int nvol, const int nex)
{
  m_nin  = nin;
  m_nvol = nvol;
  m_nex  = nex;
  m_nsize = m_nin * m_nvol * m_nex;
  if(m_nsize > 0){
    m_field = (real_t*)malloc(m_nsize*sizeof(real_t));
    vout.paranoiac("%s: data memory allocated\n", class_name.c_str());
  }else{
    m_field = NULL;
    vout.paranoiac("%s: null data size\n", class_name.c_str());
  }

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD2>::tidyup()
{
  if(m_field != NULL) free(m_field);
}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD2>::reset(const int nin, const int nvol, const int nex,
                             const element_type cmpl)
{
  if(check_size(nin,nvol,nex) && m_element_type == cmpl) return;

  int nsize_prev = m_nsize;
  m_nin  = nin;
  m_nvol = nvol;
  m_nex  = nex;
  m_element_type = cmpl;
  m_nsize = m_nin * m_nvol * m_nex;

  if(m_nsize != nsize_prev){
    if(m_field != NULL) free(m_field);
    if(m_nsize > 0){
      m_field = (real_t*)malloc(m_nsize*sizeof(real_t));
    }else{
      m_field = NULL;
    }
    vout.paranoiac("%s: data resized\n", class_name.c_str());
  }

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD2>::set(const REALTYPE a)
{
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, size());

  for (int i = is; i < ns; ++i) {
    m_field[i] = a;
  }

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD2>::copy(const AField<REALTYPE,SIMD2>& w)
{
  assert(check_size(w));

  for(int ex = 0; ex < nex(); ++ex){
    copy(ex, w, ex);
  }

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD2>::copy(const int ex,
                            const AField<REALTYPE,SIMD2>& w, const int ex_w)
{
  assert( nin() == w.nin());
  assert(nvol() == w.nvol());

  int size = m_nin * m_nvol;
  int sizev = size/VLEN;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, sizev);

  for (int i = is; i < ns; ++i) {
    for(int k = 0; k < VLEN; ++k){
      m_field[k + i*VLEN + size*ex] = w.m_field[k + i*VLEN + size*ex_w];
    }
  }

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD2>::axpy(const REALTYPE a, const AField<REALTYPE,SIMD2>& w)
{
  assert(check_size(w));

  for(int ex = 0; ex < nex(); ++ex){
    axpy(ex, a, w, ex);
  }

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD2>::axpy(const int ex, const REALTYPE a,
                            const AField<REALTYPE,SIMD2>& w, const int ex_w)
{
  assert( nin() == w.nin());
  assert(nvol() == w.nvol());
  assert(ex < nex());
  assert(ex_w < w.nex());

  int ith, nth, is, ns;
  int size2 = m_nin * m_nvol;
  int sizev = size2/VLEN;
  set_threadtask(ith, nth, is, ns, sizev);

  Vsimd_t vt, wt;

  for (int i = is; i < ns; ++i) {
    load_vec(&vt, &m_field[i*VLEN + size2*ex], 1);
    load_vec(&wt, &w.m_field[i*VLEN + size2*ex_w], 1);
    axpy_vec(&vt, a, &wt, 1);
    save_vec(&m_field[i*VLEN + size2*ex], &vt, 1);
  }

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD2>::axpy(const REALTYPE ar, const REALTYPE ai,
                                          const AField<REALTYPE,SIMD2>& w)
{
  assert(check_size(w));

  for(int ex = 0; ex < nex(); ++ex){
    axpy(ex, ar, ai, w, ex);
  }

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD2>::axpy(const int ex,
                            const REALTYPE ar, const REALTYPE ai,
                            const AField<REALTYPE,SIMD2>& w, const int ex_w)
{
  assert( nin() == w.nin());
  assert(nvol() == w.nvol());
  assert(ex < nex());
  assert(ex_w < w.nex());

  int size2 = m_nin * m_nvol;
  int sizev = size2/VLEN;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, sizev);

  Vsimd_t vt, wt;
  for (int i = is; i < ns; ++i) {
    load_vec(&vt, &m_field[i*VLEN + size2*ex], 1);
    load_vec(&wt, &w.m_field[i*VLEN + size2*ex_w], 1);
    axpy_vec(&vt, ar, ai, &wt, 1);
    save_vec(&m_field[i*VLEN + size2*ex], &vt, 1);
  }

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD2>::axpy(
            const typename ComplexTraits<REALTYPE>::complex_t a,
            const AField<REALTYPE,SIMD2>& w)
{
  REALTYPE ar = real(a);
  REALTYPE ai = imag(a);
  axpy(ar, ai, w);
}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD2>::axpy(
                const int ex,
                const typename ComplexTraits<REALTYPE>::complex_t a,
                const AField<REALTYPE,SIMD2>& w,
                const int ex_w)
{
  REALTYPE ar = real(a);
  REALTYPE ai = imag(a);
  axpy(ex, ar, ai, w, ex_w);
}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD2>::aypx(const REALTYPE a,
                                  const AField<REALTYPE,SIMD2>& w)
{
  assert(check_size(w));

  for(int ex = 0; ex < nex(); ++ex){
    aypx(ex, a, w, ex);
  }

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD2>::aypx(const int ex,
                            const REALTYPE a,
                            const AField<REALTYPE,SIMD2>& w, const int ex_w)
{
  assert( nin() == w.nin());
  assert(nvol() == w.nvol());
  assert(ex < nex());
  assert(ex_w < w.nex());

  int size2 = m_nin * m_nvol;
  int sizev = size2/VLEN;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, sizev);

  Vsimd_t vt, wt;
  for (int i = is; i < ns; ++i) {
    load_vec(&vt, &m_field[i*VLEN + size2*ex], 1);
    load_vec(&wt, &w.m_field[i*VLEN + size2*ex_w], 1);
    // aypx_vec(&vt, a, &wt, 1);
    aypx_vec(a, &vt, &wt, 1);
    save_vec(&m_field[i*VLEN + size2*ex], &vt, 1);
  }

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD2>::aypx(const int ex,
                     const REALTYPE ar, const REALTYPE ai,
                     const AField<REALTYPE,SIMD2>& w, const int ex_w)
{
  assert( nin() == w.nin());
  assert(nvol() == w.nvol());
  assert(ex < nex());
  assert(ex_w < w.nex());

  int size2 = m_nin * m_nvol;
  int sizev = size2/VLEN;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, sizev);

  Vsimd_t vt, wt;
  for (int i = is; i < ns; ++i) {
    load_vec(&vt, &m_field[i*VLEN + size2*ex], 1);
    load_vec(&wt, &w.m_field[i*VLEN + size2*ex_w], 1);
    aypx_vec(ar, ai, &vt, &wt, 1);
    save_vec(&m_field[i*VLEN + size2*ex], &vt, 1);
  }

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD2>::aypx(
                     const REALTYPE ar, const REALTYPE ai,
                     const AField<REALTYPE,SIMD2>& w)
{
  assert(check_size(w));

  for(int ex = 0; ex < nex(); ++ex){
    aypx(ex, ar, ai, w, ex);
  }

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD2>::aypx(const int ex,
                 const typename ComplexTraits<REALTYPE>::complex_t a,
                 const AField<REALTYPE,SIMD2>& w, const int ex_w)
{
  REALTYPE ar = real(a);
  REALTYPE ai = imag(a);

  aypx(ex, ar, ai, w, ex_w);

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD2>::aypx(
                 const typename ComplexTraits<REALTYPE>::complex_t a,
                 const AField<REALTYPE,SIMD2>& w)
{
  REALTYPE ar = real(a);
  REALTYPE ai = imag(a);

  for(int ex = 0; ex < nex(); ++ex){
    aypx(ar, ai, w);
  }

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD2>::scal(const REALTYPE a)
{
  for(int ex = 0; ex < nex(); ++ex){
    scal(a, ex);
  }
}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD2>::scal(const REALTYPE a, const int ex)
{
  int size2 = m_nin * m_nvol;
  int sizev = size2/VLEN;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, sizev);

  Vsimd_t vt;

  for(int i = is; i < ns; ++i){
    load_vec(&vt, &m_field[i*VLEN + size2*ex], 1);
    scal_vec(&vt, a, 1);
    save_vec(&m_field[i*VLEN + size2*ex], &vt, 1);
  }

}

//====================================================================
//qqqqq
template<typename REALTYPE>
void AField<REALTYPE,SIMD2>::scal(
                 const typename ComplexTraits<REALTYPE>::complex_t a)
{
  for(int ex = 0; ex < nex(); ++ex){
    scal(a, ex);
  }
}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD2>::scal(
               const typename ComplexTraits<REALTYPE>::complex_t a,
               const int ex)
{
  int size2 = m_nin * m_nvol;
  int sizev = size2/VLEN;

  REALTYPE ar = real(a);
  REALTYPE ai = imag(a);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, sizev);

  Vsimd_t vt;

  for(int i = is; i < ns; ++i){
    load_vec(&vt, &m_field[i*VLEN + size2*ex], 1);
    scal_vec(&vt, ar, ai, 1);
    save_vec(&m_field[i*VLEN + size2*ex], &vt, 1);
  }

}

//====================================================================
template<typename REALTYPE>
REALTYPE AField<REALTYPE,SIMD2>::dot(const AField<REALTYPE,SIMD2>& w)
{
  assert(check_size(w));

  int size2 = m_nin * m_nvol;
  int sizev = size2/VLEN;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, sizev);

  Vsimd_t av, vt, wt;
  real_t a2[VLEN];

  clear_vec(&av, 1);
  for(int ex = 0; ex < nex(); ++ex){
    for (int i = is; i < ns; ++i) {
      load_vec(&vt, &m_field[i*VLEN + size2*ex], 1);
      load_vec(&wt, &w.m_field[i*VLEN + size2*ex], 1);
      multadd_vec(av, vt, wt);
    }
  }
  save_vec(a2, &av, 1);

  real_t a = 0.0;
  for(int k = 0; k < VLEN; ++k){
    a += a2[k];
  }

  ThreadManager_OpenMP::reduce_sum_global(a, ith, nth);

  return a;

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD2>::dotc(REALTYPE& ar, REALTYPE& ai,
                                  const AField<REALTYPE,SIMD2>& w) const
{
  assert(check_size(w));

  /* the following implementation should be better if confirmed.

  int size2 = m_nin * m_nvol;
  int sizev = size2/VLEN;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, sizev);

  Vsimd_t av, vt, wt;
  real_t a2[VLEN];

  clear_vec(&av, 1);
  for(int ex = 0; ex < nex(); ++ex){
    for (int i = is; i < ns; ++i) {
      load_vec(&vt, &m_field[i*VLEN + size2*ex], 1);
      load_vec(&wt, &w.m_field[i*VLEN + size2*ex], 1);
      multadd_dn_vec(av, vt, wt);
    }
  }
  save_vec(a2, &av, 1);

  real_t ar = 0.0;
  real_t ai = 0.0;
  for(int k = 0; k < VLEN2; ++k){
    ar += a2[2*k];
    ai += a2[2*k+1];
  }
  */

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, w.size());

  ar = 0.0;
  ai = 0.0;

  // this is very naive implementation: to be improved !!
  int is2 = is/2;
  int ns2 = ns/2;
  for (int k = is2; k < ns2; ++k) {
    int kr = 2*k;
    int ki = 2*k + 1;
    ar += m_field[kr] * w.m_field[kr] + m_field[ki] * w.m_field[ki];
    ai += m_field[kr] * w.m_field[ki] - m_field[ki] * w.m_field[kr];
  }

  ThreadManager_OpenMP::reduce_sum_global(ar, ith, nth);
  ThreadManager_OpenMP::reduce_sum_global(ai, ith, nth);

}

//====================================================================
template<typename REALTYPE>
REALTYPE AField<REALTYPE,SIMD2>::norm2() const
{
  int size2 = m_nin * m_nvol;
  int sizev = size2/VLEN;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, sizev);

  Vsimd_t av, vt;
  real_t a2[VLEN];

  clear_vec(&av, 1);
  for(int ex = 0; ex < nex(); ++ex){
    for (int i = is; i < ns; ++i) {
      load_vec(&vt, &m_field[i*VLEN + size2*ex], 1);
      multadd_vec(av, vt, vt);
    }
  }
  save_vec(a2, &av, 1);

  real_t a = 0.0;
  for(int k = 0; k < VLEN; ++k){
    a += a2[k];
  }

  /*
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, size());

  REALTYPE a = 0.0;

  for (int k = is; k < ns; ++k) {
    a += m_field[k] * m_field[k];
  }
  */

  ThreadManager_OpenMP::reduce_sum_global(a, ith, nth);

  return a;

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD2>::xI()
{
  if(field_element_type() != Element_type::COMPLEX){
    vout.general("%s: xI is not relevant opearation\n",
                 class_name.c_str());
    return;
  }

  int Nin2 = nin()/2;
  int Nstv = nvol()/VLEN2;
  int Nex  = nex();

  real_t* vp = ptr(0);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstv);

  Vsimd_t vt;

  for(int ex = 0; ex < Nex; ++ex){

    for(int site = is; site < ns; ++site){
      int iv = VLEN * Nin2 * (site + Nstv * ex);

      for(int in2 = 0; in2 < Nin2; ++in2){
        int iv1 = iv + VLEN * in2;
        load_vec(&vt, &vp[iv1], 1);
        xI_vec<real_t>(&vt, 1);
        save_vec(&vp[iv1], &vt, 1);
      }

    }

  }

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD2>::conjg()
{
  if(field_element_type() == Element_type::REAL) return;

  int Nin2 = nin()/2;
  int Nstv = nvol()/VLEN2;
  int Nex  = nex();

  real_t* vp = ptr(0);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstv);

  Vsimd_t vt;

  for(int ex = 0; ex < Nex; ++ex){

    for(int site = is; site < ns; ++site){
      int iv = VLEN * Nin2 * (site + Nstv * ex);

      for(int in2 = 0; in2 < Nin2; ++in2){
        int iv1 = iv + VLEN * in2;
        load_vec(&vt, &vp[iv1], 1);
        conjg_vec(&vt, 1);
        save_vec(&vp[iv1], &vt, 1);
      }

    }

  }

}

//============================================================END=====
