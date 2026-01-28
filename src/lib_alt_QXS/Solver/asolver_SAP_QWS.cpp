/*!
        @file    asolver_SAP_QWS.cpp
        @brief   SAP solver (qws version)
        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate: 2023-03-20 10:52:44 +0900 (Mon, 20 Mar 2023) $
        @version $LastChangedRevision: 2499 $
*/

//====================================================================
#include "asolver_SAP_QWS.h"
#include "lib_alt_QXS/Field/afield-inc.h"
namespace {
#ifndef AFILED_HAS_SUB
  template<typename AFIELD>
  inline void sub(AFIELD& v, const AFIELD& w)
  {
    typedef typename AFIELD::real_t real_t;
    axpy(v, real_t(-1.0), w);
  }
#endif

#ifndef AFILED_HAS_ADD
  template<typename AFIELD>
  inline void add(AFIELD& v, const AFIELD& w)
  {
    typedef typename AFIELD::real_t real_t;
    axpy(v, real_t(1.0), w);
  }
#endif
}

template<typename AFIELD>
const std::string ASolver_SAP_QWS<AFIELD>::class_name = "ASolver_SAP_QWS";

#ifdef USE_QWSLIB
//====================================================================
template<typename AFIELD>
void ASolver_SAP_QWS<AFIELD>::init(void)
{
  ThreadManager::assert_single_thread(class_name);

  int nin  = m_fopr->field_nin();
  int nvol = m_fopr->field_nvol();
  int nex  = m_fopr->field_nex();

  m_x.reset(nin, nvol, nex);

#ifdef DEBUG
  m_r.reset(nin, nvol, nex);
  m_b.reset(nin, nvol, nex);
#endif

  { // working are for qws
    int  ret = 0;
    void *p;
    ret = posix_memalign(&p, CLS, sizeof(scs_t) * vols * 2);
    m_s = (scs_t *)p;
    ret = posix_memalign(&p, CLS, sizeof(scs_t) * vols * 2);
    m_q = (scs_t *)p;
    if (ret) {
      vout.crucial(m_vl, "%s: allocation failed\n", class_name.c_str());
      exit(EXIT_FAILURE);
    }
  }

  m_fopr->set_mode("D");

  m_nconv = -1;
}


//====================================================================
template<typename AFIELD>
void ASolver_SAP_QWS<AFIELD>::tidyup(void)
{
  // free working area for qws
  if (m_s) {
    free(m_s);
  }
  if (m_q) {
    free(m_q);
  }
  // ThreadManager::assert_single_thread(class_name);
  //delete m_sap_minres;
}


//====================================================================
template<typename AFIELD>
void ASolver_SAP_QWS<AFIELD>::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  int    Niter, Nrestart;
  double Stop_cond;

  int err = 0;
  err += params.fetch_int("maximum_number_of_iteration", Niter);
  err += params.fetch_int("maximum_number_of_restart", Nrestart);
  err += params.fetch_double("convergence_criterion_squared", Stop_cond);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  int Niter2 = Niter * Nrestart;
  set_parameters(Niter2, Stop_cond);
}


//====================================================================
template<typename AFIELD>
void ASolver_SAP_QWS<AFIELD>::set_parameters(const int Niter,
                                             const real_t Stop_cond)
{
  ThreadManager::assert_single_thread(class_name);

  m_Niter     = Niter;
  m_Stop_cond = Stop_cond;
  std::string prec = "double";
  if (sizeof(real_t) == 4) prec = "float";

  // N.B. Stop_cond is irrelevant
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Precision: %s\n", prec.c_str());
  vout.general(m_vl, "  Niter     = %d\n", m_Niter);
  vout.general(m_vl, "  nm        = %d\n", m_nm);
  vout.general(m_vl, "  Stop_cond is not used\n");

  // flop coung
  int    NPE      = CommonParameters::NPE();
  int    Nx       = CommonParameters::Nx();
  int    Ny       = CommonParameters::Ny();
  int    Nz       = CommonParameters::Nz();
  int    Nt       = CommonParameters::Nt();
  int    Nx2      = Nx / 2;
  int    Nvol2    = Nx * Ny * Nz * Nt / 2;
  double flop_clv = 576;
  double flop_sap_hops
    = 2 * ((Nx2 - 1) * Ny * Nz * Nt * 168 + Nx2 * (Ny - 1) * Nz * Nt * 168
           + Nx2 * Ny * (Nz - 1) * Nt * 168 + Nx2 * Ny * Nz * (Nt - 1) * 156);
  double flop_DEE_proc = (48 + flop_clv) * Nvol2 + flop_sap_hops;        // per domain
  double flop_DEO_proc                                                   // contains accum
    = (2 * Nvol2 - 2 * (Nx2 - 2) * (Ny - 2) * (Nz - 2) * (Nt - 2)) * 624 // clov
      + 2 * (Nx2 * Nz * Nt * 168 + Nx2 * Nz * Nt * 168 + Nx2 * Ny * Nt * 168 + Nx2 * Ny * Nz * 156);

  double flop_AEE_proc = m_nm * (flop_DEE_proc + 48 * Nvol2); // jinv_ddd
  double flop
    = m_Niter * (2 * flop_AEE_proc                            // jinv_ddd
                 + 2 * flop_DEE_proc                          // ddd_in
                 + 2 * flop_DEO_proc                          // ddd_out_pre, ddd_out_pos
                 + 2 * 2 * 12 * Nvol2);                       // accum_addsub
  flop  += 2 * flop_AEE_proc + flop_DEO_proc;                 // jinv_dd, ddd_out_pre, ddd_out_pos
  flop  += flop_clv * 2 * Nvol2;                              // clov_inv on the source
  m_flop = flop * NPE;
}


//====================================================================
template<typename AFIELD>
void ASolver_SAP_QWS<AFIELD>::solve(AFIELD& x, const AFIELD& b,
                                    int& Nconv, real_t& diff)
{
  // vout.paranoiac(m_vl, "asolver_SAP_QWS: solve, started in %s\n", __func__);
  vout.detailed(m_vl, "asolver_SAP_QWS: solve, started in %s\n", __func__);
  using real_t = typename AFIELD::real_t;

#ifdef USE_QWSLIB
  m_fopr->convert(x, b);

#ifdef DEBUG
  {
    double x2 = x.norm2();
    vout.general("hoge: %s:       after convert, |in|^2 = %15.7e\n", class_name.c_str(), x2);
  }
#endif
  int nsap = m_Niter;
  int nm   = m_nm;
  m_fopr->mult_clv_inv(x);
#pragma omp barrier
#ifdef DEBUG
  {
    double x2 = x.norm2();
    vout.general("hoge: %s:  after mult_clv_inv, |in|^2 = %15.7e\n", class_name.c_str(), x2);
  }
#endif

#pragma omp barrier
  prec_s_noprl_((scs_t *)m_x.ptr(0), (scs_t *)x.ptr(0), &nsap, &nm, m_s, m_q);

  diff = -1.0;
#ifdef DEBUG
  {
    double x2 = m_x.norm2();
    m_fopr->convert(m_b, b);

    m_fopr->mult(m_r, m_x);
    axpy(m_r, (real_t)-1.0, m_b);
    double r2 = m_r.norm2();
    diff = r2;
    vout.general("hoge: %s: before reverse, |res|^2 = %15.7e\n", class_name.c_str(), r2);
    vout.general("hoge: %s:                 |out|^2 = %15.7e\n", class_name.c_str(), x2);
  }
#endif

  m_fopr->reverse(x, m_x);

#ifdef DEBUG
  {
    double x2 = x.norm2();
    vout.general("hoge: %s:  after reverse, |out|^2 = %15.7e\n", class_name.c_str(), x2);
  }
#endif


  m_Nconv = m_Niter;

#pragma omp barrier
#else
  vout.crucial(m_vl, "%s: in solve(), USE_QWSLIB is not defined\n", class_name.c_str());
#endif
}


//====================================================================
template<typename AFIELD>
double ASolver_SAP_QWS<AFIELD>::flop_count()
{
  return m_flop;
}


template<>
const std::string ASolver_SAP_QWS<AField<float, QXS> >::class_name
  = "ASolver_SAP_QWS<AField<float,QXS> >";


template class ASolver_SAP_QWS<AField<float, QXS> >;

#endif // USE_QWSLIB
//============================================================END=====
