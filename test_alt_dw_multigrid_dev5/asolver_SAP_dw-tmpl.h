/*!
        @file    asolver_SAP_dw-tmpl.h
        @brief   SAP solver (trial version for Domainwall fermion)
        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate: 2024-05-07 13:24:28 +0900 (2024年05月07日 (火)) $
        @version $LastChangedRevision: 2602 $
*/

//====================================================================
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
const std::string ASolver_SAP_dw<AFIELD>::class_name = "ASolver_SAP_dw";

//====================================================================
template<typename AFIELD>
void ASolver_SAP_dw<AFIELD>::init(void)
{
  ThreadManager::assert_single_thread(class_name);

  int nin  = m_fopr->field_nin();
  int nvol = m_fopr->field_nvol();
  int nex  = m_fopr->field_nex();

  m_r.reset(nin, nvol, nex);
  m_p.reset(nin, nvol, nex);
  m_b.reset(nin, nvol, nex);
  //  m_sap_minres.reset(new ASolver_SAP_dw_MINRES<AFIELD>(m_fopr, m_block_index));
  //  m_sap_minres->set_parameters(m_min_res_iter, 0.0);


  m_nconv = -1;
}


//====================================================================
template<typename AFIELD>
void ASolver_SAP_dw<AFIELD>::tidyup(void)
{
  // ThreadManager::assert_single_thread(class_name);
  //delete m_sap_minres;
}


//====================================================================
template<typename AFIELD>
void ASolver_SAP_dw<AFIELD>::set_parameters(const Parameters& params)
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
void ASolver_SAP_dw<AFIELD>::set_parameters(const int Niter,
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
  //  vout.general(m_vl, "  Stop_cond = %16.8e\n", m_Stop_cond);
}


//====================================================================
template<typename AFIELD>
void ASolver_SAP_dw<AFIELD>::solve(AFIELD& x, const AFIELD& b,
                                   int& Nconv, real_t& diff)
{
  // vout.paranoiac(m_vl, "asolver_SAP: solve, started in %s\n", __func__);
  vout.detailed(m_vl, "asolver_SAP: solve, started in %s, m_Niter=%d\n",
                __func__, m_Niter);
  using real_t = typename AFIELD::real_t;

  x.set(real_t(0.0));
  //  m_fopr->set_mode("Ddag");
  //  m_fopr->mult(m_b, b);
  copy(m_b, b);
  copy(m_r, m_b);
  m_r.scal(real_t(-1.0));  // -|r> = D|x> - |b>
  int    nconv = -1;
  real_t diff0 = 0.0;

#ifdef DEBUG
  {
    real_t r2 = norm2(m_r);
    vout.general(m_vl, "initial, r2=%e\n", r2);
  }
#endif
  Nconv = m_Niter;
  for (int iter = 0; iter < m_Niter; ++iter) {
    int eo = 0;
    m_solver_in_block->solve(m_p, m_r, nconv, diff0, eo);

#ifdef DEBUG
    {
      real_t p2 = norm2(m_p);
      vout.general(m_vl, "   after solver_in_block: |m_p|^2 = %e  [nconv = %d ]\n",
                   p2, nconv);
    }
#endif
    sub(x, m_p);
    m_fopr->set_mode("D");
    m_fopr->mult(m_p, x);
    sub(m_p, m_b); // |p> =  -|r> = D|x> - |b>
#ifdef DEBUG
    {
      double r2 = norm2(m_p);
      double x2 = norm2(x);
      vout.general(m_vl, "iter = %d, p2 = %e, x2 = %e\n", iter, r2, x2);
    }
#endif

    eo = 1;
    m_solver_in_block->solve(m_r, m_p, nconv, diff0, eo);
#ifdef DEBUG
    {
      double r2 = norm2(m_r);
      vout.general(m_vl, "   after solver_in_block: |m_r|^2 = %e [ nconv = %d ]\n",
                   r2, nconv);
    }
#endif
    sub(x, m_r);
    m_fopr->set_mode("D");
    m_fopr->mult(m_r, x);
    sub(m_r, m_b);

#ifdef DEBUG
    {
      double r2 = norm2(m_r);
      double x2 = norm2(x);
      vout.general(m_vl, "iter = %d, r2 = %e, x2 = %e\n", iter, r2, x2);
    }
#endif
  } // iter loop

#ifdef DEBUG
  real_t r2 = norm2(m_r);
  diff = sqrt(r2);
#else
  diff = -1.0;
#endif

  m_Nconv = Nconv;

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
double ASolver_SAP_dw<AFIELD>::flop_count()
{
  int Nin  = m_fopr->field_nin();
  int Nvol = m_fopr->field_nvol();
  int Nex  = m_fopr->field_nex();
  int NPE  = CommonParameters::NPE();

  double flop_minres = m_solver_in_block->flop_count();
  double flop_mult   = m_fopr->flop_count();
  double flop_sub    = 2.0 * Nin * Nvol * NPE;
  double flop
    = m_Nconv * (2 * flop_minres + 2 * flop_mult + 2 * flop_sub);
  return flop;
}


//============================================================END=====
