/*!
        @file    asolver_SAP_dw2-tmpl.h
        @brief   SAP solver (trial version for Domainwall fermion)
        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate: 2024-03-17 22:33:12 +0900 (2024年03月17日 (日)) $
        @version $LastChangedRevision: 2586 $
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
const std::string ASolver_SAP_dw2<AFIELD>::class_name = "ASolver_SAP_dw2";

//====================================================================
template<typename AFIELD>
void ASolver_SAP_dw2<AFIELD>::init(void)
{
  ThreadManager::assert_single_thread(class_name);

  int nin  = m_fopr->field_nin();
  int nvol = m_fopr->field_nvol();
  int nex  = m_fopr->field_nex();

  m_r.reset(nin, nvol, nex);
  m_p.reset(nin, nvol, nex);
  m_b.reset(nin, nvol, nex);
  m_v.reset(nin, nvol, nex);
  m_x.reset(nin, nvol, nex);
  //  m_sap_minres.reset(new ASolver_SAP_dw2_MINRES<AFIELD>(m_fopr, m_block_index));
  //  m_sap_minres->set_parameters(m_min_res_iter, 0.0);

  ASolver_SAP_MINRES_dw2<AFIELD>
    *solver_in_block = dynamic_cast<ASolver_SAP_MINRES_dw2<AFIELD> *>(m_solver_in_block);
  if(solver_in_block == nullptr){
    vout.crucial(m_vl, "%s: m_solver_in_block must be ASolver_SAP_MINRES_dw2\n",
                 class_name.c_str());
    abort();
  }
  solver_in_block->set_PV_solver(m_solver_PV);
  m_fopr_PV=dynamic_cast<AFopr_dd<AFIELD> *>(m_solver_PV->get_fopr());
  if(m_fopr_PV == nullptr){
    vout.crucial(m_vl, "%s: m_solver_PV->get_fopr() = %p\n", class_name.c_str(), m_solver_PV->get_fopr());
    vout.crucial(m_vl, "%s: m_fopr_PV is NOT Afopr_dd\n", class_name.c_str());
    abort();
  }

  m_nconv = -1;
}


//====================================================================
template<typename AFIELD>
void ASolver_SAP_dw2<AFIELD>::tidyup(void)
{
  // ThreadManager::assert_single_thread(class_name);
  //delete m_sap_minres;
}


//====================================================================
template<typename AFIELD>
void ASolver_SAP_dw2<AFIELD>::set_parameters(const Parameters& params)
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
void ASolver_SAP_dw2<AFIELD>::set_parameters(const int Niter,
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
void ASolver_SAP_dw2<AFIELD>::solve(AFIELD& x, const AFIELD& b,
                                   int& Nconv, real_t& diff)
{
  // vout.paranoiac(m_vl, "asolver_SAP: solve, started in %s\n", __func__);
  vout.detailed(m_vl, "asolver_SAP: solve, started in %s, m_Niter=%d\n",
                __func__, m_Niter);
  using real_t = typename AFIELD::real_t;

  m_x.set(real_t(0.0));
  //m_fopr->set_mode("Ddag");
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
    sub(m_x, m_p);
    m_fopr->set_mode("D");
    m_solver_PV->solve(m_v, m_x, nconv, diff0);
    m_fopr->mult(m_p, m_v);
    sub(m_p, m_b); // |p> =  -|r> = D|x> - |b>
#ifdef DEBUG
    {
      double r2 = norm2(m_p);
      double x2 = norm2(m_x);
      vout.general(m_vl, "iter = %d, p2 = %e, x2 = %e\n", iter, r2, x2);
    }
#endif

    eo = 1;
    m_solver_in_block->solve(m_r, m_p, nconv, diff0, eo);
    sub(m_x, m_r);
    m_fopr->set_mode("D");
    m_solver_PV->solve(m_v, m_x, nconv, diff0);
    m_fopr->mult(m_r, m_v);

    sub(m_r, m_b);

#ifdef DEBUG
    {
      double r2 = norm2(m_r);
      double x2 = norm2(m_x);
      vout.general(m_vl, "iter = %d, r2 = %e, x2 = %e\n", iter, r2, x2);
    }
#endif
  } // iter loop

  copy(x, m_x);
#ifdef DEBUG
  {
    m_fopr->mult(m_v, x);
    m_solver_PV->solve(m_r, m_v, nconv, diff0);
    sub(m_r, b);
    real_t b2 = norm2(b);
    real_t r2 = norm2(m_r);
    vout.general(m_vl, "Nconv=%d,  r2 = %e  [relative: %e]\n",
                 Nconv, r2, r2/b2);
    diff = sqrt(r2/b2);
  }
#else
  diff = -1.0;
#endif

  m_Nconv = Nconv;

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
double ASolver_SAP_dw2<AFIELD>::flop_count()
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
