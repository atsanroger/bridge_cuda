template<typename AFIELD>
const std::string ASolver_CG<AFIELD>::class_name = "ASolver_CG";
//====================================================================
template<typename AFIELD>
void ASolver_CG<AFIELD>::init(void)
{
  ThreadManager::assert_single_thread(class_name);

  int nin  = m_fopr->field_nin();
  int nvol = m_fopr->field_nvol();
  int nex  = m_fopr->field_nex();

  m_x.reset(nin, nvol, nex);
  m_r.reset(nin, nvol, nex);
  m_p.reset(nin, nvol, nex);
  m_s.reset(nin, nvol, nex);

  //  m_vl = Bridge::DETAILED;

  m_nconv = -1;
}


//====================================================================
template<typename AFIELD>
void ASolver_CG<AFIELD>::tidyup(void)
{
  // ThreadManager::assert_single_thread(class_name);
  // nothing is to be deleted.
}


//====================================================================
template<typename AFIELD>
void ASolver_CG<AFIELD>::set_parameters(const Parameters& params)
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

  InitialGuess init_guess_mode = InitialGuess::RHS;
  if (params.find_string("initial_guess_mode")) {
    const string initial_guess_mode = params.get_string("initial_guess_mode");
    vout.detailed(m_vl, "  initila_guess_mode %s\n", initial_guess_mode.c_str());
    if (initial_guess_mode == "RHS") {
      init_guess_mode = InitialGuess::RHS;
    } else if (initial_guess_mode == "GIVEN") {
      init_guess_mode = InitialGuess::GIVEN;
    } else if (initial_guess_mode == "ZERO") {
      init_guess_mode = InitialGuess::ZERO;
    } else {
      vout.crucial(m_vl, "Error at %s: unknown initial guess mode, %s\n", class_name.c_str(), initial_guess_mode.c_str());
      exit(EXIT_FAILURE);
    }
  }

  int Niter2 = Niter * Nrestart;
  set_parameters(Niter2, Stop_cond, init_guess_mode);

  m_use_QDW = false;
  m_use_QTW = false;
  params.fetch_bool("use_QDW", m_use_QDW);
  params.fetch_bool("use_QTW", m_use_QTW);

  if (m_use_QDW) {
    m_use_qdw_int = 1;
  } else if (m_use_QTW) {
    m_use_qdw_int = 2;
  } else {
    m_use_qdw_int = 0;
  }

  if (m_use_QDW) {
    int nin = 2 * m_fopr->field_nin();
    int nvol = m_fopr->field_nvol();
    int nex = m_fopr->field_nex();
    m_x.reset(nin, nvol, nex);
    m_r.reset(nin, nvol, nex);
    m_p.reset(nin, nvol, nex);
    m_s.reset(nin, nvol, nex);
  }

  if (m_fopr) {
    m_fopr->set_use_QDW(m_use_QDW);
  }
}


//====================================================================
template<typename AFIELD>
void ASolver_CG<AFIELD>::set_parameters(const int Niter,
                                        const real_t Stop_cond)
{
  set_parameters(Niter, Stop_cond, InitialGuess::RHS);
}


//====================================================================
template<typename AFIELD>
void ASolver_CG<AFIELD>::set_parameters(const int Niter,
                                        const real_t Stop_cond,
                                        const InitialGuess init_guess_mode)
{
  ThreadManager::assert_single_thread(class_name);

  m_Niter        = Niter;
  m_Stop_cond    = Stop_cond;
  m_initial_mode = init_guess_mode;
  std::string prec = "double";
  if (sizeof(real_t) == 4) prec = "float";

  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Precision: %s\n", prec.c_str());
  vout.general(m_vl, "  Niter     = %d\n", m_Niter);
  vout.general(m_vl, "  Stop_cond = %16.8e\n", m_Stop_cond);
  vout.general(m_vl, "  init_guess_mode: %d\n", m_initial_mode);
}


//====================================================================
template<typename AFIELD>
void ASolver_CG<AFIELD>::solve(AFIELD& xq, const AFIELD& b,
                               int& Nconv, real_t& diff)
{
  copy(m_s, b);

  real_t sr    = m_s.norm2(m_use_qdw_int);
  real_t snorm = 1.0 / sr;
  vout.detailed(m_vl, "  snorm = %22.15e\n", snorm);

  real_t rr, rrp;
  int    nconv = -1;

  solve_CG_init(rrp, rr);
  vout.detailed(m_vl, "  init: %22.15e\n", rr * snorm);

  for (int iter = 0; iter < m_Niter; ++iter) {
    solve_CG_step(rrp, rr);
    vout.detailed(m_vl, "%6d  %22.15e\n", iter, rr * snorm);

    if (rr * snorm < m_Stop_cond) {
      nconv = iter;
#pragma omp master
      {
        m_nconv = nconv + 1;
      }
      break;
    }
  }

  if (nconv == -1)  {
    vout.crucial(m_vl, "Error at %s: not converged\n",
                 class_name.c_str());
    vout.crucial(m_vl, "  iter(final): %8d  %22.15e\n",
                 m_Niter, rr * snorm);
#pragma omp master
    {
      m_nconv = nconv;
    }
#pragma omp barrier
    //exit(EXIT_FAILURE);
  }

  vout.detailed(m_vl, "converged:\n");
  vout.detailed(m_vl, "  nconv = %d\n", nconv);

  copy(xq, m_x);

  AFIELD s_temp(b.nin(), b.nvol(), b.nex());
  m_fopr->mult(s_temp, xq);

  s_temp.axpy(real_t(-1.0), b, false);
  real_t diff2 = s_temp.norm2(false);

#pragma omp master
  {
    diff  = diff2;
    Nconv = m_nconv;
  }
#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void ASolver_CG<AFIELD>::solve_CG_init(real_t& rrp, real_t& rr)
{
  if (m_initial_mode == InitialGuess::RHS) {
#ifdef DEBUG
    vout.general(m_vl, "%s: using InitialGuess::RHS\n", class_name.c_str());
#endif
    copy(m_r, m_s);
    copy(m_x, m_s);
    m_fopr->set_mode("DdagD");
    vout.general(m_vl, "ASolver_CG solve_CG_init: m_x.nin() = %d, m_s.nin() = %d\n", m_x.nin(), m_s.nin());
    m_fopr->mult(m_s, m_x);
    m_s.normalize(m_use_qdw_int);
    m_r.axpy(real_t(-1.0), m_s, m_use_qdw_int);
    m_r.normalize(m_use_qdw_int);
    copy(m_p, m_r);
    rr  = m_r.norm2(m_use_qdw_int);
    rrp = rr;
  } else if (m_initial_mode == InitialGuess::GIVEN) {
    vout.crucial("%s: InitialGuess::GIVEN is not yet ready\n", class_name.c_str());
    exit(EXIT_FAILURE);
  } else if (m_initial_mode == InitialGuess::ZERO) {
#ifdef DEBUG
    vout.general(m_vl, "%s: using InitialGuess::ZERO\n", class_name.c_str());
#endif
    copy(m_r, m_s);
    m_s.set(0.0);
    m_x.set(0.0);
    copy(m_p, m_r);
    rr  = m_r.norm2(m_use_qdw_int);
    rrp = rr;
  } else {
    vout.crucial("%s: unkown init guess mode\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void ASolver_CG<AFIELD>::solve_CG_step(real_t& rrp, real_t& rr)
{
  using complex_t = typename AFIELD::complex_t;

  m_fopr->mult(m_s, m_p);
  m_s.normalize(m_use_qdw_int);

  real_t pap = m_s.dot(m_p, m_use_qdw_int);
  //    m_fopr->mult_normA_dev(pap, m_s, m_p);
  real_t cr = rrp / pap;

  m_x.axpy(cr, m_p, m_use_qdw_int);
  m_x.normalize(m_use_qdw_int);

  m_r.axpy(-cr, m_s, m_use_qdw_int);
  m_r.normalize(m_use_qdw_int);
  rr = m_r.norm2(m_use_qdw_int);

  real_t bk = rr / rrp;

  aypx(bk, m_p, m_r, m_use_qdw_int);
  m_p.normalize(m_use_qdw_int);

  rrp = rr;
}


//====================================================================
template<typename AFIELD>
double ASolver_CG<AFIELD>::flop_count()
{
  int Nin  = m_fopr->field_nin();
  int Nvol = m_fopr->field_nvol();
  int Nex  = m_fopr->field_nex();
  int NPE  = CommonParameters::NPE();

  int ninit = 1;

  double flop_field  = static_cast<double>(Nin * Nvol * Nex) * NPE;
  double flop_vector = (6 + ninit * 4 + m_nconv * 11) * flop_field;
  double flop_fopr   = (1 + ninit + m_nconv) * m_fopr->flop_count();

  double flop = flop_vector + flop_fopr;

  return flop;
}


//====================================================================
//============================================================END=====
