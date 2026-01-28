/*!
        @file    asolver_CG-tmpl.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2#$
        @version $LastChangedRevision: 2569 $
*/

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

  if(Stop_cond < 5.96e-8){  // limitation of FP16
    vout.crucial(m_vl, "Too small Stop_cond for FP16: %e\n", Stop_cond);
    exit(EXIT_FAILURE);
  }

  int Niter2 = Niter * Nrestart;
  set_parameters(Niter2, real_t(Stop_cond), init_guess_mode);
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
  if (sizeof(real_t) == 2) prec = "half";

  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Precision: %s\n", prec.c_str());
  vout.general(m_vl, "  Niter     = %d\n", m_Niter);
  vout.general(m_vl, "  Stop_cond = %12.4e\n", double(m_Stop_cond));
  vout.general(m_vl, "  init_guess_mode: %d\n", m_initial_mode);
}


//====================================================================
template<typename AFIELD>
void ASolver_CG<AFIELD>::solve(AFIELD& xq, const AFIELD& b,
                               int& Nconv, real_t& diff)
{
  copy(m_s, b);

  double sr;
  m_s.norm2_double(sr);
  double snorm = 1.0 / sr;
  vout.detailed(m_vl, "  snorm = %22.15e\n", snorm);

  double rr, rrp;
  int    nconv = -1;

  solve_CG_init_double(rrp, rr);
  vout.detailed(m_vl, "  init: %22.15e\n", rr * snorm);

  for (int iter = 0; iter < m_Niter; ++iter) {
    solve_CG_step_double(rrp, rr);
    vout.detailed(m_vl, "%6d  %22.15e\n", iter, rr * snorm);

    if (rr * snorm < m_Stop_cond) {
      nconv = iter;

      int ith = ThreadManager::get_thread_id();
      if (ith == 0) {
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
#pragma omp barrier
    //exit(EXIT_FAILURE);
  }

  vout.detailed(m_vl, "converged:\n");
  vout.detailed(m_vl, "  nconv = %d\n", nconv);

  copy(xq, m_x);

  m_fopr->mult(m_s, xq);

  axpy(m_s, real_t(-1.0), b);
  double diff2;
  m_s.norm2_double(diff2);

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) {
    diff  = real_t(diff2);
    Nconv = m_nconv;
  }
#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void ASolver_CG<AFIELD>::solve_CG_init_double(double& rrp, double& rr)
{
  if (m_initial_mode == InitialGuess::RHS) {
#ifdef DEBUG
    vout.general(m_vl, "%s: using InitialGuess::RHS\n", class_name.c_str());
#endif
    copy(m_r, m_s);
    copy(m_x, m_s);
    m_fopr->mult(m_s, m_x);
    axpy(m_r, real_t(-1.0), m_s);
    copy(m_p, m_r);
    //rr  = norm2(m_r);
    m_r.norm2_double(rr);
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
    //rr  = norm2(m_r);
    m_r.norm2_double(rr);
    rrp = rr;
  } else {
    vout.crucial("%s: unkown init guess mode\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void ASolver_CG<AFIELD>::solve_CG_step_double(double& rrp, double& rr)
{
  using complex_t = typename AFIELD::complex_t;

  m_fopr->mult(m_s, m_p);

  //real_t pap = dot(m_s, m_p);
  double pap;
  m_s.dot_double(pap, m_p);

  double cr = rrp / pap;

  axpy(m_x, real_t(cr), m_p);

  axpy(m_r, real_t(-cr), m_s);
  //rr = norm2(m_r);
  m_r.norm2_double(rr);

  double bk = rr / rrp;

  aypx(real_t(bk), m_p, m_r);

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


//============================================================END=====
