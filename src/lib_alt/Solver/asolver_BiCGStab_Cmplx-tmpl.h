/*!
        @file    solver_BiCGStab_Cmplx-tmpl.h
        @brief
        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
        @version $LastChangedRevision: 2499 $
*/

template<typename AFIELD>
const std::string ASolver_BiCGStab_Cmplx<AFIELD>::class_name
  = "ASolver_BiCGStab_Cmplx";

namespace {
#ifndef AFIELD_HAS_SUB
  template<typename AFIELD>
  void sub(AFIELD& v, const AFIELD& w)
  {
    axpy(v, -1.0, w);
  }
#endif

#ifndef AFIELD_HAS_DOTC_AND_NORM2
  template<typename AFIELD>
  void dotc_and_norm2(typename AFIELD::complex_t& dotc_vw,
                      typename AFIELD::real_t& v_norm2,
                      typename AFIELD::real_t& w_norm2,
                      const AFIELD& v, const AFIELD& w)
  {
    dotc_vw = dotc(v, w);
    v_norm2 = v.norm2();
    w_norm2 = w.norm2();
  }
#endif
}

//====================================================================
template<typename AFIELD>
void ASolver_BiCGStab_Cmplx<AFIELD>::init(void)
{
  ThreadManager::assert_single_thread(class_name);

  vout.detailed(m_vl, "%s: being setup.\n", class_name.c_str());

  int nin  = m_fopr->field_nin();
  int nvol = m_fopr->field_nvol();
  int nex  = m_fopr->field_nex();

  m_r.reset(nin, nvol, nex);
  m_p.reset(nin, nvol, nex);
  m_v.reset(nin, nvol, nex);
  m_t.reset(nin, nvol, nex);
  m_rh.reset(nin, nvol, nex);

  m_ecrit           = 1.e-32;
  m_nconv           = -1;
  m_Omega_tolerance = 0.7;

  vout.detailed(m_vl, "%s: setup finished.\n", class_name.c_str());
}


//====================================================================
template<typename AFIELD>
void ASolver_BiCGStab_Cmplx<AFIELD>::tidyup(void)
{
  // ThreadManager::assert_single_thread(class_name);
  // nothing is to be deleted.
}


//====================================================================
template<typename AFIELD>
void ASolver_BiCGStab_Cmplx<AFIELD>::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  int          Niter, Nrestart;
  double       Stop_cond;
  double       Omega_tolerance;
  InitialGuess initial_guess;
  int          err = 0;
  err += params.fetch_int("maximum_number_of_iteration", Niter);
  err += params.fetch_int("maximum_number_of_restart", Nrestart);
  err += params.fetch_double("convergence_criterion_squared", Stop_cond);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
  err = params.fetch_double("Omega_tolerance", Omega_tolerance);
  if (err) {
    Omega_tolerance = 0.7;
  }
  std::string initial_guess_str;
  err           = params.fetch_string("initial_guess_mode", initial_guess_str);
  initial_guess = InitialGuess::RHS;
  if (!err) {
    vout.detailed(m_vl, "%s: initial_guess_str=%s\n",
                  class_name.c_str(), initial_guess_str.c_str());

    if (initial_guess_str == "ZERO") {
      initial_guess = InitialGuess::ZERO;
    } else if (initial_guess_str == "GIVEN") {
      initial_guess = InitialGuess::GIVEN;
    } else {
      vout.crucial(m_vl, "%s: unknown initial guess mode string: %d\n",
                   class_name.c_str(), initial_guess_str.c_str());
      exit(EXIT_FAILURE);
    }
  }

  // "restart" is not yet supported
  int Niter2 = Niter * Nrestart;

  // N.B.
  // - default initial mode is RHS
  // - negative Stop_cond is for fixed nubmer of iterations
  set_parameters(Niter2, Stop_cond, initial_guess);
  set_parameters_BiCGStab_series(Omega_tolerance);
}


//====================================================================
template<typename AFIELD>
void ASolver_BiCGStab_Cmplx<AFIELD>::set_parameters(const int Niter,
                                                    const real_t Stop_cond)
{
  set_parameters(Niter, Stop_cond, InitialGuess::RHS);
}


//====================================================================
template<typename AFIELD>
void ASolver_BiCGStab_Cmplx<AFIELD>::set_parameters(const int Niter,
                                                    const real_t Stop_cond,
                                                    const InitialGuess init_guess_mode)
{
  ThreadManager::assert_single_thread(class_name);

  m_Niter     = Niter;
  m_Stop_cond = Stop_cond;
  std::string prec = "double";
  if (sizeof(real_t) == 4) prec = "float";
  m_initial_mode = init_guess_mode;

  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Precision: %s\n", prec.c_str());
  vout.general(m_vl, "  Niter     = %d\n", m_Niter);
  vout.general(m_vl, "  Stop_cond = %16.8e\n", m_Stop_cond);
  vout.general(m_vl, "  init_guess_mode: %d\n", m_initial_mode);
}


//====================================================================
template<typename AFIELD>
void ASolver_BiCGStab_Cmplx<AFIELD>::set_parameters_BiCGStab_series(const real_t Omega_tolerance)
{
  ThreadManager::assert_single_thread(class_name);

  //- print input parameters
  vout.general(m_vl, "  Omega_tolerance = %8.2e\n", Omega_tolerance);

  //- range check
  // NB. Omega_tolerance == 0.0 is allowed.

  //- store values
  m_Omega_tolerance = Omega_tolerance;
}


//====================================================================
template<typename AFIELD>
void ASolver_BiCGStab_Cmplx<AFIELD>::solve(AFIELD& xq, const AFIELD& b,
                                           int& Nconv, real_t& diff)
{
  vout.paranoiac(m_vl, "%s: solver start:\n", class_name.c_str());

#pragma omp master
  {
    m_init_mult = 0;
  }

  copy(m_t, b);
  real_t snorm2    = norm2(b);
  real_t snorm     = sqrt(snorm2);
  real_t snorm_inv = real_t(1.0) / snorm;
  vout.detailed(m_vl, "  |b| = snorm = %22.15e\n", snorm);
  scal(m_t, snorm_inv);

  // solutoin vector
  m_px = &xq;
  AFIELD& m_x = *m_px;

  // coefficents
  coeff_t coeff;

  // initialication
  int    nconv = -1;
  real_t rr;
  solve_init(m_t, rr, coeff, snorm_inv, m_initial_mode);
  vout.detailed(m_vl, "  init: %22.15e\n", rr);

  for (int iter = 0; iter < m_Niter; ++iter) {
    int iflg = 0;
    solve_step(rr, iflg, coeff);
    vout.detailed(m_vl, "%6d  %22.15e\n", iter, rr);
    if (iflg != 0) { // restart
      copy(m_t, b);
      scal(m_t, snorm_inv);
      solve_init_GIVEN(m_t, rr, coeff, real_t(1.0));
      vout.detailed(m_vl, "%s: restarted.\n", class_name.c_str());
      vout.detailed(m_vl, "%6d  %22.15e\n", iter, rr);
    }

    if (rr < m_Stop_cond) {
      nconv = 2 * (iter + 1);  // counting fermion mult
      break;
    }
  }
  if (nconv == -1) {
    if ((m_Stop_cond < 0) || (Nconv < 0)) {
      vout.detailed(m_vl, "  truncating with iter= %d\n", 2 * m_Niter);
      nconv = 2 * m_Niter;
    } else {
      vout.crucial(m_vl, "Error at %s: not converged\n",
                   class_name.c_str());
      vout.crucial(m_vl, "  iter(final): %8d  %22.15e\n",
                   m_Niter, rr * snorm2);
      exit(EXIT_FAILURE);
    }
  } else {
    vout.detailed(m_vl, "converged:\n");
  }
  vout.detailed(m_vl, "  nconv = %d\n", nconv);

  scal(xq, snorm);

  m_fopr->mult(m_t, xq);
  sub(m_t, b);
  real_t diff2 = norm2(m_t);

#pragma omp master
  {
    m_nconv = nconv;           // include solve_init();
  }
  Nconv = nconv + m_init_mult; // include solve_init();
  diff  = sqrt(diff2 / snorm2);

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
inline void ASolver_BiCGStab_Cmplx<AFIELD>::solve_init(const AFIELD& b,
                                                       real_t& rr,
                                                       coeff_t& prev,
                                                       const real_t scale,
                                                       const InitialGuess init_mode)
{
  if (init_mode == InitialGuess::RHS) {
    solve_init_RHS(b, rr, prev);
  } else if (init_mode == InitialGuess::GIVEN) {
    solve_init_GIVEN(b, rr, prev, scale);
  } else if (init_mode == InitialGuess::ZERO) {
    solve_init_ZERO(b, rr, prev);
  } else {
    vout.crucial("%s: unkown init guess mode\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void ASolver_BiCGStab_Cmplx<AFIELD>::solve_init_RHS(const AFIELD& b,
                                                    real_t& rr,
                                                    coeff_t& prev)
{
  AFIELD& m_x = *m_px;
  copy(m_x, b);
  copy(m_r, b);

  {
    real_t temp = norm2(m_x);
    vout.detailed(m_vl, "  |m_x|^2 = %23.16e [init_RHS]\n", temp);
  }

  m_fopr->mult(m_v, m_x);
  axpy(m_r, real_t(-1.0), m_v);

  {
    real_t temp = norm2(m_v);
    vout.detailed(m_vl, "  |m_v|^2 = %23.16e [init_RHS]\n", temp);
  }

  copy(m_rh, m_r);
  rr = norm2(m_r);
  m_p.set(real_t(0.0));
  m_v.set(real_t(0.0));

  prev.rho   = real_t(1.0);
  prev.alpha = real_t(1.0);
  prev.omega = real_t(1.0);

  vout.detailed(m_vl, "ASolver_BiCGStab_Cmplx   rr      = %23.16e [init_RHS]\n", rr);


#pragma omp master
  {
    m_init_mult++;
  }
}


//====================================================================
template<typename AFIELD>
void ASolver_BiCGStab_Cmplx<AFIELD>::solve_init_GIVEN(const AFIELD& b,
                                                      real_t& rr,
                                                      coeff_t& prev,
                                                      const real_t scale)
{
  AFIELD& m_x = *m_px;
  if (scale != real_t(1.0)) {
    scal(m_x, scale);
  }
  copy(m_r, b);
  m_fopr->mult(m_v, m_x);
  axpy(m_r, real_t(-1.0), m_v);

  copy(m_rh, m_r);
  rr = norm2(m_r);
  m_p.set(real_t(0.0));
  m_v.set(real_t(0.0));

  prev.rho   = real_t(1.0);
  prev.alpha = real_t(1.0);
  prev.omega = real_t(1.0);

#pragma omp master
  {
    m_init_mult++;
  }
  vout.detailed(m_vl, "ASolver_BiCGStab_Cmplx  rr      = %23.16e [init_GIVEN]\n", rr);
}


//====================================================================
template<typename AFIELD>
void ASolver_BiCGStab_Cmplx<AFIELD>::solve_init_ZERO(const AFIELD& b,
                                                     real_t& rr,
                                                     coeff_t& prev)
{
  // x0 = 0
  AFIELD& m_x = *m_px;
  m_x.set(0.0);
  copy(m_r, b);
  copy(m_rh, m_r);
  rr = norm2(m_r);
  m_p.set(real_t(0.0));
  m_v.set(real_t(0.0));

  prev.rho   = real_t(1.0);
  prev.alpha = real_t(1.0);
  prev.omega = real_t(1.0);
  vout.detailed(m_vl, "ASolver_BiCGStab_Cmplx  rr      = %23.16e [init_ZERO]\n", rr);
}


//====================================================================
template<typename AFIELD>
void ASolver_BiCGStab_Cmplx<AFIELD>::solve_step(real_t& rr, int& iflg, coeff_t& prev)
{
  using complex_t = typename AFIELD::complex_t;
  const real_t Omega_tolerance = m_Omega_tolerance;
  AFIELD&      m_x             = *m_px;

  complex_t rho = dotc(m_rh, m_r);
#ifdef DEBUG_BICGSTAB_CMPLX
  vout.detailed(m_vl, "ASolver_BiCGStab_Cmplx  rho     = %23.16e %23.16e\n", real(rho), imag(rho));
#endif
  if (abs(rho) < m_ecrit) {
    iflg = 1;
    vout.detailed(m_vl, "too small value of rho = %23.16e %23.16e\n", real(rho), imag(rho));
#ifdef DEBUG
    double r2  = norm2(m_r);
    double rh2 = norm2(m_rh);
    vout.detailed(m_vl, "   |m_r|^2  = %23.16e\n", r2);
    vout.detailed(m_vl, "   |m_rh|^2 = %23.16e\n", rh2);
#endif
    return;
  }
  complex_t bet = rho * prev.alpha / (prev.rho * prev.omega);

  // p = r + bet * (p - omega_prev * v);
  axpy(m_p, -prev.omega, m_v);
  aypx(bet, m_p, m_r);

  m_fopr->mult(m_v, m_p);

  complex_t aden = dotc(m_rh, m_v);
#ifdef DEBUG_BICGSTAB_CMPLX
  vout.detailed(m_vl, "ASolver_BiCGStab_Cmplx  bet     = %23.16e %23.16e\n", real(bet), imag(bet));
  vout.detailed(m_vl, "ASolver_BiCGStab_Cmplx  aden    = %23.16e %23.16e\n", real(aden), imag(aden));
#endif
  if (abs(aden) < m_ecrit) {
    iflg = 1;
    vout.detailed(m_vl, "too small denominator  aden = %23.16e %23.16e\n", real(aden), imag(aden));
    return;
  }
  complex_t alpha = rho / aden;
#ifdef DEBUG_BICGSTAB_CMPLX
  vout.detailed(m_vl, "ASolver_BiCGStab_Cmplx  alpha   = %23.16e %23.16e\n", real(alpha), imag(alpha));
#endif

  // "r" from here is an alias of "s"
  axpy(m_r, -alpha, m_v);
  m_fopr->mult(m_t, m_r);

  // todo: use something like dotc_and_norm2()
  //  complex_t omega_n = dotc(m_t, m_r);
  //  real_t omega_d = m_t.norm2();
  complex_t omega_n;
  real_t    omega_d, s_norm2;
  dotc_and_norm2(omega_n, omega_d, s_norm2, m_t, m_r);
  complex_t omega = omega_n / omega_d;

  // for omega prescription
  //  const real_t s_norm2=m_r.norm2();
  const double abs_rho = abs(omega_n) / sqrt(omega_d * s_norm2);

#ifdef DEBUG_BICGSTAB_CMPLX
  vout.detailed(m_vl, "ASolver_BiCGStab_Cmplx  omega_n = %23.16e %23.16e\n", real(omega_n), imag(omega_n));
  vout.detailed(m_vl, "ASolver_BiCGStab_Cmplx  omega_d = %23.16e\n", omega_d);
  vout.detailed(m_vl, "ASolver_BiCGStab_Cmplx  omega   = %23.16e %23.16e\n", real(omega), imag(omega));
  vout.detailed(m_vl, "ASolver_BiCGStab_Cmplx  s_norm2 = %23.16e\n", s_norm2);
  vout.detailed(m_vl, "ASolver_BiCGStab_Cmplx  abs_rho = %23.16e\n", abs_rho);
#endif


  //- a prescription to improve stability of BiCGStab_Cmplx
  if (abs_rho < Omega_tolerance) {
    vout.detailed(m_vl, "ASolver_BiCGStab_Cmplx  abs_rho = %23.16e,  omega   = %23.16e %23.16e  [Omega presciption applied]\n", abs_rho, real(omega), imag(omega));
    omega *= Omega_tolerance / abs_rho;
  }

  // x = x + omega "s" + alpha p
  axpy(m_x, omega, m_r);
  axpy(m_x, alpha, m_p);

  // "r" below is the residual vector
  axpy(m_r, -omega, m_t);

  prev.rho   = rho;
  prev.alpha = alpha;
  prev.omega = omega;
  rr         = norm2(m_r);
}


//====================================================================
template<typename AFIELD>
double ASolver_BiCGStab_Cmplx<AFIELD>::flop_count()
{
  int Nin  = m_fopr->field_nin();
  int Nvol = m_fopr->field_nvol();
  int Nex  = m_fopr->field_nex();
  int NPE  = CommonParameters::NPE();

  int ninit = 1;
  int iter  = m_nconv;
  //  if(m_nconv<0){ iter = m_Niter;}
  double flop_field  = static_cast<double>(Nin * Nvol * Nex) * NPE;
  double flop_vector = (6 + 4 * ninit + iter * 20) * flop_field;

  double flop_fopr = (1 + ninit + iter) * m_fopr->flop_count();

  if (m_initial_mode == InitialGuess::ZERO) {
    flop_vector -= (2 * ninit) * flop_field;
    flop_fopr   -= m_fopr->flop_count();
  }

  double flop = flop_vector + flop_fopr;

  return flop;
}


//====================================================================
//============================================================END=====
