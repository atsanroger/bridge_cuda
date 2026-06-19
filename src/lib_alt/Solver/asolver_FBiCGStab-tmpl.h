/*!
        @file    $Id:: asolver_FBiCGStab-tmpl.h #$

        @brief

        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2#$

        @version $LastChangedRevision: 2499 $
*/

#include "asolver_FBiCGStab.h"

#include  "lib/ResourceManager/threadManager.h"


template<typename AFIELD>
const std::string ASolver_FBiCGStab<AFIELD>::class_name
  = "ASolver_FBiCGStab";
//====================================================================
template<typename AFIELD>
void ASolver_FBiCGStab<AFIELD>::init()
{
  ThreadManager::assert_single_thread(class_name);

  int nin  = m_fopr->field_nin();
  int nvol = m_fopr->field_nvol();
  int nex  = m_fopr->field_nex();

  m_r.reset(nin, nvol, nex);
  m_x.reset(nin, nvol, nex);
  m_p.reset(nin, nvol, nex);
  m_s.reset(nin, nvol, nex);
  m_v.reset(nin, nvol, nex);
  m_t.reset(nin, nvol, nex);
  m_rh.reset(nin, nvol, nex);
  m_u.reset(nin, nvol, nex);

  m_nconv           = -1;
  m_Omega_tolerance = 0.7;
}


//====================================================================
template<typename AFIELD>
void ASolver_FBiCGStab<AFIELD>::tidyup()
{
  //  ThreadManager::assert_single_thread(class_name);
  // nothing is to be deleted.
}


//====================================================================
template<typename AFIELD>
void ASolver_FBiCGStab<AFIELD>::set_parameters(
  const Parameters& params)
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

  double Omega_tolerance;
  err = params.fetch_double("Omega_tolerance", Omega_tolerance);
  if (err) {
    Omega_tolerance = 0.7;
  }

  // "restart" is not yet supported
  int Niter2 = Niter * Nrestart;
  set_parameters(Niter2, Stop_cond);

  set_parameters_BiCGStab_series(Omega_tolerance);
}


//====================================================================
template<typename AFIELD>
void ASolver_FBiCGStab<AFIELD>::set_parameters(const int Niter,
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
  vout.general(m_vl, "  InitialGuess = %d\n", m_initial_mode);
}


//====================================================================
template<typename AFIELD>
void ASolver_FBiCGStab<AFIELD>::set_parameters(const int Niter,
                                               const real_t Stop_cond)
{
  set_parameters(Niter, Stop_cond, InitialGuess::RHS);
}


//====================================================================
template<typename AFIELD>
void ASolver_FBiCGStab<AFIELD>::set_parameters(const int Niter,
                                               const real_t Stop_cond,
                                               const bool use_init_guess)
{
  // for backward compatibility
  if (use_init_guess) {
    set_parameters(Niter, Stop_cond, InitialGuess::GIVEN);
  } else {
    set_parameters(Niter, Stop_cond, InitialGuess::RHS);
  }
}


//====================================================================
template<typename AFIELD>
void ASolver_FBiCGStab<AFIELD>::set_parameters_BiCGStab_series(const real_t Omega_tolerance)
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
void ASolver_FBiCGStab<AFIELD>::solve(AFIELD& xq,
                                      const AFIELD& b,
                                      int& Nconv, real_t& diff)
{
  vout.paranoiac(m_vl, "%s: solver start.\n", class_name.c_str());


  real_t  bnorm2, bnorm, bnorm_inv, bnorm2_inv;
  real_t  rr;
  coeff_t prev;

  // initial barrier
#pragma omp barrier

  bnorm2 = norm2(b);
  if (!(bnorm2 > 0.0)) {
    xq.set(0.0);
    return;
  }
  bnorm      = sqrt(bnorm2);
  bnorm_inv  = real_t(1.0) / bnorm;
  bnorm2_inv = real_t(1.0) / bnorm2;

  m_prec->reset_flop_count();

  int nconv = -1;

  solve_init(b, xq, rr, prev);
  vout.detailed(m_vl, " FBiCGStab init: %22.15e\n", rr * bnorm2_inv);
  assert(rr > 0.0);

  int iter2 = 0;
  for (int iter = 0; iter < m_Niter; ++iter) {
    solve_step1(rr, prev);
    iter2++;
    double flop = 0.0;
#pragma omp master
    {
      flop = flop_count_intermediate(iter2);
    }
    vout.detailed(m_vl, "FBiCGStab %d  %22.15e   Flop(double+float) %e\n", iter2, rr * bnorm2_inv, flop);
    if (rr * bnorm2_inv < m_Stop_cond) {
      nconv = iter2;  // counting fermion mult
      break;
    }
    solve_step2(rr, prev);
    iter2++;
#pragma omp master
    {
      flop = flop_count_intermediate(iter2);
    }
    vout.detailed(m_vl, "FBiCGStab %d  %22.15e   Flop(double+float) %e\n", iter2, rr * bnorm2_inv, flop);
    if (rr * bnorm2_inv < m_Stop_cond) {
      nconv = iter2;  // counting fermion mult
      break;
    }
  }

  if (nconv == -1) {
    vout.detailed(m_vl, "FBiCGStab NOT converged:\n");
    //    vout.crucial(m_vl, "Error at %s: not converged\n",
    //                 class_name.c_str());
    //    vout.crucial(m_vl, "  iter(final): %8d  %22.15e\n",
    //                 m_Niter, rr*snorm);
    //    exit(EXIT_FAILURE);
  }

  vout.detailed(m_vl, "FBiCGStab converged:\n");

  copy(xq, m_x);

  // check the solution
  m_fopr->mult(m_s, xq);
  axpy(m_s, real_t(-1.0), b);
  real_t diff2 = norm2(m_s);
  vout.general(m_vl, "FBiCGStab  nconv = %d, diff2 = %22.15e\n", nconv, diff2 * bnorm2_inv);


#pragma omp master
  {
    m_nconv = nconv;
    Nconv   = m_nconv + 1; // include solve_init();
    diff    = sqrt(diff2 * bnorm2_inv);
  }
#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void ASolver_FBiCGStab<AFIELD>::solve_init(const AFIELD& b,
                                           const AFIELD& xq,
                                           real_t& rr,
                                           coeff_t& prev)
{
  // initial guess and its residual
  if (m_initial_mode == InitialGuess::RHS) {
    copy(m_r, b);
    copy(m_x, b);
    m_fopr->mult(m_v, m_x);
    axpy(m_r, real_t(-1.0), m_v);
  } else if (m_initial_mode == InitialGuess::GIVEN) {
    copy(m_r, b);
    copy(m_x, xq);
    m_fopr->mult(m_v, m_x);
    axpy(m_r, real_t(-1.0), m_v);
  } else if (m_initial_mode == InitialGuess::ZERO) {
    copy(m_r, b);
    m_x.set(0.0);
  } else {
    vout.crucial("%s: unkown init guess mode\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  copy(m_rh, m_r);

  rr = norm2(m_r);
  assert(rr > 0.0);

  copy(m_p, m_r);

  // initial parameters
  prev.rho   = rr;
  prev.alpha = 0.0;

  vout.detailed(m_vl, "ASolver_FBiCGStab  rr      = %23.16e [init]\n", rr);
}


//====================================================================
template<typename AFIELD>
void ASolver_FBiCGStab<AFIELD>::solve_step1(real_t& rr, coeff_t& prev)
{
  using complex_t = std::complex<real_t>;

  m_prec->mult(m_u, m_p);
  m_fopr->mult(m_v, m_u); // m_v = D M p

  complex_t alpha_d = dotc(m_rh, m_v);
  complex_t alpha   = prev.rho / alpha_d;
  prev.alpha = alpha;

  axpy(m_r, -alpha, m_v);
  axpy(m_x, alpha, m_u);
  rr = norm2(m_r);

#ifdef DEBUG_FBICGSTAB
  vout.detailed(m_vl, "ASolver_FBiCGStab  alpha_d = %23.16e %23.16e\n", real(alpha_d), std::imag(alpha_d));
  vout.detailed(m_vl, "ASolver_FBiCGStab  alpha   = %23.16e %23.16e\n", real(alpha), std::imag(alpha));
  vout.detailed(m_vl, "ASolver_FBiCGStab  rr      = %23.16e\n", rr);
#endif
}


//====================================================================
template<typename AFIELD>
void ASolver_FBiCGStab<AFIELD>::solve_step2(real_t& rr, coeff_t& prev)
{
  using complex_t = std::complex<real_t>;
  m_prec->mult(m_u, m_r);
  m_fopr->mult(m_t, m_u);

  complex_t omega_n = dotc(m_t, m_r); // omega_n = t * s;
  real_t    omega_d = norm2(m_t);     // omega_d = t * t;
  real_t    r_norm2 = norm2(m_r);
  complex_t omega   = omega_n / omega_d;

  // for omega prescrition
  const real_t Omega_tolerance = m_Omega_tolerance;
  const real_t abs_rho         = abs(omega_n) / sqrt(omega_d * r_norm2);

#ifdef DEBUG_FBICGSTAB
  vout.detailed(m_vl, "ASolver_FBiCGStab  omega_n = %23.16e %23.16e\n", real(omega_n), std::imag(omega_n));
  vout.detailed(m_vl, "ASolver_FBiCGStab  omega_d = %23.16e\n", omega_d);
  vout.detailed(m_vl, "ASolver_FBiCGStab  omega   = %23.16e %23.16e\n", real(omega), std::imag(omega));
  vout.detailed(m_vl, "ASolver_FBiCGStab  r_norm2 = %23.16e\n", r_norm2);
  vout.detailed(m_vl, "ASolver_FBiCGStab  abs_rho = %23.16e\n", abs_rho);
#endif

  //- a prescription to improve stability of BiCGStab_Cmplx
  if (abs_rho < Omega_tolerance) {
    omega *= Omega_tolerance / abs_rho;
    vout.detailed(m_vl, "ASolver_FBiCGStab  omega   = %23.16e %23.16e  [Omega presciption applied]\n", real(omega), std::imag(omega));
  }

  axpy(m_x, omega, m_u);
  axpy(m_r, -omega, m_t);

  complex_t rho  = dotc(m_rh, m_r);
  complex_t beta = (rho / prev.rho) * (prev.alpha * omega);
  prev.rho = rho;

  // p = r + bet * (p - omega * v);
  axpy(m_p, -omega, m_v);  // p := p - omega v
  aypx(beta, m_p, m_r);

  rr = norm2(m_r);      // rr = r * r;
#ifdef DEBUG_FBICGSTAB
  vout.detailed(m_vl, "ASolver_FBiCGStab  rho     = %23.16e %23.16e\n", real(rho), std::imag(rho));
  vout.detailed(m_vl, "ASolver_FBiCGStab  beta    = %23.16e %23.16e\n", real(beta), std::imag(beta));
  vout.detailed(m_vl, "ASolver_FBiCGStab  rr      = %23.16e\n", rr);
#endif
}


//====================================================================
template<typename AFIELD>
void ASolver_FBiCGStab<AFIELD>::prec(AFIELD& v, AFIELD& w)
{
  copy(v, w);
}


//====================================================================
template<typename AFIELD>
double ASolver_FBiCGStab<AFIELD>::flop_count()
{
  return flop_count_intermediate(m_nconv);
}


//====================================================================
template<typename AFIELD>
double ASolver_FBiCGStab<AFIELD>::flop_count_intermediate(const int iter)
{
  int Nin  = m_fopr->field_nin();
  int Nvol = m_fopr->field_nvol();
  int Nex  = m_fopr->field_nex();
  int NPE  = CommonParameters::NPE();

  int ninit = 1;

  double flop_field = static_cast<double>(Nin * Nvol * Nex) * NPE;
  // step1:  4+4+4+2 = 14
  // step2: 4+2+4+4+4+4+4+2 = 4*6 + 2*2= 28
  double flop_vector = (6 + ninit * 4 + (iter / 2) * 42) * flop_field;
  if (iter % 2 == 1) {
    flop_vector += 14 * flop_field;
  }

  double flop_fopr = (1 + ninit + iter) * m_fopr->flop_count();
  double flop_prec = m_prec->flop_count();

  double flop = flop_vector + flop_fopr + flop_prec;

  return flop;
}


//====================================================================
