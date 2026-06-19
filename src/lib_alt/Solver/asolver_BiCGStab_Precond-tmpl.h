/*!
      @file    asolver_BiCGStab_Precond-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#include "lib_alt/Solver/asolver_BiCGStab_Precond.h"


template<typename AFIELD>
const std::string ASolver_BiCGStab_Precond<AFIELD>::class_name
  = "ASolver_BiCGStab_Precond";
//====================================================================
template<typename AFIELD>
void ASolver_BiCGStab_Precond<AFIELD>::init()
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

  m_nconv = -1;
}


//====================================================================
template<typename AFIELD>
void ASolver_BiCGStab_Precond<AFIELD>::tidyup()
{
  //  ThreadManager::assert_single_thread(class_name);
  // nothing is to be deleted.
}


//====================================================================
template<typename AFIELD>
void ASolver_BiCGStab_Precond<AFIELD>::set_parameters(
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

  int Niter2 = Niter * Nrestart;
  set_parameters(Niter2, Stop_cond);
}


//====================================================================
template<typename AFIELD>
void ASolver_BiCGStab_Precond<AFIELD>::set_parameters(const int Niter,
                                                      const real_t Stop_cond)
{
  ThreadManager::assert_single_thread(class_name);

  m_Niter     = Niter;
  m_Stop_cond = Stop_cond;
  std::string prec = "double";
  if (sizeof(real_t) == 4) prec = "float";

  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Precision: %s\n", prec.c_str());
  vout.general(m_vl, "  Niter     = %d\n", m_Niter);
  vout.general(m_vl, "  Stop_cond = %16.8e\n", m_Stop_cond);
}


//====================================================================
template<typename AFIELD>
void ASolver_BiCGStab_Precond<AFIELD>::solve(AFIELD& xq,
                                             const AFIELD& b,
                                             int& Nconv, real_t& diff)
{
  vout.paranoiac(m_vl, "%s: solver start.\n", class_name.c_str());

  copy(m_s, b);

  real_t snorm, sr;
  real_t rr;

  sr    = norm2(m_s);
  snorm = real_t(1.0) / sr;

  m_prec->reset_flop_count();

  int nconv = -1;

  solve_init(b, rr);
  vout.detailed(m_vl, "  init: %22.15e\n", rr * snorm);

  for (int iter = 0; iter < m_Niter; ++iter) {
    solve_step(rr);
    vout.detailed(m_vl, "%6d  %22.15e\n", iter, rr * snorm);

    if (rr * snorm < m_Stop_cond) {
      nconv = 2 * (iter + 1);  // counting fermion mult
      break;
    }
  }

  if (nconv == -1) {
    vout.crucial(m_vl, "Error at %s: not converged\n",
                 class_name.c_str());
    vout.crucial(m_vl, "  iter(final): %8d  %22.15e\n",
                 m_Niter, rr * snorm);
    exit(EXIT_FAILURE);
  }

  vout.detailed(m_vl, "converged:\n");
  vout.detailed(m_vl, "  nconv = %d\n", nconv);

  copy(xq, m_x);

  m_fopr->mult(m_s, xq);

  axpy(m_s, real_t(-1.0), b);

  real_t diff2 = norm2(m_s);

#pragma omp master
  {
    m_nconv = nconv;
    Nconv   = m_nconv + 1; // include solve_init();
    diff    = diff2;
  }
#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void ASolver_BiCGStab_Precond<AFIELD>::solve_init(const AFIELD& b,
                                                  real_t& rr)
{
  copy(m_x, m_s);

  m_fopr->mult(m_v, m_x);
  copy(m_r, b);
  axpy(m_r, real_t(-1.0), m_v);
  copy(m_rh, m_r);

  rr = norm2(m_r);

  m_p.set(real_t(0.0));
  m_v.set(real_t(0.0));

  m_rho_prev   = real_t(1.0);
  m_alpha_prev = real_t(1.0);
  m_omega_prev = real_t(1.0);
}


//====================================================================
template<typename AFIELD>
void ASolver_BiCGStab_Precond<AFIELD>::solve_step(real_t& rr)
{
  real_t rho = dot(m_rh, m_r);
  real_t bet = rho * m_alpha_prev / (m_rho_prev * m_omega_prev);

  // p = r + bet * (p - omega_prev * v);
  axpy(m_p, -m_omega_prev, m_v);
  aypx(bet, m_p, m_r);

  m_prec->mult(m_u, m_p);
  m_fopr->mult(m_v, m_u);

  real_t aden  = dot(m_rh, m_v);   // dcomplex aden = rh * v;
  real_t alpha = rho / aden;

  copy(m_s, m_r);                   // s  = r
  axpy(m_s, -alpha, m_v);           // s += - alpha * v;

  m_prec->mult(m_v, m_s);           // t  = m_fopr->mult(s);
  m_fopr->mult(m_t, m_v);           // t  = m_fopr->mult(s);

  real_t omega_n = dot(m_t, m_s);   // omega_n = t * s;
  real_t omega_d = dot(m_t, m_t);   // omega_d = t * t;
  real_t omega   = omega_n / omega_d;

  //    axpy(m_x, omega, m_s);
  axpy(m_x, omega, m_v);
  axpy(m_x, alpha, m_u);

  copy(m_r, m_s);         // r  = s
  axpy(m_r, -omega, m_t); // r += - omega * t;

  rr = norm2(m_r);        // rr = r * r;

  m_rho_prev   = rho;
  m_alpha_prev = alpha;
  m_omega_prev = omega;
}


//====================================================================
template<typename AFIELD>
void ASolver_BiCGStab_Precond<AFIELD>::prec(AFIELD& v, AFIELD& w)
{
  copy(v, w);
}


//====================================================================
template<typename AFIELD>
double ASolver_BiCGStab_Precond<AFIELD>::flop_count()
{
  int Nin  = m_fopr->field_nin();
  int Nvol = m_fopr->field_nvol();
  int Nex  = m_fopr->field_nex();
  int NPE  = CommonParameters::NPE();

  int ninit = 1;

  double flop_field = static_cast<double>(Nin * Nvol * Nex) * NPE;

  double flop_vector = (6 + ninit * 4 + m_nconv * 11) * flop_field;
  double flop_fopr   = (1 + ninit + m_nconv) * m_fopr->flop_count();
  double flop_prec   = m_prec->flop_count();

  double flop = flop_vector + flop_fopr + flop_prec;

  return flop;
}


//============================================================END=====
