/*!
        @file    $Id:: asolver_Richardson-tmpl.h #$

        @brief

        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2#$

        @version $LastChangedRevision: 2499 $
*/

#include "asolver_Richardson.h"

#include  "lib/ResourceManager/threadManager.h"


template<typename AFIELD>
const std::string ASolver_Richardson<AFIELD>::class_name
  = "ASolver_Richardson";
//====================================================================
template<typename AFIELD>
void ASolver_Richardson<AFIELD>::init()
{
  ThreadManager::assert_single_thread(class_name);

  int nin  = m_fopr->field_nin();
  int nvol = m_fopr->field_nvol();
  int nex  = m_fopr->field_nex();

  m_r.reset(nin, nvol, nex);
  m_s.reset(nin, nvol, nex);
  m_nconv = -1;
}


//====================================================================
template<typename AFIELD>
void ASolver_Richardson<AFIELD>::tidyup()
{
  //  ThreadManager::assert_single_thread(class_name);
  // nothing is to be deleted.
}


//====================================================================
template<typename AFIELD>
void ASolver_Richardson<AFIELD>::set_parameters(
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


  // "restart" is not yet supported
  int Niter2 = Niter * Nrestart;
  set_parameters(Niter2, Stop_cond);
}


//====================================================================
template<typename AFIELD>
void ASolver_Richardson<AFIELD>::set_parameters(const int Niter,
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

  if (m_initial_mode != InitialGuess::ZERO) {
    vout.crucial(m_vl, "%s: WARNING: initial guess is not InitialGuess::ZERO\n",
                 class_name.c_str());
  }
}


//====================================================================
template<typename AFIELD>
void ASolver_Richardson<AFIELD>::set_parameters(const int Niter,
                                                const real_t Stop_cond)
{
  set_parameters(Niter, Stop_cond, InitialGuess::ZERO);
}


//====================================================================
template<typename AFIELD>
void ASolver_Richardson<AFIELD>::set_parameters(const int Niter,
                                                const real_t Stop_cond,
                                                const bool use_init_guess)
{
  // for backward compatibility
  if (use_init_guess) {
    set_parameters(Niter, Stop_cond, InitialGuess::GIVEN);
  } else {
    set_parameters(Niter, Stop_cond, InitialGuess::ZERO);
  }
}


//====================================================================
template<typename AFIELD>
void ASolver_Richardson<AFIELD>::solve(AFIELD& xq,
                                       const AFIELD& b,
                                       int& Nconv, real_t& diff)
{
  vout.paranoiac(m_vl, "%s: solver start.\n", class_name.c_str());


  real_t bnorm2, bnorm, bnorm_inv, bnorm2_inv;

  // initial barrier
#pragma omp barrier

  bnorm2 = norm2(b);
  if (!(bnorm2 > 0.0)) {
    xq.set(0.0);
    Nconv = 0;
    diff  = 0.0;
    return;
  }

  bnorm      = sqrt(bnorm2);
  bnorm_inv  = real_t(1.0) / bnorm;
  bnorm2_inv = real_t(1.0) / bnorm2;

  copy(m_r, b);

  //  solve_init(b, xq, rr);
  m_r.scal(-1.0);  // assumes the initial guess is zero
  xq.set(0.0);

  m_prec->reset_flop_count();


  int nconv = -1;


  double r2 = bnorm2; // m_s.norm2();
  vout.detailed(m_vl, " Richardson init: %22.15e\n", r2);
#pragma omp barrier

  for (int iter = 0; iter < m_Niter; ++iter) {
    if (r2 < bnorm2 * m_Stop_cond) {
      diff  = r2 * bnorm2_inv;
      nconv = iter;  // counting fermion mult
      break;
    }
    double normalization = sqrt(r2);
    m_r.scal(1.0 / normalization);
#pragma omp barrier
    m_prec->mult(m_s, m_r);

    axpy(xq, -normalization, m_s);
    m_fopr->mult(m_r, xq);
    axpy(m_r, -1.0, b);

    r2 = m_r.norm2();
  }
  if (nconv == -1) {
    vout.detailed(m_vl, "Richardson NOT converged:\n");
    //    vout.crucial(m_vl, "Error at %s: not converged\n",
    //                 class_name.c_str());
    //    vout.crucial(m_vl, "  iter(final): %8d  %22.15e\n",
    //                 m_Niter, rr*snorm);
    //    exit(EXIT_FAILURE);
  } else {
    vout.detailed(m_vl, "Richardson converged:\n");
  }

  // check the solution
  m_fopr->mult(m_r, xq);
  axpy(m_r, real_t(-1.0), b);
  real_t diff2 = norm2(m_r);
  vout.general(m_vl, "Richardson  nconv = %d, diff2 = %22.15e\n", nconv, diff2 * bnorm2_inv);

#pragma omp master
  {
    m_nconv = nconv;
    Nconv   = m_nconv; // include solve_init();
    diff    = sqrt(diff2 * bnorm2_inv);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void ASolver_Richardson<AFIELD>::solve_init(const AFIELD& b,
                                            const AFIELD& xq,
                                            real_t& rr)
{
}


//====================================================================
template<typename AFIELD>
double ASolver_Richardson<AFIELD>::flop_count()
{
  return flop_count_intermediate(m_nconv);
}


//====================================================================
template<typename AFIELD>
double ASolver_Richardson<AFIELD>::flop_count_intermediate(const int iter)
{
  int Nin  = m_fopr->field_nin();
  int Nvol = m_fopr->field_nvol();
  int Nex  = m_fopr->field_nex();
  int NPE  = CommonParameters::NPE();

  int ninit = 1;

  double flop_field = static_cast<double>(Nin * Nvol * Nex) * NPE;
  // norm: 2
  // scal: 1
  // axpy: 2
  // sub:  1
  double flop_vector = (1 + (2 + 1 + 2 + 1) * iter) * flop_field;
  double flop_fopr   = iter * m_fopr->flop_count();
  double flop_prec   = m_prec->flop_count();

  double flop = flop_vector + flop_fopr + flop_prec;

  return flop;
}


//====================================================================
