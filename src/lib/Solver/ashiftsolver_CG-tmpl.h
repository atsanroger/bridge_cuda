/*!
        @file    ashiftsolver_CG-tmpl.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-07-03 23:32:09 #$

        @version $LastChangedRevision: 2529 $
*/

#include "lib/Solver/ashiftsolver_CG.h"
#include "lib/ResourceManager/threadManager.h"

template<typename AFIELD, typename AFOPR>
const std::string AShiftsolver_CG<AFIELD, AFOPR>::class_name
  = "AShiftsolver_CG";

//====================================================================
template<typename AFIELD, typename AFOPR>
void AShiftsolver_CG<AFIELD, AFOPR>::set_parameters(
  const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  int    Niter;
  double Stop_cond;

  int err = 0;
  err += params.fetch_int("maximum_number_of_iteration", Niter);
  err += params.fetch_double("convergence_criterion_squared", Stop_cond);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(Niter, real_t(Stop_cond));
}


//====================================================================
template<typename AFIELD, typename AFOPR>
void AShiftsolver_CG<AFIELD, AFOPR>::get_parameters(Parameters& params) const
{
  params.set_int("maximum_number_of_iteration", m_Niter);
  params.set_double("convergence_criterion_squared", double(m_Stop_cond));

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
template<typename AFIELD, typename AFOPR>
void AShiftsolver_CG<AFIELD, AFOPR>::set_parameters(
                                              const int Niter,
                                              const real_t Stop_cond)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Niter     = %d\n", Niter);
  vout.general(m_vl, "  Stop_cond = %8.2e\n", Stop_cond);

  //- range check
  int err = 0;
  err += ParameterCheck::non_negative(Niter);
  err += ParameterCheck::square_non_zero(Stop_cond);

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  //- store values
  m_Niter     = Niter;
  m_Stop_cond = Stop_cond;
}


//====================================================================
template<typename AFIELD, typename AFOPR>
void AShiftsolver_CG<AFIELD, AFOPR>::solve(
                                 std::vector<AFIELD>& xq,
                                 const std::vector<real_t>& sigma,
                                 const AFIELD& b,
                                 int& Nconv,
                                 real_t& diff)
{
#pragma omp barrier

  int Nshift = sigma.size();

  vout.paranoiac(m_vl, "  Shift CG solver start.\n");
  vout.paranoiac(m_vl, "    number of shift = %d\n", Nshift);
  vout.paranoiac(m_vl, "    values of shift:\n");
  for (int i = 0; i < Nshift; ++i) {
    vout.paranoiac(m_vl, "    %d  %12.8f\n", i, sigma[i]);
  }

  m_snorm = 1.0 / b.norm2();

  int Nconv2 = -1;

  reset_field(b, sigma, Nshift);

  copy(m_s, b);
  copy(m_r, b);

  real_t rr = 0.0;

  solve_init(rr);

  vout.detailed(m_vl, "    iter: %8d  %22.15e\n", 0, rr * m_snorm);

  bool is_converged = false;

  for (int iter = 0; iter < m_Niter; iter++) {
    solve_step(rr);

    Nconv2 += 1;

    vout.detailed(m_vl, "    iter: %8d  %22.15e  %4d\n",
                  (iter + 1), rr * m_snorm, m_Nshift2);

    if (rr * m_snorm < m_Stop_cond) {
      is_converged = true;
      break;
    }
  }

  if (!is_converged) {
    vout.crucial(m_vl, "Error at %s: not converged.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }


  std::vector<real_t> diffs(Nshift);
  for (int i = 0; i < Nshift; ++i) {
    diffs[i] = 0.0;
  }

  for (int i = 0; i < Nshift; ++i) {
    m_fopr->mult(m_s, m_x[i]);
    axpy(m_s, sigma[i], m_x[i]);
    axpy(m_s, -1.0, b);

    real_t diff1 = sqrt(m_s.norm2() * m_snorm);

    vout.paranoiac(m_vl, "    %4d  %22.15e\n", i, diff1);

    // if (diff1 > diff2) diff2 = diff1;
    diffs[i] = diff1;
  }

#pragma omp barrier
#pragma omp master
  {
    real_t diff2 = -1.0;

    for (int i = 0; i < Nshift; ++i) {
      if (diffs[i] > diff2) diff2 = diffs[i];
    }

    diff = diff2;

    Nconv = Nconv2;
  }
#pragma omp barrier

  for (int i = 0; i < Nshift; ++i) {
    copy(xq[i], m_x[i]);
  }

  vout.paranoiac(m_vl, "   diff(max) = %22.15e  \n", diff);
#pragma omp barrier

}


//====================================================================
template<typename AFIELD, typename AFOPR>
void AShiftsolver_CG<AFIELD, AFOPR>::solve_init(real_t& rr)
{
  int Nshift = m_p.size();

  vout.paranoiac(m_vl, "number of shift = %d\n", Nshift);

  for (int i = 0; i < Nshift; ++i) {
    copy(m_p[i], m_s);
    scal(m_x[i], 0.0);
  }

  copy(m_r, m_s);
  rr = m_r.norm2();

#pragma omp barrier
#pragma omp master
  {
    m_alpha_p = 0.0;
    m_beta_p  = 1.0;
  }
#pragma omp barrier
}


//====================================================================
template<typename AFIELD, typename AFOPR>
void AShiftsolver_CG<AFIELD, AFOPR>::solve_step(real_t& rr)
{
  m_fopr->mult(m_s, m_p[0]);
  axpy(m_s, m_sigma0, m_p[0]);

  real_t rr_p = rr;
  real_t pa_p = dot(m_s, m_p[0]);
  real_t beta = -rr_p / pa_p;

  axpy(m_x[0], -beta, m_p[0]);
  axpy(m_r, beta, m_s);
  rr = m_r.norm2();

  real_t alpha = rr / rr_p;

  aypx(alpha, m_p[0], m_r);

#pragma omp barrier
#pragma omp master
  {
    m_pp[0] = rr;
  }
#pragma omp barrier

  real_t alpha_h = 1.0 + m_alpha_p * beta / m_beta_p;

  for (int ish = 1; ish < m_Nshift2; ++ish) {
    real_t zeta = (alpha_h - m_csh2[ish] * beta) / m_zeta1[ish]
                  + (1.0 - alpha_h) / m_zeta2[ish];
    zeta = 1.0 / zeta;
    real_t zr      = zeta / m_zeta1[ish];
    real_t beta_s  = beta * zr;
    real_t alpha_s = alpha * zr * zr;

    axpy(m_x[ish], -beta_s, m_p[ish]);
    scal(m_p[ish], alpha_s);
    axpy(m_p[ish], zeta, m_r);

    real_t ppr = m_p[ish].norm2();

#pragma omp barrier
#pragma omp master
    {
      m_pp[ish] = ppr * m_snorm;

      m_zeta2[ish] = m_zeta1[ish];
      m_zeta1[ish] = zeta;
    }
#pragma omp barrier
  }

  int ish1 = m_Nshift2;

  for (int ish = m_Nshift2 - 1; ish >= 0; --ish) {
    vout.paranoiac(m_vl, "%4d %16.8e\n", ish, m_pp[ish]);
    if (m_pp[ish] > m_Stop_cond) {
      ish1 = ish + 1;
      break;
    }
  }

#pragma omp barrier
#pragma omp master
  {
    m_Nshift2 = ish1;

    m_alpha_p = alpha;
    m_beta_p  = beta;
  }
#pragma omp barrier
}


//====================================================================
template<typename AFIELD, typename AFOPR>
void AShiftsolver_CG<AFIELD, AFOPR>::reset_field(
                                    const AFIELD& b,
                                    const std::vector<real_t>& sigma,
                                    const int Nshift)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if(ith == 0){
    int Nin  = b.nin();
    int Nvol = b.nvol();
    int Nex  = b.nex();

    m_p.resize(Nshift);
    m_x.resize(Nshift);
    m_zeta1.resize(Nshift);
    m_zeta2.resize(Nshift);
    m_csh2.resize(Nshift);
    m_pp.resize(Nshift);

    for (int i = 0; i < Nshift; ++i) {
      m_p[i].reset(Nin, Nvol, Nex);
      m_x[i].reset(Nin, Nvol, Nex);
      m_zeta1[i] = 1.0;
      m_zeta2[i] = 1.0;
      m_csh2[i]  = sigma[i] - sigma[0];

      m_pp[i] = 0.0;
    }

    m_s.reset(Nin, Nvol, Nex);
    m_r.reset(Nin, Nvol, Nex);

    m_sigma0 = sigma[0];

    m_Nshift2 = Nshift;
  }
#pragma omp barrier
}


//====================================================================
template<typename AFIELD, typename AFOPR>
double AShiftsolver_CG<AFIELD, AFOPR>::flop_count()
{
  vout.general(m_vl, "Warning at %s: flop_count() not yet implemented.\n",
               class_name.c_str());
  return 0.0;
}


//============================================================END=====
