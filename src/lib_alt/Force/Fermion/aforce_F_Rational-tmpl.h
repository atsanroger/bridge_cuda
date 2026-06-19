/*!
        @file    aforce_F_Rational-tmpl.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-10-17 13:51:32 #$

        @version $LastChangedRevision: 2550 $
*/

#include "aforce_F_Rational.h"

template<typename AFIELD>
const std::string AForce_F_Rational<AFIELD>::class_name
                                             = "AForce_F_Rational";
//====================================================================
template<typename AFIELD>
void AForce_F_Rational<AFIELD>::init(const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  } else {
    m_vl = CommonParameters::Vlevel();
  }

  vout.general(m_vl, "%s: construction\n", class_name.c_str());
  vout.increase_indent();

  set_parameters(params);

  int NinF  = m_fopr->field_nin();
  int NvolF = m_fopr->field_nvol();
  int NexF  = m_fopr->field_nex();

  /*
  Parameters params_solver;
  params_solver.set_int("maximum_number_of_iteration", m_Niter);
  params_solver.set_double("convergence_criterion_squared", double(m_Stop_cond));
  params_solver.set_string("verbose_level", vout.get_verbose_level(m_vl));
  */

  m_solver =
    //new AShiftsolver_CG<AFIELD>(m_fopr, params_solver);
    new AShiftsolver_CG<AFIELD, AFopr<AFIELD> >(m_fopr, m_Niter, m_Stop_cond);

  m_psi.resize(m_Np);
  m_psi2.resize(m_Np);
  for (int i = 0; i < m_Np; ++i) {
    m_psi[i].reset(NinF, NvolF, NexF);
    m_psi2[i].reset(NinF, NvolF, NexF);
  }

  m_v1.reset(NinF, NvolF, NexF);
  m_v2.reset(NinF, NvolF, NexF);


  int Nc   = CommonParameters::Nc();
  int NinG = 2 * Nc * Nc;
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  m_force1.reset(NinG, Nvol, Ndim);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}

//====================================================================
template<typename AFIELD>
void AForce_F_Rational<AFIELD>::tidyup()
{
  delete m_solver;

}

//====================================================================
template<typename AFIELD>
void AForce_F_Rational<AFIELD>::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  int    Np, n_exp, d_exp;
  double x_min, x_max;
  int    Niter;
  double Stop_cond;

  int err = 0;
  err += params.fetch_int("number_of_poles", Np);
  err += params.fetch_int("exponent_numerator", n_exp);
  err += params.fetch_int("exponent_denominator", d_exp);
  err += params.fetch_double("lower_bound", x_min);
  err += params.fetch_double("upper_bound", x_max);
  err += params.fetch_int("maximum_number_of_iteration", Niter);
  err += params.fetch_double("convergence_criterion_squared", Stop_cond);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(Np, n_exp, d_exp, real_t(x_min), real_t(x_max),
                 Niter, real_t(Stop_cond));
}


//====================================================================
template<typename AFIELD>
void AForce_F_Rational<AFIELD>::set_parameters(const int Np,
                                  const int n_exp, const int d_exp,
                                  const real_t x_min, const real_t x_max,
                                  const int Niter, const real_t Stop_cond)
{
  //- range check
  int err = 0;
  err += ParameterCheck::non_zero(Np);
  err += ParameterCheck::non_zero(n_exp);
  err += ParameterCheck::non_zero(d_exp);
  // NB. x_min,x_max=0 is allowed.
  err += ParameterCheck::non_zero(Niter);
  err += ParameterCheck::square_non_zero(Stop_cond);

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) {
    m_Np        = Np;
    m_n_exp     = n_exp;
    m_d_exp     = d_exp;
    m_x_min     = x_min;
    m_x_max     = x_max;
    m_Niter     = Niter;
    m_Stop_cond = Stop_cond;

    m_cl.resize(m_Np);
    m_bl.resize(m_Np);

    //- Rational approximation
    const real_t x_min2 = m_x_min * m_x_min;
    const real_t x_max2 = m_x_max * m_x_max;

    Math_Rational rational;
    rational.set_parameters(m_Np, m_n_exp, m_d_exp, x_min2, x_max2);
    rational.get_parameters(m_a0, m_bl, m_cl);
  }
#pragma omp barrier


  //- print input parameters
  vout.general(m_vl, "%s: parameters:\n", class_name.c_str());
  vout.general(m_vl, "  Np        = %d\n", m_Np);
  vout.general(m_vl, "  n_exp     = %d\n", m_n_exp);
  vout.general(m_vl, "  d_exp     = %d\n", m_d_exp);
  vout.general(m_vl, "  x_min     = %12.8f\n", m_x_min);
  vout.general(m_vl, "  x_max     = %12.8f\n", m_x_max);
  vout.general(m_vl, "  Niter     = %d\n", m_Niter);
  vout.general(m_vl, "  Stop_cond = %8.2e\n", m_Stop_cond);

  vout.general(m_vl, " a0 = %18.14f\n", m_a0);
  for (int i = 0; i < m_Np; i++) {
    vout.general(m_vl, " bl[%d] = %18.14f  cl[%d] = %18.14f\n",
                 i, m_bl[i], i, m_cl[i]);
  }
}


//====================================================================
template<typename AFIELD>
void AForce_F_Rational<AFIELD>::get_parameters(Parameters& params) const
{
  params.set_int("number_of_poles", m_Np);
  params.set_int("exponent_numerator", m_n_exp);
  params.set_int("exponent_denominator", m_d_exp);
  params.set_double("lower_bound", double(m_x_min));
  params.set_double("upper_bound", double(m_x_max));
  params.set_int("maximum_number_of_iteration", m_Niter);
  params.set_double("convergence_criterion_squared", double(m_Stop_cond));

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}

//====================================================================
template<typename AFIELD>
void AForce_F_Rational<AFIELD>::set_config(Field *U)
{
  m_U = (Field_G *)U;
  m_fopr->set_config(U);
  m_forceF->set_config(U);

  AIndex_lex<real_t,AFIELD::IMPL> index_lex;
  convert_gauge(index_lex, m_Ucp, *U);
}

//====================================================================
template<typename AFIELD>
void AForce_F_Rational<AFIELD>::force_udiv(AFIELD& force,
                                           const AFIELD& eta)
{
  force_udiv_impl(force, eta);
}


//====================================================================
template<typename AFIELD>
void AForce_F_Rational<AFIELD>::force_udiv1(AFIELD& force,
                                            const AFIELD& zeta,
                                            const AFIELD& eta)
{
  force_udiv1_impl(force, zeta, eta);
}


//====================================================================
template<typename AFIELD>
void AForce_F_Rational<AFIELD>::force_udiv_impl(AFIELD& force,
                                                const AFIELD& eta)
{
#pragma omp barrier

  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  int NinF  = eta.nin();
  int NvolF = eta.nvol();
  int NexF  = eta.nex();

  int Nc   = CommonParameters::Nc();
  int NinG = 2 * Nc * Nc;

  m_fopr->set_mode("DdagD");

  int   Nconv;
  real_t diff;
  m_solver->solve(m_psi, m_cl, eta, Nconv, diff);
  vout.general(m_vl, "      diff(max) = %22.15e  \n", diff);

  force.set(0.0);
#pragma omp barrier

  for (int i = 0; i < m_Np; ++i) {
    m_forceF->force_udiv(m_force1, m_psi[i]);
    axpy(force, m_bl[i], m_force1);
#pragma omp barrier
  }

}


//====================================================================
template<typename AFIELD>
void AForce_F_Rational<AFIELD>::force_udiv1_impl(AFIELD& force,
                                                const AFIELD& zeta,
                                                const AFIELD& eta)
{
#pragma omp barrier

  vout.general(m_vl, "%s: force_udiv1 start.\n", class_name.c_str());
  vout.increase_indent();

  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  int NinF  = eta.nin();
  int NvolF = eta.nvol();
  int NexF  = eta.nex();

  int Nc   = CommonParameters::Nc();
  int NinG = 2 * Nc * Nc;

  m_fopr->set_mode("DdagD");

  int   Nconv;
  real_t diff;
  m_solver->solve(m_psi, m_cl, eta, Nconv, diff);
  vout.general(m_vl, "shiftsolver: nconv = %d  diff(max) = %22.15e\n",
               Nconv, diff);
  m_solver->solve(m_psi2, m_cl, zeta, Nconv, diff);
  vout.general(m_vl, "shiftsolver: nconv = %d  diff(max) = %22.15e\n",
               Nconv, diff);

  force.set(0.0);
#pragma omp barrier

  m_fopr->set_mode("H");
  for (int i = 0; i < m_Np; ++i) {

    m_fopr->mult(m_v1, m_psi[i]);
    m_forceF->set_mode("Hdag");
    m_forceF->force_udiv1(m_force1, m_psi2[i], m_v1);
    axpy(force, m_bl[i], m_force1);
#pragma omp barrier

    m_fopr->mult(m_v1, m_psi2[i]);
    m_forceF->set_mode("H");
    m_forceF->force_udiv1(m_force1, m_v1, m_psi[i]);
    axpy(force, m_bl[i], m_force1);
#pragma omp barrier

  }

  vout.decrease_indent();

}


//====================================================================
template<typename AFIELD>
void AForce_F_Rational<AFIELD>::force_core1(
                             AFIELD&, const AFIELD&, const AFIELD&)
{
  vout.crucial(m_vl, "Error at %s: not implemented.\n", __func__);
  exit(EXIT_FAILURE);
}


//===========================================================END======
