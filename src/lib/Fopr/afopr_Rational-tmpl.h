/*!
        @file    afopr_Rational-tmpl.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-10-17 13:51:32 #$

        @version $LastChangedRevision: 2550 $
*/

#include "Fopr/afopr_Rational.h"

#include "lib/Parameters/commonParameters.h"
#include "lib/Communicator/communicator.h"
#include "lib/ResourceManager/threadManager.h"

//#ifdef USE_FACTORY_AUTOREGISTER
//namespace {
//  bool init = Fopr_Rational::register_factory();
//}
//#endif

template<typename AFIELD>
const std::string AFopr_Rational<AFIELD>::class_name = "AFopr_Rational";

//====================================================================
template<typename AFIELD>
void AFopr_Rational<AFIELD>::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

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
void AFopr_Rational<AFIELD>::get_parameters(Parameters& params) const
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
void AFopr_Rational<AFIELD>::set_parameters(
                          const int Np,
                          const int n_exp,    const int d_exp,
                          const real_t x_min, const real_t x_max,
                          const int Niter,    const real_t Stop_cond)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Np        = %d\n", Np);
  vout.general(m_vl, "  n_exp     = %d\n", n_exp);
  vout.general(m_vl, "  d_exp     = %d\n", d_exp);
  vout.general(m_vl, "  x_min     = %12.8f\n", x_min);
  vout.general(m_vl, "  x_max     = %12.8f\n", x_max);
  vout.general(m_vl, "  Niter     = %d\n", Niter);
  vout.general(m_vl, "  Stop_cond = %8.2e\n", Stop_cond);

  //- range check
  int err = 0;
  err += ParameterCheck::non_zero(Np);
  err += ParameterCheck::non_zero(n_exp);
  err += ParameterCheck::non_zero(d_exp);
  // NB. x_min,x_max=0 is allowed.
  err += ParameterCheck::non_zero(Niter);
  err += ParameterCheck::square_non_zero(Stop_cond);

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  //- store values
  m_Np        = Np;
  m_n_exp     = n_exp;
  m_d_exp     = d_exp;
  m_x_min     = x_min;
  m_x_max     = x_max;
  m_Niter     = Niter;
  m_Stop_cond = Stop_cond;

  //- post-process
  init_parameters();
}



//====================================================================
template<typename AFIELD>
void AFopr_Rational<AFIELD>::set_config(Field *U)
{
    m_fopr->set_config(U);
}

//====================================================================
template<typename AFIELD>
void AFopr_Rational<AFIELD>::init_parameters()
{
  const int Nin  = m_fopr->field_nin();
  const int Nvol = m_fopr->field_nvol();
  const int Nex  = m_fopr->field_nex();

  const int Nshift = m_Np;

  const double x_min2 = m_x_min * m_x_min;
  const double x_max2 = m_x_max * m_x_max;

  m_cl.resize(m_Np);
  m_bl.resize(m_Np);

  m_xq.resize(m_Np);
  for (int i = 0; i < Nshift; ++i) {
    m_xq[i].reset(Nin, Nvol, Nex);
  }

  m_solver = new AShiftsolver_CG<AFIELD, AFopr<AFIELD> >(
  m_fopr, m_Niter, m_Stop_cond);

  std::vector<double> bl(m_Np), cl(m_Np);
  double a0;

  Math_Rational rational;
  rational.set_parameters(m_Np, m_n_exp, m_d_exp,
                          double(x_min2), double(x_max2));
  rational.get_parameters(a0, bl, cl);

  m_a0 = a0;

  vout.general(m_vl, "  a0     = %18.14f\n", m_a0);
  for (int i = 0; i < m_Np; i++) {
    m_bl[i] = real_t(bl[i]);    
    m_cl[i] = real_t(cl[i]);    
    vout.general(m_vl, "  bl[%d] = %18.14f  cl[%d] = %18.14f\n",
                 i, m_bl[i], i, m_cl[i]);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Rational<AFIELD>::mult(AFIELD& v, const AFIELD& b)
{
#pragma omp barrier

  assert(v.nin() == b.nin());
  assert(v.nvol() == b.nvol());
  assert(v.nex() == b.nex());

  vout.general(m_vl, "Shift solver in rational function\n");
  vout.increase_indent();

  vout.detailed(m_vl, "Number of shift values = %d\n", m_Np);

  int    Nconv;
  real_t diff;

  m_fopr->set_mode("DdagD");
  m_solver->solve(m_xq, m_cl, b, Nconv, diff);

  vout.general(m_vl, "diff(max) = %22.15e  \n", diff);

  copy(v, b);
  scal(v, real_t(m_a0));
  for (int i = 0; i < m_Np; i++) {
    axpy(v, real_t(m_bl[i]), m_xq[i]);
  }
  vout.decrease_indent();

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
typename AFIELD::real_t AFopr_Rational<AFIELD>::func(const real_t x)
{
  real_t y = m_a0;

  for (int k = 0; k < m_Np; ++k) {
    y += m_bl[k] / (x + m_cl[k]);
  }

  return y;
}


//============================================================END=====
