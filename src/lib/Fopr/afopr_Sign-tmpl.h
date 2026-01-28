/*!
        @file    afopr_Sign-tmpl.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "Fopr/afopr_Sign.h"
#include "lib/ResourceManager/threadManager.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = AFopr_Sign<AFIELD>::register_factory();
}
#endif

template<typename AFIELD>
const std::string AFopr_Sign<AFIELD>::class_name = "AFopr_Sign";

//====================================================================
template<typename AFIELD>
void AFopr_Sign<AFIELD>::init(const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: construction\n", class_name.c_str());
  vout.general(m_vl, "  -- Sign function with Zolotarev approximation\n");
  vout.increase_indent();

  m_Nin  = m_fopr->field_nin();
  m_Nvol = m_fopr->field_nvol();
  m_Nex  = m_fopr->field_nex();

  m_Nsbt = 0;
  m_ev   = 0;
  m_vk   = 0;

  set_parameters(params);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());
}


//====================================================================
template<typename AFIELD>
void AFopr_Sign<AFIELD>::init()
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: construction (obsolete)\n", class_name.c_str());
  vout.general(m_vl, "  -- Sign function with Zolotarev approximation\n");

  m_Nin  = m_fopr->field_nin();
  m_Nvol = m_fopr->field_nvol();
  m_Nex  = m_fopr->field_nex();

  m_Nsbt = 0;
  m_ev   = 0;
  m_vk   = 0;

  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());
}


//====================================================================
template<typename AFIELD>
void AFopr_Sign<AFIELD>::tidyup()
{
  delete m_solver;
}


//====================================================================
template<typename AFIELD>
void AFopr_Sign<AFIELD>::set_parameters(const Parameters& params)
{
#pragma omp barrier
  int ith = ThreadManager::get_thread_id();

  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    if (ith == 0) m_vl = vout.set_verbose_level(vlevel);
  }
#pragma omp barrier

  //- fetch and check input parameters
  int    Np;
  double x_min, x_max;
  int    Niter;
  double Stop_cond;

  int err = 0;
  err += params.fetch_int("number_of_poles", Np);
  err += params.fetch_double("lower_bound", x_min);
  err += params.fetch_double("upper_bound", x_max);
  err += params.fetch_int("maximum_number_of_iteration", Niter);
  err += params.fetch_double("convergence_criterion_squared", Stop_cond);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(Np, real_t(x_min), real_t(x_max),
                 Niter, real_t(Stop_cond));
}


//====================================================================
template<typename AFIELD>
void AFopr_Sign<AFIELD>::set_parameters(const int Np,
                                        const real_t x_min,
                                        const real_t x_max,
                                        const int Niter,
                                        const real_t Stop_cond)
{
#pragma omp barrier
  int ith = ThreadManager::get_thread_id();

  //- range check
  int err = 0;
  err += ParameterCheck::non_zero(Np);
  // NB. x_min,x_max == 0 is allowed.
  err += ParameterCheck::non_zero(Niter);
  err += ParameterCheck::square_non_zero(Stop_cond);

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  if (ith == 0) {
    m_Np        = Np;
    m_x_min     = x_min;
    m_x_max     = x_max;
    m_Niter     = Niter;
    m_Stop_cond = Stop_cond;

    m_sigma.resize(m_Np);
    m_cl.resize(2 * m_Np);
    m_bl.resize(m_Np);
  }
#pragma omp barrier

  //- print input parameters
  vout.general(m_vl, "%s: input parameters\n", class_name.c_str());
  vout.general(m_vl, "  Np        = %4d\n", m_Np);
  vout.general(m_vl, "  x_min     = %12.8f\n", m_x_min);
  vout.general(m_vl, "  x_max     = %12.8f\n", m_x_max);
  vout.general(m_vl, "  Niter     = %6d\n", m_Niter);
  vout.general(m_vl, "  Stop_cond = %8.2e\n", m_Stop_cond);

  init_parameters();
}


//====================================================================
template<typename AFIELD>
void AFopr_Sign<AFIELD>::init_parameters()
{
#pragma omp barrier
  int ith = ThreadManager::get_thread_id();

  vout.increase_indent();

  if (ith == 0) {
    m_sigma.resize(m_Np);
    m_cl.resize(2 * m_Np);
    m_bl.resize(m_Np);

    // Zolotarev coefficient defined
    const real_t bmax = m_x_max / m_x_min;

    Math_Sign_Zolotarev sign_func(m_Np, bmax);
    sign_func.get_sign_parameters(m_cl, m_bl);

    for (int i = 0; i < m_Np; i++) {
      m_sigma[i] = m_cl[2 * i] * m_x_min * m_x_min;
    }

    for (int i = 0; i < m_Np; i++) {
      vout.general(m_vl, " %3d %12.4e %12.4e %12.4e\n",
                   i, m_cl[i], m_cl[i + m_Np], m_bl[i]);
    }

    m_xq.resize(m_Np);
    for (int i = 0; i < m_Np; ++i) {
      m_xq[i].reset(m_Nin, m_Nvol, m_Nex);
    }

    m_solver = new AShiftsolver_CG<AFIELD, AFOPR>(m_fopr,
                                                  m_Niter,
                                                  m_Stop_cond);
    // m_solver->set_parameter_verboselevel(m_vl);

    m_w1.reset(m_Nin, m_Nvol, m_Nex);
  }

  vout.decrease_indent();

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Sign<AFIELD>::get_parameters(Parameters& params) const
{
  params.set_int("number_of_poles", m_Np);
  params.set_double("lower_bound", double(m_x_min));
  params.set_double("upper_bound", double(m_x_max));
  params.set_int("maximum_number_of_iteration", m_Niter);
  params.set_double("convergence_criterion_squared", double(m_Stop_cond));

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
template<typename AFIELD>
void AFopr_Sign<AFIELD>::set_config(Field *U)
{
  m_fopr->set_config(U);
}


//====================================================================
template<typename AFIELD>
void AFopr_Sign<AFIELD>::set_mode(const std::string mode)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_mode = mode;

  m_fopr->set_mode(mode);    // this operation is irrlevant since
                             // reset in mult.

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Sign<AFIELD>::set_lowmodes(const int Nsbt,
                                      std::vector<real_t> *ev,
                                      std::vector<AFIELD> *vk)
{
  if ((Nsbt > ev->size()) || (Nsbt > vk->size())) {
    vout.crucial(m_vl, "Error at %s: Nsbt is larger than array size\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  m_Nsbt = Nsbt;
  m_ev   = ev;
  m_vk   = vk;
}


//====================================================================
template<typename AFIELD>
void AFopr_Sign<AFIELD>::mult(AFIELD& v, const AFIELD& b)
{
  assert(b.check_size(m_Nin, m_Nvol, m_Nex));
  assert(v.check_size(m_Nin, m_Nvol, m_Nex));

#pragma omp barrier

  copy(m_w1, b);
  if (m_Nsbt > 0) subtract_lowmodes(m_w1);

  m_fopr->set_mode("DdagD");

  int    Nconv;
  real_t diff;
  m_solver->solve(m_xq, m_sigma, m_w1, Nconv, diff);

  //  scal(v, 0.0);
  v.set(0.0);

  for (int i = 0; i < m_Np; i++) {
    axpy(v, m_bl[i], m_xq[i]);
  }
#pragma omp barrier

  m_fopr->mult(m_w1, v);

  const real_t coeff = m_cl[2 * m_Np - 1] * m_x_min * m_x_min;
  axpy(m_w1, coeff, v);

  m_fopr->set_mode("H");
  m_fopr->mult(v, m_w1);

  scal(v, 1.0 / m_x_min);     // v *= 1/m_x_min;
#pragma omp barrier

  if (m_Nsbt > 0) evaluate_lowmodes(v, b);
}


//====================================================================
template<typename AFIELD>
void AFopr_Sign<AFIELD>::subtract_lowmodes(AFIELD& w)
{
#pragma omp barrier

  for (int k = 0; k < m_Nsbt; ++k) {
    complex_t prod = dotc((*m_vk)[k], w);
    axpy(w, -prod, (*m_vk)[k]);
  }
#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Sign<AFIELD>::evaluate_lowmodes(AFIELD& x, const AFIELD& w)
{
#pragma omp barrier

  for (int k = 0; k < m_Nsbt; ++k) {
    complex_t prod = dotc((*m_vk)[k], w);

    real_t ev  = (*m_ev)[k];
    real_t sgn = ev / fabs(ev);
    prod *= sgn;
    axpy(x, prod, (*m_vk)[k]);
  }
#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
typename AFIELD::real_t AFopr_Sign<AFIELD>::sign_zolotarev(const real_t x)
{
//  cl[2*Np], bl[Np]: coefficients of rational approx.

  real_t x2R = 0.0;

  for (int l = 0; l < m_Np; l++) {
    x2R += m_bl[l] / (x * x + m_cl[2 * l]);
  }
  x2R = x2R * (x * x + m_cl[2 * m_Np - 1]);

  return x * x2R;
}


//====================================================================
template<typename AFIELD>
double AFopr_Sign<AFIELD>::flop_count()
{
  flop_count(m_mode);
}


//====================================================================
template<typename AFIELD>
double AFopr_Sign<AFIELD>::flop_count(const std::string mode)
{
  int NPE = CommonParameters::NPE();

  double gflop_solver = m_solver->flop_count();

  double gflop_fopr = m_fopr->flop_count("DdagD")
                      + m_fopr->flop_count("H");

  double flop_blas = m_Nin * m_Nex * ((m_Np + 1) * 2 + 1);
  // (Np + 1) axpy + 1 scal

  int flop_subt = m_Nsbt * m_Nin * m_Nex * 4 * 2 * 2;
  // for each subt vector, (dotc + complex axpy) * 2 (subt and eval)
  // dotc and axpy respectively amount Nin * 4 * Nex.
  // Note that this counting is for a complex Field.

  double flop_site = double(flop_blas) + double(flop_subt);

  double gflop = flop_site * double(m_Nvol) * double(NPE) * 1.0e-9;

  gflop += gflop_solver + gflop_fopr;

  return gflop;
}


//============================================================END=====
