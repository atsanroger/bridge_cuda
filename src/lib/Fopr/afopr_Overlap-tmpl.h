/*!
        @file    afopr_Overlap-tmpl.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "Fopr/afopr_Overlap.h"

#include "ResourceManager/threadManager.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = AFopr_Overlap<AFIELD>::register_factory();
}
#endif

template<typename AFIELD>
const std::string AFopr_Overlap<AFIELD>::class_name = "Fopr_Overlap";

//====================================================================
template<typename AFIELD>
void AFopr_Overlap<AFIELD>::init(const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: construction\n", class_name.c_str());
  vout.increase_indent();

  m_is_initial_step = true;

  std::string kernel_type;
  int         err = params.fetch_string("kernel_type", kernel_type);
  if (err > 0) {
    vout.crucial(m_vl, "Error at %s: kernel_type is not specified.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
  m_kernel_type = kernel_type;

  double M0;
  params.fetch_double("domain_wall_height", M0);

  Parameters params_kernel = params;
  double     kappa         = 1.0 / (8.0 - 2.0 * M0);
  params_kernel.set_double("hopping_parameter", kappa);

  m_fopr_w         = AFOPR::New(kernel_type, params_kernel);
  m_kernel_created = true;

  m_Nin  = m_fopr_w->field_nin();
  m_Nvol = m_fopr_w->field_nvol();
  m_Nex  = m_fopr_w->field_nex();

  m_sign = new AFopr_Sign<AFIELD>(m_fopr_w, params);

  const int Ndim = CommonParameters::Ndim();
  m_boundary.resize(Ndim);

  m_w0.reset(m_Nin, m_Nvol, m_Nex);
  m_w1.reset(m_Nin, m_Nvol, m_Nex);
  m_w2.reset(m_Nin, m_Nvol, m_Nex);

  set_parameters(params);

  m_is_initial_step = false;

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());
}


//====================================================================
template<typename AFIELD>
void AFopr_Overlap<AFIELD>::init(AFOPR *fopr, const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: construction\n", class_name.c_str());
  vout.increase_indent();

  m_is_initial_step = true;

  m_fopr_w         = fopr;
  m_kernel_created = false;

  m_kernel_type = "unknown";
  m_repr        = "unknown";

  m_Nin  = m_fopr_w->field_nin();
  m_Nvol = m_fopr_w->field_nvol();
  m_Nex  = m_fopr_w->field_nex();

  m_sign = new AFopr_Sign<AFIELD>(m_fopr_w, params);

  const int Ndim = CommonParameters::Ndim();
  m_boundary.resize(Ndim);

  m_w0.reset(m_Nin, m_Nvol, m_Nex);
  m_w1.reset(m_Nin, m_Nvol, m_Nex);
  m_w2.reset(m_Nin, m_Nvol, m_Nex);

  set_parameters(params);

  m_is_initial_step = false;

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());
}


//====================================================================
template<typename AFIELD>
void AFopr_Overlap<AFIELD>::init()
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: construction (obsolete)\n",
               class_name.c_str());
  vout.increase_indent();

  m_is_initial_step = true;
  m_kernel_created  = false;

  m_kernel_type = "unknown";
  m_repr        = "unknown";

  m_Nin  = m_fopr_w->field_nin();
  m_Nvol = m_fopr_w->field_nvol();
  m_Nex  = m_fopr_w->field_nex();

  m_sign = new AFopr_Sign<AFIELD>(m_fopr_w);

  const int Ndim = CommonParameters::Ndim();
  m_boundary.resize(Ndim);

  m_w0.reset(m_Nin, m_Nvol, m_Nex);
  m_w1.reset(m_Nin, m_Nvol, m_Nex);
  m_w2.reset(m_Nin, m_Nvol, m_Nex);

  m_is_initial_step = false;

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());
}


//====================================================================
template<typename AFIELD>
void AFopr_Overlap<AFIELD>::tidyup()
{
  delete m_sign;
  if (m_kernel_created == true) delete m_fopr_w;
}


//====================================================================
template<typename AFIELD>
void AFopr_Overlap<AFIELD>::set_parameters(const Parameters& params)
{
#pragma omp barrier
  int ith = ThreadManager::get_thread_id();

  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    if (ith == 0) m_vl = vout.set_verbose_level(vlevel);
  }
#pragma omp barrier

  //- fetch and check input parameters
  double           mq, M0, x_min, x_max;
  int              Np;
  std::vector<int> bc;

  int    Niter;
  double Stop_cond;

  int err = 0;
  err += params.fetch_double("quark_mass", mq);
  err += params.fetch_double("domain_wall_height", M0);
  err += params.fetch_int("number_of_poles", Np);
  err += params.fetch_double("lower_bound", x_min);
  err += params.fetch_double("upper_bound", x_max);
  err += params.fetch_int("maximum_number_of_iteration", Niter);
  err += params.fetch_double("convergence_criterion_squared", Stop_cond);
  err += params.fetch_int_vector("boundary_condition", bc);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  std::string repr, kernel;
  if (!params.fetch_string("gamma_matrix_type", repr)) m_repr = repr;
  if (!params.fetch_string("kernel_type", kernel)) m_kernel_type = kernel;

  set_parameters(real_t(mq), real_t(M0),
                 Np, real_t(x_min), real_t(x_max),
                 Niter, real_t(Stop_cond), bc);

  if (!m_is_initial_step) {
    Parameters params_kernel = params;
    double     kappa         = 1.0 / (8.0 - 2.0 * M0);
    params_kernel.set_double("hopping_parameter", kappa);
    m_fopr_w->set_parameters(params_kernel);

    m_sign->set_parameters(params);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Overlap<AFIELD>::set_parameters(const real_t mq,
                                           const real_t M0,
                                           const int Np,
                                           const real_t x_min,
                                           const real_t x_max,
                                           const int Niter,
                                           const real_t Stop_cond,
                                           const std::vector<int> bc)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  //- range check
  int err = 0;
  err += ParameterCheck::non_zero(mq);
  err += ParameterCheck::non_zero(M0);
  err += ParameterCheck::non_zero(Np);
  // NB. x_min,x_max == 0 is allowed.
  err += ParameterCheck::non_zero(Niter);
  err += ParameterCheck::square_non_zero(Stop_cond);

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  const int Ndim = CommonParameters::Ndim();
  assert(bc.size() == Ndim);

  if (ith == 0) {
    m_mq        = mq;
    m_M0        = M0;
    m_Np        = Np;
    m_x_min     = x_min;
    m_x_max     = x_max;
    m_Niter     = Niter;
    m_Stop_cond = Stop_cond;
    m_boundary  = bc;
  }
#pragma omp barrier

  //- propagate parameters
  // m_sign->set_parameters(m_Np, m_x_min, m_x_max, m_Niter, m_Stop_cond);

  //- print input parameters
  vout.general(m_vl, "%s parameters:\n", class_name.c_str());
  vout.general(m_vl, "  mq        = %12.8f\n", m_mq);
  vout.general(m_vl, "  M0        = %12.8f\n", m_M0);
  vout.general(m_vl, "  Np        = %4d\n", m_Np);
  vout.general(m_vl, "  x_min     = %12.8f\n", m_x_min);
  vout.general(m_vl, "  x_max     = %12.8f\n", m_x_max);
  vout.general(m_vl, "  Niter     = %6d\n", m_Niter);
  vout.general(m_vl, "  Stop_cond = %8.2e\n", m_Stop_cond);
  for (int mu = 0; mu < Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Overlap<AFIELD>::get_parameters(Parameters& params) const
{
  params.set_double("quark_mass", double(m_mq));
  params.set_double("domain_wall_height", double(m_M0));
  params.set_int("number_of_poles", m_Np);
  params.set_double("lower_bound", double(m_x_min));
  params.set_double("upper_bound", double(m_x_max));
  params.set_int("maximum_number_of_iteration", m_Niter);
  params.set_double("convergence_criterion_squared", double(m_Stop_cond));
  params.set_int_vector("boundary_condition", m_boundary);
  params.set_string("kernel_type", m_kernel_type);
  params.set_string("gamma_matrix_type", m_repr);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
template<typename AFIELD>
void AFopr_Overlap<AFIELD>::set_config(Field *U)
{
  m_sign->set_config(U);
}


//====================================================================
template<typename AFIELD>
void AFopr_Overlap<AFIELD>::set_lowmodes(const int Nsbt,
                                         std::vector<real_t> *ev,
                                         std::vector<AFIELD> *vk)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  if (ith == 0) {
    m_Nsbt = Nsbt;
    m_ev   = ev;
    m_vk   = vk;
  }

#pragma omp barrier

  m_sign->set_lowmodes(m_Nsbt, m_ev, m_vk);
}


//====================================================================
template<typename AFIELD>
void AFopr_Overlap<AFIELD>::set_mode(const std::string mode)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_mode = mode;

  m_sign->set_mode(mode);

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Overlap<AFIELD>::mult(AFIELD& v, const AFIELD& w)
{
  if (m_mode == "D") {
    D(v, w);
  } else if (m_mode == "Ddag") {
    Ddag(v, w);
  } else if (m_mode == "DdagD") {
    DdagD(v, w);
  } else if (m_mode == "H") {
    H(v, w);
  } else {
    vout.crucial(m_vl, "Error at %s: mode undefined.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Overlap<AFIELD>::mult_dag(AFIELD& v, const AFIELD& w)
{
  if (m_mode == "D") {
    Ddag(v, w);
  } else if (m_mode == "Ddag") {
    D(v, w);
  } else if (m_mode == "DdagD") {
    DdagD(v, w);
  } else if (m_mode == "H") {
    H(v, w);
  } else {
    vout.crucial(m_vl, "Error at %s: mode undefined.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Overlap<AFIELD>::mult_gm5(AFIELD& v, const AFIELD& w)
{
  m_fopr_w->mult_gm5(v, w);
}


//====================================================================
template<typename AFIELD>
void AFopr_Overlap<AFIELD>::DdagD(AFIELD& v, const AFIELD& b)
{
  D(m_w1, b);
  m_fopr_w->mult_gm5(m_w2, m_w1);

  D(m_w1, m_w2);
  m_fopr_w->mult_gm5(v, m_w1);
}


//====================================================================
template<typename AFIELD>
void AFopr_Overlap<AFIELD>::H(AFIELD& v, const AFIELD& b)
{
  D(m_w1, b);
  m_fopr_w->mult_gm5(v, m_w1);
}


//====================================================================
template<typename AFIELD>
void AFopr_Overlap<AFIELD>::Ddag(AFIELD& v, const AFIELD& b)
{
  m_fopr_w->mult_gm5(m_w1, b);
  D(m_w2, m_w1);
  m_fopr_w->mult_gm5(v, m_w2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Overlap<AFIELD>::D(AFIELD& v, const AFIELD& b)
{
  assert(b.check_size(m_Nin, m_Nvol, m_Nex));
  assert(v.check_size(m_Nin, m_Nvol, m_Nex));

#pragma omp barrier

  real_t f1 = m_M0 + 0.5 * m_mq;
  real_t f2 = m_M0 - 0.5 * m_mq;

  m_sign->mult(m_w0, b);

  m_fopr_w->mult_gm5(v, m_w0);

  scal(v, f2);
  axpy(v, f1, b);

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
double AFopr_Overlap<AFIELD>::flop_count()
{
  flop_count(m_mode);
}


//====================================================================
template<typename AFIELD>
double AFopr_Overlap<AFIELD>::flop_count(const std::string mode)
{
  int NPE = CommonParameters::NPE();

  double gflop_sign = m_sign->flop_count();

  double flop_blas = double(m_Nin * m_Nex * 3);
  // scal(1 flop per component) + axpy(2 flops per component)

  double gflop = flop_blas * double(m_Nvol) * double(NPE) * 1.0e-9;

  gflop += gflop_sign;

  if (mode == "DdagD") gflop *= 2.0;

  return gflop;
}


//============================================================END=====
