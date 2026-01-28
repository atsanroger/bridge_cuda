/*!
        @file    afopr_Clover_Chemical-tmpl.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "Fopr/afopr_Clover_Chemical.h"
#include "ResourceManager/threadManager.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = AFopr_Clover_Chemical<AFIELD>::register_factory();
}
#endif

template<typename AFIELD>
const std::string AFopr_Clover_Chemical<AFIELD>::class_name
  = "AFopr_Clover_Chemical<AFIELD>";

//====================================================================
template<typename AFIELD>
void AFopr_Clover_Chemical<AFIELD>::init(const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();
  m_is_initial_step = true;

  vout.general(m_vl, "%s: construction\n", class_name.c_str());
  vout.increase_indent();

  int Nc = CommonParameters::Nc();
  int Nd = CommonParameters::Nd();
  m_NinF = 2 * Nc * Nd;
  m_Nvol = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();

  m_boundary.resize(m_Ndim);

  std::string kernel_type;
  if (!params.fetch_string("kernel_type", kernel_type)) {
    m_kernel_type = kernel_type;
  } else {
    m_kernel_type = "Clover";  // default
    vout.general("kernel_type is set to default: %s\n",
                 m_kernel_type.c_str());
  }

  std::string repr;
  if (!params.fetch_string("gamma_matrix_type", repr)) {
    m_repr = repr;
  } else {
    m_repr = "Dirac";  // default
  }
  if ((m_repr != "Dirac") && (m_repr != "Chiral")) {
    vout.crucial("Error in %s: irrelevant mult mode = %s\n",
                 class_name.c_str(), m_repr.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(params);

  m_U = 0;

  m_fopr_w = AFopr<AFIELD>::New(m_kernel_type, params);
  m_fopr_w->set_mode("D");

  vout.decrease_indent();

  m_v1.reset(m_NinF, m_Nvol, 1);
  m_v2.reset(m_NinF, m_Nvol, 1);

  m_is_initial_step = false;

  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_Chemical<AFIELD>::init(const std::string repr)
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();
  m_is_initial_step = true;

  vout.general(m_vl, "%s: construction (obsolete)\n",
               class_name.c_str());
  vout.increase_indent();

  int Nc = CommonParameters::Nc();
  int Nd = CommonParameters::Nd();
  m_NinF = 2 * Nc * Nd;
  m_Nvol = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();

  m_boundary.resize(m_Ndim);

  std::string kernel_type;
  m_kernel_type = "Clover";

  m_repr = repr;
  if ((m_repr != "Dirac") && (m_repr != "Chiral")) {
    vout.crucial("Error in %s: irrelevant mult mode = %s\n",
                 class_name.c_str(), m_repr.c_str());
    exit(EXIT_FAILURE);
  }

  m_U = 0;

  m_fopr_w = AFopr<AFIELD>::New(m_kernel_type, m_repr);
  m_fopr_w->set_mode("D");

  vout.decrease_indent();

  m_v1.reset(m_NinF, m_Nvol, 1);
  m_v2.reset(m_NinF, m_Nvol, 1);

  m_is_initial_step = false;

  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_Chemical<AFIELD>::tidyup()
{
  delete m_fopr_w;
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_Chemical<AFIELD>::set_parameters(
  const Parameters& params)
{
#pragma omp barrier
  int         ith = ThreadManager::get_thread_id();
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    if (ith == 0) m_vl = vout.set_verbose_level(vlevel);
  }
#pragma omp barrier

  //- fetch and check input parameters
  double           kappa, cSW, mu;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_double("hopping_parameter", kappa);
  err += params.fetch_double("clover_coefficient", cSW);
  err += params.fetch_double("chemical_potential", mu);
  err += params.fetch_int_vector("boundary_condition", bc);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters_impl(real_t(kappa), real_t(cSW), real_t(mu), bc);

  if (!m_is_initial_step) {
    vout.increase_indent();
    m_fopr_w->set_parameters(params);
    vout.decrease_indent();
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_Chemical<AFIELD>::set_parameters_impl(
  const real_t kappa,
  const real_t cSW,
  const real_t mu,
  const std::vector<int> bc)
{
  assert(bc.size() == m_Ndim);

#pragma omp barrier
  int ith = ThreadManager::get_thread_id();
  if (ith == 0) {
    m_kappa    = kappa;
    m_cSW      = cSW;
    m_mu       = mu;
    m_boundary = bc;
    m_exp_mu   = exp(mu);
  }
#pragma omp barrier

  vout.general(m_vl, "%s: input parameters\n", class_name.c_str());
  vout.general(m_vl, "  kernel_type = %s\n", m_kernel_type.c_str());
  vout.general(m_vl, "  gamma matrix type = %s\n", m_repr.c_str());
  vout.general(m_vl, "  kappa = %12.8f\n", m_kappa);
  vout.general(m_vl, "  cSW   = %12.8f\n", m_cSW);
  vout.general(m_vl, "  mu    = %12.8f\n", m_mu);
  for (int dir = 0; dir < m_Ndim; ++dir) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", dir, m_boundary[dir]);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_Chemical<AFIELD>::get_parameters(
  Parameters& params) const
{
  params.set_double("hopping_parameter", m_kappa);
  params.set_double("clover_coefficient", m_cSW);
  params.set_double("chemical_potential", m_mu);
  params.set_int_vector("boundary_condition", m_boundary);

  params.set_string("kernel_type", m_kernel_type);
  params.set_string("gamma_matrix_type", m_repr);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_Chemical<AFIELD>::set_config(Field *U)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_U = (Field_G *)U;
#pragma omp barrier

  m_fopr_w->set_config(U);
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_Chemical<AFIELD>::set_mode(const std::string mode)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_mode = mode;
#pragma omp barrier

  m_fopr_w->set_mode("D");
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_Chemical<AFIELD>::mult(AFIELD& v, const AFIELD& w)
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
    vout.crucial(m_vl, "Error at %s: irrelevant input mode:\n",
                 class_name.c_str(), m_mode.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_Chemical<AFIELD>::mult_dag(AFIELD& v, const AFIELD& w)
{
  if (m_mode == "D") {
    Ddag(v, w);
  } else if (m_mode == "Ddag") {
    D(v, w);
  } else if (m_mode == "DdagD") {
    DdagD(v, w);
  } else if (m_mode == "H") {
    Hdag(v, w);
  } else {
    vout.crucial(m_vl, "Error at %s: irrelevant input mode:\n",
                 class_name.c_str(), m_mode.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_Chemical<AFIELD>::D(AFIELD& v, const AFIELD& w)
{
  m_fopr_w->mult(v, w);

  m_v1.set(0.0);
#pragma omp barrier

  m_fopr_w->mult_up(3, m_v1, w);
  axpy(v, -m_kappa * (m_exp_mu - 1.0), m_v1);

  m_v1.set(0.0);
#pragma omp barrier

  m_fopr_w->mult_dn(3, m_v1, w);
  axpy(v, -m_kappa * (1.0 / m_exp_mu - 1.0), m_v1);

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_Chemical<AFIELD>::Dminmu(AFIELD& v, const AFIELD& w)
{
  m_fopr_w->mult(v, w);

  m_v1.set(0.0);
#pragma omp barrier

  m_fopr_w->mult_up(3, m_v1, w);
  axpy(v, -m_kappa * (1.0 / m_exp_mu - 1.0), m_v1);

  m_v1.set(0.0);
#pragma omp barrier

  m_fopr_w->mult_dn(3, m_v1, w);
  axpy(v, -m_kappa * (m_exp_mu - 1.0), m_v1);

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_Chemical<AFIELD>::mult_gm5(AFIELD& v, const AFIELD& w)
{
  m_fopr_w->mult_gm5(v, w);
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_Chemical<AFIELD>::mult_gm5p(const int mu,
                                              AFIELD& v, const AFIELD& w)
{
  m_v1.set(0.0);
#pragma omp barrier

  m_fopr_w->mult_up(mu, m_v1, w);
  m_fopr_w->mult_gm5(v, m_v1);
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_Chemical<AFIELD>::DdagD(AFIELD& v, const AFIELD& w)
{
  D(m_v2, w);
  mult_gm5(v, m_v2);
  Dminmu(m_v2, v);
  mult_gm5(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_Chemical<AFIELD>::Ddag(AFIELD& v, const AFIELD& w)
{
  mult_gm5(v, w);
  Dminmu(m_v2, v);
  mult_gm5(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_Chemical<AFIELD>::H(AFIELD& v, const AFIELD& w)
{
  D(m_v2, w);
  mult_gm5(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_Chemical<AFIELD>::Hdag(AFIELD& v, const AFIELD& w)
{
  Dminmu(m_v2, w);
  mult_gm5(v, m_v2);
}


//====================================================================
template<typename AFIELD>
double AFopr_Clover_Chemical<AFIELD>::flop_count()
{
  int NPE = CommonParameters::NPE();

  double gflop_w = m_fopr_w->flop_count();

  double gflop_tm = 2.0 * double(m_NinF) * double(m_Nvol)
                    * double(NPE) * 1.0e-9;

  if ((m_mode == "DdagD") || (m_mode == "DDdag")) gflop_tm *= 2;

  double gflop = gflop_w + gflop_tm;

  return gflop;
}


//============================================================END=====
