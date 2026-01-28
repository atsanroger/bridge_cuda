/*!
        @file    afopr_Wilson_TwistedMass-tmpl.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "Fopr/afopr_Wilson_TwistedMass.h"
#include "ResourceManager/threadManager.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = AFopr_Wilson_TwistedMass<AFIELD>::register_factory();
}
#endif

template<typename AFIELD>
const std::string AFopr_Wilson_TwistedMass<AFIELD>::class_name
  = "AFopr_Wilson_TwistedMass";

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_TwistedMass<AFIELD>::init(const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();
  m_is_initial_step = true;

  vout.general(m_vl, "%s: construction\n", class_name.c_str());
  vout.increase_indent();

  m_Nc   = CommonParameters::Nc();
  m_Nd   = CommonParameters::Nd();
  m_NinF = 2 * m_Nc * m_Nd;
  m_Nvol = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();
  m_boundary.resize(m_Ndim);

  std::string kernel_type;
  if (!params.fetch_string("kernel_type", kernel_type)) {
    m_kernel_type = kernel_type;
  } else {
    m_kernel_type = "Wilson";  // default
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

  m_fopr_w = AFopr<AFIELD>::New(m_kernel_type, params);

  vout.decrease_indent();

  m_U = 0;

  m_v2.reset(m_NinF, m_Nvol, 1);

  m_is_initial_step = false;

  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson_TwistedMass<AFIELD>::init(const std::string& repr)
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();
  m_is_initial_step = true;

  vout.general(m_vl, "%s: construction (obsolete)\n",
               class_name.c_str());
  vout.increase_indent();

  m_Nc   = CommonParameters::Nc();
  m_Nd   = CommonParameters::Nd();
  m_NinF = 2 * m_Nc * m_Nd;
  m_Nvol = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();
  m_boundary.resize(m_Ndim);

  m_repr = repr;

  m_U = 0;

  std::string kernel_type = "Wilson";
  m_fopr_w = AFopr<AFIELD>::New(kernel_type, repr);

  m_v2.reset(m_NinF, m_Nvol, 1);

  m_is_initial_step = false;

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson_TwistedMass<AFIELD>::tidyup()
{
  delete m_fopr_w;
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson_TwistedMass<AFIELD>::set_parameters(const Parameters& params)
{
  std::string vlevel;
  int         ith = ThreadManager::get_thread_id();
  if (!params.fetch_string("verbose_level", vlevel)) {
    if (ith == 0) m_vl = vout.set_verbose_level(vlevel);
  }
#pragma omp barrier

  //- fetch and check input parameters
  double           kappa, tw_mass;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_double("hopping_parameter", kappa);
  err += params.fetch_double("twisted_mass", tw_mass);
  err += params.fetch_int_vector("boundary_condition", bc);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(kappa, tw_mass, bc);

  if (!m_is_initial_step) {
    m_fopr_w->set_parameters(params);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson_TwistedMass<AFIELD>::set_parameters(const double kappa,
                                                      const double tw_mass,
                                                      const std::vector<int> bc)
{
  //- range check
  // NB. kappa,tw_mass = 0 is allowed.
  assert(bc.size() == m_Ndim);

  int ith = ThreadManager::get_thread_id();
#pragma omp barrier
  if (ith == 0) {
    m_kappa   = kappa;
    m_tw_mass = tw_mass;
    m_boundary.resize(m_Ndim);
    m_boundary = bc;
  }
#pragma omp barrier

  //- print input parameters
  vout.general(m_vl, "%s: input parameters\n", class_name.c_str());
  vout.general(m_vl, "  kernel_type = %s\n", m_kernel_type.c_str());
  vout.general(m_vl, "  gamma matrix type = %s\n", m_repr.c_str());
  vout.general(m_vl, "  kappa   = %12.8f\n", m_kappa);
  vout.general(m_vl, "  tw_mass = %12.8f\n", m_tw_mass);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson_TwistedMass<AFIELD>::get_parameters(Parameters& params) const
{
  params.set_double("hopping_parameter", m_kappa);
  params.set_double("twisted_mass", m_tw_mass);
  params.set_int_vector("boundary_condition", m_boundary);

  params.set_string("kernel_type", m_kernel_type);
  params.set_string("gamma_matrix_type", m_repr);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson_TwistedMass<AFIELD>::set_config(Field *U)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_U = (Field_G *)U;
#pragma omp barrier

  m_fopr_w->set_config(U);
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson_TwistedMass<AFIELD>::mult(AFIELD& v, const AFIELD& w)
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
    vout.crucial(m_vl, "Error at %s: irrelevant mult mode.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson_TwistedMass<AFIELD>::mult_dag(AFIELD& v, const AFIELD& w)
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
    vout.crucial(m_vl, "Error at %s: irrelevant mult mode.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson_TwistedMass<AFIELD>::mult_gm5(AFIELD& v, const AFIELD& w)
{
  m_fopr_w->mult_gm5(v, w);
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson_TwistedMass<AFIELD>::set_mode(const std::string mode)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_mode = mode;
#pragma omp barrier

  m_fopr_w->set_mode(mode);
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson_TwistedMass<AFIELD>::mult_gm5p(const int mu, AFIELD& v,
                                                 const AFIELD& w)
{
  //  m_fopr_w->mult_gm5p(mu, v, w);

#pragma omp barrier
  m_v2.set(0.0);
#pragma omp barrier

  m_fopr_w->mult_up(mu, m_v2, w);
  m_fopr_w->mult_gm5(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson_TwistedMass<AFIELD>::D(AFIELD& v, const AFIELD& w)
{
  //m_fopr_w->D(v, w);
  m_fopr_w->mult(v, w, "D");

  m_fopr_w->mult_gm5(m_v2, w);
  m_v2.xI();
#pragma omp barrier

  const real_t tm_2K = 2.0 * m_kappa * m_tw_mass;

  axpy(v, tm_2K, m_v2);  // v += tm_2K * (AFIELD)w2;

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson_TwistedMass<AFIELD>::Ddag(AFIELD& v, const AFIELD& w)
{
  // m_fopr_w->Ddag(v, w);
  m_fopr_w->mult(v, w, "Ddag");

  m_fopr_w->mult_gm5(m_v2, w);
  m_v2.xI();
#pragma omp barrier

  const real_t tm_2K = 2.0 * m_kappa * m_tw_mass;

  axpy(v, -tm_2K, m_v2);

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson_TwistedMass<AFIELD>::H(AFIELD& v, const AFIELD& w)
{
  // m_fopr_w->H(v, w);
  m_fopr_w->mult(v, w, "H");

  copy(m_v2, w);
  m_v2.xI();
#pragma omp barrier

  const real_t tm_2K = 2.0 * m_kappa * m_tw_mass;

  axpy(v, tm_2K, m_v2);

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson_TwistedMass<AFIELD>::Hdag(AFIELD& v, const AFIELD& w)
{
  // m_fopr_w->H(v, w);
  m_fopr_w->mult(v, w, "H");

  copy(m_v2, w);
  m_v2.xI();
#pragma omp barrier

  const real_t tm_2K = 2.0 * m_kappa * m_tw_mass;

  axpy(v, -tm_2K, m_v2);

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Wilson_TwistedMass<AFIELD>::DdagD(AFIELD& v, const AFIELD& w)
{
  //m_fopr_w->DdagD(v, w);
  m_fopr_w->mult(v, w, "DdagD");

  const real_t tm_2K = 2.0 * m_kappa * m_tw_mass;

  axpy(v, tm_2K * tm_2K, w);

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
double AFopr_Wilson_TwistedMass<AFIELD>::flop_count()
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
