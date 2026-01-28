/*!
        @file    fopr_Clover.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "fopr_Clover.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Fopr_Clover::register_factory();
}
#endif

const std::string Fopr_Clover::class_name = "Fopr_Clover";

//====================================================================
void Fopr_Clover::init(const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();
  m_is_initial_step = true;

  vout.general(m_vl, "%s: construction\n", class_name.c_str());
  vout.increase_indent();

  setup();

  std::string repr;
  if (!params.fetch_string("gamma_matrix_type", repr)) {
    m_repr = repr;
  } else {
    m_repr = "Dirac";  // default
    vout.general(m_vl, "gamma_matrix_type is not given: defalt = %s\n",
                 m_repr.c_str());
  }
  if ((m_repr != "Dirac") && (m_repr != "Chiral")) {
    vout.crucial("Error in %s: irrelevant mult mode = %s\n",
                 class_name.c_str(), m_repr.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(params);

  m_fopr_w   = new Fopr_Wilson(params);
  m_fopr_csw = new Fopr_CloverTerm(params);

  m_is_initial_step = false;

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());
}


//====================================================================
void Fopr_Clover::init(const std::string repr)
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();
  m_is_initial_step = true;

  vout.general(m_vl, "%s: construction (obsolete)\n",
               class_name.c_str());
  vout.increase_indent();

  setup();

  m_repr = repr;

  m_fopr_w   = new Fopr_Wilson(repr);
  m_fopr_csw = new Fopr_CloverTerm(repr);

  m_is_initial_step = false;

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());
}


//====================================================================
void Fopr_Clover::setup()
{
  m_Nvol = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();
  m_Nc   = CommonParameters::Nc();
  m_Nd   = CommonParameters::Nd();
  m_NinF = 2 * m_Nc * m_Nd;

  m_boundary.resize(m_Ndim);

  m_U = 0;

  m_v1.reset(m_NinF, m_Nvol, 1);
  m_v2.reset(m_NinF, m_Nvol, 1);
}


//====================================================================
void Fopr_Clover::tidyup()
{
  delete m_fopr_w;
  delete m_fopr_csw;
}


//====================================================================
void Fopr_Clover::set_parameters(const Parameters& params)
{
#pragma omp barrier
  int ith = ThreadManager::get_thread_id();

  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    if (ith == 0) m_vl = vout.set_verbose_level(vlevel);
  }
#pragma omp barrier

  //- fetch and check input parameters
  double           kappa, cSW;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_double("hopping_parameter", kappa);
  err += params.fetch_double("clover_coefficient", cSW);
  err += params.fetch_int_vector("boundary_condition", bc);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(kappa, cSW, bc);
}


//====================================================================
void Fopr_Clover::set_parameters(const double kappa, const double cSW,
                                 const std::vector<int> bc)
{
  assert(bc.size() == m_Ndim);

#pragma omp barrier
  int ith = ThreadManager::get_thread_id();
  if (ith == 0) {
    m_kappa    = kappa;
    m_cSW      = cSW;
    m_boundary = bc;
  }
#pragma omp barrier

  //- print input parameters
  vout.general(m_vl, "%s: input parameters\n", class_name.c_str());
  vout.general(m_vl, "  gamma matrix type = %s\n", m_repr.c_str());
  vout.general(m_vl, "  kappa = %12.8f\n", m_kappa);
  vout.general(m_vl, "  cSW   = %12.8f\n", m_cSW);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }

  //- propagate parameters to components
  if (!m_is_initial_step) {
    Parameters params;
    get_parameters(params);
    m_fopr_w->set_parameters(params);
    m_fopr_csw->set_parameters(params);
  }
}


//====================================================================
void Fopr_Clover::get_parameters(Parameters& params) const
{
  params.set_double("hopping_parameter", m_kappa);
  params.set_double("clover_coefficient", m_cSW);
  params.set_int_vector("boundary_condition", m_boundary);
  params.set_string("gamma_matrix_type", m_repr);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Fopr_Clover::set_config(Field *U)
{
#pragma omp barrier
  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_U = (Field_G *)U;
#pragma omp barrier

  m_fopr_w->set_config(U);
  m_fopr_csw->set_config(U);
}


//====================================================================
void Fopr_Clover::set_mode(const std::string mode)
{
#pragma omp barrier
  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_mode = mode;
#pragma omp barrier
}


//====================================================================
void Fopr_Clover::mult(Field& v, const Field& w)
{
  if (m_mode == "D") {
    D(v, w);
  } else if (m_mode == "Ddag") {
    Ddag(v, w);
  } else if (m_mode == "DdagD") {
    DdagD(v, w);
  } else if (m_mode == "DDdag") {
    DDdag(v, w);
  } else if (m_mode == "H") {
    H(v, w);
  } else {
    vout.crucial(m_vl, "Error at %s: undefined mode = %s.\n",
                 class_name.c_str(), m_mode.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void Fopr_Clover::mult_dag(Field& v, const Field& w)
{
  if (m_mode == "D") {
    Ddag(v, w);
  } else if (m_mode == "Ddag") {
    D(v, w);
  } else if (m_mode == "DdagD") {
    DdagD(v, w);
  } else if (m_mode == "DDdag") {
    DDdag(v, w);
  } else if (m_mode == "H") {
    H(v, w);
  } else {
    vout.crucial(m_vl, "Error at %s: undefined mode = %s.\n",
                 class_name.c_str(), m_mode.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void Fopr_Clover::mult_gm5(Field& v, const Field& w)
{
  m_fopr_w->mult_gm5(v, w);
}


//====================================================================
void Fopr_Clover::mult_up(const int mu, Field& v, const Field& w)
{
  m_fopr_w->mult_up(mu, v, w);
}


//====================================================================
void Fopr_Clover::mult_dn(const int mu, Field& v, const Field& w)
{
  m_fopr_w->mult_dn(mu, v, w);
}


//====================================================================
void Fopr_Clover::D(Field& v, const Field& w)
{
  assert(w.nex() == 1);

  m_fopr_w->D(v, w);
  m_fopr_csw->mult_sigmaF(m_v1, w);

  axpy(v, -1.0, m_v1);  //  w -= m_v1;
#pragma omp barrier
}


//====================================================================
void Fopr_Clover::Ddag(Field& v, const Field& w)
{
  mult_gm5(v, w);
  D(m_v2, v);
  mult_gm5(v, m_v2);
}


//====================================================================
void Fopr_Clover::DdagD(Field& v, const Field& w)
{
  D(m_v2, w);
  mult_gm5(v, m_v2);
  D(m_v2, v);
  mult_gm5(v, m_v2);
}


//====================================================================
void Fopr_Clover::DDdag(Field& v, const Field& w)
{
  mult_gm5(m_v2, w);
  D(v, m_v2);
  mult_gm5(m_v2, v);
  D(v, m_v2);
}


//====================================================================
void Fopr_Clover::H(Field& v, const Field& w)
{
  D(m_v2, w);
  mult_gm5(v, m_v2);
}


//====================================================================
void Fopr_Clover::mult_isigma(Field_F& v, const Field_F& w,
                              const int mu, const int nu)
{
  m_fopr_csw->mult_isigma(v, w, mu, nu);
}


//====================================================================
double Fopr_Clover::flop_count()
{
  // Counting of floating point operations in giga unit.
  // defined only for D, Dag, H, DDdag, DdagD which can be called
  // from the solver algorithms.
  // Since the flop_count() of Fopr_Wilson_eo defines flop of
  // (1 - Meo*Moe), flop of clover term is twice added together with
  // contribution of addition.

  const int Nvol = CommonParameters::Nvol();
  const int NPE  = CommonParameters::NPE();

  const double gflop_w = m_fopr_w->flop_count();

  double gflop_csw = m_fopr_csw->flop_count();

  gflop_csw += 2 * m_Nc * m_Nd / Nvol / NPE / 1.0e+9;

  double gflop = gflop_w + gflop_csw;

  //- additional twice mult of clover term
  if ((m_mode == "DdagD") || (m_mode == "DDdag")) gflop += gflop_csw;

  return gflop;
}


//============================================================END=====
