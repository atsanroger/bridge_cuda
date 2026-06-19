/*!
        @file    fopr_Clover_eo.cpp

        @brief

        @author  Satoru Ueda (maintained by h.Matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "fopr_Clover_eo.h"

#include "Fopr/fopr_thread-inc.h"


#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Fopr_Clover_eo::register_factory();
}
#endif


const std::string Fopr_Clover_eo::class_name = "Fopr_Clover_eo";

//====================================================================
void Fopr_Clover_eo::init(const Parameters& params)
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
  }
  if ((m_repr != "Dirac") && (m_repr != "Chiral")) {
    vout.crucial("Error in %s: irrelevant mult mode = %s\n",
                 class_name.c_str(), m_repr.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(params);

  m_fopr_w   = new Fopr_Wilson_eo(params);
  m_fopr_csw = new Fopr_CloverTerm_eo(params);

  m_is_initial_step = false;

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());
}


//====================================================================
void Fopr_Clover_eo::init(const std::string repr)
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();
  m_is_initial_step = true;

  vout.general(m_vl, "%s: construction (obsolete)\n",
               class_name.c_str());
  vout.increase_indent();

  setup();

  m_repr = repr;

  m_fopr_w   = new Fopr_Wilson_eo(m_repr);
  m_fopr_csw = new Fopr_CloverTerm_eo(m_repr);

  m_is_initial_step = false;

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());
}


//====================================================================
void Fopr_Clover_eo::setup()
{
  m_Nc    = CommonParameters::Nc();
  m_Nd    = CommonParameters::Nd();
  m_Ndim  = CommonParameters::Ndim();
  m_NinF  = 2 * m_Nc * m_Nd;
  m_Nvol  = CommonParameters::Nvol();
  m_Nvol2 = m_Nvol / 2;

  m_boundary.resize(m_Ndim);

  m_w1.reset(m_NinF, m_Nvol2, 1);
  m_w2.reset(m_NinF, m_Nvol2, 1);

  m_v1.reset(m_NinF, m_Nvol2, 1);
}


//====================================================================
void Fopr_Clover_eo::tidyup()
{
  delete m_fopr_w;
  delete m_fopr_csw;
}


//====================================================================
void Fopr_Clover_eo::set_parameters(const Parameters& params)
{
#pragma omp barrier
  int         ith = ThreadManager::get_thread_id();
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

  //- propagate parameters to components
  if (!m_is_initial_step) {
    m_fopr_w->set_parameters(params);
    m_fopr_csw->set_parameters(params);
  }
}


//====================================================================
void Fopr_Clover_eo::set_parameters(const double kappa,
                                    const double cSW,
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

  vout.general(m_vl, "%s: input parameters\n", class_name.c_str());
  vout.general(m_vl, "  gamma matrix type = %s\n", m_repr.c_str());
  vout.general(m_vl, "  kappa = %12.8f\n", kappa);
  vout.general(m_vl, "  cSW   = %12.8f\n", cSW);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, bc[mu]);
  }
}


//====================================================================
void Fopr_Clover_eo::get_parameters(Parameters& params) const
{
  params.set_double("hopping_parameter", m_kappa);
  params.set_double("clover_coefficient", m_cSW);
  params.set_int_vector("boundary_condition", m_boundary);
  params.set_string("gamma_matrix_type", m_repr);
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Fopr_Clover_eo::set_config(Field *U)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_U = (Field_G *)U;
#pragma omp barrier

  m_fopr_w->set_config(U);
  m_fopr_csw->set_config(U);
}


//====================================================================
void Fopr_Clover_eo::set_mode(const std::string mode)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_mode = mode;

#pragma omp barrier
}


//====================================================================
void Fopr_Clover_eo::mult(Field& v, const Field& w)
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
void Fopr_Clover_eo::mult_dag(Field& v, const Field& w)
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
void Fopr_Clover_eo::mult_gm5(Field& v, const Field& w)
{
  m_fopr_w->mult_gm5(v, w);
}


//====================================================================
void Fopr_Clover_eo::preProp(Field& Be, Field& bo, const Field& b)
{
  m_idx.convertField(m_w1, b, 1);
  m_fopr_csw->mult_csw_inv(bo, m_w1, 1);

  m_idx.convertField(m_w1, b, 0);
  m_fopr_csw->mult_csw_inv(Be, m_w1, 0);

  if (m_mode == "D") {
    Meo(m_w1, bo, 0);
  } else {
    Mdageo(m_w1, bo, 0);
  }
  axpy(Be, -1.0, m_w1);
#pragma omp barrier
}


//====================================================================
void Fopr_Clover_eo::postProp(Field& x,
                              const Field& xe, const Field& bo)
{
  if (m_mode == "D") {
    Meo(m_w1, xe, 1);
  } else {
    Mdageo(m_w1, xe, 1);
  }

  aypx(-1.0, m_w1, bo);
#pragma omp barrier

  m_idx.reverseField(x, xe, 0);
  m_idx.reverseField(x, m_w1, 1);

#pragma omp barrier
}


//====================================================================
void Fopr_Clover_eo::DdagD(Field& v, const Field& w)
{
  D(m_w1, w);
  Ddag(v, m_w1);
}


//====================================================================
void Fopr_Clover_eo::DDdag(Field& v, const Field& w)
{
  Ddag(m_w1, w);
  D(v, m_w1);
}


//====================================================================
void Fopr_Clover_eo::H(Field& v, const Field& w)
{
  D(m_w1, w);
  mult_gm5(v, m_w1);
}


//====================================================================
void Fopr_Clover_eo::D(Field& v, const Field& w)
{
  Meo(m_w2, w, 1);
  Meo(v, m_w2, 0);
  aypx(-1.0, v, w);

#pragma omp barrier
}


//====================================================================
void Fopr_Clover_eo::Ddag(Field& v, const Field& w)
{
  Mdageo(m_w2, w, 1);
  Mdageo(v, m_w2, 0);
  aypx(-1.0, v, w);

#pragma omp barrier
}


//====================================================================
void Fopr_Clover_eo::Meo(Field& v, const Field& f, const int ieo)
{
  // ieo=0: even <-- odd
  // ieo=1: odd  <-- even

  m_fopr_w->Meo(m_v1, f, ieo);
  m_fopr_csw->mult_csw_inv(v, m_v1, ieo);
}


//====================================================================
void Fopr_Clover_eo::Mdageo(Field& v, const Field& w, const int ieo)
{
  m_fopr_w->mult_gm5(m_v1, w);
  m_fopr_w->Meo(v, m_v1, ieo);
  m_fopr_w->mult_gm5(m_v1, v);

  m_fopr_csw->mult_csw_inv(v, m_v1, ieo);
}


//====================================================================
void Fopr_Clover_eo::mult_isigma(Field_F& w, const Field_F& f,
                                 const int mu, const int nu)
{
  m_fopr_csw->mult_isigma(w, f, mu, nu);
}


//====================================================================
double Fopr_Clover_eo::flop_count()
{
  // Counting of floating point operations in giga unit.
  // defined only for D, Dag, H, DDdag, DdagD which can be called
  // from the solver algorithms.
  // Since the flop_count() of Fopr_Wilson_eo defines flop of
  // (1 - Meo*Moe), flop of clover term is twice added together with
  // contribution of addition.

  // this is for aypx + Meo * 2 with Wilson_eo
  const double gflop_w = m_fopr_w->flop_count();

  // this is for inversion of (1 - clover term)
  const double gflop_csw = m_fopr_csw->flop_count();

  double gflop = gflop_w + 2 * gflop_csw;

  // Additional twice mult of clover term
  if ((m_mode == "DdagD") || (m_mode == "DDdag")) gflop += 2 * gflop_csw;

  return gflop;
}


//====================================================================
//============================================================END=====
