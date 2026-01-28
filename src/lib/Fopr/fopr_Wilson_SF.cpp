/*!
        @file    fopr_Wilson_SF.cpp

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "fopr_Wilson_SF.h"

#include "lib/Field/field_SF.h"
#include "lib/Fopr/fopr_thread-inc.h"


#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Fopr_Wilson_SF::register_factory();
}
#endif

const std::string Fopr_Wilson_SF::class_name = "Fopr_Wilson_SF";

//====================================================================
void Fopr_Wilson_SF::init()
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();

  m_repr = "Dirac";

  vout.general(m_vl, "%s: construction\n", class_name.c_str());

  m_Nvol = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();
  m_Nc   = CommonParameters::Nc();
  m_Nd   = CommonParameters::Nd();
  m_NinF = 2 * m_Nc * m_Nd;

  m_boundary.resize(m_Ndim);
  m_U      = 0;
  m_fopr_w = new Fopr_Wilson;

  m_w1.reset(m_NinF, m_Nvol, 1);
  m_w2.reset(m_NinF, m_Nvol, 1);

  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());
}


//====================================================================
void Fopr_Wilson_SF::tidyup()
{
  delete m_fopr_w;
}


//====================================================================
void Fopr_Wilson_SF::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  double           kappa;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_double("hopping_parameter", kappa);
  err += params.fetch_int_vector("boundary_condition", bc);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(kappa, bc);
}


//====================================================================
void Fopr_Wilson_SF::set_parameters(const double kappa,
                                    const std::vector<int> bc)
{
#pragma omp barrier

  assert(bc.size() == m_Ndim);

  int ith = ThreadManager::get_thread_id();

  //- store values
  if (ith == 0) {
    m_kappa = kappa;

    m_boundary.resize(m_Ndim);
    m_boundary = bc;
  }

#pragma omp barrier

  //- print input parameters
  vout.general(m_vl, "%s: set parameters\n", class_name.c_str());
  vout.general(m_vl, "  kappa = %12.8f\n", m_kappa);
  for (int dir = 0; dir < m_Ndim; ++dir) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", dir, m_boundary[dir]);
  }

  //- propagate parameters

  Parameters params;
  get_parameters(params);
  m_fopr_w->set_parameters(params);
  // m_fopr_w->set_parameters(m_kappa, m_boundary);

#pragma omp barrier
}


//====================================================================
void Fopr_Wilson_SF::get_parameters(Parameters& params) const
{
  params.set_double("hopping_parameter", m_kappa);
  params.set_int_vector("boundary_condition", m_boundary);
  params.set_string("gamma_matrix_type", m_repr);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Fopr_Wilson_SF::set_config(Field *U)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_U = (Field_G *)U;
#pragma omp barrier

  m_fopr_w->set_config(U);
}


//====================================================================
void Fopr_Wilson_SF::mult(Field& v, const Field& f)
{
  if (m_mode == "D") {
    D(v, f);
  } else if (m_mode == "DdagD") {
    DdagD(v, f);
  } else if (m_mode == "Ddag") {
    Ddag(v, f);
  } else if (m_mode == "H") {
    H(v, f);
  } else {
    vout.crucial(m_vl, "Error at %s: mode undefined.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void Fopr_Wilson_SF::mult_dag(Field& v, const Field& f)
{
  if (m_mode == "D") {
    Ddag(v, f);
  } else if (m_mode == "DdagD") {
    DdagD(v, f);
  } else if (m_mode == "Ddag") {
    D(v, f);
  } else if (m_mode == "H") {
    H(v, f);
  } else {
    vout.crucial(m_vl, "Error at %s: mode undefined.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void Fopr_Wilson_SF::mult(Field& v, const Field& f,
                          const std::string mode)
{
  if (mode == "D") {
    D(v, f);
  } else if (mode == "DdagD") {
    DdagD(v, f);
  } else if (mode == "Ddag") {
    Ddag(v, f);
  } else if (mode == "H") {
    H(v, f);
  } else {
    vout.crucial(m_vl, "Error at %s: mode undefined.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void Fopr_Wilson_SF::mult_dag(Field& v, const Field& f,
                              const std::string mode)
{
  if (mode == "D") {
    Ddag(v, f);
  } else if (mode == "DdagD") {
    DdagD(v, f);
  } else if (mode == "Ddag") {
    D(v, f);
  } else if (mode == "H") {
    H(v, f);
  } else {
    vout.crucial(m_vl, "Error at %s: mode undefined.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void Fopr_Wilson_SF::set_mode(const std::string mode)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_mode = mode;

#pragma omp barrier
}


//====================================================================
void Fopr_Wilson_SF::DdagD(Field& w, const Field& f)
{
  D(m_w2, f);
  mult_gm5(w, m_w2);
  D(m_w2, w);
  mult_gm5(w, m_w2);
}


//====================================================================
void Fopr_Wilson_SF::Ddag(Field& w, const Field& f)
{
  mult_gm5(w, f);
  D(m_w2, w);
  mult_gm5(w, m_w2);
}


//====================================================================
void Fopr_Wilson_SF::H(Field& w, const Field& f)
{
  D(m_w2, f);
  mult_gm5(w, m_w2);
}


//====================================================================
void Fopr_Wilson_SF::mult_gm5(Field& v, const Field& w)
{
  m_fopr_w->mult_gm5(v, w);
}


//====================================================================
void Fopr_Wilson_SF::D(Field& v, const Field& w)
{
#pragma omp barrier

  copy(m_w1, w);
  Field_SF::set_boundary_zero(m_w1);

  m_fopr_w->D(v, m_w1);
  Field_SF::set_boundary_zero(v);
}


//====================================================================
void Fopr_Wilson_SF::mult_gm5p(const int mu, Field& v, const Field& w)
{
  ThreadManager::assert_single_thread(class_name);

  //#pragma omp barrier

  copy(m_w2, w);
  //#pragma omp barrier

  Field_SF::set_boundary_zero(m_w2);

  m_fopr_w->mult_gm5p(mu, m_w1, m_w2);

  Field_SF::set_boundary_zero(m_w1);
  copy(v, m_w1);
}


//====================================================================
double Fopr_Wilson_SF::flop_count()
{
  //- Counting of floating point operations in giga unit.
  //  not implemented, yet.

  vout.general(m_vl, "Warning at %s: flop_count() has not been implemented.\n",
               class_name.c_str());

  const double gflop = 0;

  return gflop;
}


//============================================================END=====
