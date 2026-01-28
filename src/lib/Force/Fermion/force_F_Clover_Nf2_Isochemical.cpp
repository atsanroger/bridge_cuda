/*!
        @file    force_F_Clover_Nf2_Isochemical.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "force_F_Clover_Nf2_Isochemical.h"

const std::string Force_F_Clover_Nf2_Isochemical::class_name = "Force_F_Clover_Nf2_Isochemical";

//====================================================================
void Force_F_Clover_Nf2_Isochemical::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  double           kappa, cSW, mu;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_double("hopping_parameter", kappa);
  err += params.fetch_double("clover_coefficient", cSW);
  err += params.fetch_double("chemical_potential", mu);
  err += params.fetch_int_vector("boundary_condition", bc);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(kappa, cSW, mu, bc);

  m_fopr_c->set_parameters(params);
  m_force_w->set_parameters(params);
  m_force_csw->set_parameters(params);
}


//====================================================================
void Force_F_Clover_Nf2_Isochemical::get_parameters(Parameters& params) const
{
  params.set_double("hopping_parameter", m_kappa);
  params.set_double("clover_coefficient", m_cSW);
  params.set_double("chemical_potential", m_mu);
  params.set_int_vector("boundary_condition", m_boundary);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Force_F_Clover_Nf2_Isochemical::set_parameters(
  const double kappa,
  const double cSW,
  const double mu,
  const std::vector<int> bc)
{
  set_parameters_impl(kappa, cSW, mu, bc);

  Parameters params;
  get_parameters(params);
  m_fopr_c->set_parameters(params);
  m_force_w->set_parameters(params);
  m_force_csw->set_parameters(params);
}


//====================================================================
void Force_F_Clover_Nf2_Isochemical::set_parameters_impl(
  const double kappa,
  const double cSW,
  const double mu,
  const std::vector<int> bc)
{
#pragma omp barrier

  const int Ndim = CommonParameters::Ndim();
  assert(bc.size() == Ndim);

  int ith = ThreadManager::get_thread_id();

  if (ith == 0) {
    m_kappa = kappa;
    m_cSW   = cSW;
    m_mu    = mu;

    m_boundary.resize(Ndim);
    m_boundary = bc;
  }
#pragma omp barrier

  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  kappa = %12.8f\n", m_kappa);
  vout.general(m_vl, "  cSW   = %12.8f\n", m_cSW);
  vout.general(m_vl, "  mu    = %12.8f\n", m_mu);
  for (int dir = 0; dir < Ndim; ++dir) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", dir, m_boundary[dir]);
  }
}


//====================================================================
void Force_F_Clover_Nf2_Isochemical::init(std::string repr)
{
  m_repr = repr;

  m_Ndim = CommonParameters::Ndim();

  m_fopr_c    = new Fopr_Clover_Chemical(m_repr);
  m_force_w   = new Force_F_Wilson_Nf2_Isochemical(m_repr);
  m_force_csw = new Force_F_CloverTerm(m_repr);

  m_boundary.resize(CommonParameters::Ndim());
}


//====================================================================
void Force_F_Clover_Nf2_Isochemical::tidyup()
{
  delete m_force_csw;
  delete m_force_w;
  delete m_fopr_c;
}


//====================================================================
void Force_F_Clover_Nf2_Isochemical::force_udiv(Field& force_, const Field& eta_)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  Field_F zeta;
  Field_F eta(eta_);

  m_fopr_c->H(zeta, eta);

  Field_G force(Nvol, Ndim);
  set_mode("H");
  force_udiv1_impl(force, zeta, eta);
  copy(force_, force); // force_ = force;

  set_mode("Hdag");
  force_udiv1_impl(force, eta, zeta);
  axpy(force_, 1.0, force); // force_ += force;
}


//====================================================================
void Force_F_Clover_Nf2_Isochemical::force_udiv1(Field& force_, const Field& zeta_, const Field& eta_)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  Field_G force(Nvol, Ndim);
  Field_F zeta(zeta_);
  Field_F eta(eta_);

  force_udiv1_impl(force, zeta, eta);

  copy(force_, force); // force_ = force;
}


//====================================================================
void Force_F_Clover_Nf2_Isochemical::force_udiv1_impl(Field_G& force, const Field_F& zeta, const Field_F& eta)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  force.set(0.0);

  m_force_w->set_mode(m_mode);
  m_force_w->force_udiv1(force, zeta, eta);

  Field_G force2(Nvol, Ndim);
  m_force_csw->force_udiv1(force2, zeta, eta);

  axpy(force, 1.0, force2); // force += force2;
}


//====================================================================
//============================================================END=====
