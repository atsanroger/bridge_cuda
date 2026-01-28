/*!
        @file    force_F_Wilson_TwistedMass_Nf2.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "force_F_Wilson_TwistedMass_Nf2.h"

const std::string Force_F_Wilson_TwistedMass_Nf2::class_name = "Force_F_Wilson_TwistedMass_Nf2";

//====================================================================
void Force_F_Wilson_TwistedMass_Nf2::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  double           kappa, tw_mass;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_double("hopping_parameter", kappa);
  err += params.fetch_double("twisted_mass", tw_mass);
  err += params.fetch_int_vector("boundary_condition", bc);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(kappa, tw_mass, bc);
}


//====================================================================
void Force_F_Wilson_TwistedMass_Nf2::get_parameters(Parameters& params) const
{
  params.set_double("hopping_parameter", m_kappa);
  params.set_double("twisted_mass", m_tw_mass);
  params.set_int_vector("boundary_condition", m_boundary);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Force_F_Wilson_TwistedMass_Nf2::set_parameters(const double kappa, const double tw_mass,
                                                    const std::vector<int> bc)
{
  const int Ndim = CommonParameters::Ndim();

  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  kappa   = %12.8f\n", kappa);
  vout.general(m_vl, "  tw_mass = %12.8f\n", tw_mass);
  for (int mu = 0; mu < Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, bc[mu]);
  }

  //- range check
  // NB. kappa,tw_mass == 0 is allowed.
  assert(bc.size() == Ndim);

  //- store values
  m_kappa   = kappa;
  m_tw_mass = tw_mass;

  m_boundary.resize(Ndim);
  m_boundary = bc;

  //- post-process
  m_fopr_tw->set_parameters(m_kappa, m_tw_mass, m_boundary);
}


//====================================================================
void Force_F_Wilson_TwistedMass_Nf2::force_udiv(Field& force_, const Field& eta_)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  Field_F zeta;
  Field_F eta(eta_);

  m_fopr_tw->set_mode("H");
  m_fopr_tw->mult(zeta, eta);

  Field_G force1(Nvol, Ndim);
  force_udiv1(force1, zeta, eta);
  copy(force_, force1); // force_ = force1;

  force_udiv1(force1, eta, zeta);
  axpy(force_, 1.0, force1); // force_ += force1;
}


//====================================================================
void Force_F_Wilson_TwistedMass_Nf2::force_udiv1(Field& force_, const Field& zeta_, const Field& eta_)
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
void Force_F_Wilson_TwistedMass_Nf2::force_udiv1_impl(Field_G& force, const Field_F& zeta, const Field_F& eta)
{
  const int Ndim = CommonParameters::Ndim();

  force.set(0.0);

  for (int mu = 0; mu < Ndim; ++mu) {
    Field_F eta2;
    m_fopr_tw->mult_gm5p(mu, eta2, eta);

    Field_F eta3;
    mult_Field_Gd(eta3, 0, *m_U, mu, eta2, 0);

    tensorProd_Field_F(force, mu, zeta, eta3);
  }

  scal(force, -m_kappa); // force *= -m_kappa;
}


//====================================================================
//============================================================END=====
