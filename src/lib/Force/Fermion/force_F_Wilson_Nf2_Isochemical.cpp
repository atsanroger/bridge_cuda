/*!
        @file    force_F_Wilson_Nf2_Isochemical.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "force_F_Wilson_Nf2_Isochemical.h"

const std::string Force_F_Wilson_Nf2_Isochemical::class_name
  = "Force_F_Wilson_Nf2_Isochemical";

//====================================================================
void Force_F_Wilson_Nf2_Isochemical::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  double           kappa, mu;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_double("hopping_parameter", kappa);
  err += params.fetch_double("chemical_potential", mu);
  err += params.fetch_int_vector("boundary_condition", bc);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters_impl(kappa, mu, bc);

  vout.increase_indent();
  m_fopr_w->set_parameters(params);
  vout.decrease_indent();
}


//====================================================================
void Force_F_Wilson_Nf2_Isochemical::set_parameters(
  const double kappa,
  const double mu,
  const std::vector<int> bc)
{
  set_parameters_impl(kappa, mu, bc);

  Parameters params;
  get_parameters(params);
  m_fopr_w->set_parameters(params);
}


//====================================================================
void Force_F_Wilson_Nf2_Isochemical::set_parameters_impl(
  const double kappa,
  const double mu,
  const std::vector<int> bc)
{
#pragma omp barrier

  const int Ndim = CommonParameters::Ndim();
  assert(bc.size() == Ndim);

  int ith = ThreadManager::get_thread_id();

  if (ith == 0) {
    m_kappa  = kappa;
    m_mu     = mu;
    m_exp_mu = exp(mu);

    m_boundary.resize(Ndim);
    m_boundary = bc;
  }
#pragma omp barrier

  //- print input parameters
  vout.general(m_vl, "%s: set parameters\n", class_name.c_str());
  vout.general(m_vl, "  kappa = %12.8f\n", m_kappa);
  vout.general(m_vl, "  mu    = %12.8f\n", m_mu);
  for (int dir = 0; dir < Ndim; ++dir) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", dir, m_boundary[dir]);
  }
}


//====================================================================
void Force_F_Wilson_Nf2_Isochemical::get_parameters(Parameters& params) const
{
  params.set_double("hopping_parameter", m_kappa);
  params.set_double("chemical_potential", m_mu);
  params.set_int_vector("boundary_condition", m_boundary);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Force_F_Wilson_Nf2_Isochemical::force_udiv(Field& force_, const Field& eta_)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  Field_F zeta;
  Field_F eta(eta_);

  m_fopr_w->set_mode("H");
  m_fopr_w->mult(zeta, eta);

  Field_G force(Nvol, Ndim);
  set_mode("H");
  force_udiv1_impl(force, zeta, eta);
  copy(force_, force); // force_ = force;

  set_mode("Hdag");
  force_udiv1_impl(force, eta, zeta);
  axpy(force_, 1.0, force); // force_ += force;
}


//====================================================================
void Force_F_Wilson_Nf2_Isochemical::force_udiv1(Field& force_, const Field& zeta_, const Field& eta_)
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
void Force_F_Wilson_Nf2_Isochemical::force_udiv1_impl(Field_G& force, const Field_F& zeta, const Field_F& eta)
{
  const int Ndim = CommonParameters::Ndim();

  force.set(0.0);

  for (int dir = 0; dir < Ndim - 1; ++dir) {
    Field_F eta2;
    m_fopr_w->mult_gm5p(dir, eta2, eta);

    Field_F eta3;
    mult_Field_Gd(eta3, 0, *m_U, dir, eta2, 0);
    scal(eta3, -m_kappa); // eta3 *= -m_kappa;

    tensorProd_Field_F(force, dir, zeta, eta3);
  }

  {
    const int dir = Ndim - 1;

    Field_F eta2;
    m_fopr_w->mult_gm5p(dir, eta2, eta);

    Field_F eta3;
    mult_Field_Gd(eta3, 0, *m_U, dir, eta2, 0);

    if (m_mode == "H") {
      scal(eta3, -(m_kappa * m_exp_mu)); // eta3 *= -(m_kappa * m_exp_mu);
    } else if (m_mode == "Hdag") {
      scal(eta3, -(m_kappa / m_exp_mu)); // eta3 *= -(m_kappa / m_exp_mu);
    } else {
      vout.crucial(m_vl, "Error at %s: illegal mode.\n", class_name.c_str());
      exit(EXIT_FAILURE);
    }

    tensorProd_Field_F(force, dir, zeta, eta3);
  }
}


//====================================================================
//============================================================END=====
