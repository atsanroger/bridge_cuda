/*!
        @file    smear_APE_SF.cpp

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "smear_APE_SF.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Smear_APE_SF::register_factory();
}
#endif

const std::string Smear_APE_SF::class_name = "Smear_APE_SF";

//====================================================================
void Smear_APE_SF::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  double              rho1;
  std::vector<double> phi, phipr;

  int err = 0;
  err += params.fetch_double("rho_uniform", rho1);
  err += params.fetch_double_vector("phi", phi);
  err += params.fetch_double_vector("phipr", phipr);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(rho1, phi, phipr);
}


//====================================================================
void Smear_APE_SF::get_parameters(Parameters& params) const
{
  std::vector<double> rho(m_rho.size());
  for (size_t i = 0; i < m_rho.size(); ++i) {
    rho[i] = m_rho[i];
  }
  params.set_double_vector("rho", rho);
  params.set_double_vector("phi", m_phi);
  params.set_double_vector("phipr", m_phipr);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Smear_APE_SF::set_parameters(const double rho1,
                                  const std::vector<double>& phi,
                                  const std::vector<double>& phipr)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  rho = %8.4f\n", rho1);

  vout.general(m_vl, "  phi1  = %12.6f\n", phi[0]);
  vout.general(m_vl, "  phi2  = %12.6f\n", phi[1]);
  vout.general(m_vl, "  phi3  = %12.6f\n", phi[2]);
  vout.general(m_vl, "  phipr1= %12.6f\n", phipr[0]);
  vout.general(m_vl, "  phipr2= %12.6f\n", phipr[1]);
  vout.general(m_vl, "  phipr3= %12.6f\n", phipr[2]);

  //- range check
  // NB. rho == 0 is allowed.
  // NB. phi,phipr == 0 is allowed.

  //- store values
  // m_rho.resize(m_Ndim * m_Ndim);  // already resized in init.
  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      m_rho[mu + nu * m_Ndim] = rho1;
    }
  }

  m_phi.resize(3);
  m_phipr.resize(3);
  for (int i = 0; i < 3; ++i) {
    m_phi[i]   = phi[i];
    m_phipr[i] = phipr[i];
  }
  Field_SF::set_boundary_matrix(m_wk, m_phi);
  Field_SF::set_boundary_matrix(m_wkpr, m_phipr);
}


//====================================================================
void Smear_APE_SF::set_parameters(const std::vector<double>& rho,
                                  const std::vector<double>& phi,
                                  const std::vector<double>& phipr)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  rho[%d] = %8.4f\n", mu, rho[mu]);
  }

  vout.general(m_vl, "  phi1  = %12.6f\n", phi[0]);
  vout.general(m_vl, "  phi2  = %12.6f\n", phi[1]);
  vout.general(m_vl, "  phi3  = %12.6f\n", phi[2]);
  vout.general(m_vl, "  phipr1= %12.6f\n", phipr[0]);
  vout.general(m_vl, "  phipr2= %12.6f\n", phipr[1]);
  vout.general(m_vl, "  phipr3= %12.6f\n", phipr[2]);

  //- range check
  // NB. rho == 0 is allowed.
  // NB. phi,phipr == 0 is allowed.
  assert(rho.size() == m_Ndim * m_Ndim);

  // store values
  // m_rho.resize(m_Ndim * m_Ndim);  // already resized in init.
  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      m_rho[mu + nu * m_Ndim] = rho[mu + nu * m_Ndim];
    }
  }

  m_phi.resize(3);
  m_phipr.resize(3);
  for (int i = 0; i < 3; ++i) {
    m_phi[i]   = phi[i];
    m_phipr[i] = phipr[i];
  }

  Field_SF::set_boundary_matrix(m_wk, m_phi);
  Field_SF::set_boundary_matrix(m_wkpr, m_phipr);
}


//====================================================================
void Smear_APE_SF::smear(Field_G& Usmear, const Field_G& U)
{
  const int Nvol = CommonParameters::Nvol();

  assert(U.nvol() == Nvol);
  assert(U.nex() == m_Ndim);
  assert(Usmear.nvol() == Nvol);
  assert(Usmear.nex() == m_Ndim);

  Staple_SF staple;
  staple.set_parameters(m_phi, m_phipr);

  Usmear.set(0.0);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    Field_G c_tmp;
    c_tmp.set(0.0);

    Field_G u_tmp;
    copy(u_tmp, 0, U, mu);

    Field_G u_tmp2;

    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu != mu) {
        double rho = m_rho[mu + m_Ndim * nu];
        staple.upper(u_tmp2, U, mu, nu);
        axpy(c_tmp, 0, rho, u_tmp2, 0);

        staple.lower(u_tmp2, U, mu, nu);
        axpy(c_tmp, 0, rho, u_tmp2, 0);
      }
    }

    double rho0 = m_rho[mu + m_Ndim * mu];
    m_proj->project(u_tmp2, rho0, c_tmp, u_tmp);

    if (mu != 3) Field_SF::set_boundary_wk(u_tmp2, m_wk);
    copy(Usmear, mu, u_tmp2, 0);
  }
}


//============================================================END=====
