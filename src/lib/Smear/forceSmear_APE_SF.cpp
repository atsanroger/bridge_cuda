/*!
        @file    forceSmear_APE_SF.cpp

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "forceSmear_APE_SF.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = ForceSmear_APE_SF::register_factory();
}
#endif

const std::string ForceSmear_APE_SF::class_name = "ForceSmear_APE_SF";

//====================================================================
void ForceSmear_APE_SF::set_parameters(const Parameters& params)
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
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(rho1, phi, phipr);
}


//====================================================================
void ForceSmear_APE_SF::get_parameters(Parameters& params) const
{
  params.set_double_vector("rho", m_rho);
  params.set_double_vector("phi", m_phi);
  params.set_double_vector("phipr", m_phipr);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void ForceSmear_APE_SF::set_parameters(const double rho1,
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

  //- propagate parameters
  //  set_wk.set_parameters(m_phi, m_phipr);

  const int Lx     = CommonParameters::Lx();
  double    Lx_inv = 1.0 / double(Lx);
  double    c0r    = cos(phi[0] * Lx_inv);
  double    c0i    = sin(phi[0] * Lx_inv);
  double    c1r    = cos(phi[1] * Lx_inv);
  double    c1i    = sin(phi[1] * Lx_inv);
  double    c2r    = cos(phi[2] * Lx_inv);
  double    c2i    = sin(phi[2] * Lx_inv);
  m_wk.zero();
  m_wk.set(0, 0, c0r, c0i);
  m_wk.set(1, 1, c1r, c1i);
  m_wk.set(2, 2, c2r, c2i);

  c0r = cos(phipr[0] * Lx_inv);
  c0i = sin(phipr[0] * Lx_inv);
  c1r = cos(phipr[1] * Lx_inv);
  c1i = sin(phipr[1] * Lx_inv);
  c2r = cos(phipr[2] * Lx_inv);
  c2i = sin(phipr[2] * Lx_inv);

  m_wkpr.zero();
  m_wkpr.set(0, 0, c0r, c0i);
  m_wkpr.set(1, 1, c1r, c1i);
  m_wkpr.set(2, 2, c2r, c2i);
}


//====================================================================
void ForceSmear_APE_SF::set_parameters(const std::vector<double>& rho,
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

  // propagate parameters
  //  set_wk.set_parameters(m_phi, m_phipr);

  const int Lx     = CommonParameters::Lx();
  double    Lx_inv = 1.0 / double(Lx);
  double    c0r    = cos(phi[0] * Lx_inv);
  double    c0i    = sin(phi[0] * Lx_inv);
  double    c1r    = cos(phi[1] * Lx_inv);
  double    c1i    = sin(phi[1] * Lx_inv);
  double    c2r    = cos(phi[2] * Lx_inv);
  double    c2i    = sin(phi[2] * Lx_inv);
  m_wk.zero();
  m_wk.set(0, 0, c0r, c0i);
  m_wk.set(1, 1, c1r, c1i);
  m_wk.set(2, 2, c2r, c2i);

  c0r = cos(phipr[0] * Lx_inv);
  c0i = sin(phipr[0] * Lx_inv);
  c1r = cos(phipr[1] * Lx_inv);
  c1i = sin(phipr[1] * Lx_inv);
  c2r = cos(phipr[2] * Lx_inv);
  c2i = sin(phipr[2] * Lx_inv);

  m_wkpr.zero();
  m_wkpr.set(0, 0, c0r, c0i);
  m_wkpr.set(1, 1, c1r, c1i);
  m_wkpr.set(2, 2, c2r, c2i);
}


//====================================================================
void ForceSmear_APE_SF::init()
{
  m_Ndim = CommonParameters::Ndim();
  m_Nvol = CommonParameters::Nvol();

  m_rho.resize(m_Ndim * m_Ndim);
  m_U.resize(m_Ndim);
  m_iTheta.resize(m_Ndim);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      m_rho[mu + nu * m_Ndim] = 0.0;
    }
  }
}


//====================================================================

/*!
<ul>
<li>See the implementation note "note_cloverHMC.pdf" (21 Mar 2012) by H.Matsufuru.
<li>Sigmap \f$=\Sigma_\mu'(x)\f$ in eq.(93) in terms of (k)-th smearing.
<li>U \f$=U_\mu(x)\f$ in (k-1)-th smearing.
<li>Xi is used as \f$\Sigma_\mu'(x)\exp(iQ_\mu(x))+C_\mu^\dagger i\Lambda_\mu(x)\f$.
<li>m_iTheta[mu] is \f$i\Lambda_\mu(x)U_\mu(x)\f$.
<li>The pre-force field Sigma \f$=\Sigma_\mu(x)\f$ is set to zero for the boundary spatial link.
<li>[Y.Taniguchi 2012.04.16]
</ul>
*/

/*
Field ForceSmear_APE_SF::force_udiv(const Field_G& Sigmap, const Field_G& U)
{
  Field_G Sigma(Sigmap.nvol(), Sigmap.nex());

  force_udiv(Sigma, Sigmap, U);

  return Sigma;
}
*/

//====================================================================
void ForceSmear_APE_SF::force_udiv(Field_G& Sigma, const Field_G& Sigmap, const Field_G& U)
{
  const int Nc = CommonParameters::Nc();

  assert(Sigmap.nin() == (2 * Nc * Nc));
  assert(Sigmap.nvol() == m_Nvol);
  assert(Sigmap.nex() == m_Ndim);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    copy(m_U[mu], 0, U, mu);
    if (mu != 3) Field_SF::set_boundary_wk(m_U[mu], m_wk);
  }

  for (int mu = 0; mu < m_Ndim; ++mu) {
    Field_G C;
    C.set(0.0);

    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;
      double rho = m_rho[mu + m_Ndim * nu];

      Field_G c_tmp;
      staple(c_tmp, m_U[mu], m_U[nu], mu, nu);
      axpy(C, 0, rho, c_tmp, 0);
    }

    Field_G sigmap_tmp;
    copy(sigmap_tmp, 0, Sigmap, mu);

    double alpha = m_rho[mu + m_Ndim * mu];

    Field_G Xi;
    m_proj->force_recursive(Xi, m_iTheta[mu],
                            alpha, sigmap_tmp, C, m_U[mu]);
    copy(Sigma, mu, Xi, 0);
  }

  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;

      double rho = m_rho[mu + m_Ndim * nu];

      Field_G sigma_tmp;
      force_each(sigma_tmp, m_U[mu], m_U[nu],
                 m_iTheta[mu], m_iTheta[nu], mu, nu);
      axpy(Sigma, mu, rho, sigma_tmp, 0);
    }
  }
  Field_SF::set_boundary_spatial_link_zero(Sigma);
}


//====================================================================

/*!
  <ul>
  <li>Each block evaluates:
  <li>The 1st block: \f$-U_\nu(x+\mu)U_\mu^\dagger(x+\nu)U_\nu^\dagger(x)\Lambda_\nu(x)\f$
  <li>The 2nd block: \f$-U_\nu^\dagger(x+\mu-\nu)U_\mu^\dagger(x-\nu)\Lambda_\mu(x-\nu)U_\nu(x-\nu)\f$
  <li>The 3rd block: \f$U_\nu^\dagger(x+\mu-\nu)U_\mu^\dagger(x-\nu)\Lambda_\nu(x-\nu)U_\nu(x-\nu)\f$
  <li>The 4th block: \f$\Lambda_\nu(x+\mu)U_\nu(x+\mu)U_\mu^\dagger(x+\nu)U_\nu^\dagger(x)\f$
  <li>The 5th block: \f$-U_\nu^\dagger(x+\mu-\nu)\Lambda_\mu(x+\mu-\nu)U_\mu^\dagger(x-\nu)U_\nu(x-\nu)\f$
  <li>The 6th block: \f$-U_\nu(x+\mu)U_\mu^\dagger(x+\nu)\Lambda_\mu(x+\nu)U_\nu^\dagger(x)\f$
  <li>
  </ul>
*/
void ForceSmear_APE_SF::force_each(Field_G& Sigma_mu,
                                   const Field_G& V_mu, const Field_G& V_nu,
                                   const Field_G& iTheta_mu,
                                   const Field_G& iTheta_nu,
                                   const int mu, const int nu)
{
  Sigma_mu.set(0.0);

  //- The 1st block
  Field_G vt1;
  m_shift.backward(vt1, V_nu, mu);
  if (mu == 3) Field_SF::set_boundary_wkpr(vt1, m_wkpr);

  Field_G vt2;
  m_shift.backward(vt2, V_mu, nu);
  if (nu == 3) Field_SF::set_boundary_wkpr(vt2, m_wkpr);

  Field_G vt3;
  mult_Field_Gnd(vt3, 0, vt1, 0, vt2, 0);
  multadd_Field_Gnd(Sigma_mu, 0, vt3, 0, iTheta_nu, 0, 1.0);

  //- The 2nd block
  mult_Field_Gdn(vt3, 0, iTheta_mu, 0, V_nu, 0);
  mult_Field_Gdn(vt2, 0, vt1, 0, vt3, 0);
  m_shift.forward(vt3, vt2, nu);
  axpy(Sigma_mu, 1.0, vt3); // Sigma_mu += vt3;

  //- The 3rd block
  mult_Field_Gdn(vt3, 0, V_mu, 0, iTheta_nu, 0);
  mult_Field_Gdn(vt2, 0, vt1, 0, vt3, 0);
  m_shift.forward(vt3, vt2, nu);
  axpy(Sigma_mu, 1.0, vt3); // Sigma_mu += vt3;

  //- The 4th block
  m_shift.backward(vt1, iTheta_nu, mu);
  m_shift.backward(vt2, V_mu, nu);
  if (nu == 3) Field_SF::set_boundary_wkpr(vt2, m_wkpr);
  mult_Field_Gnd(vt3, 0, vt1, 0, vt2, 0);
  multadd_Field_Gnd(Sigma_mu, 0, vt3, 0, V_nu, 0, 1.0);

  //- The 5th block
  mult_Field_Gdd(vt2, 0, vt1, 0, V_mu, 0);
  mult_Field_Gnn(vt3, 0, vt2, 0, V_nu, 0);
  m_shift.forward(vt2, vt3, nu);
  axpy(Sigma_mu, 1.0, vt2); // Sigma_mu += vt2;

  //- The 6th block
  m_shift.backward(vt1, V_nu, mu);
  if (mu == 3) Field_SF::set_boundary_wkpr(vt1, m_wkpr);
  m_shift.backward(vt2, iTheta_mu, nu);
  mult_Field_Gnd(vt3, 0, vt1, 0, vt2, 0);
  multadd_Field_Gnd(Sigma_mu, 0, vt3, 0, V_nu, 0, 1.0);
}


//====================================================================
void ForceSmear_APE_SF::staple(Field_G& c,
                               const Field_G& u_mu, const Field_G& u_nu,
                               const int mu, const int nu)
{
  //- upper direction
  Field_G v1;

  m_shift.backward(v1, u_mu, nu);
  if (nu == 3) Field_SF::set_boundary_wkpr(v1, m_wkpr);

  Field_G v2;
  mult_Field_Gnn(v2, 0, u_nu, 0, v1, 0);

  m_shift.backward(v1, u_nu, mu);
  if (mu == 3) Field_SF::set_boundary_wkpr(v1, m_wkpr);

  mult_Field_Gnd(c, 0, v2, 0, v1, 0);

  //- lower direction
  m_shift.backward(v2, u_nu, mu);
  if (mu == 3) Field_SF::set_boundary_wkpr(v2, m_wkpr);

  mult_Field_Gnn(v1, 0, u_mu, 0, v2, 0);
  mult_Field_Gdn(v2, 0, u_nu, 0, v1, 0);
  m_shift.forward(v1, v2, nu);
  axpy(c, 0, 1.0, v1, 0);
}


//============================================================END=====
