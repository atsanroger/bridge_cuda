/*!
        @file    force_F_Clover_SF.cpp

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "force_F_Clover_SF.h"

const std::string Force_F_Clover_SF::class_name = "Force_F_Clover_SF";

//====================================================================
void Force_F_Clover_SF::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  double              kappa, cSW;
  std::vector<int>    bc;
  std::vector<double> phi, phipr;

  int err = 0;
  err += params.fetch_double("hopping_parameter", kappa);
  err += params.fetch_double("clover_coefficient", cSW);
  err += params.fetch_int_vector("boundary_condition", bc);
  err += params.fetch_double_vector("phi", phi);
  err += params.fetch_double_vector("phipr", phipr);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(kappa, cSW, bc, &phi[0], &phipr[0]);
}


//====================================================================
void Force_F_Clover_SF::get_parameters(Parameters& params) const
{
  params.set_double("hopping_parameter", m_kappa);
  params.set_double("clover_coefficient", m_cSW);
  params.set_int_vector("boundary_condition", m_boundary);
  params.set_double_vector("phi", m_phi);
  params.set_double_vector("phipr", m_phipr);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Force_F_Clover_SF::set_parameters(const double kappa, const double cSW, const std::vector<int> bc,
                                       double *phi, double *phipr)
{
  const int Ndim = CommonParameters::Ndim();

  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  kappa = %12.8f\n", kappa);
  vout.general(m_vl, "  cSW   = %12.8f\n", cSW);
  for (int mu = 0; mu < Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, bc[mu]);
  }

  vout.general(m_vl, "  phi1  = %12.6f\n", phi[0]);
  vout.general(m_vl, "  phi2  = %12.6f\n", phi[1]);
  vout.general(m_vl, "  phi3  = %12.6f\n", phi[2]);
  vout.general(m_vl, "  phipr1= %12.6f\n", phipr[0]);
  vout.general(m_vl, "  phipr2= %12.6f\n", phipr[1]);
  vout.general(m_vl, "  phipr3= %12.6f\n", phipr[2]);

  //- range check
  // NB. kappa,cSW == 0 is allowed.
  assert(bc.size() == Ndim);
  // NB. phi,phipr == 0 is allowed.

  //- store values
  m_kappa = kappa;
  m_cSW   = cSW;

  m_boundary.resize(Ndim);
  m_boundary = bc;

  m_phi.resize(3);
  m_phipr.resize(3);
  for (int i = 0; i < 3; ++i) {
    m_phi[i]   = phi[i];
    m_phipr[i] = phipr[i];
  }

  //- propagate parameters
  m_fopr_c->set_parameters(m_kappa, m_cSW, m_boundary, m_phi, m_phipr);
  m_force_w->set_parameters(m_kappa, m_boundary);

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
  //  set_wk.set_parameters(m_phi, m_phipr);
}


//====================================================================
void Force_F_Clover_SF::force_udiv(Field& force_, const Field& eta_)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  Field_F zeta;
  Field_F eta(eta_);

  m_fopr_c->set_mode("H");
  m_fopr_c->mult(zeta, eta);

  Field_G force(Nvol, Ndim);
  force_udiv1_impl(force, eta, zeta);

  Field_G force1(Nvol, Ndim);
  force_udiv1_impl(force1, zeta, eta);

  axpy(force, 1.0, force1); // force += force1;

  Field_SF::set_boundary_spatial_link_zero(force);

  copy(force_, force); // force_ = force;
}


//====================================================================
void Force_F_Clover_SF::force_udiv1(Field& force_, const Field& zeta_, const Field& eta_)
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
void Force_F_Clover_SF::force_udiv1_impl(Field_G& force, const Field_F& zeta, const Field_F& eta)
{
  const int Nc   = CommonParameters::Nc();
  const int Nd   = CommonParameters::Nd();
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  force.set(0.0);

  Field_F zeta1(zeta);
  set_boundary_zero(zeta1);

  m_force_w->force_udiv1(force, zeta, eta);

  Field_F eta2;
  m_fopr_c->mult_gm5(eta2, eta);
  set_boundary_zero(eta2);

  for (int mu = 0; mu < Ndim; ++mu) {
    for (int nu = 0; nu < Ndim; ++nu) {
      if (nu == mu) continue;

      Field_F eta3;
      m_fopr_c->mult_isigma(eta3, eta2, mu, nu);

      Field_G Umu;
      Umu.setpart_ex(0, *m_U, mu);

      Field_G Unu;
      Unu.setpart_ex(0, *m_U, nu);
      if (mu != 3) Field_SF::set_boundary_wk(Umu, m_wk);
      if (nu != 3) Field_SF::set_boundary_wk(Unu, m_wk);

      const int ex = 0;

      // R(1) and R(5)
      Field_F vt1;
      mult_Field_Gd(vt1, 0, *m_Cud, index_dir(mu, nu), eta3, ex);

      Field_G force1;
      tensorProd_Field_F(force1, zeta1, vt1);

      Field_G force2;
      copy(force2, force1); // force2 = force1;

      // R(2)
      Field_F vt3;
      mult_Field_Gd(vt3, 0, Umu, 0, eta3, ex);

      ShiftField_lex shift;
      shift.backward(vt1, vt3, nu);

      Field_F vt2;
      shift.backward(vt2, zeta1, nu);

      Field_G Utmp;
      shift.backward(Utmp, Unu, mu);
      if (mu == 3) Field_SF::set_boundary_wkpr(Utmp, m_wkpr);

      mult_Field_Gn(vt3, 0, Utmp, 0, vt1, ex);

      Field_F vt4;
      mult_Field_Gn(vt4, 0, Unu, 0, vt2, ex);
      tensorProd_Field_F(force1, vt4, vt3);
      axpy(force2, 1.0, force1); // force2 += force1;

      // R(4) and R(8)
      shift.backward(vt1, eta3, mu);

      Field_F zeta_mu;
      shift.backward(zeta_mu, zeta1, mu);
      mult_Field_Gn(vt4, 0, *m_Cud, index_dir(mu, nu), zeta_mu, ex);
      tensorProd_Field_F(force1, vt4, vt1);
      axpy(force2, 1.0, force1); // force2 += force1;

      // R(3)
      shift.backward(vt1, eta3, nu);
      mult_Field_Gn(vt3, 0, Unu, 0, vt1, ex);
      mult_Field_Gn(vt4, 0, Umu, 0, zeta_mu, ex);
      shift.backward(vt1, vt3, mu);
      shift.backward(vt2, vt4, nu);
      mult_Field_Gn(vt4, 0, Unu, 0, vt2, ex);
      tensorProd_Field_F(force1, vt4, vt1);
      axpy(force2, 1.0, force1); // force2 += force1;

      // R(6)
      shift.backward(Utmp, Unu, mu);
      if (mu == 3) Field_SF::set_boundary_wkpr(Utmp, m_wkpr);

      Field_G Utmp2;
      mult_Field_Gdd(Utmp2, 0, Utmp, 0, Umu, 0);
      mult_Field_Gn(vt1, 0, Utmp2, 0, eta3, ex);
      mult_Field_Gd(vt2, 0, Unu, 0, zeta1, ex);
      shift.forward(vt3, vt1, nu);
      shift.forward(vt4, vt2, nu);
      tensorProd_Field_F(force1, vt4, vt3);
      axpy(force2, -1.0, force1); // force2 -= force1;

      // R(7)
      mult_Field_Gd(vt1, 0, Unu, 0, eta3, ex);
      mult_Field_Gn(vt2, 0, Umu, 0, zeta_mu, ex);
      shift.backward(vt3, vt1, mu);
      shift.forward(vt1, vt3, nu);
      mult_Field_Gd(vt4, 0, Unu, 0, vt2, ex);
      shift.forward(vt2, vt4, nu);
      tensorProd_Field_F(force1, vt2, vt1);
      axpy(force2, -1.0, force1);           // force2 -= force1;

      scal(force2, -m_kappa * m_cSW / 8.0); // force2 *= -m_kappa * m_cSW / 8.0;
      force.addpart_ex(mu, force2, 0);
    }
  }
}


//====================================================================
void Force_F_Clover_SF::set_component()
{
  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;

      Staple_SF staple;
      staple.set_parameters(m_phi, m_phipr);

      Field_G Cmu_ud1;
      staple.upper(Cmu_ud1, *m_U, mu, nu);

      Field_G Cmu_ud2;
      staple.lower(Cmu_ud2, *m_U, mu, nu);

      axpy(Cmu_ud1, -1.0, Cmu_ud2);
      m_Cud->setpart_ex(index_dir(mu, nu), Cmu_ud1, 0);
    }
  }
}


//====================================================================
void Force_F_Clover_SF::set_boundary_zero(Field_F& f)
{
#pragma omp barrier

  int Nc   = CommonParameters::Nc();
  int Nd   = CommonParameters::Nd();
  int Nvol = CommonParameters::Nvol();
  int Svol = Nvol / CommonParameters::Nt();
  int Nc2  = 2 * Nc;

  if (Communicator::ipe(3) == 0) {
    for (int site = 0; site < Svol; ++site) {
      for (int s = 0; s < Nd; ++s) {
        for (int cc = 0; cc < Nc2; ++cc) {
          f.set(cc + Nc2 * s, site, 0, 0.0);
        }
      }
    }
  }

#pragma omp barrier
}


//============================================================END=====
