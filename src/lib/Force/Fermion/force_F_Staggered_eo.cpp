/*!
        @file    force_F_Staggered_eo.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-04-17 11:32:37 #$

        @version $LastChangedRevision: 2513 $
*/

#include "force_F_Staggered_eo.h"

const std::string Force_F_Staggered_eo::class_name = "Force_F_Staggered_eo";

//====================================================================
void Force_F_Staggered_eo::init()
{
  m_vl = CommonParameters::Vlevel();

  m_Nvol = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();

  m_boundary.resize(m_Ndim);

  m_fopr_ks = new Fopr_Staggered_eo();

  m_Ueo = new Field_G(m_Nvol, m_Ndim);
}


//====================================================================
void Force_F_Staggered_eo::tidyup()
{
  delete m_fopr_ks;
  delete m_Ueo;
}


//====================================================================
void Force_F_Staggered_eo::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  double           mq;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_double("quark_mass", mq);
  err += params.fetch_int_vector("boundary_condition", bc);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(mq, bc);
}


//====================================================================
void Force_F_Staggered_eo::get_parameters(Parameters& params) const
{
  params.set_double("quark_mass", m_mq);
  params.set_int_vector("boundary_condition", m_boundary);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Force_F_Staggered_eo::set_parameters(const double mq, const std::vector<int> bc)
{
  const int Ndim = CommonParameters::Ndim();

  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  mq   = %12.8f\n", mq);
  for (int mu = 0; mu < Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, bc[mu]);
  }

  //- range check
  int err = 0;
  err += ParameterCheck::non_zero(mq);

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  assert(bc.size() == Ndim);

  //- store values
  m_mq = mq;

  // m_boundary.resize(Ndim);  // already resized in the constructor.
  m_boundary = bc;

  //- post-process
  m_fopr_ks->set_parameters(m_mq, m_boundary);
}


//====================================================================
void Force_F_Staggered_eo::set_config(Field *U)
{
  m_index_eo.convertField(*m_Ueo, *U);
  m_fopr_ks->set_config(U);
}


//====================================================================
void Force_F_Staggered_eo::force_core(Field& force_, const Field& eta_)
{
  const int Nc    = CommonParameters::Nc();
  const int Nvol  = CommonParameters::Nvol();
  const int Ndim  = CommonParameters::Ndim();
  const int Nvol2 = Nvol / 2;

  assert(eta_.nin() == 2 * Nc);
  assert(eta_.nvol() == Nvol2);
  assert(eta_.nex() == 1);

  Field_F_1spinor zeta(Nvol2, 1);
  Field_F_1spinor eta(eta_);
  //m_fopr_ks->Meo(zeta, eta, 1);
  m_fopr_ks->mult(zeta, eta, "Doe");

  Field_G force_e(Nvol2, Ndim);
  force_udiv1(force_e, eta, zeta, 0);

  Field_G force_o(Nvol2, Ndim);
  force_udiv1(force_o, zeta, eta, 1);
  scal(force_o, -1.0); // force_o *= -1.0;

  Field_G force_eo(Nvol, Ndim);
  for (int mu = 0; mu < Ndim; ++mu) {
    for (int site = 0; site < Nvol2; ++site) {
      int gx = m_index_eo.site(site, 0);

      Mat_SU_N ut(Nc);
      ut = m_Ueo->mat(gx, mu) * force_e.mat(site, mu);
      ut.at();
      ut *= -2.0;
      force_eo.set_mat(gx, mu, ut);

      gx = m_index_eo.site(site, 1);
      ut = m_Ueo->mat(gx, mu) * force_o.mat(site, mu);
      ut.at();
      ut *= -2.0;
      force_eo.set_mat(gx, mu, ut);
    }
  }

  Field_G force(Nvol, Ndim);
  m_index_eo.reverseField(force, force_eo);

  copy(force_, force); // force_ = force;
}


//====================================================================
void Force_F_Staggered_eo::force_udiv(Field& force_, const Field& eta_)
{
  const int Nc    = CommonParameters::Nc();
  const int Nvol  = CommonParameters::Nvol();
  const int Ndim  = CommonParameters::Ndim();
  const int Nvol2 = Nvol / 2;

  assert(eta_.nin() == 2 * Nc);
  assert(eta_.nvol() == Nvol2);
  assert(eta_.nex() == 1);

  Field_F_1spinor zeta(Nvol2, 1);
  Field_F_1spinor eta(eta_);
  // m_fopr_ks->Meo(zeta, eta, 1);
  m_fopr_ks->mult(zeta, eta, "Doe");

  Field_G force_e(Nvol2, Ndim);
  force_udiv1(force_e, eta, zeta, 0);

  Field_G force_o(Nvol2, Ndim);
  force_udiv1(force_o, zeta, eta, 1);
  scal(force_o, -1.0); // force_o *= -1.0;

  Field_G force(Nvol, Ndim);
  m_index_eo.reverseField(force, force_e, 0);
  m_index_eo.reverseField(force, force_o, 1);

  copy(force_, force); // force_ = force;
}


//====================================================================
void Force_F_Staggered_eo::force_udiv1(Field_G& force_h,
                                       const Field_F_1spinor& zeta, const Field_F_1spinor& eta, int ieo)
{
  const int Nc    = CommonParameters::Nc();
  const int Nvol  = CommonParameters::Nvol();
  const int Ndim  = CommonParameters::Ndim();
  const int Nvol2 = Nvol / 2;



  for (int mu = 0; mu < Ndim; ++mu) {
    Field_F_1spinor eta2(Nvol2, 1);
    m_shift_eo.backward_h(eta2, eta, m_boundary[mu], mu, ieo);
    m_fopr_ks->mult_staggered_phase(eta2, mu, ieo);

    for (int site = 0; site < Nvol2; ++site) {
      Mat_SU_N ut(Nc);

      for (int c1 = 0; c1 < Nc; ++c1) {
        for (int c2 = 0; c2 < Nc; ++c2) {
          double ut_r = zeta.cmp_r(c2, site) * eta2.cmp_r(c1, site)
                        + zeta.cmp_i(c2, site) * eta2.cmp_i(c1, site);
          double ut_i = zeta.cmp_r(c2, site) * eta2.cmp_i(c1, site)
                        - zeta.cmp_i(c2, site) * eta2.cmp_r(c1, site);
          // ut_r *= -0.5 / m_mq;  // hopping normalization
          // ut_i *= -0.5 / m_mq;
          ut_r *= -0.5;   // mass normalization
          ut_i *= -0.5;

          ut.set(c1, c2, ut_r, ut_i);
        }
      }

      force_h.set_mat(site, mu, ut);
    }
  }
}


//====================================================================
void Force_F_Staggered_eo::force_core1(Field&, const Field&, const Field&)
{
  vout.crucial(m_vl, "Error at %s: unimplemented.\n", __func__);
  exit(EXIT_FAILURE);
}


//====================================================================
void Force_F_Staggered_eo::force_udiv1(Field&, const Field&, const Field&)
{
  vout.crucial(m_vl, "Error at %s: unimplemented.\n", __func__);
  exit(EXIT_FAILURE);
}


//====================================================================
//============================================================END=====
