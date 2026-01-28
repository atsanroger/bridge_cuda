/*!
        @file    force_F_Staggered.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "Force/Fermion/force_F_Staggered.h"

using Bridge::vout;

const std::string Force_F_Staggered::class_name = "Force_F_Staggered";

//====================================================================
void Force_F_Staggered::init()
{
  m_boundary.resize(CommonParameters::Ndim());

  m_fopr_ks = new Fopr_Staggered();
}


//====================================================================
void Force_F_Staggered::tidyup()
{
  delete m_fopr_ks;
}


//====================================================================
void Force_F_Staggered::set_parameters(const Parameters& params)
{
  const string vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(vlevel);

  //- fetch and check input parameters
  double           mq;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_double("quark_mass", mq);
  err += params.fetch_int_vector("boundary_condition", bc);

  if (err) {
    vout.crucial(m_vl, "%s: fetch error, input parameter not found.\n", class_name.c_str());
    abort();
  }


  set_parameters(mq, bc);
}


//====================================================================
void Force_F_Staggered::set_parameters(double mq, const std::vector<int> bc)
{
  int Ndim = CommonParameters::Ndim();

  //- print input parameters
  vout.general(m_vl, "Parameters of %s:\n", class_name.c_str());
  vout.general(m_vl, "  mq   = %8.4f\n", mq);
  for (int mu = 0; mu < Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, bc[mu]);
  }

  //- range check
  int err = 0;
  err += ParameterCheck::non_zero(mq);

  if (err) {
    vout.crucial(m_vl, "%s: parameter range check failed.\n",
                 class_name.c_str());
    abort();
  }

  assert(bc.size() == Ndim);

  //- store values
  m_mq = mq;

  // m_boundary.resize(Ndim);  // already resized in the constructor.
  for (int mu = 0; mu < Ndim; ++mu) {
    m_boundary[mu] = bc[mu];
  }

  //- post-process
  m_fopr_ks->set_parameters(m_mq, m_boundary);
}


//====================================================================
void Force_F_Staggered::get_parameters(Parameters& params) const
{
  params.set_double("quark_mass", m_mq);
  params.set_int_vector("boundary_condition", m_boundary);
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}

//====================================================================
void Force_F_Staggered::force_udiv(Field& force_, const Field& eta_)
{
  int Nc   = CommonParameters::Nc();
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  assert(eta_.nin() == 2 * Nc);
  assert(eta_.nvol() == Nvol);
  assert(eta_.nex() == 1);

  Field_G force1(Nvol, Ndim);

  Field_F_1spinor zeta(Nvol, 1);
  Field_F_1spinor eta(eta_);

  m_fopr_ks->H(zeta, eta);

  force_udiv1(force1, eta, zeta);
  copy(force_, force1);

  force_udiv1(force1, zeta, eta);
  axpy(force_, 1.0, force1);
}


//====================================================================
void Force_F_Staggered::force_udiv1(Field& force_,
                                    const Field& eta_,
                                    const Field& zeta_)
{
  int Nc   = CommonParameters::Nc();
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  assert(eta_.nin() == 2 * Nc);
  assert(eta_.nvol() == Nvol);
  assert(eta_.nex() == 1);
  assert(zeta_.nin() == 2 * Nc);
  assert(zeta_.nvol() == Nvol);
  assert(zeta_.nex() == 1);

  Field_G         force(Nvol, Ndim);
  Field_F_1spinor zeta(zeta_);
  Field_F_1spinor eta(eta_);

  force_udiv1(force, zeta, eta);

  copy(force_, force);
}


//====================================================================
void Force_F_Staggered::force_udiv1(Field_G& force,
                                    const Field_F_1spinor& zeta,
                                    const Field_F_1spinor& eta)
{
  int Nc   = CommonParameters::Nc();
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  Field_F_1spinor eta2(Nvol, 1);
  Mat_SU_N        ut(Nc);
  double          utr, uti;

  for (int mu = 0; mu < Ndim; ++mu) {
    m_shift.backward(eta2, eta, m_boundary[mu], mu);
    m_fopr_ks->mult_staggered_phase(eta2, mu);
    m_fopr_ks->mult_gm5(eta2);

    for (int site = 0; site < Nvol; ++site) {
      for (int c1 = 0; c1 < Nc; ++c1) {
        for (int c2 = 0; c2 < Nc; ++c2) {
          utr = zeta.cmp_r(c2, site) * eta2.cmp_r(c1, site)
                + zeta.cmp_i(c2, site) * eta2.cmp_i(c1, site);
          uti = zeta.cmp_r(c2, site) * eta2.cmp_i(c1, site)
                - zeta.cmp_i(c2, site) * eta2.cmp_r(c1, site);
          // hopping normalization
          // utr *= 0.5/m_mq;
          // uti *= 0.5/m_mq;
          // mass normalization
          utr *= 0.5;
          uti *= 0.5;
          ut.set(c1, c2, utr, uti);
        }
      }

      force.set_mat(site, mu, ut);
    }
  }
}


//====================================================================
//============================================================END=====
