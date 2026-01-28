/*!
        @file    force_F_Wilson_TwistedMass_Nf2.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef FORCE_F_TMFWILSON_NF2_INCLUDED
#define FORCE_F_TMFWILSON_NF2_INCLUDED

#include "force_F.h"
#include "tensorProd.h"

#include "Fopr/fopr_Wilson_TwistedMass.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Force for the twisted-mass Wilson fermion operator

/*!
    This class calculates the force of the twisted-mass Wilson
    fermion.
    The gamma matrix representation is given as control string
    "Dirac"(default) or "Chiral" at the construction, which is
    used to construct the Fopr_Wilson instance.
                                     [23 Dec 2011 H.Matusfuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.             [14 Nov 2012 Y.Namekawa]
 */


class Force_F_Wilson_TwistedMass_Nf2 : public Force
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  double m_kappa;
  double m_tw_mass;
  std::vector<int> m_boundary;
  Fopr_Wilson_TwistedMass *m_fopr_tw;
  Field_F m_psf;
  std::string m_repr;

 public:

  DEPRECATED
  Force_F_Wilson_TwistedMass_Nf2()
    : m_vl(CommonParameters::Vlevel())
  {
    m_repr    = "Dirac";
    m_fopr_tw = new Fopr_Wilson_TwistedMass(m_repr);
    m_boundary.resize(CommonParameters::Ndim());
  }

  DEPRECATED
  Force_F_Wilson_TwistedMass_Nf2(const std::string repr)
    : m_vl(CommonParameters::Vlevel())
  {
    m_repr    = repr;
    m_fopr_tw = new Fopr_Wilson_TwistedMass(m_repr);
    m_boundary.resize(CommonParameters::Ndim());
  }

  Force_F_Wilson_TwistedMass_Nf2(const Parameters& params)
    : m_vl(CommonParameters::Vlevel())
  {
    std::string repr = params.get_string("gamma_matrix_type");

    // m_repr    = repr;
    // m_fopr_tw = new Fopr_Wilson_TwistedMass(m_repr);
    // m_boundary.resize(CommonParameters::Ndim());

    m_repr    = repr;
    m_fopr_tw = new Fopr_Wilson_TwistedMass(params);
    m_boundary.resize(CommonParameters::Ndim());

    set_parameters(params);
  }

  ~Force_F_Wilson_TwistedMass_Nf2()
  {
    delete m_fopr_tw;
  }

  void set_parameters(const Parameters& params);

  void set_parameters(const double kappa, const double tw_mass,
                      const std::vector<int> bc);

  void get_parameters(Parameters& params) const;

  void set_config(Field *U)
  {
    m_U = (Field_G *)U;
    m_fopr_tw->set_config(U);
  }

  void force_udiv(Field& force, const Field& eta);
  void force_udiv1(Field& force, const Field& zeta, const Field& eta);

 private:
  void force_udiv1_impl(Field_G& force, const Field_F& zeta, const Field_F& eta);
};
#endif
