/*!
        @file    force_F_Wilson_SF.h

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-01-26 21:43:53 #$

        @version $LastChangedRevision: 2457 $
*/

#ifndef FORCE_F_WILSON_SF_INCLUDED
#define FORCE_F_WILSON_SF_INCLUDED

#include "force_F.h"
#include "tensorProd.h"

#include "Field/field_SF.h"
#include "Fopr/fopr_Wilson_SF.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Force for the standard Wilson fermion operator

/*!
    This class calculates the force of the standard Wilson fermion.
    The gamma matrix representation is given as control string
    "Dirac"(default) or "Chiral" at the construction, which is
    used to construct the Fopr_Wilson instance.
                                     [23 Dec 2011 H.Matusfuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.             [14 Nov 2012 Y.Namekawa]
 */


class Force_F_Wilson_SF : public Force
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  // Field_G *m_U;

  double m_kappa;
  std::vector<int> m_boundary;
  Fopr_Wilson_SF *m_fopr_w;
  Field_F m_psf;
  std::string m_repr;

 public:

  DEPRECATED
  Force_F_Wilson_SF() : m_vl(CommonParameters::Vlevel())
  {
    //    m_repr = "Dirac";
    //    m_fopr_w = new Fopr_Wilson_SF(m_repr);
    m_fopr_w = new Fopr_Wilson_SF();
    m_boundary.resize(CommonParameters::Ndim());
  }

  Force_F_Wilson_SF(const Parameters& params)
    : m_vl(CommonParameters::Vlevel())
  {
    std::string repr = params.get_string("gamma_matrix_type");

    // //    m_repr = "Dirac";
    // //    m_fopr_w = new Fopr_Wilson_SF(m_repr);
    // m_fopr_w = new Fopr_Wilson_SF();
    // m_boundary.resize(CommonParameters::Ndim());

    //    m_repr = "Dirac";
    //    m_fopr_w = new Fopr_Wilson_SF(m_repr);
    m_fopr_w = new Fopr_Wilson_SF(params);
    m_boundary.resize(CommonParameters::Ndim());

    set_parameters(params);
  }

  ~Force_F_Wilson_SF()
  {
    delete m_fopr_w;
  }

  void set_parameters(const Parameters& params);
  void set_parameters(const double kappa, const std::vector<int> bc);

  void get_parameters(Parameters& params) const;

  void set_config(Field *U)
  {
    m_U = (Field_G *)U;
    m_fopr_w->set_config(U);
  }

  void force_udiv(Field& force, const Field& eta);
  void force_udiv1(Field& force, const Field& zeta, const Field& eta);

 private:
  void force_udiv1_impl(Field_G& force, const Field_F& zeta, const Field_F& eta);
};
#endif
