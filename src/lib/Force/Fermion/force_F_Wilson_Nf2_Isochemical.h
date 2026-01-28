/*!
        @file    force_F_Wilson_Nf2_Isochemical.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef FORCE_F_WILSON_NF2_ISOCHEMICAL_INCLUDED
#define FORCE_F_WILSON_NF2_ISOCHEMICAL_INCLUDED

#include "force_F.h"
#include "tensorProd.h"

#include "Fopr/fopr_Wilson_Chemical.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Force for the Wilson fermion operator with isospin chemical potential.

/*!
    This class calculates the force of the standard Wilson fermion with
    isospin chemical potential with two flavors.
                                     [24 Aug 2011 H.Matusfuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.             [14 Nov 2012 Y.Namekawa]
 */


class Force_F_Wilson_Nf2_Isochemical : public Force
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  double m_kappa;                      //!< hopping parameter
  double m_mu;                         //!< isospin chemical potential
  double m_exp_mu;                     //!< exp(mu)
  std::vector<int> m_boundary;
  Fopr_Wilson_Chemical *m_fopr_w;
  Field_F m_psf;

  std::string m_repr;
  std::string m_mode;

 public:
  DEPRECATED
  Force_F_Wilson_Nf2_Isochemical()
    : m_vl(CommonParameters::Vlevel())
  {
    m_repr   = "Dirac";
    m_fopr_w = new Fopr_Wilson_Chemical(m_repr);
    m_boundary.resize(CommonParameters::Ndim());
  }

  DEPRECATED
  Force_F_Wilson_Nf2_Isochemical(const std::string repr)
    : m_vl(CommonParameters::Vlevel())
  {
    m_repr   = repr;
    m_fopr_w = new Fopr_Wilson_Chemical(m_repr);
    m_boundary.resize(CommonParameters::Ndim());
  }

  Force_F_Wilson_Nf2_Isochemical(const Parameters& params)
    : m_vl(CommonParameters::Vlevel())
  {
    std::string repr = params.get_string("gamma_matrix_type");

    // m_repr   = repr;
    // m_fopr_w = new Fopr_Wilson_Isochemical(m_repr);
    // m_boundary.resize(CommonParameters::Ndim());
    m_repr   = repr;
    m_fopr_w = new Fopr_Wilson_Chemical(params);
    m_boundary.resize(CommonParameters::Ndim());

    set_parameters(params);
  }

  ~Force_F_Wilson_Nf2_Isochemical()
  {
    delete m_fopr_w;
  }

  void set_parameters(const Parameters& params);

  void set_parameters(const double kappa,
                      const double mu,
                      const std::vector<int> bc);

  void get_parameters(Parameters& params) const;

  void set_config(Field *U)
  {
    m_U = (Field_G *)U;
    m_fopr_w->set_config(U);
  }

  void set_mode(const std::string& mode)
  {
    m_mode = mode;
  }

  void force_udiv(Field& force, const Field& eta);
  void force_udiv1(Field& force, const Field& zeta, const Field& eta);

 private:
  void set_parameters_impl(const double kappa,
                           const double mu,
                           const std::vector<int> bc);

  void force_udiv1_impl(Field_G& force,
                        const Field_F& zeta, const Field_F& eta);
};
#endif
