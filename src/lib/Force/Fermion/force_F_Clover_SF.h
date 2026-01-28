/*!
        @file    force_F_Clover_SF.h

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef FORCE_F_CLOVER_SF_INCLUDED
#define FORCE_F_CLOVER_SF_INCLUDED

#include "force_F_Wilson_SF.h"
#include "tensorProd.h"

#include "Field/field_SF.h"
#include "Fopr/fopr_Clover_SF.h"
#include "Measurements/Gauge/staple_SF.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Force calculation for clover quark action with SF BC.

/*!
    At present, only the Dirac representation for gamma-matrix
    is available.
                                [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.        [14 Nov 2012 Y.Namekawa]
 */


class Force_F_Clover_SF : public Force
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  double m_kappa;                //!< hopping parameter
  double m_cSW;                  //!< clover coefficient
  std::vector<int> m_boundary;   //!< boundary conditions

  int m_Ndim;
  Field_G *m_Cud;                  //!< for force calculation

  Fopr_Clover_SF *m_fopr_c;
  Force_F_Wilson_SF *m_force_w;
  Force_F_Clover_SF *m_imp;

  //! SF boundary condition at t=0
  std::vector<double> m_phi;
  //! SF boundary condition at t=Nt
  std::vector<double> m_phipr;

  //! In order to set the boundary field.
  //  Field_F_SF set_zero;
  //  Field_G_SF set_wk;

  Mat_SU_N m_wk;    //!< SF boundary condition at t=0
  Mat_SU_N m_wkpr;  //!< SF boundary condition at t=Nt

 public:
  Force_F_Clover_SF() : m_vl(CommonParameters::Vlevel()),
    m_wk(CommonParameters::Nc()), m_wkpr(CommonParameters::Nc())
  {
    init();
  }

  Force_F_Clover_SF(const Parameters& params)
    : m_vl(CommonParameters::Vlevel()),
    m_wk(CommonParameters::Nc()), m_wkpr(CommonParameters::Nc())
  {
    init();
    set_parameters(params);
  }

  ~Force_F_Clover_SF()
  {
    tidyup();
  }

 private:
  void init()
  {
    m_fopr_c  = new Fopr_Clover_SF;
    m_force_w = new Force_F_Wilson_SF;

    const int Nvol = CommonParameters::Nvol();
    m_Ndim = CommonParameters::Ndim();

    m_boundary.resize(m_Ndim);
    m_Cud = new Field_G(Nvol, m_Ndim * m_Ndim);
  }

  void tidyup()
  {
    delete m_Cud;
    delete m_force_w;
    delete m_fopr_c;
  }

 public:

  void set_parameters(const Parameters& params);

  //! Setting parameters of clover fermion.
  void set_parameters(const double kappa, const double cSW, const std::vector<int> bc,
                      double *phi, double *phipr);

  void get_parameters(Parameters& params) const;

  //! Setting gauge configuration
  void set_config(Field *U)
  {
    m_U = (Field_G *)U;
    m_fopr_c->set_config(m_U);
    m_force_w->set_config(U);
    set_component();
  }

  //! For recursive calculation of smeared force.
  void force_udiv(Field& force, const Field& eta);

  //! For recursive calculation of smeared force.
  void force_udiv1(Field& force, const Field& zeta, const Field& eta);


 private:
  //! Core implemetation of clover force calculation.
  void force_udiv1_impl(Field_G& force, const Field_F& zeta, const Field_F& eta);

  //! Set building components for force calculation.
  void set_component();

  void set_boundary_zero(Field_F& f);

  int index_dir(const int mu, const int nu)
  {
    return mu + m_Ndim * nu;
  }
};

#endif
