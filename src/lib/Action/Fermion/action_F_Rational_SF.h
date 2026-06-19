/*!
        @file    action_F_Rational_SF.h

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-01-26 21:43:53 #$

        @version $LastChangedRevision: 2457 $
*/

#ifndef ACTION_FRATIONAL_SF_INCLUDED
#define ACTION_FRATIONAL_SF_INCLUDED

#include "Action/action.h"

#include "Fopr/fopr_Rational.h"
#include "Force/Fermion/force_F_Rational.h"

#include "Field/field_SF.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! action class for RHMC, with externally constructed Fopr_Rational.

/*!
    For the class, Fopr and Force objects are instantiated outside
    the class and specified at the construction.
    This class just provides the framework of rational actions.
                                         [28 Dec 2011 H.Matsufuru]
    unique_ptr is introduced to avoid memory leaks
                                         [21 Mar 2015 Y.Namekawa]
 */


class Action_F_Rational_SF : public Action {
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  std::string m_label;   // label of action

  Fopr *m_fopr_langev;
  Fopr *m_fopr_H;
  Force *m_fopr_force_MD;

  Field *m_U;

  Field m_psf;

 public:
  //! constructor requires pointers to Fopr and Force instances.
  Action_F_Rational_SF(Fopr *fopr_langev, Fopr *fopr_H,
                       Force *fopr_force_MD)
    : m_vl(CommonParameters::Vlevel())
  {
    m_fopr_langev   = fopr_langev;
    m_fopr_H        = fopr_H;
    m_fopr_force_MD = fopr_force_MD;
  }

  Action_F_Rational_SF(Fopr *fopr_langev, Fopr *fopr_H,
                       Force *fopr_force_MD,
                       const Parameters& params)
    : m_vl(CommonParameters::Vlevel())
  {
    m_fopr_langev   = fopr_langev;
    m_fopr_H        = fopr_H;
    m_fopr_force_MD = fopr_force_MD;

    set_parameters(params);
  }

  //! destructor. constructed instances are deconstructed in tydyup().
  ~Action_F_Rational_SF()
  {}

  //! setting parameters and creating class instances.
  void set_parameters(const Parameters& params);

  //! getting parameters
  void get_parameters(Parameters& params) const;

  //! set the label of action.
  void set_label(const std::string label)
  {
    m_label = label;
    vout.detailed(m_vl, "  label: %s\n", m_label.c_str());
  }

  //! returns the label of action.
  std::string get_label()
  {
    return m_label;
  }

  //! setting gauge configuration.
  void set_config(Field *U)
  {
    m_U = U;
    m_fopr_langev->set_config(U);
    m_fopr_H->set_config(U);
    m_fopr_force_MD->set_config(U);
  }

  //! Langevin step called at the beginning of HMC.
  double langevin(RandomNumbers *);

  //! calculation of Hamiltonian.
  double calcH();

  //! returns the force for updating conjugate momentum.
  //const Field force();
  void force(Field&);

 private:
};
#endif
