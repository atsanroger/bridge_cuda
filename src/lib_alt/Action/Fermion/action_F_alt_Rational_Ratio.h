/*!
      @file    action_F_alt_Rational_Ratio.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-06-29 22:43:06 #$
      @version $LastChangedRevision: 2527 $
*/

#ifndef ACTION_F_ALT_RATIONAL_RATIO_INCLUDED
#define ACTION_F_ALT_RATIONAL_RATIO_INCLUDED

#include "lib/Action/action.h"

#include "lib/Fopr/afopr_Rational.h"
#include "lib_alt/Force/Fermion/aforce_F_Rational.h"

#include "lib/IO/bridgeIO.h"
using Bridge::vout;

// include files in alt-code
//#include "lib_alt/Field/afield.h"
//#include "lib_alt/Measurements/Fermion/afprop.h"

//! action class for RHMC, with externally constructed AFopr_Rational.

/*!
    Action class for RHMC that is an alternative to Action_F_Rational
    in the core library.
                                            [05 Feb 2019 H.Matsufuru]
 */

template<typename AFIELD>
class Action_F_alt_Rational_Ratio : public Action
{
 public:
  //typedef AField<double> AFIELD;
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

 private:
  std::string m_label;    // label of action

  AFopr<AFIELD>  *m_fopr1_langev;
  AFopr<AFIELD>  *m_fopr1_H;
  AForce_F<AFIELD> *m_forceF1_MD;

  AFopr<AFIELD>  *m_fopr2_langev;
  AFopr<AFIELD>  *m_fopr2_H;
  AForce_F<AFIELD> *m_forceF2_MD;

  Field *m_U;

  AFIELD m_psf;

  AFIELD m_v1, m_v2, m_v3;

  Bridge::VerboseLevel m_vl;   //!< verbose level

 public:
  //! constructor.
  Action_F_alt_Rational_Ratio(AFopr<AFIELD> *fopr1_langev,
                              AFopr<AFIELD> *fopr1_H,
                              AForce_F<AFIELD> *forceF1_MD,
                              AFopr<AFIELD> *fopr2_langev,
                              AFopr<AFIELD> *fopr2_H,
                              AForce_F<AFIELD> *forceF2_MD)
    : Action(),
      m_fopr1_langev(fopr1_langev), m_fopr1_H(fopr1_H),
      m_forceF1_MD(forceF1_MD),
      m_fopr2_langev(fopr2_langev), m_fopr2_H(fopr2_H),
      m_forceF2_MD(forceF2_MD)
      { init(); }

  //! destructor.
  ~Action_F_alt_Rational_Ratio() { tidyup(); }

  void set_parameters(const Parameters& params);

  void get_parameters(Parameters&) const;

  //! set the label of action.
  void set_label(const std::string label)
  {
    m_label = label;
    vout.detailed(m_vl, "  label: %s\n", m_label.c_str());
  }

  //! returns the label of action.
  std::string get_label()
  { return m_label; }

  //! setting gauge configuration.
  void set_config(Field *U);

  //! Langevin step called at the beginning of HMC.
  double langevin(RandomNumbers *);

  //! calculation of Hamiltonian.
  double calcH();

  //! returns the force for updating conjugate momentum.
  //const Field force();
  void force(Field&);

 private:
  void init();
  void tidyup();

};
#endif
