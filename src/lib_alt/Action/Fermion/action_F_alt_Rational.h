/*!
      @file    action_F_alt_Rational.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-10-17 13:51:32 #$
      @version $LastChangedRevision: 2550 $
*/

#ifndef ACTION_F_ALT_RATIONAL_INCLUDED
#define ACTION_F_ALT_RATIONAL_INCLUDED

#include "lib/Action/action.h"

#include "lib/Fopr/afopr_Rational.h"
#include "lib_alt/Force/Fermion/aforce_F_Rational.h"

#include "lib/Fopr/afopr.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

//! action class for RHMC, with externally constructed AFopr_Rational.

/*!
    Action class for RHMC that is an alternative to Action_F_Rational
    in the core library.
                                            [05 Feb 2019 H.Matsufuru]
 */

template<typename AFIELD>
class Action_F_alt_Rational : public Action
{
 public:
  //typedef AField<double> AFIELD;
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

 private:
  std::string m_label;    // label of action

  AFopr<AFIELD>  *m_fopr_langev;
  AFopr<AFIELD>  *m_fopr_H;
  AForce_F<AFIELD> *m_forceF_MD;

  Field *m_U;

  AFIELD m_psf;

  Bridge::VerboseLevel m_vl;   //!< verbose level

 public:
  //! constructor.
  Action_F_alt_Rational(AFopr<AFIELD> *fopr_langev,
                        AFopr<AFIELD> *fopr_H,
                        AForce_F<AFIELD> *forceF_MD)
    : Action(),
      m_fopr_langev(fopr_langev), m_fopr_H(fopr_H),
      m_forceF_MD(forceF_MD)
      { init(); }

  //! destructor.
  ~Action_F_alt_Rational() { tidyup(); }

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
