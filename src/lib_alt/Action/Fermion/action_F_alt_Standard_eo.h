/*!
      @file    action_F_alt_Standard_eo.h
      @brief
      @author  Hideo Matsufuru
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-05-23 16:09:22 #$
      @version $LastChangedRevision: 2517 $
*/

#ifndef ACTION_F_ALT_STANDARD_EO_INCLUDED
#define ACTION_F_ALT_STANDARD_EO_INCLUDED

#include "lib/Action/action.h"

#include "lib/Force/Fermion/force_F.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib_alt/Measurements/Fermion/fprop_alt.h"

//! Standard fermion action for HMC.

/*!
    Standard fermion action with alternative implementation.
    [04 Oct 2018 H.Matsufuru]
 */

template<typename AFIELD>
class Action_F_alt_Standard_eo : public Action
{
 public:
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

 private:
  AFopr<AFIELD> *m_fopr;
  AForce_F<AFIELD> *m_fopr_force;
  AFIELD m_psf;
  std::string m_label;
  Bridge::VerboseLevel m_vl;   //!< verbose level

  Fprop_alt<AFIELD> *m_fprop_MD;
  Fprop_alt<AFIELD> *m_fprop_H;

  Field *m_U;

 public:
  Action_F_alt_Standard_eo(
    AFopr<AFIELD> *fopr, AForce_F<AFIELD> *fopr_force,
    Fprop_alt<AFIELD> *fprop_MD, Fprop_alt<AFIELD> *fprop_H)
    : Action(),
      m_fopr(fopr), m_fopr_force(fopr_force),
      m_fprop_MD(fprop_MD), m_fprop_H(fprop_H)
   { init(); }

  ~Action_F_alt_Standard_eo()
    { tidyup(); }

  void set_parameters(const Parameters&);

  void get_parameters(Parameters&) const;

  void set_label(const std::string label)
  {
    m_label = label;
    vout.detailed(m_vl, "  label: %s\n", m_label.c_str());
  }

  std::string get_label()
  { return m_label; }

  void set_config(Field *U);

  double langevin(RandomNumbers *);

  double calcH();

  void force(Field&);

 private:
  void init();
  void tidyup();

};
#endif
