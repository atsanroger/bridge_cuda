/*!
        @file    action_F_Staggered_lex.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef ACTION_F_STAGGERED_LEX_INCLUDED
#define ACTION_F_STAGGERED_LEX_INCLUDED

#include "Action/action.h"
#include "Measurements/Fermion/fprop.h"
#include "Fopr/fopr_Staggered.h"
#include "Force/Fermion/force_F.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Standard fermion action for HMC.

/*!
    This class is used to define an ataggered fermion action
    used in HMC. The pseudo-fermion field is set on only even
    stes so as to halve the number of flavors.
    The implementation is based on the Action_F_Standard_lex
    class. Added code discards the pseudo-fermion field on odd
    sites.
                                        [03 Oct 2015 H.Matsufuru]
 */

class Action_F_Staggered_lex : public Action
{
 public:
  static const std::string class_name;

 private:
  Fopr_Staggered *m_fopr;
  Force *m_fopr_force;
  Field m_psf;
  std::string m_label;

  Fprop *m_fprop_MD;
  Fprop *m_fprop_H;

  Field *m_U;

  Bridge::VerboseLevel m_vl;


 public:
  Action_F_Staggered_lex(
    Fopr_Staggered *fopr, Force *fopr_force,
    Fprop *fprop_MD, Fprop *fprop_H)
    : Action(),
    m_fopr(fopr), m_fopr_force(fopr_force),
    m_fprop_MD(fprop_MD), m_fprop_H(fprop_H)
  {
    set_parameters();
  }

  Action_F_Staggered_lex(
    unique_ptr<Fopr>& fopr, unique_ptr<Force>& fopr_force,
    unique_ptr<Fprop>& fprop_MD, unique_ptr<Fprop>& fprop_H)
    : Action(),
    m_fopr((Fopr_Staggered *)fopr.get()),
    m_fopr_force(fopr_force.get()),
    m_fprop_MD(fprop_MD.get()), m_fprop_H(fprop_H.get())
  {
    set_parameters();
  }

  ~Action_F_Staggered_lex()
  {
    // delete m_fprop;
  }

  void set_parameters(const Parameters&);
  void set_parameters();

  void get_parameters(Parameters& params) const;

  void set_label(std::string label)
  {
    m_label = label;
    vout.detailed(m_vl, "  label: %s\n", m_label.c_str());
  }

  std::string get_label()
  {
    return m_label;
  }

  void set_config(Field *U);

  double langevin(RandomNumbers *);

  double calcH();

  void force(Field&);
};
#endif
