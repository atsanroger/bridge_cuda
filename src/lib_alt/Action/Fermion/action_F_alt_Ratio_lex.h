/*!
      @file    action_F_alt_Ratio_lex.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-10-17 13:51:32 #$
      @version $LastChangedRevision: 2550 $
*/

#ifndef ACTION_F_ALT_RATIO_LEX_INCLUDED
#define ACTION_F_ALT_RATIO_LEX_INCLUDED

#include "lib/Action/action.h"

#include "lib/Force/Fermion/force_F.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib_alt/Measurements/Fermion/fprop_alt.h"

//! HMC action for Hasenbusch preconditioned fermions.

/*!
    Hasenbusch preconditioned fermion action in alternative
    implementation.
                                    [04 Feb 2019 H.Matsufuru]
 */

template<typename AFIELD>
class Action_F_alt_Ratio_lex : public Action
{
 public:
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;        //!< verbose level

  AFopr<AFIELD>    *m_fopr_prec;    //!< preconditioner
  AForce_F<AFIELD> *m_forceF_prec;  //!< force of preconditioner
  AFopr<AFIELD>    *m_fopr;         //!< dynamical fermion
  AForce_F<AFIELD> *m_forceF;       //!< force of dynamical fermion
  AFIELD m_psf;                     //!< pseudofermion field
  std::string m_label;              //!< label of action

  Fprop_alt<AFIELD> *m_fprop_H_prec;
  Fprop_alt<AFIELD> *m_fprop_MD;
  Fprop_alt<AFIELD> *m_fprop_H;

  Field *m_U;

 public:
  Action_F_alt_Ratio_lex(
    AFopr<AFIELD> *fopr_prec, AForce_F<AFIELD> *forceF_prec,
    AFopr<AFIELD> *fopr,      AForce_F<AFIELD> *forceF,
    Fprop_alt<AFIELD> *fprop_H_prec,
    Fprop_alt<AFIELD> *fprop_MD, Fprop_alt<AFIELD> *fprop_H)
    : Action(),
      m_fopr_prec(fopr_prec), m_forceF_prec(forceF_prec),
      m_fopr(fopr), m_forceF(forceF),
      m_fprop_H_prec(fprop_H_prec),
      m_fprop_MD(fprop_MD), m_fprop_H(fprop_H)
  { init(); }

  ~Action_F_alt_Ratio_lex()
  { tidyup();}

  void set_parameters(const Parameters&);
  void set_parameters();

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
