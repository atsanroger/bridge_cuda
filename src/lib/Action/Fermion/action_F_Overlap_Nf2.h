/*!
        @file    action_F_Overlap_Nf2.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-09-08 06:55:16 #$

        @version $LastChangedRevision: 2314 $
*/

#ifndef ACTION_F_OVERLAP_NF2_INCLUDED
#define ACTION_F_OVERLAP_NF2_INCLUDED

#include "Action/action.h"

#include "Force/Fermion/force_F.h"
#include "Measurements/Fermion/fprop.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Action for overlap fermion: temporary implementation.

/*!
    This is Action class dedicated to overlap fermion, as a
    temporary implementation.
    It may be not well organized and to be improved.
                                 [28 Dec 2011 H.Matsufuru]
    unique_ptr is introduced to avoid memory leaks
                                 [21 Mar 2015 Y.Namekawa]
 */

class Action_F_Overlap_Nf2 : public Action
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  std::string m_label;

  Fopr *m_fopr;
  Force *m_fopr_force;
  Fprop *m_fprop_MD;
  Fprop *m_fprop_H;

  Field_G *m_U;

  Field m_psf;

 public:
  Action_F_Overlap_Nf2(Fopr *fopr, Force *fopr_force,
                       Fprop *fprop_MD, Fprop *fprop_H)
    : m_vl(CommonParameters::Vlevel()),
    m_fopr(fopr),
    m_fopr_force(fopr_force),
    m_fprop_MD(fprop_MD),
    m_fprop_H(fprop_H)
  {}

  Action_F_Overlap_Nf2(Fopr *fopr, Force *fopr_force,
                       Fprop *fprop_MD, Fprop *fprop_H,
                       const Parameters& params)
    : m_vl(CommonParameters::Vlevel()),
    m_fopr(fopr),
    m_fopr_force(fopr_force),
    m_fprop_MD(fprop_MD),
    m_fprop_H(fprop_H)
  {
    set_parameters(params);
  }

  ~Action_F_Overlap_Nf2() {}

  void set_parameters(const Parameters& param);

  void get_parameters(Parameters& param) const;

  void set_label(const std::string label)
  {
    m_label = label;
    vout.detailed(m_vl, "  label: %s\n", m_label.c_str());
  }

  std::string get_label()
  {
    return m_label;
  }

  void set_config(Field *U)
  {
    m_U = (Field_G *)U;

    m_fopr->set_config(U);
    m_fopr_force->set_config(U);
  }

  double langevin(RandomNumbers *);

  double calcH();

  void force(Field& force);

  // not implemented function
  // Field force_core(Field_F&, Field_F&);
  // Field force_Wilson_core(Field_F&, Field_F&);
};
#endif
