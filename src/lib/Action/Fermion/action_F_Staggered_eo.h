/*!
        @file    action_F_Staggered_eo.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-09-08 06:55:16 #$

        @version $LastChangedRevision: 2314 $
*/


#ifndef ACTION_F_STAGGERED_EO_INCLUDED
#define ACTION_F_STAGGERED_EO_INCLUDED

#include "Action/action.h"

#include "Fopr/fopr_Staggered_eo.h"
#include "Force/Fermion/force_F.h"

#include "Smear/projection.h"
#include "Smear/smear.h"
//#include "Force/Fermion/forceSmear.h"
#include "Smear/forceSmear.h"
#include "Solver/solver.h"

#include "Tools/mat_SU_N.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Action for staggered fermion: temporary implementation.

/*!
   This is Action class dedicated to staggered fermion in even-odd
   site index, as a temporary implementation.
   It may be not well organized and to be improved.
                                 [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.         [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                 [21 Mar 2015 Y.Namekawa]
 */

class Action_F_Staggered_eo : public Action {
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  Field_G *m_U;

  Fopr *m_fopr;
  Force *m_fopr_force;
  Field m_psf;

  Index_eo m_index_eo;

  Projection *m_proj;
  Smear *m_smear;
  ForceSmear *m_force_smear;
  Solver *m_solver;

  int m_Nsmear;
  double m_rho;
  std::vector<Field> m_Usmear;

 public:
  Action_F_Staggered_eo(Fopr *fopr, Force *fopr_force)
    : m_vl(CommonParameters::Vlevel())
  {
    m_fopr       = fopr;
    m_fopr_force = fopr_force;

    init();
  }

  Action_F_Staggered_eo(Fopr *fopr, Force *fopr_force, const Parameters& params)
    : m_vl(CommonParameters::Vlevel())
  {
    m_fopr       = fopr;
    m_fopr_force = fopr_force;

    init();

    set_parameters(params);
  }

  ~Action_F_Staggered_eo()
  {
    tidyup();
  }

  void set_parameters(const Parameters& params);
  void set_parameters(const int Nsmear, const double rho);

  void get_parameters(Parameters& params) const;

  void set_config(Field *U)
  {
    m_U = (Field_G *)U;
    m_fopr->set_config(U);
  }

  double langevin(RandomNumbers *);

  double calcH();

  void force(Field&);

 private:
  void init();
  void tidyup();
};
#endif
