/*!
        @file    force_F_Wilson_eo_Nf2.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-07-03 23:32:09 #$

        @version $LastChangedRevision: 2529 $
*/

#ifndef AFORCE_F_WILSON_EO_NF2_INCLUDED
#define AFORCE_F_WILSON_EO_NF2_INCLUDED

#include "lib/Force/Fermion/aforce_F.h"

#include "lib/Fopr/afopr.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib_alt/Force/Fermion/aforce_F_Wilson_Nf2.h"


//! Force for the even-odd Wilson fermion operator

/*!
    This class calculates the force of the even-odd Wilson fermion.
    The gamma matrix representation is given as control string
    "Dirac"(default) or "Chiral" at the construction, which is
    used to construct the Fopr_Wilson instance.
                                     [23 May 2023 H.Matusfuru]
 */

template<typename AFIELD>
class AForce_F_Wilson_eo_Nf2 : public AForce_F<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

 private:
  real_t m_kappa;               //!< hopping parameter
  std::vector<int> m_boundary;  //!< pointer to boundary condition
  std::string m_repr;           //!< gamma matrix representation
  Bridge::VerboseLevel m_vl;

  using AForce_F<AFIELD>::m_Ucp;
 
  AForce_F_Wilson_Nf2<AFIELD> *m_force_w;
  AFopr<AFIELD> *m_fopr_w;

  bool m_at_construct;

  //  AFIELD  m_U, m_force1;
  AFIELD  m_force1; //, m_force2;
  AFIELD  m_zeta; //, m_eta2, m_eta3;
  AFIELD  m_v1, m_eta_lex, m_zeta_lex;

 public:
  //! constructor.
  AForce_F_Wilson_eo_Nf2(const Parameters& params)
    :AForce_F<AFIELD>()
  { init(params); }

  ~AForce_F_Wilson_eo_Nf2() { tidyup(); }

  void set_parameters(const Parameters& params);
  void set_parameters(const double kappa, const std::vector<int> bc);

  void get_parameters(Parameters& params) const;

  void set_config(Field *U);

  void force_udiv(AFIELD& force, const AFIELD& eta);
  void force_udiv1(AFIELD& force, const AFIELD& zeta, const AFIELD& eta);

 private:
  //! initial setup.
  void init(const Parameters& params);

  void tidyup();

  void force_udiv1_impl(AFIELD& force, const AFIELD& zeta, const AFIELD& eta);

};
#endif
