/*!
        @file    force_F_Wilson_Nf2.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-04-17 11:32:37 #$

        @version $LastChangedRevision: 2513 $
*/

#ifndef AFORCE_F_WILSON_NF2_INCLUDED
#define AFORCE_F_WILSON_NF2_INCLUDED

#include "lib/Force/Fermion/aforce_F.h"

#include "lib/Fopr/afopr.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Force for the standard Wilson fermion operator

/*!
    This class calculates the force of the standard Wilson fermion.
    The gamma matrix representation is given as control string
    "Dirac"(default) or "Chiral" at the construction, which is
    used to construct the Fopr_Wilson instance.
                                     [23 Dec 2011 H.Matusfuru]
 */

template<typename AFIELD>
class AForce_F_Wilson_Nf2 : public AForce_F<AFIELD>
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
 
  AFopr<AFIELD> *m_fopr_w;

  AFIELD  m_U, m_force1;
  AFIELD  m_zeta, m_eta2, m_eta3;

 public:
  //! constructor.
  AForce_F_Wilson_Nf2(const Parameters& params) { init(params); }

  ~AForce_F_Wilson_Nf2() { tidyup(); }

  void set_parameters(const Parameters& params);
  void set_parameters(const double kappa, const std::vector<int> bc);

  void get_parameters(Parameters& params) const;

  void set_config(Field *U);

  //void force_udiv(Field& force, const Field& eta);
  //void force_udiv1(Field& force, const Field& zeta, const Field& eta);

  void force_udiv(AFIELD& force, const AFIELD& eta);
  void force_udiv1(AFIELD& force, const AFIELD& zeta, const AFIELD& eta);

 private:
  //! initial setup.
  void init(const Parameters& params);

  void tidyup();

  void force_udiv1_impl(AFIELD& force, const AFIELD& zeta, const AFIELD& eta);

};
#endif
