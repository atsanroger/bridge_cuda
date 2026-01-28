/*!
        @file    force_F_Clover_Nf2.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-10-17 13:51:32 #$

        @version $LastChangedRevision: 2550 $
*/

#ifndef AFORCE_F_CLOVER_NF2_INCLUDED
#define AFORCE_F_CLOVER_NF2_INCLUDED

#include "lib/Force/Fermion/aforce_F.h"
//#include "lib/Force/Fermion/force_F_Clover_Nf2.h"
//#include "lib/Force/Fermion/force_F_CloverTerm.h"
#include "lib_alt/Force/Fermion/aforce_F_Wilson_Nf2.h"
#include "lib_alt/Force/Fermion/aforce_F_CloverTerm.h"

#include "lib/Fopr/afopr.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Force for the Clover fermion operator

/*!
    This class calculates the force of the Clover fermion.
    This is a temporary implementation using corelib code.
                                     [10 Oct 2023 H.Matusfuru]
 */

template<typename AFIELD>
class AForce_F_Clover_Nf2 : public AForce_F<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

 private:
  real_t m_kappa;               //!< hopping parameter
  real_t m_cSW;                 //!< clover coefficient
  std::vector<int> m_boundary;  //!< pointer to boundary condition
  std::string m_repr;           //!< gamma matrix representation
  Bridge::VerboseLevel m_vl;

  using AForce_F<AFIELD>::m_Ucp;
 
  AFopr<AFIELD> *m_fopr_w;

  //Force_F_Clover_Nf2 *m_forceF;
  //Force_F_Wilson_Nf2 *m_forceF_w;
  AForce_F_Wilson_Nf2<AFIELD> *m_forceF_w;
  AForce_F_CloverTerm<AFIELD> *m_forceF_ct;
  //  Force_F_CloverTerm *m_forceF_ct;

  AFIELD  m_U, m_force1, m_force2;
  AFIELD  m_zeta, m_eta2, m_eta3;

 public:
  //! constructor.
  AForce_F_Clover_Nf2(const Parameters& params) { init(params); }

  ~AForce_F_Clover_Nf2() { tidyup(); }

  void set_parameters(const Parameters& params);

  void set_parameters(const double kappa,
                      const double cSW,
                      const std::vector<int> bc);

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
