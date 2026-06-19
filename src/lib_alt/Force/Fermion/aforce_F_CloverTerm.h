/*!
        @file    aforce_F_CloverTerm.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-10-17 13:51:32 #$

        @version $LastChangedRevision: 2550 $
*/

#ifndef AFORCE_F_CLOVERTERM_INCLUDED
#define AFORCE_F_CLOVERTERM_INCLUDED

#include "lib/Force/Fermion/aforce_F.h"
#include "lib/Force/Fermion/force_F_CloverTerm.h"

#include "lib/Fopr/afopr.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

template<typename AFIELD>
class ShiftAField_lex;

template<typename AFIELD>
class AFopr_CloverTerm;

template<typename AFIELD>
class AStaple_lex;


//! Force for the Clover fermion operator

/*!
    This class calculates the force of the Clover fermion.
    This is a temporary implementation using corelib code.
                                     [10 Oct 2023 H.Matusfuru]
 */

template<typename AFIELD>
class AForce_F_CloverTerm : public AForce_F<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

 private:
  int m_Ndim;
  real_t m_kappa;               //!< hopping parameter
  real_t m_cSW;                 //!< clover coefficient
  std::vector<int> m_boundary;  //!< pointer to boundary condition
  std::string m_repr;           //!< gamma matrix representation
  Bridge::VerboseLevel m_vl;

  using AForce_F<AFIELD>::m_Ucp;
 
  AFopr_CloverTerm<AFIELD> *m_fopr_csw;

  AStaple_lex<AFIELD> *m_staple;

  ShiftAField_lex<AFIELD>* m_shift;

  //Force_F_CloverTerm *m_forceF_ct;

  AFIELD  m_zeta, m_eta2, m_eta3;
  AFIELD  m_zeta_mu;
  AFIELD  m_vt1, m_vt2, m_vt3;

  AFIELD  m_force1;
  AFIELD  m_Cup, m_Cdn;
  AFIELD  m_Utmp1, m_Utmp2;

 public:
  //! constructor.
  AForce_F_CloverTerm(const Parameters& params)
  { init(params); }

  ~AForce_F_CloverTerm() { tidyup(); }

  void set_parameters(const Parameters& params);

  void set_parameters(const double kappa,
                      const double cSW,
                      const std::vector<int> bc);

  void get_parameters(Parameters& params) const;

  void set_config(Field *U);
  void set_config_omp(Field *U);
  void set_config_impl(Field *U);

  void force_udiv(AFIELD& force, const AFIELD& eta);
  void force_udiv1(AFIELD& force, const AFIELD& zeta, const AFIELD& eta);

 private:
  //! initial setup.
  void init(const Parameters& params);

  void tidyup();

  void set_component();

  int index_dir(const int mu, const int nu);

  void force_udiv1_impl(AFIELD& force, const AFIELD& zeta, const AFIELD& eta);

  void force_udiv1_impl1(AFIELD& force, const AFIELD& zeta, const AFIELD& eta);

  void force_udiv1_impl2(AFIELD& force, const AFIELD& zeta, const AFIELD& eta);

};
#endif
