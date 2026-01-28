/*!
        @file    aforce_F_Rational.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-10-17 13:51:32 #$

        @version $LastChangedRevision: 2550 $
*/

#ifndef AFORCE_F_RATIONAL_INCLUDED
#define AFORCE_F_RATIONAL_INCLUDED

#include "Force/Fermion/aforce_F.h"

#include "Field/field_F.h"
#include "Fopr/afopr.h"
//#include "Fopr/fopr_Rational.h"
#include "Solver/ashiftsolver_CG.h"
#include "Tools/math_Rational.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Force calculation for smeared fermion operators.

/*!
    This class determines the force of a rational approximation
    for a given fermion operator.
                                    [28 Dec 2011 H.Matsufuru]

    - Template version              [30 Jun 2023 H.Matsufuru]
 */

template<typename AFIELD>
class AForce_F_Rational : public AForce_F<AFIELD>
{
 public:
  static const std::string class_name;
  typedef typename AFIELD::real_t real_t;

 private:
  int    m_Np;              // number of poles in rational approx.
  int    m_n_exp, m_d_exp;  // numerator and denominator of the exponent
  real_t m_x_min, m_x_max;  // valid range of approximate sign function
  int    m_Niter;           // max iteration of shiftsolver
  real_t m_Stop_cond;       // stopping condition of shift solver
  Bridge::VerboseLevel m_vl;

  using AForce_F<AFIELD>::m_Ucp;

  Field_G *m_U;
  AFopr<AFIELD>    *m_fopr;
  AForce_F<AFIELD> *m_forceF;

  // rational approx. coefficients
  real_t  m_a0;
  std::vector<real_t> m_bl;
  std::vector<real_t> m_cl;

  std::vector<AFIELD> m_psi;
  std::vector<AFIELD> m_psi2;
  AFIELD m_v1, m_v2;

  AShiftsolver_CG<AFIELD, AFopr<AFIELD> >* m_solver;

  AFIELD m_force1;

  public:

  AForce_F_Rational(AFopr<AFIELD> *fopr, AForce_F<AFIELD> *forceF,
                    const Parameters& params)
    : AForce_F<AFIELD>(), m_fopr(fopr), m_forceF(forceF)
  { init(params); }

  ~AForce_F_Rational() { tidyup(); }

  void set_parameters(const Parameters& params);
  void set_parameters(const int Np, const int n_exp, const int d_exp,
                      const real_t x_min, const real_t x_max,
                      const int Niter, const real_t Stop_cond);

  void get_parameters(Parameters& params) const;

  void set_config(Field *U);

  void force_udiv(AFIELD&, const AFIELD&);

  void force_core1(AFIELD&, const AFIELD&, const AFIELD&);
  void force_udiv1(AFIELD&, const AFIELD&, const AFIELD&);

 private:

  void init(const Parameters& params);

  void tidyup();

  void force_udiv_impl(AFIELD&, const AFIELD&);
  void force_udiv1_impl(AFIELD&, const AFIELD&, const AFIELD&);

#ifdef USE_FACTORY
 private:
  //static AForce_F<AFIELD> *create_object(AFopr<AFIELD> *fopr,
  //                                       AForce_F<AFIELD> *force_F)
  //{ return new AForce_F_Rational(fopr, force_F); }

 public:
  //static bool register_factory()
  //{
  //  bool init1 = AForce_F<AFIELD>::Factory_fopr_force::Register(
  //                                         "Rational", create_object);
  //  return init1;
  // }
#endif


};
#endif
