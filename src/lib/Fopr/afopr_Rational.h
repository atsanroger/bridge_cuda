/*!
        @file    afopr_Rational.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-06-29 22:43:06 #$

        @version $LastChangedRevision: 2527 $
*/


#ifndef AFOPR_RATIONAL_INCLUDED
#define AFOPR_RATIONAL_INCLUDED

#include "Fopr/afopr.h"

#include "Field/field_G.h"
#include "Solver/ashiftsolver_CG.h"
#include "Tools/math_Rational.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Fermion operator for rational approximation.

/*!
    This class generates fermion operator with a rational approximation
    for a given fermion operator (given to the constructer).
    Shift-solver is used which is at present set to the CG solver
    explicitly.
                                     [05 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.             [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                     [21 Mar 2015 Y.Namekawa]
 */

template<typename AFIELD>
class AFopr_Rational : public AFopr<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  int m_Np;                 // number of poles in rational approx.
  int m_n_exp, m_d_exp;     // numerator and denominator of the exponent
  real_t m_x_min, m_x_max;  // valid range of approximate sign function
  int m_Niter;              // max iteration of shiftsolver
  real_t m_Stop_cond;       // stopping condition of shift solver

  AFopr<AFIELD> *m_fopr;
  AShiftsolver_CG<AFIELD, AFopr<AFIELD> > *m_solver;

  real_t m_a0;
  std::vector<real_t> m_cl;
  std::vector<real_t> m_bl;
  std::vector<AFIELD> m_xq;

 public:
  AFopr_Rational(AFopr<AFIELD> *fopr)
    : m_vl(CommonParameters::Vlevel()), m_fopr(fopr) {}

  AFopr_Rational(AFopr<AFIELD> *fopr, const Parameters& params)
    : m_vl(CommonParameters::Vlevel()), m_fopr(fopr)
  {
    set_parameters(params);
  }

  ~AFopr_Rational()
  {
    delete m_solver;
  }

  void set_parameters(const Parameters& params);
  void set_parameters(const int Np, const int n_exp, const int d_exp,
                      const real_t x_min, const real_t x_max,
                      const int Niter, const real_t Stop_cond);

  void get_parameters(Parameters& params) const;

  void set_config(Field *U);

  void mult(AFIELD& v, const AFIELD& f);

  void mult_dag(AFIELD& v, const AFIELD& f)
  {
    mult(v, f);
  }

  real_t func(const real_t x);

  int field_nvol() { return m_fopr->field_nvol(); }
  int field_nin() { return m_fopr->field_nin(); }
  int field_nex() { return m_fopr->field_nex(); }

  //! returns true if additional field conversion is needed.
  virtual bool needs_convert()
  { return m_fopr->needs_convert(); }

  //! converts a Field object into other format if necessary.
  virtual void convert(AFIELD& v, const Field& w)
  { m_fopr->convert(v, w); }

  //! reverses to a Field object from other format if necessary.
  virtual void reverse(Field& v, const AFIELD& w)
  { m_fopr->reverse(v, w); }

 private:
  void init_parameters();

#ifdef USE_FACTORY
 private:
  static AFopr<AFIELD> *create_object(AFopr<AFIELD> *fopr)
  { return new AFopr_Rational<AFIELD>(fopr); }

  static AFopr<AFIELD> *create_object_with_params(AFopr<AFIELD> *fopr, const Parameters& params)
  { return new AFopr_Rational<AFIELD>(fopr, params); }

 public:
  static bool register_factory()
  {
    bool                init = true;
    init &= AFopr<AFIELD>::Factory_fopr::Register("Rational", create_object);
    init &= AFopr<AFIELD>::Factory_fopr_params::Register("Rational", create_object_with_params);
    return init;
  }
#endif
};
#endif
