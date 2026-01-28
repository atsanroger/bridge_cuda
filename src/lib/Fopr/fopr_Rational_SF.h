/*!
        @file    fopr_Rational_SF.h

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/


#ifndef FOPR_RATIONAL_SF_INCLUDED
#define FOPR_RATIONAL_SF_INCLUDED

#include "fopr_Rational.h"

#include "Field/field_SF.h"
#include "Solver/shiftsolver_CG.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Fermion operator for rational approximation.

/*!
    This class generates fermion operator with rational approximation
    for a given fermion operator (given to the constructer).
    Shift-solver is used which is at present set to the CG solver
    explicitly.
    <ul>
    <li>Modified for a Wilson type Dirac operator with SF BC.
    <li>A modification is to set the field value at the temporal boundary to zero before and after a multiplication of the Dirac operator.
    <li>[07 Apr 2012 Y.Taniguchi]
    </ul>
    (Coding history will be recovered from trac.)
    YAML is implemented.     [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                             [21 Mar 2015 Y.Namekawa]
 */

class Fopr_Rational_SF : public Fopr
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  int m_Np;                 // number of poles in rational approx.
  int m_n_exp, m_d_exp;     // numerator and denominator of the exponent
  double m_x_min, m_x_max;  // valid range of approximate sign function
  int m_Niter;              // max iteration of shiftsolver
  double m_Stop_cond;       // stopping condition of shift solver

  Fopr *m_fopr;
  Shiftsolver_CG *m_solver;

  double m_a0;
  std::vector<double> m_cl;
  std::vector<double> m_bl;
  std::vector<Field> m_xq;

 public:

  Fopr_Rational_SF(Fopr *fopr)
    : m_vl(CommonParameters::Vlevel()), m_fopr(fopr) {}

  Fopr_Rational_SF(Fopr *fopr, const Parameters& params)
    : m_vl(CommonParameters::Vlevel()), m_fopr(fopr)
  {
    set_parameters(params);
  }

  ~Fopr_Rational_SF()
  {
    delete m_solver;
  }

  void set_parameters(const Parameters& params);
  void set_parameters(const int Np, const int n_exp, const int d_exp, const double x_min, const double x_max,
                      const int Niter, const double Stop_cond);

  void get_parameters(Parameters& params) const;

  void set_config(Field *U)
  {
    m_fopr->set_config(U);
  }

  void mult(Field& v, const Field& f);

  void mult_dag(Field& v, const Field& f)
  {
    mult(v, f);
  }

  double func(const double x);

  int field_nvol() { return m_fopr->field_nvol(); }
  int field_nin() { return m_fopr->field_nin(); }
  int field_nex() { return m_fopr->field_nex(); }

  //! this returns the number of floating point operations.
  double flop_count();

 private:
  void init_parameters();

#ifdef USE_FACTORY
 private:
  static Fopr *create_object(Fopr *fopr)
  {
    return new Fopr_Rational_SF(fopr);
  }

  static Fopr *create_object_with_params(Fopr *fopr, const Parameters& params)
  {
    return new Fopr_Rational_SF(fopr, params);
  }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= Fopr::Factory_fopr::Register("Rational_SF", create_object);
    init &= Fopr::Factory_fopr_params::Register("Rational_SF", create_object_with_params);
    return init;
  }
#endif
};
#endif
