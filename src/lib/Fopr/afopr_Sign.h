/*!
        @file    afopr_Sign.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-07-05 20:48:53 #$

        @version $LastChangedRevision: 2530 $
*/

#ifndef AFOPR_SIGN_INCLUDED
#define AFOPR_SIGN_INCLUDED

#include "Fopr/afopr.h"

#include "Field/field_G.h"
#include "Solver/ashiftsolver_CG.h"
#include "Tools/math_Sign_Zolotarev.h"
#include "bridge_complex.h"
#include "complexTraits.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Sign function of a given fermion operator.

/*!
    This class defines a sign function of a given fermion operator,
    which must be provided at the construction, and thus itself
    being a kind of fermion operator.
    The sign function is approximated with partially fractional
    series whose coefficients are given by Zolotarev approximation.
    When eigenvalues and eigenvectors of the base fermion operator
    are given, contribution of these modes to sign function are
    explicitly calculated and subtracted from the approximation
    formula.
    Number of subtracted modes are initially set to zero, and if
    the function set_lowmodes() is not called, no subtraction is
    performed.
    In the present implementation, shiftsolver is explicitly set
    to CG, which might be generalized.
                                     [20 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.             [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                     [21 Mar 2015 Y.Namekawa]
    Changed to a template class in ver.2.0.
                                     [12 Feb 2022 H.Matsufuru]
 */


template<typename AFIELD>
class AFopr_Sign : public AFopr<AFIELD>
{
 public:
  typedef AFopr<AFIELD>                               AFOPR;
  typedef typename AFIELD::real_t                     real_t;
  typedef typename ComplexTraits<real_t>::complex_t   complex_t;

  static const std::string class_name;

 private:
  // input parameters
  int m_Np;           //!< number of poles in rational approx.
  real_t m_x_min;     //!< lower range of approximate sign function
  real_t m_x_max;     //!< upper range of approximate sign function
  int m_Niter;        //!< max iteration of shiftsolver
  real_t m_Stop_cond; //!< stopping condition of shift solver

  Bridge::VerboseLevel m_vl;

  std::string m_mode;

  AFOPR *m_fopr;

  // local parameters
  int m_Nin, m_Nvol, m_Nex;

  std::vector<real_t> m_cl;
  std::vector<real_t> m_bl;
  std::vector<real_t> m_sigma;

  //- low-mode subtraction
  int m_Nsbt;
  std::vector<real_t> *m_ev;
  std::vector<AFIELD> *m_vk;

  //- shiftsolver
  AShiftsolver_CG<AFIELD, AFOPR> *m_solver;
  std::vector<AFIELD> m_xq;

  //- workarea
  AFIELD m_w1;

 public:

  AFopr_Sign(AFOPR *fopr, const Parameters& params) : m_fopr(fopr)
  { init(params); }

  DEPRECATED
  AFopr_Sign(AFOPR *fopr) : m_fopr(fopr) { init(); }

  ~AFopr_Sign() { tidyup(); }

  void set_parameters(const Parameters& params);

  void set_parameters(const int Np, const real_t x_min, const real_t x_max,
                      const int Niter, const real_t Stop_cond);

  void get_parameters(Parameters& params) const;

  void set_config(Field *U);

  void set_mode(const std::string mode);

  std::string get_mode() const { return m_mode; }

  void set_lowmodes(const int Nsbt, std::vector<real_t> *,
                    std::vector<AFIELD> *);

  void mult(AFIELD& v, const AFIELD& w);

  void mult_dag(AFIELD& v, const AFIELD& w) { mult(v, w); }

  //! returns true if additional field conversion is needed.
  virtual bool needs_convert()
  { return m_fopr->needs_convert(); }

  //! converts a Field object into other format if necessary.
  virtual void convert(AFIELD& v, const Field& w)
  { m_fopr->convert(v, w); }

  //! reverses to a Field object from other format if necessary.
  virtual void reverse(Field& v, const AFIELD& w)
  { m_fopr->reverse(v, w); }

  int field_nin() { return m_Nin; }
  int field_nvol() { return m_Nvol; }
  int field_nex() { return m_Nex; }

  //! returns the number of floating point operations.
  double flop_count();

  //! returns the flops per site for specified mode.
  double flop_count(const std::string mode);

 private:
  void init();

  void init(const Parameters& params);

  void tidyup();

  void init_parameters();

  void subtract_lowmodes(AFIELD&);
  void evaluate_lowmodes(AFIELD&, const AFIELD&);

  real_t sign_zolotarev(const real_t x);

#ifdef USE_FACTORY
 private:
  static AFOPR *create_object(AFOPR *fopr)
  {
    return new AFopr_Sign<AFIELD>(fopr);
  }

  static AFOPR *create_object_with_params(AFOPR *fopr, const Parameters& params)
  {
    return new AFopr_Sign<AFIELD>(fopr, params);
  }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= AFOPR::Factory_fopr::Register("Sign", create_object);
    init &= AFOPR::Factory_fopr_params::Register("Sign", create_object_with_params);
    return init;
  }
#endif
};
#endif
