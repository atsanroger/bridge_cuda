/*!
        @file    asolver_BiCGStab_Cmplx.h

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2#$

        @version $LastChangedRevision: 2569 $
*/

#ifndef ASOLVER_BICGSTAB_Cmplx_H
#define ASOLVER_BICGSTAB_Cmplx_H


#include <string>
using std::string;

#include "lib_alt/Solver/asolver.h"
#include "lib/Fopr/afopr.h"


/*!
    This class implements BiCGStab algorithm for nonhermitian matrix.
    The product of vectors is treated in complex.
                                   [12 Feb 2012 Y.Namekawa]
    (Coding history will be recovered from trac.)
    YAML is implemented.           [14 Nov 2012 Y.Namekawa]
    Multi-threaded.                [10 Jul 2014 H.Matsufuru]
    Introduce unique_ptr to avoid memory leaks.
                                   [21 Mar 2015 Y.Namekawa]
    Add restart.                   [22 Feb 2016 Y.Namekawa]
    Add flop_count.                [ 8 Aug 2016 Y.Namekawa]
    Add use_init_guess.            [ 7 Jul 2017 Y.Namekawa]
    Add a prescription to improve stability of BiCGStab, recommended
    by Kanamori-san. See G.L.G.Sleijpen and H.A.van der Vorst,
    Numerical Algorithms 10(1995)203-22.
                                   [26 Apr 2018 Y.Namekawa]
    Class template version.        [27 Jul 2019 H.Matsufuru]
    Add initial guess, coeff_t     [ 4 Jun 2020 I.Kanamori]
 */

/*
    Uses a prescription to improve stability of BiCGStab.
    See G.L.G.Sleijpen and H.A.van der Vorst,
    Numerical Algorithms 10(1995)203-22.
 */

template<typename AFIELD>
class ASolver_BiCGStab_Cmplx : public ASolver<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  using ASolver<AFIELD>::m_vl;
  static const std::string class_name;
  using InitialGuess = typename ASolver<AFIELD>::InitialGuess;
 private:

  AFopr<AFIELD> *m_fopr;    //!< fermion operator.

  int m_Niter;              //!< maximum iteration number.
  real_t m_Stop_cond;       //!< stopping criterion (squared).

  real_t m_Omega_tolerance; //!< tolerance for the stability

  real_t m_ecrit;           //!< to avoid too small denominator denominator.

  //! Matsufuru added: new AField_dev implementation
  //! Kanamori Modified: m_x is now a reference through pointer
  AFIELD *m_px;
  AFIELD m_r, m_p, m_rh, m_v, m_t;

  struct coeff_t
  {
    std::complex<real_t> rho;
    std::complex<real_t> alpha;
    std::complex<real_t> omega;
  };

  //! for FP16
  struct dcoeff_t
  {
    std::complex<double> rho;
    std::complex<double> alpha;
    std::complex<double> omega;
  };

 protected:

  //! to remember convergence iteration to provide flop count.
  int m_nconv;

  //! to remember convergence iteration to provide flop count.
  int m_init_mult;

  //! mode switch for initial guess
  InitialGuess m_initial_mode;

  //! calling constructor without fermion operator is forbidden.
  ASolver_BiCGStab_Cmplx() { }

 public:
  //! constructor.
  ASolver_BiCGStab_Cmplx(AFopr<AFIELD> *fopr)
    : m_fopr(fopr), m_Niter(0), m_Stop_cond(0.0L), m_initial_mode(InitialGuess::RHS)
  {
    vout.paranoiac("setting fopr to ASolver_BiCGStab_Cmplx: fopr=%p\n", fopr);
    this->init();
  }

  //! destructor.
  ~ASolver_BiCGStab_Cmplx() { this->tidyup(); }

  //! setting parameters by a Parameter object.
  void set_parameters(const Parameters& params);

  //! setting parameters.
  void set_parameters(const int Niter, const real_t Stop_cond);

  //! setting parameters.
  void set_parameters(const int Niter, const real_t Stop_cond, const InitialGuess init_guess);

  void set_init_mode(const InitialGuess init_guess)
  {
    m_initial_mode = init_guess;
  }

  //! setting BiCGStab specific parameters.
  void set_parameters_BiCGStab_series(const real_t Omega_tolerance);

  //! solver main.
  void solve(AFIELD& xq, const AFIELD& b, int& nconv, real_t& diff);

  //! returns the pointer to the fermion operator.
  AFopr<AFIELD> *get_fopr() { return m_fopr; }

  //! returns the floating point operation count.
  double flop_count();


 protected:

  void init(void);

  void tidyup(void);

  void solve_init(const AFIELD& b, real_t& rr, coeff_t&, const real_t scale, const InitialGuess);
  void solve_init_RHS(const AFIELD& b, real_t& rr, coeff_t&);
  void solve_init_ZERO(const AFIELD& b, real_t& rr, coeff_t&);
  void solve_init_GIVEN(const AFIELD& b, real_t& rr, coeff_t&, const real_t);

  void solve_step(real_t& rr, int& iflg, coeff_t&);

  //! for FP16
  void solve_init_double(const AFIELD& b, double& rr, dcoeff_t&, const double scale, const InitialGuess);
  void solve_init_double_RHS(const AFIELD& b, double& rr, dcoeff_t&);
  void solve_init_double_ZERO(const AFIELD& b, double& rr, dcoeff_t&);
  void solve_init_double_GIVEN(const AFIELD& b, double& rr, dcoeff_t&, const double);

  void solve_step_double(double& rr, int& iflg, dcoeff_t&);

#ifdef USE_FACTORY
 private:
  static ASolver<AFIELD> *create_object_with_fopr(AFopr<AFIELD> *fopr)
  {
    return new ASolver_BiCGStab_Cmplx<AFIELD>(fopr);
  }

 public:
  static bool register_factory()
  {
    return ASolver<AFIELD>::Factory_fopr::Register("BiCGStab_Cmplx", create_object_with_fopr);
  }
#endif
};

#endif // ASOLVER_BICGSTAB_Cmplx_H
