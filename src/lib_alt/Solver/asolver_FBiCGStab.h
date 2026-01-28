/*!
        @file    $Id:: asolver_FBiCGStab.h #$

        @brief

        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2#$

        @version $LastChangedRevision: 2499 $
*/


#ifndef ASOLVER_FBICGSTAB_H
#define ASOLVER_FBICGSTAB_H

#include <string>

#include "asolver.h"
#include "lib/Fopr/afopr.h"
#include "aprecond.h"

/*
  Flexible BiCGStab:
    a different impelemntation from asolver_BiCGStab_Precond
    this version allows exiting after each mult (not each 2 mults)

 */

template<typename AFIELD>
class ASolver_FBiCGStab : public ASolver<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  using ASolver<AFIELD>::m_vl;
  static const std::string class_name;
  using InitialGuess = typename ASolver<AFIELD>::InitialGuess;

  struct coeff_t
  {
    std::complex<real_t> rho;
    std::complex<real_t> alpha;
  }
  coeff;

 private:

  AFopr<AFIELD> *m_fopr;    //!< fermion operator.
  APrecond<AFIELD> *m_prec; //!< preconditioner.

  int m_Niter;              //!< maximum iteration number.
  real_t m_Stop_cond;       //!< stopping criterion (squared).
  real_t m_Stop_cond2;      //!< stopping criterion for inner solver.

  //! Matsufuru added: new AField_dev implementation
  AFIELD m_x, m_r, m_p, m_s, m_rh, m_v, m_t, m_u;

  //  std::complex<real_t>  m_rho_prev, m_alpha;

 protected:

  //! to remember convergence iteration to provide flop count.
  int m_nconv;

  InitialGuess m_initial_mode;

  double m_Omega_tolerance = 0.7;
  //! calling constructor without fermion operator is forbidden.
  ASolver_FBiCGStab() { }

 public:
  //! constructor.
  ASolver_FBiCGStab(AFopr<AFIELD> *fopr,
                    APrecond<AFIELD> *prec)
    : m_Niter(0), m_Stop_cond(0.0L), m_initial_mode(InitialGuess::RHS)
  {
    m_fopr = fopr;
    m_prec = prec;
    init();
  }

  //! destructor.
  ~ASolver_FBiCGStab() { tidyup(); }

  //! setting parameters by a Parameter object.
  void set_parameters(const Parameters& params);

  //! setting parameters.
  void set_parameters(const int Niter, const real_t Stop_cond);

  //! setting parameters.
  void set_parameters(const int Niter, const real_t Stop_cond, const InitialGuess init_guess_mode);

  //! setting parameters.
  void set_parameters(const int Niter, const real_t Stop_cond, const bool use_init_guess);

  //! setting BiCGStab specific parameters.
  void set_parameters_BiCGStab_series(const real_t Omega_tolerance);

  //! solver main.
  void solve(AFIELD& xq, const AFIELD& b, int& nconv, real_t& diff);

  //! returns the pointer to the fermion operator.
  AFopr<AFIELD> *get_afopr() { return m_fopr; }

  //! returns the floating point operation count.
  double flop_count();

 protected:

  //! initial setup.
  void init(void);

  //! final tidy-up.
  void tidyup(void);

  void solve_init(const AFIELD& b, const AFIELD& xq, real_t& rr, coeff_t& prev);

  void solve_step1(real_t& rr, coeff_t& prev);
  void solve_step2(real_t& rr, coeff_t& prev);

  void prec(AFIELD&, AFIELD&);

  double flop_count_intermediate(const int iter);

  /*
#ifdef USE_FACTORY
 private:
  static ASolver<AFIELD> *create_object_with_fopr(AFopr<AFIELD>* fopr)
  { return new ASolver_FBiCGStab<AFIELD>(fopr); }

 public:
  static bool register_factory()
  {
    bool init = ASolver<AFIELD>::Factory_fopr::Register("ASolver_FBiCGStab",
                                                        create_object_with_fopr);
    return init;
  }
#endif
  */
};

#endif // ASOLVER_FBICGSTAB_H
