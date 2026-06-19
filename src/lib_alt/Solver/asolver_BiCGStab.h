/*!
        @file    asolver_BiCGStab.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2#$
        @version $LastChangedRevision: 2569 $
*/

#ifndef ASOLVER_BICGSTAB_H
#define ASOLVER_BICGSTAB_H

#include <cstdio>
#include <cstdlib>
#include <vector>
using std::vector;
#include <string>
using std::string;

#include "asolver.h"
#include "lib/Fopr/afopr.h"

template<typename AFIELD>
class ASolver_BiCGStab : public ASolver<AFIELD>
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
  AFIELD m_x, m_r, m_p, m_s, m_rh, m_v, m_t;
  //  real_t  m_rho_prev, m_alpha_prev, m_omega_prev;

 protected:

  //! to remember convergence iteration to provide flop count.
  int m_nconv;

  //! to remember convergence iteration to provide flop count.
  int m_init_mult;

  //! mode switch for initial guess
  InitialGuess m_initial_mode;

  //! calling constructor without fermion operator is forbidden.
  ASolver_BiCGStab() { }

 public:
  //! constructor.
  ASolver_BiCGStab(AFopr<AFIELD> *fopr)
    : m_Niter(0), m_Stop_cond(0.0L)
  {
    m_fopr = fopr;
    this->init();
  }

  //! destructor.
  ~ASolver_BiCGStab() { this->tidyup(); }

  //! setting parameters by a Parameter object.
  void set_parameters(const Parameters& params);

  //! setting parameters.
  void set_parameters(const int Niter, const real_t Stop_cond);

  //! setting parameters.
  void set_parameters(const int Niter, const real_t Stop_cond, const InitialGuess init_mode);

  //! setting the initial guess mode
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

  void solve_init(const AFIELD& b, real_t& rr,
                  real_t& alpha_prev, real_t& rho_prev,
                  real_t& omega_prev, const InitialGuess init_mode);

  void solve_step(real_t& rr, int& iflg,
                  real_t& alpha_prev, real_t& rho_prev,
                  real_t& omega_prev);

  //! for FP16 implementation
  void solve_init_double(const AFIELD& b, double& rr,
                         double& alpha_prev, double& rho_prev,
                         double& omega_prev, const InitialGuess init_mode);

  //! for FP16 implementation
  void solve_step_double(double& rr, int& iflg,
                         double& alpha_prev, double& rho_prev,
                         double& omega_prev);


#ifdef USE_FACTORY
 private:
  static ASolver<AFIELD> *create_object_with_fopr(AFopr<AFIELD> *fopr)
  { return new ASolver_BiCGStab<AFIELD>(fopr); }

 public:
  static bool register_factory()
  {
    bool init = ASolver<AFIELD>::Factory_fopr::Register("BiCGStab",
                                                        create_object_with_fopr);
    return init;
  }
#endif
};

#endif // ASOLVER_BICGSTAB_H
