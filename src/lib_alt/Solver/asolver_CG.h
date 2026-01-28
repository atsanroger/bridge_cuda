/*!
        @file    asolver_CG.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2#$
        @version $LastChangedRevision: 2562 $
*/

#ifndef ASOLVER_CG_H
#define ASOLVER_CG_H

#include <cstdio>
#include <cstdlib>
#include <vector>
using std::vector;
#include <string>
using std::string;

#include "lib/Fopr/afopr.h"

#include "lib_alt/Solver/asolver.h"

template<typename AFIELD>
class ASolver_CG : public ASolver<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  using ASolver<AFIELD>::m_vl;
  static const std::string class_name;
  using InitialGuess = typename ASolver<AFIELD>::InitialGuess;

 protected:

  AFopr<AFIELD> *m_fopr; //!< fermion operator.

  int m_Niter;           //!< maximum iteration number.
  real_t m_Stop_cond;    //!< stopping criterion (squared).

  //! to remember convergence iteration to provide flop count.
  int m_nconv;

  //! mode switch for initial guess
  InitialGuess m_initial_mode;

  //! calling constructor without fermion operator is forbidden.
  ASolver_CG() { }

  //! working vectors.
  AFIELD m_x, m_r, m_p, m_s;

 public:
  //! constructor.
  ASolver_CG(AFopr<AFIELD> *fopr)
    : m_Niter(0), m_Stop_cond(0.0L)
  {
    m_fopr = fopr;
    this->init();
  }

  //! destructor.
  ~ASolver_CG() { this->tidyup(); }

  //! setting parameters by a Parameter object.
  void set_parameters(const Parameters& params);

  //! setting parameters.
  void set_parameters(const int Niter, const real_t Stop_cond);

  //! setting parameters.
  void set_parameters(const int Niter, const real_t Stop_cond, const InitialGuess init_guess_mode);

  void set_init_mode(const InitialGuess init_guess)
  {
    m_initial_mode = init_guess;
  }

  //! solver main.
  void solve(AFIELD& xq, const AFIELD& b, int& nconv, real_t& diff);

  //! returns the pointer to the fermion operator.
  AFopr<AFIELD> *get_fopr() { return m_fopr; }

  //! returns the floating point operation count.
  double flop_count();


 protected:

  void init(void);

  void tidyup(void);

  void solve_CG_init(real_t& rrp, real_t& rr);

  void solve_CG_step(real_t& rrp, real_t& rr);

  //! for FP16 implementation
  void solve_CG_init_double(double& rrp, double& rr);

  //! for FP16 implementation
  void solve_CG_step_double(double& rrp, double& rr);


#ifdef USE_FACTORY
 private:
  static ASolver<AFIELD> *create_object_with_fopr(AFopr<AFIELD> *fopr)
  { return new ASolver_CG<AFIELD>(fopr); }

 public:
  static bool register_factory()
  {
    bool init = ASolver<AFIELD>::Factory_fopr::Register("CG",
                                                        create_object_with_fopr);
    return init;
  }
#endif
};

#endif // ASOLVER_CG_H
