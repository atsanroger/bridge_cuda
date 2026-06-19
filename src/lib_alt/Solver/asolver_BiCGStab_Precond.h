/*!
      @file    asolver_BiCGStab_Precond.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#ifndef ASOLVER_BICGSTAB_PRECOND_H
#define ASOLVER_BICGSTAB_PRECOND_H

#include <cstdio>
#include <cstdlib>
#include <vector>
using std::vector;
#include <string>
using std::string;

#include "lib_alt/Solver/asolver.h"
#include "lib/Fopr/afopr.h"
#include "lib_alt/Solver/aprecond.h"

template<typename AFIELD>
class ASolver_BiCGStab_Precond : public ASolver<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  using ASolver<AFIELD>::m_vl;
  static const std::string class_name;

 private:

  AFopr<AFIELD> *m_fopr;    //!< fermion operator.
  APrecond<AFIELD> *m_prec; //!< preconditioner.

  int m_Niter;              //!< maximum iteration number.
  real_t m_Stop_cond;       //!< stopping criterion (squared).
  real_t m_Stop_cond2;      //!< stopping criterion for inner solver.

  //! Matsufuru added: new AField_dev implementation
  AFIELD m_x, m_r, m_p, m_s, m_rh, m_v, m_t, m_u;
  real_t m_rho_prev, m_alpha_prev, m_omega_prev;

 protected:

  //! to remember convergence iteration to provide flop count.
  int m_nconv;

  //! calling constructor without fermion operator is forbidden.
  ASolver_BiCGStab_Precond() { }

 public:
  //! constructor.
  ASolver_BiCGStab_Precond(AFopr<AFIELD> *fopr,
                           APrecond<AFIELD> *prec)
    : m_Niter(0), m_Stop_cond(0.0L)
  {
    m_fopr = fopr;
    m_prec = prec;
    init();
  }

  //! destructor.
  ~ASolver_BiCGStab_Precond() { tidyup(); }

  //! setting parameters by a Parameter object.
  void set_parameters(const Parameters& params);

  //! setting parameters.
  void set_parameters(const int Niter, const real_t Stop_cond);

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

  void solve_init(const AFIELD& b, real_t& rr);

  void solve_step(real_t& rr);

  void prec(AFIELD&, AFIELD&);
};

#endif // ASOLVER_BICGSTAB_PRECOND_H
