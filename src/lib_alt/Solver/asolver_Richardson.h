/*!
        @file    $Id:: asolver_Richardson.h#$

        @brief

        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2#$

        @version $LastChangedRevision: 2499 $
*/


#ifndef ASOLVER_RICHARDSON_H
#define ASOLVER_RICHARDSON_H

#include <string>

#include "asolver.h"
#include "lib/Fopr/afopr.h"
#include "aprecond.h"

/*
  preconditionted Rihcradson iteration

 */

template<typename AFIELD>
class ASolver_Richardson : public ASolver<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  using ASolver<AFIELD>::m_vl;
  static const std::string class_name;
  using InitialGuess = typename ASolver<AFIELD>::InitialGuess;


 private:

  AFopr<AFIELD> *m_fopr;    //!< fermion operator.
  APrecond<AFIELD> *m_prec; //!< preconditioner.

  int m_Niter;              //!< maximum iteration number.
  real_t m_Stop_cond;       //!< stopping criterion (squared).
  real_t m_Stop_cond2;      //!< stopping criterion for inner solver.

  //! Matsufuru added: new AField_dev implementation
  AFIELD m_r, m_s;


 protected:

  //! to remember convergence iteration to provide flop count.
  int m_nconv;

  InitialGuess m_initial_mode;

  //! calling constructor without fermion operator is forbidden.
  ASolver_Richardson() { }

 public:
  //! constructor.
  ASolver_Richardson(AFopr<AFIELD> *fopr,
                     APrecond<AFIELD> *prec)
    : m_Niter(0), m_Stop_cond(0.0L), m_initial_mode(InitialGuess::ZERO)
  {
    m_fopr = fopr;
    m_prec = prec;
    init();
  }

  //! destructor.
  ~ASolver_Richardson() { tidyup(); }

  //! setting parameters by a Parameter object.
  void set_parameters(const Parameters& params);

  //! setting parameters.
  void set_parameters(const int Niter, const real_t Stop_cond);

  //! setting parameters.
  void set_parameters(const int Niter, const real_t Stop_cond, const InitialGuess init_guess_mode);

  //! setting parameters.
  void set_parameters(const int Niter, const real_t Stop_cond, const bool use_init_guess);

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

  void solve_init(const AFIELD& b, const AFIELD& xq, real_t& rr);

  //  void solve(real_t& rr, coeff_t &prev);

  void prec(AFIELD&, AFIELD&);

  double flop_count_intermediate(const int iter);

  /*
#ifdef USE_FACTORY
 private:
  static ASolver<AFIELD> *create_object_with_fopr(AFopr<AFIELD>* fopr)
  { return new ASolver_Richardson<AFIELD>(fopr); }

 public:
  static bool register_factory()
  {
    bool init = ASolver<AFIELD>::Factory_fopr::Register("ASolver_Richardson",
                                                        create_object_with_fopr);
    return init;
  }
#endif
  */
};

#endif // ASOLVER_RICHARDSON_H
