/*!
        @file    $Id:: asolver_FGMRES.h 258#$

        @brief

        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2#$

        @version $LastChangedRevision: 2585 $
*/


#ifndef ASOLVER_FGMRES_H
#define ASOLVER_FGMRES_H

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
class ASolver_FGMRES : public ASolver<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  typedef typename AFIELD::complex_t complex_t;
  using ASolver<AFIELD>::m_vl;
  static const std::string class_name;
  using InitialGuess = typename ASolver<AFIELD>::InitialGuess;

  struct coeff_t
  {
    std::vector<double> ht;   //  h(i+1,i)
    std::vector<dcomplex> h;  // 2-dim
    std::vector<dcomplex> r;  // 2-dim
    std::vector<dcomplex> y;
    std::vector<double>   c;  // "cos"
    std::vector<dcomplex> s;  // "sin"
    std::vector<dcomplex> b;
  };


 private:

  AFopr<AFIELD> *m_fopr;    //!< fermion operator.
  APrecond<AFIELD> *m_prec; //!< preconditioner.

  int m_Niter;              //!< maximum iteration number.
  int m_N_M;
  real_t m_Stop_cond;       //!< stopping criterion (squared).
  real_t m_Stop_cond2;      //!< stopping criterion for inner solver.

  //!
  AFIELD m_x, m_r, m_w;
  std::vector<AFIELD> m_z;
  std::vector<AFIELD> m_v;

  //  std::complex<real_t>  m_rho_prev, m_alpha;

 protected:

  //! to remember convergence iteration to provide flop count.
  int m_nconv;

  InitialGuess m_initial_mode;

  //! calling constructor without fermion operator is forbidden.
  ASolver_FGMRES() { }

 public:
  //! constructor.
  ASolver_FGMRES(AFopr<AFIELD> *fopr,
                    APrecond<AFIELD> *prec)
    : m_Niter(0), m_Stop_cond(0.0L), m_initial_mode(InitialGuess::RHS)
  {
    m_fopr = fopr;
    m_prec = prec;
    init();
  }

  //! destructor.
  ~ASolver_FGMRES() { tidyup(); }

  //! setting parameters by a Parameter object.
  void set_parameters(const Parameters& params);

  //! setting parameters.
  void set_parameters(const int Niter, const real_t Stop_cond);

  //! setting parameters.
  void set_parameters(const int Niter, const real_t Stop_cond, const InitialGuess init_guess_mode);

  //! setting parameters.
  void set_parameters(const int Niter, const real_t Stop_cond, const InitialGuess init_guess_mode, const int N_M);

  //! setting parameters.
  void set_parameters_GMRES_m(const int N_M);

  void get_parameters(Parameters& params) const;

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

  void solve_init(const AFIELD& b, const AFIELD& xq, real_t& rr, coeff_t& coeff, const InitialGuess);

  void init_coeff(coeff_t &);

  int index_ij(const int i, const int j)
  {
    return i + m_N_M * j;
  }


  void solve_step1(real_t& rr, coeff_t& coeff, const int i);
  void solve_step2(real_t& rr, coeff_t& coeff, const int nr, const AFIELD &b);

  void prec(AFIELD&, AFIELD&);

  double flop_count_intermediate(const int iter);

  /*
#ifdef USE_FACTORY
 private:
  static ASolver<AFIELD> *create_object_with_fopr(AFopr<AFIELD>* fopr)
  { return new ASolver_FGMRES<AFIELD>(fopr); }

 public:
  static bool register_factory()
  {
    bool init = ASolver<AFIELD>::Factory_fopr::Register("ASolver_FGMRES",
                                                        create_object_with_fopr);
    return init;
  }
#endif
  */
};

#endif // ASOLVER_FBICGSTAB_H
