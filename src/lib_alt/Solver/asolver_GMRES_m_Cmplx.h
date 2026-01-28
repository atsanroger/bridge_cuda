/*!
        @file    asolver_GMRES_m_Cmplx.h

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2024-03-14 11:21:44 #$

        @version $LastChangedRevision: 2578 $
*/

#ifndef ASOLVER_GMRES_m_CMPLX_INCLUDED
#define ASOLVER_GMRES_m_CMPLX_INCLUDED

#include <vector>
#include <string>

#include "asolver.h"
#include "lib/Fopr/afopr.h"

//! GMRES(m) algorithm with complex variables.

/*!
   This class implements GMRES(m) algorithm for nonhermitian matrix.
   The product of vectors is treated in complex.
   See Y.Saad and M.H.Schultz, SIAM J.Sci.Stat.Comput. 7 (1986) 856.
                                   [8 Aug 2012  Y.Namekawa]
    (Coding history will be recovered from trac.)
    YAML is implemented.           [14 Nov 2012 Y.Namekawa]
    Multi-threaded.                [17 Jul 2014 Y.Namekawa]
    Introduce unique_ptr to avoid memory leaks.
                                   [21 Mar 2015 Y.Namekawa]
    Add restart.                   [22 Feb 2016 Y.Namekawa]
    Add flop_count.                [ 8 Aug 2016 Y.Namekawa]
    Add use_init_guess.            [ 7 Jul 2017 Y.Namekawa]
    Fix a bug of h[ij], which is kindly reported by Ji-Chong Yang.
                                   [ 7 Mar 2019 Y.Namekawa]

    AField version                 [20 Feb 2024 I.Kanamori]
*/

template<typename AFIELD>
class ASolver_GMRES_m_Cmplx : public ASolver<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  typedef typename AFIELD::complex_t complex_t;
  using ASolver<AFIELD>::m_vl;
  static const std::string class_name;
  using InitialGuess = typename ASolver<AFIELD>::InitialGuess;

 private:

  AFopr<AFIELD> *m_fopr; //!< fermion operator.

  int m_Niter;           //!< maximum iteration number.
  int m_Nrestart;
  real_t m_Stop_cond;    //!< stopping criterion (squared).

  int m_N_M;

  int m_Nrestart_count;
  int m_Nconv_count;

  //! mode switch for initial guess
  InitialGuess m_initial_mode;
  //- working area
  struct coeff_t
  {
    real_t beta_prev;
    std::vector<dcomplex> y;
    std::vector<dcomplex> h;
    std::vector<dcomplex> g;
  }
  coeff;
  //  double m_beta_prev;

  std::vector<AFIELD> m_v;
  AFIELD m_s, m_r, m_x, m_v_tmp;


  //! calling constructor without fermion operator is forbidden.
  ASolver_GMRES_m_Cmplx() { }

public:
  ASolver_GMRES_m_Cmplx(AFopr<AFIELD> *fopr)
    : m_fopr(fopr)
  {
    init();
  }

  ASolver_GMRES_m_Cmplx(AFopr<AFIELD> *fopr, const Parameters& params)
    : m_fopr(fopr)
  {
    init();
    set_parameters(params);
  }

  ~ASolver_GMRES_m_Cmplx() {}

  void set_parameters(const Parameters& params);

  void get_parameters(Parameters& params) const;

  void solve(AFIELD& solution, const AFIELD& source, int& Nconv, real_t& diff);

  AFopr<AFIELD> *get_fopr() { return m_fopr; }

  double flop_count();

 private:
  void init();

  void tidyup();

  void set_parameters(const int Niter, const int Nrestart, const real_t Stop_cond);

  void set_parameters(const int Niter, const int Nrestart, const real_t Stop_cond, const InitialGuess init_guess_mode);

  void set_parameters_GMRES_m(const int N_M);

  void set_parameters(const int Niter, const int Nrestart, const real_t Stop_cond, const InitialGuess init_guess_mode, const int N_M);


  void set_init_mode(const InitialGuess init_guess)
  {
    m_initial_mode = init_guess;
  }


  //  int fetch_initial_guess(const Parameters &, const std::string, InitialGuess &);

  //  void reset_field(const AFIELD&);

  void solve_init(const AFIELD&, const AFIELD&, real_t&, coeff_t&, const InitialGuess);
  void solve_step(const AFIELD&, real_t&, coeff_t&);

  //  void innerprod_c(double& prod_r, double& prod_i,
  //                   const Field& v, const Field& w);

  void min_J(coeff_t& coeff);

  int index_ij(const int i, const int j)
  {
    return i + (m_N_M + 1) * j;
  }

#ifdef USE_FACTORY
 private:
  static ASolver<AFIELD> *create_object_with_fopr(AFopr<AFIELD> *fopr)
  {
    return new ASolver_GMRES_m_Cmplx<AFIELD>(fopr);
  }

 public:
  static bool register_factory()
  {
    bool init = ASolver<AFIELD>::Factory_fopr::Register("GMRES_m_Cmplx",
                                                        create_object_with_fopr);
    return init;
  }
#endif
};
#endif // include guard
