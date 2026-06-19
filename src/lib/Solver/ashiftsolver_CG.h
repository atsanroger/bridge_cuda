/*!
        @file    ashiftsolver_CG.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-06-29 22:43:06 #$

        @version $LastChangedRevision: 2527 $
*/

#ifndef ASHIFTSOLVER_CG_INCLUDED
#define ASHIFTSOLVER_CG_INCLUDED

#include "Solver/ashiftsolver.h"

#include "IO/bridgeIO.h"
using Bridge::vout;


//! Multishift Conjugate Gradient solver.

/*!
                                    [23 Dec 2011  H.Matsufuru]
    unique_ptr is introduced to avoid memory leaks
                                    [21 Mar 2015 Y.Namekawa]
 */

template<typename FIELD, typename FOPR>
class AShiftsolver_CG : public AShiftsolver<FIELD>
{
 public:
  typedef typename FIELD::real_t real_t;

  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  FOPR *m_fopr;

  int m_Niter;
  real_t m_Stop_cond;

  std::vector<FIELD> m_x, m_p;
  FIELD m_r, m_s;
  std::vector<real_t> m_zeta1, m_zeta2, m_csh2, m_pp;

  real_t m_snorm, m_alpha_p, m_beta_p;
  int m_Nshift2;

  real_t m_sigma0;

 public:

  AShiftsolver_CG(FOPR *fopr)
    : m_vl(CommonParameters::Vlevel()),
    m_fopr(fopr) {}

  AShiftsolver_CG(FOPR *fopr, int niter, real_t stop_cond)
    : m_vl(CommonParameters::Vlevel()),
    m_fopr(fopr)
  { set_parameters(niter, stop_cond); }

  ~AShiftsolver_CG() {}

  void set_parameters(const Parameters& params);

  void set_parameters(const int niter, const real_t stop_cond);

  void get_parameters(Parameters& params) const;

  void solve(
    std::vector<FIELD>& solution,
    const std::vector<real_t>& shift,
    const FIELD& source,
    int& Nconv,
    real_t& diff);

  double flop_count();

 private:

  void solve_init(real_t&);

  void solve_step(real_t&);

  void reset_field(const FIELD& b,
                   const std::vector<real_t>& sigma,
                   const int Nshift);
};
#endif
