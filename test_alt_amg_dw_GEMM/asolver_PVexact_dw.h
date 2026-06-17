/*!
      @file    asolver_PVexact_dw.h
      @brief   Exact site-diagonal C_PV^{-1} wrapped as an ASolver.
      @author  Wei-Lun Chen

  Applies the Pauli-Villars operator's "Prec" mode (the exact 5D site-diagonal
  LU inverse, L_inv o U_inv) in one shot, with no iteration. This replaces the
  SAP loose inverse of D_PV (the paper's "B_PV" regime, which develops negative
  real parts) by the exact C_PV^{-1} that the paper recommends for the best-
  conditioned multigrid target operator:

      A = (D_PV C_PV^{-1})^dagger  C^{-1}  D

  Ref: I. Kanamori, W.-L. Chen, H. Matsufuru, "Spectrum of Preconditioned
  Moebius Domain-wall Operators", PoS(LATTICE2024)414.

  C_PV^{-1} is 4d-site-local (it couples only the fifth dimension at a fixed 4d
  site), so it commutes with the SAP even/odd block masking used during the
  coarse-operator construction; the `ieo` argument of the 5-arg solve() is
  therefore irrelevant and ignored.
*/
#ifndef ASOLVER_PVEXACT_DW_INCLUDED
#define ASOLVER_PVEXACT_DW_INCLUDED

#include "lib_alt/Solver/asolver.h"
#include "lib/Fopr/afopr.h"

template<typename AFIELD>
class ASolver_PVexact_dw : public ASolver<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  using InitialGuess = typename ASolver<AFIELD>::InitialGuess;
  static const std::string class_name;

 private:
  AFopr<AFIELD> *m_fopr;  //!< Pauli-Villars fermion operator (provides "Prec").

  ASolver_PVexact_dw() {}

 public:
  explicit ASolver_PVexact_dw(AFopr<AFIELD> *fopr) : m_fopr(fopr) {}
  ~ASolver_PVexact_dw() {}

  // nothing to configure: the inverse is exact and direct.
  void set_parameters(const Parameters&) {}
  void get_parameters(Parameters&) const {}

  void set_init_mode(const InitialGuess) {}

  AFopr<AFIELD> *get_fopr() { return m_fopr; }

  // x = C_PV^{-1} b  (exact, one application of "Prec").
  void solve(AFIELD& x, const AFIELD& b, int& nconv, real_t& diff)
  {
    m_fopr->mult(x, b, "Prec");
    nconv = 1;
    diff  = real_t(0);
  }

  // ieo ignored: C_PV^{-1} is site-local and commutes with eo block masking.
  void solve(AFIELD& x, const AFIELD& b, int& nconv, real_t& diff, int /*ieo*/)
  {
    m_fopr->mult(x, b, "Prec");
    nconv = 1;
    diff  = real_t(0);
  }

  double flop_count() { return 0.0; }
};

template<typename AFIELD>
const std::string ASolver_PVexact_dw<AFIELD>::class_name = "ASolver_PVexact_dw";

#endif // ASOLVER_PVEXACT_DW_INCLUDED
