/*!
    @file    spectrum_Domainwall_PVprec_alt.h
    @brief   Linear solves on the PV-preconditioned Moebius domain-wall
             operator  A = (D_PV C_PV^{-1})^dagger C^{-1} D.
    @author  Wei-Lun Chen

    @details
      Separated from eigenvalue_alt (which is for the spectrum/eigenvalue
      study). This class hosts the *solver* side built on the same A operator:

        - solve_A          : solve  A x = b  directly (diagnostic).
        - solve_propagator : solve the physical propagator  D psi = eta  through
                             A, i.e. form  b' = P_L eta  and GMRES-solve
                             A psi = b'  so that  psi = D^{-1} eta  exactly.

      A is the best-conditioned multigrid target of
        I. Kanamori, W.-L. Chen, H. Matsufuru,
        "Spectrum of Preconditioned Moebius Domain-wall Operators",
        PoS(LATTICE2024)414,
      so this file is the natural place to grow the multigrid-on-A solver next.
*/

#ifndef SPECTRUM_DOMAINWALL_PVPREC_ALT_INCLUDED
#define SPECTRUM_DOMAINWALL_PVPREC_ALT_INCLUDED

#include <vector>
#include <string>

#include "spectrum_alt.h"

#include "lib_alt/alt_impl.h"

#include "lib/Parameters/commonParameters.h"
#include "lib/Parameters/parameters.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "job_Utils.h"

template<typename AFIELD>
class Spectrum_Domainwall_PVprec_alt : public Spectrum_alt
{
 public:
  static const std::string class_name;

 protected:
  using Spectrum_alt::m_vl;

 public:
  //! constructor
  Spectrum_Domainwall_PVprec_alt() : Spectrum_alt() { init(); }

  //! destructor
  ~Spectrum_Domainwall_PVprec_alt() {}

  //! solve A x = b with GMRES(m) (factory fopr from the "Fopr" section).
  //! Diagnostic: exercises A directly with a localized source.
  int solve_A(const std::string file_params);

  //! solve the physical propagator  D psi = eta  via the preconditioned
  //! operator A: form  b' = P_L eta, GMRES-solve  A psi = b', and certify
  //! ||D psi - eta||/||eta||.  psi == D^{-1} eta == the CGNE solution, so all
  //! observables (e.g. 2pt functions) coincide with a direct CGNE solve.
  int solve_propagator(const std::string file_params);

  //! full hadron 2pt function through the preconditioned operator A.
  //! For each of the 12 color-spin point sources: build the 5d source, solve
  //! the physical propagator D psi = eta via  A psi = P_L eta  (psi = D^{-1}eta
  //! exactly), reconstruct the 4d quark propagator, then form the meson
  //! correlators with Corr2pt_4spinor -- the SAME 2pt a direct CGNE solve gives.
  //! Mirrors Spectrum_Domainwall_alt::hadron_2ptFunction, swapping only the
  //! inversion (CGNE -> preconditioned-A GMRES).
  int solve_2pt(const std::string file_params);

 private:
  //! initial setup
  void init() {}
};
#endif // SPECTRUM_DOMAINWALL_PVPREC_ALT_INCLUDED
