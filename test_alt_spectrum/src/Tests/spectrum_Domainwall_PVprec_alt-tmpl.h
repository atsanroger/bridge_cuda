/*!
      @file    spectrum_Domainwall_PVprec_alt-tmpl.h
      @brief   Definitions for Spectrum_Domainwall_PVprec_alt (linear solves on A).
      @author  Wei-Lun Chen
*/

#include "spectrum_Domainwall_PVprec_alt.h"
#include "test.h"

#include "lib/Field/field_G.h"
#include "lib/Field/field_F.h"
#include "lib/Fopr/afopr.h"
#include "lib/Fopr/fopr.h"
#include "lib/Tools/gammaMatrixSet.h"
#include "lib/Measurements/Fermion/source.h"
#include "lib/Measurements/Fermion/corr2pt_4spinor.h"

#include "lib_alt/Solver/asolver_GMRES_m_Cmplx.h"

#include <cmath>

template<typename AFIELD>
const std::string Spectrum_Domainwall_PVprec_alt<AFIELD>::class_name = "Spectrum_Domainwall_PVprec_alt";

//====================================================================
// Solve  A x = b  directly (diagnostic) on the highly non-normal target
// operator  A = (D_PV C_PV^{-1})^dag C^{-1} D.  GMRES(m) is used because it is
// residual-minimizing and does not break down the way BiCGStab does on a
// non-normal operator.
template<typename AFIELD>
int Spectrum_Domainwall_PVprec_alt<AFIELD>::solve_A(std::string file_params)
{
  typedef AFopr<AFIELD>     AFOPR;
  typedef typename AFIELD::real_t real_t;

  vout.general("\n");
  vout.general("---------------------------------------------------\n");
  vout.general("test: solve A x = b (GMRES(m)) on the MG target operator\n");
  vout.general("parameter file: %s\n", file_params.c_str());

  int Ndim = CommonParameters::Ndim();
  int Nvol = CommonParameters::Nvol();

  params_all = ParameterManager::read(file_params);
  Parameters params_test = params_all.lookup("Test_Eigensolver");
  Parameters params_fopr = params_all.lookup("Fopr");

  string str_vlevel = params_test.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);

  U.reset(new Field_G(Nvol, Ndim));

  string fopr_type = params_fopr.get_string("fermion_type");
  vout.general(m_vl, "fermion_type = %s\n", fopr_type.c_str());
  unique_ptr<AFOPR> fopr(AFOPR::New(fopr_type, params_fopr));

  setup_config(U, params_test);
  fopr->set_config(U.get());
  fopr->set_mode("D");

  int NFin  = fopr->field_nin();
  int NFvol = fopr->field_nvol();
  int NFex  = fopr->field_nex();

  AFIELD b(NFin, NFvol, NFex), x(NFin, NFvol, NFex), r(NFin, NFvol, NFex);
  b.set(0.0);
  b.set(0, 1.0);            // simple localized source
  real_t bnorm2 = b.norm2();

  ASolver_GMRES_m_Cmplx<AFIELD> solver(fopr.get());
  Parameters params_solver;
  params_solver.set_int("maximum_number_of_iteration", 1);
  params_solver.set_int("maximum_number_of_restart", 100);
  params_solver.set_int("number_of_orthonormal_vectors", 30);
  // ~3e-7 residual; kept safely ABOVE the noisy recursive-residual floor (~1e-15
  // rr/bnorm2) to avoid the restarted-GMRES overshoot-to-nan (see solve_propagator).
  params_solver.set_double("convergence_criterion_squared", 1.0e-13);
  params_solver.set_string("verbose_level", "Detailed");
  solver.set_parameters(params_solver);

  int    nconv = -1;
  real_t diff  = -1;
  unique_ptr<Timer> timer(new Timer("solve_A"));
  timer->start();
  solver.solve(x, b, nconv, diff);
  timer->stop();

  // true residual ||A x - b||^2
  fopr->mult(r, x);
  axpy(r, real_t(-1), b);
  real_t trueres2 = r.norm2();

  vout.general(m_vl, "\n");
  vout.general(m_vl, "Solve A x = b (GMRES(m)):\n");
  vout.general(m_vl, "  |b|^2            = %.6e\n", (double)bnorm2);
  vout.general(m_vl, "  nconv (iters)    = %d\n", nconv);
  vout.general(m_vl, "  recursive diff   = %.6e\n", (double)diff);
  vout.general(m_vl, "  true ||Ax-b||^2  = %.6e\n", (double)trueres2);
  timer->report();

  return EXIT_SUCCESS;
}

//====================================================================
// Solve the physical propagator  D psi = eta  through the preconditioned
// operator  A = (D_PV C_PV^{-1})^dag C^{-1} D  =  P_L D  with  P_L =
// (D_PV C_PV^{-1})^dag C^{-1}.  Left-preconditioning D psi = eta by P_L gives
//   A psi = b' ,  b' = P_L eta ,
// and because A = P_L D the GMRES solution is psi = D^{-1} eta EXACTLY -- the
// same propagator a direct CGNE solve of D psi = eta produces, so every derived
// observable (2pt functions, ...) is identical to solver tolerance. We certify
// this by the true residual on the ORIGINAL system, ||D psi - eta|| / ||eta||.
template<typename AFIELD>
int Spectrum_Domainwall_PVprec_alt<AFIELD>::solve_propagator(std::string file_params)
{
  typedef AFopr<AFIELD>     AFOPR;
  typedef typename AFIELD::real_t real_t;

  vout.general("\n");
  vout.general("---------------------------------------------------\n");
  vout.general("test: propagator D psi = eta via preconditioned A (GMRES(m))\n");
  vout.general("parameter file: %s\n", file_params.c_str());

  int Ndim = CommonParameters::Ndim();
  int Nvol = CommonParameters::Nvol();

  params_all = ParameterManager::read(file_params);
  Parameters params_test = params_all.lookup("Test_Eigensolver");
  Parameters params_fopr = params_all.lookup("Fopr");

  string str_vlevel = params_test.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);

  U.reset(new Field_G(Nvol, Ndim));

  // A = (D_PV C_PV^{-1})^dag C^{-1} D  (factory: "Domainwall_PVprec")
  string fopr_type = params_fopr.get_string("fermion_type");
  vout.general(m_vl, "fermion_type = %s\n", fopr_type.c_str());
  unique_ptr<AFOPR> fopr(AFOPR::New(fopr_type, params_fopr));

  setup_config(U, params_test);
  fopr->set_config(U.get());
  fopr->set_mode("D");                 // mult() applies the full A

  int NFin  = fopr->field_nin();
  int NFvol = fopr->field_nvol();
  int NFex  = fopr->field_nex();

  AFIELD eta(NFin, NFvol, NFex), bprime(NFin, NFvol, NFex);
  AFIELD psi(NFin, NFvol, NFex), r(NFin, NFvol, NFex);

  eta.set(0.0);
  eta.set(0, 1.0);                     // simple localized (point) source
  real_t eta2 = eta.norm2();

  // RHS of the preconditioned system: b' = P_L eta
  fopr->mult(bprime, eta, "leftprec");
  real_t bprime2 = bprime.norm2();

  // GMRES(m): robust for the highly non-normal A (no BiCGStab breakdown).
  ASolver_GMRES_m_Cmplx<AFIELD> solver(fopr.get());
  Parameters params_solver;
  params_solver.set_int("maximum_number_of_iteration", 1);
  params_solver.set_int("maximum_number_of_restart", 100);
  params_solver.set_int("number_of_orthonormal_vectors", 30);
  // ||r||^2 < 1e-13  (||r|| < ~3e-7): propagator-grade tolerance kept ABOVE the
  // achievable recursive-residual floor (~1e-15 rr/bnorm2 WITH reorthogonalization).
  // The floor is noisy run to run (GPU reductions are non-deterministic) and a
  // single GMRES(N_M) cycle can drop the residual by many orders, so a criterion
  // AT the floor (e.g. 1e-16) sometimes overshoots into it and the next Arnoldi
  // step divides by ~0 -> nan. 1e-13 leaves ~2 orders of margin (crosses near iter
  // ~1900, safely before the nan onset ~2250) and is more than enough for 2pt.
  params_solver.set_double("convergence_criterion_squared", 1.0e-13);
  params_solver.set_string("verbose_level", "Detailed");
  solver.set_parameters(params_solver);

  int    nconv = -1;
  real_t diff  = -1;
  unique_ptr<Timer> timer(new Timer("solve_propagator"));
  timer->start();
  solver.solve(psi, bprime, nconv, diff);   // A psi = b'  ->  psi = D^{-1} eta
  timer->stop();

  // certify on the ORIGINAL system: r = D psi - eta
  fopr->mult(r, psi, "physicalD");
  axpy(r, real_t(-1), eta);
  real_t res2 = r.norm2();

  vout.general(m_vl, "\n");
  vout.general(m_vl, "Propagator via preconditioned A (GMRES(m)):\n");
  vout.general(m_vl, "  |eta|^2              = %.6e\n", (double)eta2);
  vout.general(m_vl, "  |b' = P_L eta|^2     = %.6e\n", (double)bprime2);
  vout.general(m_vl, "  nconv (iters)        = %d\n", nconv);
  vout.general(m_vl, "  recursive |Apsi-b'|  = %.6e\n", (double)diff);
  vout.general(m_vl, "  true ||D psi-eta||^2 = %.6e\n", (double)res2);
  vout.general(m_vl, "  true ||D psi-eta||/||eta|| = %.6e\n",
               (double)std::sqrt(res2 / eta2));
  vout.general(m_vl, "  -> psi = D^{-1} eta (same propagator/2pt as CGNE) if "
               "this relative residual is at the solver tolerance.\n");
  timer->report();

  return EXIT_SUCCESS;
}

//====================================================================
// Full hadron 2pt function through the preconditioned operator A.
//
// This is a clone of Spectrum_Domainwall_alt::hadron_2ptFunction, with the
// quark inversion replaced by the preconditioned-A GMRES solve. For each of the
// 12 color-spin point sources eta:
//   1. build the 5d source b5 = D_- (P_+ 0 .. 0 P_-)^T b4   (Wilson kernel),
//   2. convert b5 -> AFIELD eta, form b' = P_L eta, GMRES-solve A psi = b'
//      so that psi = D^{-1} eta EXACTLY (A = P_L D), certify ||D psi - eta||,
//   3. reverse psi -> x5, project the 5d solution to the 4d propagator sq[idx].
// Then Corr2pt_4spinor::meson_all(sq, sq) gives the meson correlators. Because
// each psi equals the CGNE propagator to solver tolerance, the resulting 2pt
// function is identical (to tolerance) to a direct CGNE measurement.
template<typename AFIELD>
int Spectrum_Domainwall_PVprec_alt<AFIELD>::solve_2pt(std::string file_params)
{
  typedef AFopr<AFIELD>     AFOPR;
  typedef typename AFIELD::real_t real_t;

  const std::string test_name = class_name + ".solve_2pt";

  vout.general("\n");
  vout.general("---------------------------------------------------\n");
  vout.general("test: hadron 2pt function via preconditioned A (GMRES(m))\n");
  vout.general("parameter file: %s\n", file_params.c_str());

  int Nc   = CommonParameters::Nc();
  int Nd   = CommonParameters::Nd();
  int Ndim = CommonParameters::Ndim();
  int Nvol = CommonParameters::Nvol();

  params_all = ParameterManager::read(file_params);
  Parameters params_test   = params_all.lookup("Test_Eigensolver");
  Parameters params_fopr   = params_all.lookup("Fopr");
  Parameters params_source = params_all.lookup("Source");

  string str_vlevel = params_test.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);

  const bool   do_check        = params_test.is_set("expected_result");
  const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

  // gauge configuration
  U.reset(new Field_G(Nvol, Ndim));
  setup_config(U, params_test);

  // gamma matrix set (for the meson correlators)
  string gmset_type = params_fopr.get_string("gamma_matrix_type");
  unique_ptr<GammaMatrixSet> gmset(GammaMatrixSet::New(gmset_type));

  // Wilson kernel operator, for the 4d <-> 5d source/sink construction.
  string kernel_type = params_fopr.get_string("kernel_type");
  double M0          = params_fopr.get_double("domain_wall_height");
  double coeff_c     = params_fopr.get_double("coefficient_c");
  double kappa           = 1.0 / (8.0 - 2.0 * M0);
  double one_over_2kappa = 4.0 - M0;
  Parameters params_kernel = params_fopr;
  params_kernel.set_double("hopping_parameter", kappa);
  unique_ptr<Fopr> foprw(Fopr::New(kernel_type, params_kernel));
  foprw->set_mode("D");
  foprw->set_config(U.get());

  // A = (D_PV C_PV^{-1})^dag C^{-1} D  (factory: "Domainwall_PVprec")
  string fopr_type = params_fopr.get_string("fermion_type");
  vout.general(m_vl, "fermion_type = %s\n", fopr_type.c_str());
  unique_ptr<AFOPR> fopr(AFOPR::New(fopr_type, params_fopr));
  fopr->set_config(U.get());
  fopr->set_mode("D");

  int Ns = params_fopr.get_int("extent_of_5th_dimension");

  // alt-field working vectors
  int NFin  = fopr->field_nin();
  int NFvol = fopr->field_nvol();
  int NFex  = fopr->field_nex();
  AFIELD eta(NFin, NFvol, NFex), bprime(NFin, NFvol, NFex);
  AFIELD psi(NFin, NFvol, NFex), rchk(NFin, NFvol, NFex);
  AFIELD delta(NFin, NFvol, NFex), rphys(NFin, NFvol, NFex);

  // outer iterative-refinement controls (see the per-source loop below)
  const int    max_refine    = 8;
  const double refine_target2 = 1.0e-18;   // ||D psi - eta||^2/||eta||^2 (rel ~1e-9)

  // GMRES(m): identical settings to solve_propagator (robust on non-normal A).
  ASolver_GMRES_m_Cmplx<AFIELD> solver(fopr.get());
  Parameters params_solver;
  params_solver.set_int("maximum_number_of_iteration", 1);
  params_solver.set_int("maximum_number_of_restart", 100);
  params_solver.set_int("number_of_orthonormal_vectors", 30);
  params_solver.set_double("convergence_criterion_squared", 1.0e-13);
  params_solver.set_string("verbose_level", "General");
  solver.set_parameters(params_solver);

  // source generator
  string str_source_type = params_source.get_string("source_type");
  vout.general(m_vl, "source_type = %s\n", str_source_type.c_str());
  unique_ptr<Source> source(Source::New(str_source_type));
  source->set_parameters(params_source);

  // 4d quark propagators (12 color-spin columns): boundary-projected sq (the
  // physical surface field) and MIDPOINT-projected sq_mid (the 5th-dim midpoint
  // field, used for the residual mass). For Ns even the midpoint cut is between
  // s = Ns/2-1 and s = Ns/2.
  std::vector<Field_F> sq(Nc * Nd), sq_mid(Nc * Nd);
  for (int i = 0; i < Nc * Nd; ++i) { sq[i].set(0.0); sq_mid[i].set(0.0); }
  const int s_mid_lo = Ns / 2 - 1;   // lower-half top slice
  const int s_mid_hi = Ns / 2;       // upper-half bottom slice

  Field_F b(Nvol, 1), vt1(Nvol, 1), vt2(Nvol, 1);
  Field_F b5(Nvol, Ns), x5(Nvol, Ns);

  unique_ptr<Timer> timer(new Timer(test_name));
  timer->start();

  vout.general(m_vl, "\n");
  vout.general(m_vl, "Solving quark propagator through A (with refinement):\n");
  vout.general(m_vl, "  color spin   Nconv   nref   ||D psi-eta||/||eta||\n");

  for (int ispin = 0; ispin < Nd; ++ispin) {
    for (int icolor = 0; icolor < Nc; ++icolor) {
      int idx = icolor + Nc * ispin;

      b.set(0.0);
      source->set(b, idx);

      // 5d source: b5 = D_- ( P_+ 0 .... 0 P_- )^T b4   (same as the reference)
#pragma omp parallel
      {
        b5.set(0.0);

        // s5 = 0 : 0.5 (cD-1)(1+gm5) b
        foprw->mult_gm5(vt1, b);
        axpy(vt1, +1.0, b);
        foprw->mult(vt2, vt1);
        aypx(-one_over_2kappa * coeff_c, vt2, vt1);
        scal(vt2, -0.5);
        copy(b5, 0, vt2, 0);

        // s5 = Ns-1 : 0.5 (cD-1)(1-gm5) b
        foprw->mult_gm5(vt1, b);
        axpy(vt1, -1.0, b);
        foprw->mult(vt2, vt1);
        aypx(-one_over_2kappa * coeff_c, vt2, vt1);
        scal(vt2, 0.5);
        copy(b5, Ns - 1, vt2, 0);
      }

      // invert through A with iterative refinement on the PHYSICAL residual.
      //
      // GMRES on A = P_L D minimizes ||P_L(D psi - eta)||, NOT ||D psi - eta||;
      // for modes where P_L (cond ~2e4) amplifies, the physical residual lags
      // the preconditioned one by the conditioning factor (some 2pt sources
      // stalled at ~1e-3 even after the inner solve "converged"). A few steps of
      // outer Richardson refinement on the physical system D psi = eta drive the
      // true residual down geometrically: each step solves A delta = P_L r_phys
      // (delta ~ D^{-1} r_phys) and updates psi += delta. This is exactly the
      // left-preconditioned-defect-correction wrapping of the A solve.
      fopr->convert(eta, (Field&)b5);
      real_t eta2 = eta.norm2();

      psi.set(0.0);
      copy(rphys, eta);                 // r_phys = eta - D psi  (psi = 0)
      real_t res2     = eta2;
      int    nconv_tot = 0;
      int    nref      = 0;
      for (nref = 0; nref < max_refine; ++nref) {
        fopr->mult(bprime, rphys, "leftprec");   // b' = P_L r_phys

        int    nconv = -1;
        real_t diff  = -1;
        delta.set(0.0);
        solver.solve(delta, bprime, nconv, diff);  // A delta = b' -> delta ~ D^{-1} r_phys
        nconv_tot += nconv;

        axpy(psi, real_t(1), delta);               // psi += delta

        fopr->mult(rchk, psi, "physicalD");        // D psi
        copy(rphys, eta);
        axpy(rphys, real_t(-1), rchk);             // r_phys = eta - D psi
        res2 = rphys.norm2();

        if (res2 < refine_target2 * eta2) { ++nref; break; }
      }

      vout.general(m_vl, "   %2d   %2d   %6d  %3d   %12.4e\n",
                   icolor, ispin, nconv_tot, nref,
                   (double)std::sqrt(res2 / eta2));

      // reverse psi -> x5, then project the 5d solution to the 4d propagator
      fopr->reverse((Field&)x5, psi);

#pragma omp parallel
      {
        copy(vt1, 0, x5, 0);
        copy(vt2, 0, x5, Ns - 1);
        copy(sq[idx], vt1);
#pragma omp barrier
        axpy(sq[idx], 1.0, vt2);    // (x5[0] + x5[Ns-1])
        axpy(vt1, -1.0, vt2);       // (x5[0] - x5[Ns-1])
#pragma omp barrier
        foprw->mult_gm5(vt2, vt1);  // gm5 (x5[0] - x5[Ns-1])
        axpy(sq[idx], -1.0, vt2);   // (1-gm5) x5[0] + (1+gm5) x5[Ns-1]
        scal(sq[idx], 0.5);
      }

      // midpoint-projected propagator (for the residual mass): the SAME
      // reconstruction at the s = Ns/2-1 / Ns/2 cut --
      //   sq_mid = 0.5[(1-gm5) x5[Ns/2] + (1+gm5) x5[Ns/2-1]]
      //          = P_- x5[Ns/2] + P_+ x5[Ns/2-1].
#pragma omp parallel
      {
        copy(vt1, 0, x5, s_mid_hi);
        copy(vt2, 0, x5, s_mid_lo);
        copy(sq_mid[idx], vt1);
#pragma omp barrier
        axpy(sq_mid[idx], 1.0, vt2);   // (x5[Ns/2] + x5[Ns/2-1])
        axpy(vt1, -1.0, vt2);          // (x5[Ns/2] - x5[Ns/2-1])
#pragma omp barrier
        foprw->mult_gm5(vt2, vt1);     // gm5 (x5[Ns/2] - x5[Ns/2-1])
        axpy(sq_mid[idx], -1.0, vt2);  // (1-gm5) x5[Ns/2] + (1+gm5) x5[Ns/2-1]
        scal(sq_mid[idx], 0.5);
      }
    }
  }

  // meson correlators
  Corr2pt_4spinor corr(gmset.get());
  corr.set_parameters(params_all.lookup("Corr2pt_4spinor"));

  vout.general(m_vl, "\n");
  vout.general(m_vl, "2-point correlator (propagator solved through A):\n");

  double result = corr.meson_all(sq, sq);

  vout.general(m_vl, "RESULT_PVprec_2pt: %.17e\n", result);

  // ---- residual mass --------------------------------------------------------
  // m_res(t) = <J5q(t) P(0)> / <P(t) P(0)>  =  C_PP(sq_mid,sq_mid)(t)
  //                                            / C_PP(sq,sq)(t),
  // both pseudoscalar (gamma5 x gamma5) correlators from the SAME boundary
  // point source. P uses the boundary surface field, J5q the 5th-dim midpoint
  // field; their ratio is the residual chiral-symmetry breaking from the finite
  // bulk-boundary coupling at this Ns (and, once parameter_alpha is wired in,
  // the diagnostic of how alpha retunes that coupling). m_res is the large-t
  // plateau.
  GammaMatrix gm5 = gmset->get_GM(gmset->GAMMA5);
  const int Lt = CommonParameters::Lt();
  std::vector<dcomplex> c_PP(Lt), c_J5q(Lt);   // meson_correlator needs them pre-sized
  corr.meson_correlator(c_PP,  gm5, gm5, sq,     sq);       // <P  P>   (denominator)
  corr.meson_correlator(c_J5q, gm5, gm5, sq_mid, sq_mid);   // <J5q J5q>(numerator)
  vout.general(m_vl, "\n");
  vout.general(m_vl, "Residual mass (J5q/PP ratio):\n");
  vout.general(m_vl, "    t      C_PP(t)         C_J5q(t)        m_res(t)\n");
  // plateau average over the central half of the time extent (avoid the
  // contaminated ends near the source and its periodic image).
  double mres_sum = 0.0; int mres_n = 0;
  for (int t = 0; t < Lt; ++t) {
    double cpp  = real(c_PP[t]);
    double cj5q = real(c_J5q[t]);
    double r    = (cpp != 0.0) ? cj5q / cpp : 0.0;
    vout.general(m_vl, " %4d  %15.7e  %15.7e  %15.7e\n", t, cpp, cj5q, r);
    if (t >= Lt / 4 && t < 3 * Lt / 4) { mres_sum += r; ++mres_n; }
  }
  double m_res = (mres_n > 0) ? mres_sum / mres_n : 0.0;
  vout.general(m_vl, "RESULT_PVprec_mres (plateau avg over t in [%d,%d)): %.10e\n",
               Lt / 4, 3 * Lt / 4, m_res);
  // ---------------------------------------------------------------------------

  timer->stop();
  timer->report();

  if (do_check) {
    // Verify the A-based 2pt against the expected value read from the yaml.
    // The GMRES+refinement propagator reaches ||D psi-eta||/||eta|| ~1e-9, so a
    // wrong propagator would differ by O(1); check at 1e-6 (well above the solver
    // precision, far below any physical-difference scale). The CGNE reference in
    // Spectrum_Domainwall_alt verifies the SAME expected value at the default
    // 1e-11, so PVprec-vs-CGNE agreement to ~1e-9 is the cross-check.
    return Test::verify(result, expected_result, 1.0e-6);
  } else {
    vout.detailed(m_vl, "check skipped: expected_result not set.\n\n");
    return EXIT_SKIP;
  }
}
//============================================================END=====
