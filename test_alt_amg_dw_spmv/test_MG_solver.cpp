/*!

        @file    test_MG_solver.cpp

        @brief   implementation of the test for multigrid solver (domainwall)

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: #$

        @version $LastChangedRevision: 2595 $
*/
//====================================================================

#include "test_MG_solver.h"

#include "lib/IO/gaugeConfig.h"
#include "lib/Tools/randomNumberManager.h"
#include "lib/Tools/timer.h"
#include "lib/Fopr/fopr.h"
#include "lib/Field/field_F.h"
#include "lib/Field/field_G.h"
#include "lib/Field/shiftField_lex.h"
#include "lib/Parameters/parameterManager.h"
#include "lib/Tools/gammaMatrixSet.h"
#include "lib/Measurements/Fermion/source.h"
#include "lib/Measurements/Fermion/corr2pt_4spinor.h"

#ifdef USE_ALT_QXS
#include "lib_alt_QXS/bridge_alt_qxs.h"
#include "lib_alt_QXS/Field/aindex_block_lex.h"
#include "lib_alt_QXS/Field/aindex_coarse_lex.h"
#include "lib_alt_QXS/Fopr/afopr_Domainwall_5din.h"
#include "lib_alt_QXS/Fopr/afopr_Domainwall_5din_dd.h"
#include "lib_alt_QXS/Fopr/afopr_Domainwall_coarse.h"
const Impl IMPL = QXS;
#endif

#ifdef USE_ALT_ACCEL
#include "lib_alt_Accel/bridge_alt_accel.h"
#include "lib_alt_Accel/Field/aindex_block_lex.h"
#include "lib_alt_Accel/Field/aindex_coarse_lex.h"
#include "lib_alt_Accel/Fopr/afopr_Domainwall_5din.h"
#include "lib_alt_Accel/Fopr/afopr_Domainwall_5din_dd.h"
#include "lib_alt_Accel/Fopr/afopr_Domainwall_coarse.h"
#include "lib_alt_Accel/Fopr/afopr_Domainwall_PVprec.h"
const Impl IMPL = ACCEL;
#endif

#include "lib_alt/Solver/asolver.h"
#include "lib_alt/Solver/asolver_BiCGStab.h"
#include "lib_alt/Solver/asolver_CG.h"
#include "lib_alt/Solver/asolver_CGNR.h"
#include "lib_alt/Solver/aprecond_Mixedprec.h"
#include "lib_alt/Solver/asolver_BiCGStab_Precond.h"
#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_eo_Richardson.h"

#include "asolver_MG_dw.h"

#include "asolver_SAP_dw.h"
#include "asolver_SAP_MINRES_dw.h"


//====================================================================
// Registration as a test
namespace test_MG_solver {
  const string test_name = "Test_MG_Solver";

  int test(void)
  {
    Test_MG_Solver test;
    return test.test();
  }
}

// class name
const std::string Test_MG_Solver::class_name = "Test_MG_Solver";

//====================================================================
void Test_MG_Solver::init()
{
  // do nothing.
}


//====================================================================
namespace {
  Bridge::VerboseLevel vl;


  typedef AField<double, IMPL>   AFIELD_d;
  typedef AField<float, IMPL>    AFIELD_f;


//====================================================================
  int run_bicgstab(AFIELD_d& ax, const AFIELD_d& ab,
                   Field_G *U,
                   const Parameters params_fopr,
                   const Parameters params_outer_solver)
  {
    vout.general("running BiCGstab (Cmplx) solver\n");
    return 0;
  }


//====================================================================
  int run_mixed(AFIELD_d& ax, const AFIELD_d& ab,
                Field_G *U,
                const Parameters params_fopr,
                const Parameters params_outer_solver)
  {
    vout.general("running mixed solver with FBiCGStab\n");
    return 0;
  }


//====================================================================
  int run_mixed_eo(AFIELD_d& ax, const AFIELD_d& ab,
                Field_G *U,
                const Parameters &params_fopr,
                const Parameters &params_outer_solver)
  {
    vout.general("running even-odd mixed solver with Richardson\n");

    Parameters params_fopr_eo=params_fopr;
    Parameters params_solver=params_outer_solver;
    params_solver.set_string("solver_type", "CGNR");
    Fprop_alt_Standard_eo_Richardson<AFIELD_d, AFIELD_f> fprop(
                                     params_fopr_eo, params_solver);

    fprop.set_mode("D");
    fprop.set_config(U);

    vout.general(vl, "setup finished.\n");
    int               nconv = -1;
    double            diff  = -1.0;
    unique_ptr<Timer> timer_solve(new Timer("solve [MixedPrec_eo]"));
    timer_solve->start();

    fprop.invert(ax, ab, nconv, diff);

    timer_solve->stop();
    timer_solve->report();

    fprop.report_performance();

    return EXIT_SUCCESS;

  }

//====================================================================
  int run_MG(AFIELD_d& ax, const AFIELD_d& ab,
             Field_G *U,
             const Parameters &params_fopr,
             const Parameters &params_solver,
             const Parameters &params_all)
  {
    vout.general("running MG solver\n");

    // ####  parameter setup  ####
    const int Nc   = CommonParameters::Nc();
    const int Nd   = CommonParameters::Nd();
    const int Ndim = CommonParameters::Ndim();
    const int Nvol = CommonParameters::Nvol();

    // outer solver
    const int    max_niter = params_solver.get_int("maximum_number_of_iteration");
    const double stop_cond = params_solver.get_double("convergence_criterion_squared");

    // level1
    Parameters params_coarse  = params_solver.lookup("MultiGrid_Level1");
    const std::vector<int> sap_block_size = params_coarse.get_int_vector("sap_block");
    const int  num_vectors    = params_coarse.get_int("setup_number_of_vectors");
    const int  setup_niter    = params_coarse.get_int("setup_number_of_step");

    // coarse grid solver
    const int  coarse_niter     = params_coarse.get_int("maximum_number_of_iteration");
    const float coarse_stop_cond = params_coarse.get_double("convergence_criterion_squared");

    // smoother
    const int   smoother_niter     = params_coarse.get_int("smoother_number_of_iteration");
    const float smoother_stop_cond = params_coarse.get_double("smoother_convergence_criterion_squared");
    const string vlevel             = params_coarse.get_string("verbose_level");

    vout.general(vl, "  sap_block_size = %s\n", Parameters::to_string(sap_block_size).c_str());
    vout.general(vl, "  setup_niter  = %d\n", setup_niter);
    vout.general(vl, "  coarse_niter = %d\n", coarse_niter);
    vout.general(vl, "  coarse_stop_cond = %e\n", coarse_stop_cond);
    vout.general(vl, "  smoother_niter = %d\n", smoother_niter);
    vout.general(vl, "  smoother_stop_cond = %e\n", smoother_stop_cond);
    vout.general(vl, "  num_vectors  = %d\n", num_vectors);
    vout.general(vl, "  outer_niter  = %d\n", max_niter);
    vout.general(vl, "  outer_stop_cond  = %e\n", stop_cond);

    using complexF_t = typename AFIELD_f::complex_t;

    // OUTER operator (double): the positive MG target A = (D_PV C_PV^{-1})^dag
    // C^{-1} D, so the outer solve is the validated GMRES-on-A (bare D is non-
    // normal -> diverges). Built directly from the same domain-wall params
    // (PVprec internally makes D and D_PV with quark_mass_PauliVillars).
    unique_ptr<AFopr_Domainwall_PVprec<AFIELD_d> > afopr_fineD(
                       new AFopr_Domainwall_PVprec<AFIELD_d> (params_fopr));
    afopr_fineD->set_config(U);
    afopr_fineD->set_mode("D");                 // mult() applies the full A
    vout.general(vl, "fine grid operator A=(D_PV C_PV^-1)^dag C^-1 D (double) is ready\n");

    // fine grid operator: single prec.
    Parameters params_fopr_dd = params_fopr;
    params_fopr_dd.set_int_vector("sap_block", sap_block_size);
    params_fopr_dd.set_int_vector("block_size", sap_block_size);

    unique_ptr<AFopr_Domainwall_5din_dd<AFIELD_f> > afopr_fineF(
                 new AFopr_Domainwall_5din_dd<AFIELD_f> (params_fopr_dd));

    vout.general(vl, "fine grid operator (float) constructed\n");
    afopr_fineF->set_config(U);
    vout.general(vl, "fine grid operator (float) is ready\n");

    // PV operator
    const double m_PV= params_fopr.get_double("quark_mass_PauliVillars");
    Parameters params_fopr_PV = params_fopr_dd;
    params_fopr_PV.set_double("quark_mass", m_PV);
    unique_ptr<AFopr_Domainwall_5din_dd<AFIELD_f> > afopr_PV(
                 new AFopr_Domainwall_5din_dd<AFIELD_f> (params_fopr_PV));
    vout.general(vl, "PV operator (float) constructed\n");
    afopr_PV->set_config(U);
    vout.general(vl, "PV operator (float) is ready\n");

    // float A = (D_PV C_PV^{-1})^dag C^{-1} D for the clean two-grid V-cycle
    // (smoother = GMRES-on-A, residual r = w - A x). Same params as the outer.
    unique_ptr<AFopr_Domainwall_PVprec<AFIELD_f> > afopr_A_F(
                 new AFopr_Domainwall_PVprec<AFIELD_f> (params_fopr));
    afopr_A_F->set_config(U);
    afopr_A_F->set_mode("D");
    vout.general(vl, "float A operator (PVprec) is ready\n");

    // multigrid solver
    unique_ptr<ASolver_MG_dw<AFIELD_d> > asolver_mg(
                                            new ASolver_MG_dw<AFIELD_d> );

    asolver_mg->set_parameters(params_solver);
    // set fopr and parameters
    asolver_mg->set_foprD(afopr_fineD.get());
    asolver_mg->set_foprF(afopr_fineF.get());
    asolver_mg->set_fopr_PV(afopr_PV.get());
    asolver_mg->set_fopr_A(afopr_A_F.get());

    unique_ptr<AFopr_dd<AFIELD_f> > afopr_smoother;
    if (asolver_mg->use_fopr_for_smoother()) {
      afopr_smoother.reset(asolver_mg->new_fopr_smoother(params_fopr_dd));
      afopr_smoother->set_config(U);
      asolver_mg->set_fopr_smoother(afopr_smoother.get());
    }

    asolver_mg->init_solver();

    // setup null space vectors
    asolver_mg->run_setup();

    // Left-precondition the physical system D x = b by P_L so the outer solves
    // the positive operator:  A x = b' , b' = P_L b  =>  x = D^{-1} b exactly
    // (same recipe as the spectrum GMRES-on-A). P_L = (D_PV C_PV^{-1})^dag C^{-1}.
    AFIELD_d ab_prime;
    ab_prime.reset(ab.nin(), ab.nvol(), ab.nex());
    afopr_fineD->mult(ab_prime, ab, "leftprec");
    vout.general(vl, "|b' = P_L b|^2 = %.8e\n", (double)ab_prime.norm2());

    int               nconv = 0;
    double            diff  = -1.0;
    unique_ptr<Timer> timer_solve(new Timer("solve [MG, solve only]"));
    timer_solve->start();
#pragma omp parallel
    {
      asolver_mg->solve(ax, ab_prime, nconv, diff);
    }
    timer_solve->stop();
    timer_solve->report();
    double etime = timer_solve->elapsed_msec();

    double flop_solve = asolver_mg->flop_count();
    vout.general(" Flops (double+float) : %f GFlop/sec\n",
                 1.0e-6 * flop_solve / etime);

    // ================================================================
    // Hadron 2pt function through the AMG solve.
    //
    // Reuses the already-built A operator (afopr_fineD, the PVprec A) and the
    // AMG solver (asolver_mg). The propagator is solved EXACTLY as in the
    // spectrum GMRES-on-A 2pt (spectrum_Domainwall_PVprec_alt::solve_2pt):
    // the SAME chiral 5d source/sink construction and the SAME outer
    // iterative-refinement on the PHYSICAL residual (each step solves
    // A delta = P_L r_phys, psi += delta, r_phys = eta - D psi). The ONLY
    // difference vs that reference is the inner solver: AMG V-cycle instead of
    // unpreconditioned GMRES. The resulting C_PP / m_res therefore cross-check
    // directly against the CGNE and GMRES-on-A measurements on the same setup.
    // ================================================================
    {
      using real_t = AFIELD_d::real_t;   // double
      const int Nvol = CommonParameters::Nvol();
      const int Lt   = CommonParameters::Lt();
      const int Ns   = params_fopr.get_int("extent_of_5th_dimension");

      const Parameters params_source = params_all.lookup("Source");

      // gamma matrix set (for the meson correlators)
      std::string gmset_type = params_fopr.get_string("gamma_matrix_type");
      unique_ptr<GammaMatrixSet> gmset(GammaMatrixSet::New(gmset_type));

      // Wilson kernel operator, for the 4d <-> 5d source/sink construction.
      std::string kernel_type = params_fopr.get_string("kernel_type");
      double M0          = params_fopr.get_double("domain_wall_height");
      double coeff_c     = params_fopr.get_double("coefficient_c");
      double kappa           = 1.0 / (8.0 - 2.0 * M0);
      double one_over_2kappa = 4.0 - M0;
      Parameters params_kernel = params_fopr;
      params_kernel.set_double("hopping_parameter", kappa);
      unique_ptr<Fopr> foprw(Fopr::New(kernel_type, params_kernel));
      foprw->set_mode("D");
      foprw->set_config(U);

      // source generator
      std::string str_source_type = params_source.get_string("source_type");
      unique_ptr<Source> source(Source::New(str_source_type));
      source->set_parameters(params_source);

      // 4d propagators: boundary-projected sq and midpoint-projected sq_mid.
      std::vector<Field_F> sq(Nc * Nd), sq_mid(Nc * Nd);
      for (int i = 0; i < Nc * Nd; ++i) { sq[i].set(0.0); sq_mid[i].set(0.0); }
      const int s_mid_lo = Ns / 2 - 1;
      const int s_mid_hi = Ns / 2;

      Field_F b(Nvol, 1), vt1(Nvol, 1), vt2(Nvol, 1);
      Field_F b5(Nvol, Ns), x5(Nvol, Ns);

      int NFin  = afopr_fineD->field_nin();
      int NFvol = afopr_fineD->field_nvol();
      int NFex  = afopr_fineD->field_nex();
      AFIELD_d eta(NFin, NFvol, NFex), bprime(NFin, NFvol, NFex);
      AFIELD_d psi(NFin, NFvol, NFex), rchk(NFin, NFvol, NFex);
      AFIELD_d delta(NFin, NFvol, NFex), rphys(NFin, NFvol, NFex);

      const int    max_refine     = 8;
      const double refine_target2 = 1.0e-18;   // physical ||r||^2/||eta||^2

      vout.general(vl, "\n");
      vout.general(vl, "Hadron 2pt: 12 quark propagators through A (AMG + refinement):\n");
      vout.general(vl, "  color spin   Nconv_tot  nref   ||D psi-eta||/||eta||\n");

      for (int ispin = 0; ispin < Nd; ++ispin) {
        for (int icolor = 0; icolor < Nc; ++icolor) {
          int idx = icolor + Nc * ispin;

          b.set(0.0);
          source->set(b, idx);

          // 5d source: b5 = (cD-1)( P_+ 0 .... 0 P_- )^T b4
#pragma omp parallel
          {
            b5.set(0.0);

            foprw->mult_gm5(vt1, b);
            axpy(vt1, +1.0, b);                          // (1+gm5) b
            foprw->mult(vt2, vt1);
            aypx(-one_over_2kappa * coeff_c, vt2, vt1);  // (1 - cD)(1+gm5) b
            scal(vt2, -0.5);
            copy(b5, 0, vt2, 0);

            foprw->mult_gm5(vt1, b);
            axpy(vt1, -1.0, b);                          // (-1+gm5) b
            foprw->mult(vt2, vt1);
            aypx(-one_over_2kappa * coeff_c, vt2, vt1);  // (cD-1)(1-gm5) b
            scal(vt2, 0.5);
            copy(b5, Ns - 1, vt2, 0);
          }

          afopr_fineD->convert(eta, (Field&)b5);
          real_t eta2 = eta.norm2();

          psi.set(0.0);
          copy(rphys, eta);                 // r_phys = eta - D psi  (psi = 0)
          real_t res2      = eta2;
          int    nconv_tot = 0;
          int    nref      = 0;
          for (nref = 0; nref < max_refine; ++nref) {
            afopr_fineD->mult(bprime, rphys, "leftprec");  // b' = P_L r_phys

            int    nconv1 = -1;
            double diff1  = -1.0;
            delta.set(0.0);
#pragma omp parallel
            {
              asolver_mg->solve(delta, bprime, nconv1, diff1);  // A delta = b'
            }
            nconv_tot += nconv1;

            axpy(psi, real_t(1), delta);                  // psi += delta

            afopr_fineD->mult(rchk, psi, "physicalD");     // D psi
            copy(rphys, eta);
            axpy(rphys, real_t(-1), rchk);                 // r_phys = eta - D psi
            res2 = rphys.norm2();

            if (res2 < refine_target2 * eta2) { ++nref; break; }
          }

          vout.general(vl, "   %2d   %2d   %6d  %3d   %12.4e\n",
                       icolor, ispin, nconv_tot, nref,
                       (double)std::sqrt(res2 / eta2));

          // reverse psi -> x5, then project to the 4d propagator
          afopr_fineD->reverse((Field&)x5, psi);

#pragma omp parallel
          {
            copy(vt1, 0, x5, 0);
            copy(vt2, 0, x5, Ns - 1);
            copy(sq[idx], vt1);
#pragma omp barrier
            axpy(sq[idx], 1.0, vt2);
            axpy(vt1, -1.0, vt2);
#pragma omp barrier
            foprw->mult_gm5(vt2, vt1);
            axpy(sq[idx], -1.0, vt2);     // (1-gm5) x5[0] + (1+gm5) x5[Ns-1]
            scal(sq[idx], 0.5);
          }

#pragma omp parallel
          {
            copy(vt1, 0, x5, s_mid_hi);
            copy(vt2, 0, x5, s_mid_lo);
            copy(sq_mid[idx], vt1);
#pragma omp barrier
            axpy(sq_mid[idx], 1.0, vt2);
            axpy(vt1, -1.0, vt2);
#pragma omp barrier
            foprw->mult_gm5(vt2, vt1);
            axpy(sq_mid[idx], -1.0, vt2);
            scal(sq_mid[idx], 0.5);
          }
        }
      }

      // meson correlators
      Corr2pt_4spinor corr(gmset.get());
      corr.set_parameters(params_all.lookup("Corr2pt_4spinor"));

      vout.general(vl, "\n");
      vout.general(vl, "2-point correlator (propagator solved through A by AMG):\n");
      double result = corr.meson_all(sq, sq);
      vout.general(vl, "RESULT_AMG_2pt: %.17e\n", result);

      // residual mass: m_res(t) = C_PP(sq_mid)/C_PP(sq)
      GammaMatrix gm5 = gmset->get_GM(gmset->GAMMA5);
      std::vector<dcomplex> c_PP(Lt), c_J5q(Lt);
      corr.meson_correlator(c_PP,  gm5, gm5, sq,     sq);
      corr.meson_correlator(c_J5q, gm5, gm5, sq_mid, sq_mid);
      vout.general(vl, "\n");
      vout.general(vl, "Residual mass (J5q/PP ratio):\n");
      vout.general(vl, "    t      C_PP(t)         C_J5q(t)        m_res(t)\n");
      double mres_sum = 0.0; int mres_n = 0;
      for (int t = 0; t < Lt; ++t) {
        double cpp  = real(c_PP[t]);
        double cj5q = real(c_J5q[t]);
        double r    = (cpp != 0.0) ? cj5q / cpp : 0.0;
        vout.general(vl, " %4d  %15.7e  %15.7e  %15.7e\n", t, cpp, cj5q, r);
        if (t >= Lt / 4 && t < 3 * Lt / 4) { mres_sum += r; ++mres_n; }
      }
      double m_res = (mres_n > 0) ? mres_sum / mres_n : 0.0;
      vout.general(vl, "RESULT_AMG_mres (plateau avg over t in [%d,%d)): %.10e\n",
                   Lt / 4, 3 * Lt / 4, m_res);
    }

    return EXIT_SUCCESS;
  }

} // anonymous namespace

//====================================================================
int Test_MG_Solver::test()
{
  //  const std::string filename_input = "test_alt_Multigrid_32x64.yaml";
  const std::string filename_input = "test_alt_Multigrid.yaml";
  const std::string test_name      = class_name;

  vout.general("\n");
  vout.general("test name: %s\n", test_name.c_str());

  // ####  parameter setup  ####
  const int Nc   = CommonParameters::Nc();
  const int Nd   = CommonParameters::Nd();
  const int Ndim = CommonParameters::Ndim();
  const int Nvol = CommonParameters::Nvol();

  params_all = ParameterManager::read(filename_input);

  Parameters params_gauge = params_all.lookup("Gauge");

  const string str_gconf_status = params_gauge.get_string("gauge_config_status");
  const string str_gconf_read   = params_gauge.get_string("gauge_config_type_input");
  const string readfile         = params_gauge.get_string("config_filename_input");

  Parameters   params_test = params_all.lookup("TestType");
  const string str_vlevel  = params_test.get_string("verbose_level");
  const string solver_type = params_test.get_string("solver_type"); // MG

  Parameters   params_fopr = params_all.lookup("Fopr");
  const string fopr_type   = params_fopr.get_string("fermion_type");
  //  Parameters params_outer_solver = params_all.lookup("OuterSolver");
  Parameters params_MG_solver = params_all.lookup("MGSolver");

  vl = vout.set_verbose_level(str_vlevel);

  // setup random number manager
  RandomNumberManager::initialize("Mseries", 1234567UL);

  //- print input parameters
  vout.general(vl, "  gconf_status = %s\n", str_gconf_status.c_str());
  vout.general(vl, "  gconf_read   = %s\n", str_gconf_read.c_str());
  vout.general(vl, "  readfile     = %s\n", readfile.c_str());
  vout.general(vl, "  vlevel       = %s\n", str_vlevel.c_str());
  vout.general(vl, "  fuga:  solver_type  = %s\n", solver_type.c_str());
  vout.general(vl, "  Fopr         = %s\n", fopr_type.c_str());
  if (fopr_type != "Domainwall_5din") {
    vout.crucial("only for Domainwall_5din, sorry\n");
    abort();
  }
  const double quark_mass = params_fopr.get_double("quark_mass");
  const double PV_mass    = params_fopr.get_double("quark_mass_PauliVillars");
  const double DWheight   = params_fopr.get_double("domain_wall_height");
  const double coeff_b    = params_fopr.get_double("coefficient_b");
  const double coeff_c    = params_fopr.get_double("coefficient_c");
  const int N5 = params_fopr.get_int("extent_of_5th_dimension");
  vout.general(vl, "     quark mass  = %f\n", quark_mass);
  vout.general(vl, "     PV mass     = %f\n", PV_mass);
  vout.general(vl, "     DW height   = %f\n", DWheight);
  vout.general(vl, "     coeff_b     = %f\n", coeff_b);
  vout.general(vl, "     coeff_c     = %f\n", coeff_c);
  vout.general(vl, "     N5          = %d\n", N5);

  //- input parameter check
  int err = 0;
  err += ParameterCheck::non_NULL(str_gconf_status);

  if (err) {
    vout.crucial(vl, "Error at %s: input parameters have not been set\n", test_name.c_str());
    exit(EXIT_FAILURE);
  }

  unique_ptr<Timer> timer_total(new Timer(test_name));
  timer_total->start();

  //RandomNumberManager::initialize("Mseries", 1234567UL);
  RandomNumbers *random = RandomNumberManager::getInstance();

  // ####  Setup gauge configuration  ####
  vout.general("Field, creating: U\n");
  U.reset(new Field_G(Nvol, Ndim));

  if (str_gconf_status == "Continue") {
    GaugeConfig(str_gconf_read).read(*(U.get()), readfile);
  } else if (str_gconf_status == "Cold_start") {
    GaugeConfig("Unit").read(*(U.get()));
  } else if (str_gconf_status == "Hot_start") {
    GaugeConfig("Random").read(*(U.get()));
  } else {
    vout.crucial(m_vl, "Error at %s: unsupported gconf status \"%s\"\n",
                 test_name.c_str(), str_gconf_status.c_str());
    exit(EXIT_FAILURE);
  }


  // Source vector
  Field_F b(Nvol, N5);

  Field_F y(Nvol, N5);
  double  norm_b;
  //#pragma omp parallel
  {
    //b.set(1.0);
    random->uniform_lex_global(b);
    //    b.set(0.0);
    //    b.set(0,1.0);
    //if(Communicator::nodeid()==0){
    //  b.set(12,1.0);  // fot HOP_TP
    //  b.set(14+24*(1),1.0);  // fot HOP_TP
    //}
    norm_b = sqrt(b.norm2());
  }
  vout.general(vl, "|source|   =%.8f\n", norm_b);
  vout.general(vl, "|source|^2 =%.8f\n", norm_b * norm_b);

  //
  using complex_t = typename AFIELD_f::complex_t;

  // Source vector
  int nin = 2*Nc*Nd*N5;
  AFIELD_d ab;
  vout.general(vl, "resizing soruce\n");
  ab.reset(nin, Nvol, 1);

  // Solution vector
  AFIELD_d ax;
  vout.general(vl, "resizing solution\n");
  ax.reset(nin, Nvol, 1);

  // conversion
  {
    // fine grid operator: double prec.
    AFopr_Domainwall_5din<AFIELD_d> *afopr_fineD = new AFopr_Domainwall_5din<AFIELD_d>(params_fopr);

    vout.paranoiac(vl, "index is ready\n");
    vout.general(vl, "converting soruce\n");
#pragma omp parallel
    {
      if (afopr_fineD->needs_convert()) {
        vout.detailed(m_vl, "convert in rep. required.\n");
        afopr_fineD->convert(ab, b);
      } else {
        vout.detailed(m_vl, "convert in pep. not required.\n");
        AIndex_lex<double, AFIELD_d::IMPL> index_alt;
        convert(index_alt, ab, b);
      }
    }
  }
  vout.general(vl, "converting source, done\n");
  double n2 = ab.norm2();
  ab.scal(sqrt(1.0 / n2));


  if (solver_type == "BiCGStab") {
    run_bicgstab(ax, ab, U.get(), params_fopr, params_MG_solver);
  } else if (solver_type == "Mixed_eo") {
    run_mixed_eo(ax, ab, U.get(), params_fopr, params_MG_solver);
  }  else if (solver_type == "MG") {
    run_MG(ax, ab, U.get(), params_fopr, params_MG_solver, params_all);
  } else {
    vout.crucial("unkonwn solver type: solver_type=%s\n", solver_type.c_str());
    abort();
  }

  // convergence check
  {
    // fine grid operator: double prec.
    AFopr_Domainwall_5din<AFIELD_d> *afopr_fineD = new AFopr_Domainwall_5din<AFIELD_d>(params_fopr);
    //    afopr_fineD->set_parameters(params_fopr);
    afopr_fineD->set_config(U.get());
    afopr_fineD->set_mode("D");
    vout.general(vl, "fine grid operator (double) is ready [for convergence check]\n");

    // Solution vector
    AFIELD_d ay;
    vout.general(vl, "resizing solution\n");
    int nin  = afopr_fineD->field_nin();
    int nvol = afopr_fineD->field_nvol();
    int nex  = afopr_fineD->field_nex();
    ay.reset(nin, nvol, nex);
#pragma omp parallel
    {
      afopr_fineD->mult(ay, ax);
#pragma omp barrier
      axpy(ay, -1.0, ab);
    }
    double diff2 = ay.norm2();
    double norm2 = ab.norm2();
    vout.general("solving, done:      diff2=%22.15e\n", diff2);
    vout.general("           relatvie diff2=%22.15e\n", diff2 / norm2);
    delete afopr_fineD;
  }

  timer_total->report();

  RandomNumberManager::finalize();

  return EXIT_SUCCESS;
}


//============================================================END=====
