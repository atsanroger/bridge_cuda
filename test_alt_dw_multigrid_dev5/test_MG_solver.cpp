/*!

        @file    test_MG_solver.cpp

        @brief   implementation of the test for multigrid solver (domainwall)

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: #$

        @version $LastChangedRevision: 2602 $
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
             const Parameters &params_solver)
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

    // fine grid operator: double prec.
    unique_ptr<AFopr_Domainwall_5din<AFIELD_d> > afopr_fineD(
                       new AFopr_Domainwall_5din<AFIELD_d> (params_fopr));
    //    afopr_fineD->set_parameters(params_fopr);
    afopr_fineD->set_config(U);
    vout.general(vl, "fine grid operator (double) is ready\n");

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

    // multigrid solver
    unique_ptr<ASolver_MG_dw<AFIELD_d> > asolver_mg(
                                            new ASolver_MG_dw<AFIELD_d> );

    asolver_mg->set_parameters(params_solver);
    // set fopr and parameters
    asolver_mg->set_foprD(afopr_fineD.get());
    asolver_mg->set_foprF(afopr_fineF.get());
    asolver_mg->set_fopr_PV(afopr_PV.get());

    unique_ptr<AFopr_dd<AFIELD_f> > afopr_smoother;
    if (asolver_mg->use_fopr_for_smoother()) {
      afopr_smoother.reset(asolver_mg->new_fopr_smoother(params_fopr_dd));
      afopr_smoother->set_config(U);
      asolver_mg->set_fopr_smoother(afopr_smoother.get());
    }

    asolver_mg->init_solver();

    // setup null space vectors
    asolver_mg->run_setup();

    int               nconv = 0;
    double            diff  = -1.0;
    unique_ptr<Timer> timer_solve(new Timer("solve [MG, solve only]"));
    timer_solve->start();
#pragma omp parallel
    {
      asolver_mg->solve(ax, ab, nconv, diff);
    }
    timer_solve->stop();
    timer_solve->report();
    double etime = timer_solve->elapsed_msec();

    double flop_solve = asolver_mg->flop_count();
    vout.general(" Flops (double+float) : %f GFlop/sec\n",
                 1.0e-6 * flop_solve / etime);

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
    run_MG(ax, ab, U.get(), params_fopr, params_MG_solver);
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
