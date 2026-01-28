/*!
        @file    test_GradientFlow_EnergyMomentumTensor_Fermion_2pt_Function.cpp

        @brief

        @author  Yusuke Taniguchi (tanigchi)

        @date    $LastChangedDate:: 2021-09-08 06:55:16 #$

        @version $LastChangedRevision: 2314 $
*/

#include "test.h"

#include "IO/gaugeConfig.h"

#include "Fopr/fopr_Clover_eo.h"
#include "Fopr/fopr_Clover.h"

#include "Measurements/Gauge/gradientFlow.h"
#include "Measurements/Fermion/fermionFlow.h"
#include "Measurements/Fermion/fermionFlow_2pt_Function.h"
#include "Measurements/Fermion/fprop_Standard_lex.h"
#include "Measurements/Fermion/fprop_Standard_eo.h"

#include "Solver/solver.h"

#include "Tools/filename.h"
#include "Tools/gammaMatrixSet.h"
#include "Tools/randomNumberManager.h"

#include <sstream>
#include <iostream>
#include <iomanip>  // std::setw, std::setfill

//====================================================================
//! Test of fermion contribution to the energy-momentum tensor.

/*!
 Introduced by modifying Test_EnergyMomentumTensor_Gauge.cpp
                               [1 June 2017 Y.Taniguchi]
*/

namespace Test_GradientFlow_EnergyMomentumTensor_Fermion_2pt_Function {
  // const std::string test_name = "Test_EnergyMomentumTensor_Fermion_2pt_Function";
  const std::string test_name = "GradientFlow.EnergyMomentumTensor.Fermion.2pt_Function";

  //- test-private parameters
  namespace {
    const std::string filename_input = "test_GradientFlow_EnergyMomentumTensor_Fermion_2pt_Function.yaml";
  }

  //- prototype declaration
  int update(int iconf);

  int run_update() { return update(0); }

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
#if defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
    static const bool is_registered = TestManager::RegisterTest(
      test_name,
      run_update
      );
#endif
  }
#endif


//====================================================================
  int update(int iconf)
  {
    // ####  parameter setup  ####
    const int Nvol = CommonParameters::Nvol();
    const int Ndim = CommonParameters::Ndim();
    const int Nc   = CommonParameters::Nc();

    const Parameters params_all = ParameterManager::read(filename_input);

    const Parameters params_test = params_all.lookup("Test_EnergyMomentumTensor_Fermion_2pt_Function");
    //- NB. additional parameter for action_G is set in the following object setup
    Parameters       params_action_G = params_all.lookup("Action_G");
    const Parameters params_gflow    = params_all.lookup("GradientFlow");
    const Parameters params_fflow    = params_all.lookup("FermionFlow");
    const Parameters params_f_2pt_fn = params_all.lookup("FermionFlow_2pt_Function");

    const string        str_gconf_status = params_test.get_string("gauge_config_status");
    const string        str_gconf_read   = params_test.get_string("gauge_config_type_input");
    const string        readfile         = params_test.get_string("config_filename_input");
    const string        str_gconf_write  = params_test.get_string("gauge_config_type_output");
    const string        writefile        = params_test.get_string("config_filename_output");
    const string        str_rand_type    = params_test.get_string("random_number_type");
    const unsigned long seed             = params_test.get_unsigned_long("seed_for_random_number");
    const int           N_quark          = params_test.get_int("number_of_valence_quarks");
    const int           even_odd         = params_test.get_int("even_odd");
    const string        str_vlevel       = params_test.get_string("verbose_level");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

    const string str_action_G_type = params_action_G.get_string("action_type");

    const Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

    vector<Parameters> params_quark;
    for (int iq = 0; iq < N_quark; ++iq) {
      std::string qlabel = Filename("Quark_{id}").format(iq + 1);
      params_quark.push_back(params_all.lookup(qlabel));
    }
    std::vector<Parameters> params_clover(N_quark);
    std::vector<Parameters> params_solver(N_quark);
    for (int iq = 0; iq < N_quark; ++iq) {
      params_clover[iq] = params_quark[iq].get_Parameters("Fopr");
      params_solver[iq] = params_quark[iq].get_Parameters("Solver");
    }

    //- print input parameters
    vout.general(vl, "  gconf_status  = %s\n", str_gconf_status.c_str());
    vout.general(vl, "  gconf_read    = %s\n", str_gconf_read.c_str());
    vout.general(vl, "  readfile      = %s\n", readfile.c_str());
    vout.general(vl, "  gconf_write   = %s\n", str_gconf_write.c_str());
    vout.general(vl, "  writefile     = %s\n", writefile.c_str());
    vout.general(vl, "  number_quarks = %d\n", N_quark);
    vout.general(vl, "  rand_type     = %s\n", str_rand_type.c_str());
    vout.general(vl, "  seed          = %lu\n", seed);
    vout.general(vl, "  even_odd      = %d\n", even_odd);
    vout.general(vl, "  vlevel        = %s\n", str_vlevel.c_str());
    for (int iq = 0; iq < N_quark; ++iq) {
      vout.general(vl, "Quark_%d\n", iq + 1);
      string str_gmset_type = params_clover[iq].get_string("gamma_matrix_type");
      vout.general(vl, "    gmset_type   = %s\n", str_gmset_type.c_str());
      double hopping_parameter = params_clover[iq].get_double("hopping_parameter");
      vout.general(vl, "    hopping_parameter= %lf\n", hopping_parameter);
      double clover_coefficient = params_clover[iq].get_double("clover_coefficient");
      vout.general(vl, "    clover_coefficient= %lf\n", clover_coefficient);
      string str_solver_type = params_solver[iq].get_string("solver_type");
      vout.general(vl, "    solver_type   = %s\n", str_solver_type.c_str());
    }
    vout.general(vl, "\n");

    //- input parameter check
    int err = 0;
    err += ParameterCheck::non_NULL(str_gconf_status);

    if (err) {
      vout.crucial(vl, "Error at %s: input parameters have not been set\n", test_name.c_str());
      exit(EXIT_FAILURE);
    }


    RandomNumberManager::initialize(str_rand_type, seed);


    // #### object setup #####
    Field_G U(Nvol, Ndim);

    if (str_gconf_status == "Continue") {
      GaugeConfig(str_gconf_read).read(U, readfile);
    } else if (str_gconf_status == "Cold_start") {
      GaugeConfig("Unit").read(U);
    } else if (str_gconf_status == "Hot_start") {
      GaugeConfig("Random").read(U);
    } else if (str_gconf_status == "Read") {
      std::ostringstream num;
      num << iconf;
      string readfile_iconf = readfile + num.str();
      GaugeConfig(str_gconf_read).read(U, readfile_iconf);
      vout.general(vl, "  conf_readfile = %s\n", readfile_iconf.c_str());
    } else if (str_gconf_status == "Read_gauge_heavy") {
      std::ostringstream num;
      num << "conf" << std::setw(4) << std::setfill('0') << iconf << ".bin";
      string readfile_iconf = readfile + num.str();
      GaugeConfig(str_gconf_read).read(U, readfile_iconf);
      vout.general(vl, "  conf_readfile = %s\n", readfile_iconf.c_str());
    } else if (str_gconf_status == "Read_gauge_phys") {
      std::ostringstream num;
      num << std::setw(6) << std::setfill('0') << iconf;
      string readfile_iconf = readfile + num.str();
      GaugeConfig(str_gconf_read).read(U, readfile_iconf);
      vout.general(vl, "  conf_readfile = %s\n", readfile_iconf.c_str());
    } else {
      vout.crucial(vl, "Error at %s: unsupported gconf status \"%s\"\n", test_name.c_str(), str_gconf_status.c_str());
      exit(EXIT_FAILURE);
    }

    //- beta is overwritten to be Nc in GradientFlow
    params_action_G.set_double("beta", static_cast<double>(Nc));

    unique_ptr<Action> action_G(Action::New(str_action_G_type, params_action_G));

    struct QuarkType {
      unique_ptr<Fopr>   fopr;
      unique_ptr<Solver> solver;
      unique_ptr<Fprop>  fprop;
    };
    std::vector<QuarkType> quarks(N_quark);

    for (int iq = 0; iq < N_quark; ++iq) {
      const Parameters& params_fopr = params_clover[iq];
      const Parameters& params_qsolver = params_solver[iq];

      if (even_odd == 1) {
        quarks[iq].fopr.reset(new Fopr_Clover_eo(params_fopr));
        quarks[iq].fopr->set_config(&U);

        quarks[iq].solver.reset(
          Solver::New(params_qsolver.get_string("solver_type"),
                      quarks[iq].fopr.get(),
                      params_qsolver));

        quarks[iq].fprop.reset(new Fprop_Standard_eo(quarks[iq].solver.get()));
        
      } else {
        quarks[iq].fopr.reset(new Fopr_Clover(params_fopr));
        quarks[iq].fopr->set_config(&U);

        quarks[iq].solver.reset(
          Solver::New(params_qsolver.get_string("solver_type"),
                      quarks[iq].fopr.get(),
                      params_qsolver));

        quarks[iq].fprop.reset(new Fprop_Standard_lex(quarks[iq].solver.get()));
      }
    }
    vout.general(vl, "\n");

    std::vector<Fprop*> fprop(N_quark);
    for (int iq = 0; iq < N_quark; ++iq) {
      fprop[iq] = quarks[iq].fprop.get();
    }

    FermionFlow_2pt_Function measurement(fprop, action_G.get(), &params_clover);
    measurement.set_parameters(params_f_2pt_fn, params_gflow, params_fflow);

    Timer timer(test_name);

    // ####  Execution main part  ####
    timer.start();

    //measurement.measure_meson_correlator(U);
    double result = measurement.measure_EMT_correlator(U);

    timer.report();

    RandomNumberManager::finalize();


    if (do_check) {
      return Test::verify(result, expected_result);
    } else {
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }
  }
} // namespace Test_EnergyMomentumTensor_Gauge
