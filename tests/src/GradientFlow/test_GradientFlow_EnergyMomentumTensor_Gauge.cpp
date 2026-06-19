/*!
        @file    test_GradientFlow_EnergyMomentumTensor_Gauge.cpp

        @brief

        @author  Yusuke Taniguchi (tanigchi)

        @date    $LastChangedDate:: 2023-02-28 23:24:17 #$

        @version $LastChangedRevision: 2496 $
*/

#include "test.h"

#include "IO/gaugeConfig.h"

#include "Measurements/Gauge/gradientFlow.h"
#include "Measurements/Gauge/energyMomentumTensor.h"
//#include "Measurements/Gauge/topologicalCharge.h"

#include "Tools/randomNumberManager.h"

#include <iomanip>  // std::setw, std::setfill

//====================================================================
//! Test of gauge contribution to the energy-momentum tensor.

/*!
 Introduced by modifying Test_EnergyMomentumTensor_Gauge_Nf0.cpp
                               [26 May 2017 Y.Taniguchi]
 Register Test_EnergyMomentumTensor_Gauge to TestManager
                               [ 9 Sep 2017 Y.Namekawa]
*/

namespace Test_GradientFlow_EnergyMomentumTensor_Gauge {
  const std::string test_name             = "GradientFlow.EnergyMomentumTensor.Gauge";
  const std::string test_name_RK1         = "GradientFlow.EnergyMomentumTensor.Gauge.RK1";
  const std::string test_name_RK2         = "GradientFlow.EnergyMomentumTensor.Gauge.RK2";
  const std::string test_name_RK3         = "GradientFlow.EnergyMomentumTensor.Gauge.RK3";
  const std::string test_name_RK4         = "GradientFlow.EnergyMomentumTensor.Gauge.RK4";
  const std::string test_name_RK_adaptive = "GradientFlow.EnergyMomentumTensor.Gauge.RK_adaptive";

#if 0
  const std::string test_name             = "Test_EnergyMomentumTensor_Gauge";
  const std::string test_name_RK1         = "Test_EnergyMomentumTensor_Gauge.RK1";
  const std::string test_name_RK2         = "Test_EnergyMomentumTensor_Gauge.RK2";
  const std::string test_name_RK3         = "Test_EnergyMomentumTensor_Gauge.RK3";
  const std::string test_name_RK4         = "Test_EnergyMomentumTensor_Gauge.RK4";
  const std::string test_name_RK_adaptive = "Test_EnergyMomentumTensor_Gauge.RK_adaptive";
#endif

  //- test-private parameters
  namespace {
    const std::string filename_input_RK1         = "test_GradientFlow_EnergyMomentumTensor_Gauge_RK1.yaml";
    const std::string filename_input_RK2         = "test_GradientFlow_EnergyMomentumTensor_Gauge_RK2.yaml";
    const std::string filename_input_RK3         = "test_GradientFlow_EnergyMomentumTensor_Gauge_RK3.yaml";
    const std::string filename_input_RK4         = "test_GradientFlow_EnergyMomentumTensor_Gauge_RK4.yaml";
    const std::string filename_input_RK_adaptive = "test_GradientFlow_EnergyMomentumTensor_Gauge_RK_adaptive.yaml";
    // const std::string filename_input = filename_input_RK4;
  }

  //- prototype declaration
  int update(const std::string& filename_input);

  int run_test_RK1()
  { return update(filename_input_RK1); }

  int run_test_RK2()
  { return update(filename_input_RK2); }

  int run_test_RK3()
  { return update(filename_input_RK3); }

  int run_test_RK4()
  { return update(filename_input_RK4); }

  int run_test_RK_adaptive()
  { return update(filename_input_RK_adaptive); }

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
#if defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
    // static const bool is_registered_RK1 = TestManager::RegisterTest(
    //   test_name_RK1,
    //   run_test_RK1
    //   );
    // static const bool is_registered_RK2 = TestManager::RegisterTest(
    //   test_name_RK2,
    //   run_test_RK2
    //   );
    // static const bool is_registered_RK3 = TestManager::RegisterTest(
    //   test_name_RK3,
    //   run_test_RK3
    //   );
    static const bool is_registered_RK4 = TestManager::RegisterTest(
      test_name_RK4,
      run_test_RK4
      );
    // static const bool is_registered_RK_adaptive = TestManager::RegisterTest(
    //   test_name_RK_adaptive,
    //   run_test_RK_adaptive
    //   );
#endif
  }
#endif


  //====================================================================
  int update(const std::string& filename_input)
  {
    // ####  parameter setup  ####
    const int Nvol = CommonParameters::Nvol();
    const int Ndim = CommonParameters::Ndim();
    const int Nc   = CommonParameters::Nc();

    const Parameters params_all = ParameterManager::read(filename_input);

    const Parameters params_test = params_all.lookup("Test_EnergyMomentumTensor_Gauge");
    //- NB. additional parameter for action_G is set in the following object setup
    Parameters       params_action_G = params_all.lookup("Action_G");
    const Parameters params_g_flow   = params_all.lookup("GradientFlow");
    const Parameters params_energy_momentum_tensor = params_all.lookup("EnergyMomentumTensor");
    const Parameters params_topological_charge     = params_all.lookup("TopologicalCharge");

    const string        str_gconf_status = params_test.get_string("gauge_config_status");
    const string        str_gconf_read   = params_test.get_string("gauge_config_type_input");
    const string        readfile         = params_test.get_string("config_filename_input");
    const string        str_gconf_write  = params_test.get_string("gauge_config_type_output");
    const string        writefile        = params_test.get_string("config_filename_output");
    const string        str_rand_type    = params_test.get_string("random_number_type");
    const unsigned long seed             = params_test.get_unsigned_long("seed_for_random_number");
    int                 i_conf           = params_test.get_int("trajectory_number");
    const int           Nstep            = params_test.get_int("number_of_steps");
    const double        t_flow_max       = params_test.get_double("max_flow_time");
    const string        str_vlevel       = params_test.get_string("verbose_level");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

    const string str_action_G_type = params_action_G.get_string("action_type");

    const Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

    //- print input parameters
    vout.general(vl, "  gconf_status = %s\n", str_gconf_status.c_str());
    vout.general(vl, "  gconf_read   = %s\n", str_gconf_read.c_str());
    vout.general(vl, "  readfile     = %s\n", readfile.c_str());
    vout.general(vl, "  gconf_write  = %s\n", str_gconf_write.c_str());
    vout.general(vl, "  writefile    = %s\n", writefile.c_str());
    vout.general(vl, "  rand_type    = %s\n", str_rand_type.c_str());
    vout.general(vl, "  seed         = %lu\n", seed);
    vout.general(vl, "  i_conf       = %d\n", i_conf);
    vout.general(vl, "  Nstep        = %d\n", Nstep);
    vout.general(vl, "  t_flow_max   = %10.6f\n", t_flow_max);
    vout.general(vl, "  vlevel       = %s\n", str_vlevel.c_str());
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
    } else if (str_gconf_status == "Branch") {
      GaugeConfig(str_gconf_read).read(U, readfile);
    } else if (str_gconf_status == "Read") {
      std::ostringstream num;
      num << i_conf;
      string readfile_i_conf = readfile + num.str();
      GaugeConfig(str_gconf_read).read(U, readfile_i_conf);
      vout.general(vl, "  conf_readfile = %s\n", readfile_i_conf.c_str());
    } else if (str_gconf_status == "Read_gauge_heavy") {
      std::ostringstream num;
      num << "conf" << std::setw(4) << std::setfill('0') << i_conf << ".bin";
      string readfile_i_conf = readfile + num.str();
      GaugeConfig(str_gconf_read).read(U, readfile_i_conf);
      vout.general(vl, "  conf_readfile = %s\n", readfile_i_conf.c_str());
    } else if (str_gconf_status == "Read_gauge_phys") {
      std::ostringstream num;
      num << std::setw(6) << std::setfill('0') << i_conf;
      string readfile_i_conf = readfile + num.str();
      GaugeConfig(str_gconf_read).read(U, readfile_i_conf);
      vout.general(vl, "  conf_readfile = %s\n", readfile_i_conf.c_str());
    } else {
      vout.crucial(vl, "Error at %s: unsupported gconf status \"%s\"\n", test_name.c_str(), str_gconf_status.c_str());
      exit(EXIT_FAILURE);
    }

    //- beta is overwritten to be Nc in GradientFlow
    params_action_G.set_double("beta", static_cast<double>(Nc));

    unique_ptr<Action> action_G(Action::New(str_action_G_type, params_action_G));

    unique_ptr<GradientFlow> g_flow(new GradientFlow(action_G.get(), params_g_flow));

    EnergyMomentumTensor energy_momentum_tensor(params_energy_momentum_tensor);

    // TopologicalCharge topological_charge(params_topological_charge);

    Timer timer(test_name);


    // ####  Execution main part  ####
    timer.start();

    double t_flow = 0.0;
    double result = 0.0;

    for (int i = 0; i < Nstep; ++i) {
      double result_g_flow = g_flow->evolve(t_flow, U);

      energy_momentum_tensor.set_field_strength(U);

      result += energy_momentum_tensor.measure_EMT(t_flow);
      result += energy_momentum_tensor.measure_EMT_at_t(t_flow);
      result += energy_momentum_tensor.measure_EMT_at_t_FT(t_flow);
#if 0
      result += energy_momentum_tensor.measure_EMT_at_x(t_flow);
      result += energy_momentum_tensor.measure_EMT_at_x_FT(t_flow);
      result += energy_momentum_tensor.measure_EMT_at_y(t_flow);
      result += energy_momentum_tensor.measure_EMT_at_y_FT(t_flow);
      result += energy_momentum_tensor.measure_EMT_at_z(t_flow);
      result += energy_momentum_tensor.measure_EMT_at_z_FT(t_flow);
#endif

      // topological_charge.set_field_strength(U);
      // topological_charge.measure_topological_charge(t_flow);
      // topological_charge.measure_topological_density_t(t_flow);
      // topological_charge.measure_topological_density_t_FT(t_flow);
      // topological_charge.measure_topological_density_x(t_flow);
      // topological_charge.measure_topological_density_x_FT(t_flow);
      // topological_charge.measure_topological_density_y(t_flow);
      // topological_charge.measure_topological_density_y_FT(t_flow);
      // topological_charge.measure_topological_density_z(t_flow);
      // topological_charge.measure_topological_density_z_FT(t_flow);

      if (t_flow > t_flow_max) break;
    }
    vout.general(vl, " result = %.16e\n", result);

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
