/*!
        @file    test_HMC_Staggered_Nf4_eo.cpp

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-06-29 22:43:06 #$

        @version $LastChangedRevision: 2527 $
*/

#include "test.h"

#include "Action/Fermion/action_F_Staggered_eo.h"
#include "Action/Fermion/action_F_Rational.h"

#include "Force/Fermion/force_F_Staggered_eo.h"

#include "HMC/hmc_General.h"
#include "HMC/builder_Integrator.h"

#include "IO/gaugeConfig.h"

#include "Tools/file_utils.h"
#include "Tools/randomNumberManager.h"
#include "Tools/randomNumbers_Mseries.h"

//====================================================================
//! Test of staggered fermion.

/*!
    This class tests HMC update for dynamical staggered fermions.
                                          [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    Implement YAML.                       [14 Nov 2012 Y.Namekawa]
    Introduce unique_ptr to avoid memory leaks.
                                          [21 Mar 2015 Y.Namekawa]
    Add Nc check for USE_GROUP_SU_N.      [31 May 2021 Y.Namekawa]
*/

namespace Test_HMC_Staggered {
  const std::string test_name = "HMC.Staggered.Nf4_eo";

  //- test-private parameters
  namespace {
    const std::string filename_input = "test_HMC_Staggered_Nf4.yaml";
  }

  //- prototype declaration
  int update_Nf4_eo(void);

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
#if defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
    static const bool is_registered = TestManager::RegisterTest(
      test_name,
      update_Nf4_eo
      );
#endif
  }
#endif

  //====================================================================
  int update_Nf4_eo(void)
  {
    // #####  parameter setup  #####
    const int Nc   = CommonParameters::Nc();
    const int Nvol = CommonParameters::Nvol();
    const int Ndim = CommonParameters::Ndim();

    const Parameters params_all = ParameterManager::read(filename_input);

    const Parameters params_test       = params_all.lookup("Test_HMC_Staggered");
    const Parameters params_action_G   = params_all.lookup("Action_G");
    const Parameters params_action_F   = params_all.lookup("Action_F_Staggered");
    const Parameters params_staggered  = params_all.lookup("Fopr_Staggered");
    const Parameters params_integrator = params_all.lookup("Builder_Integrator");
    const Parameters params_hmc        = params_all.lookup("HMC_General");

    const string        str_gconf_status = params_test.get_string("gauge_config_status");
    const string        str_gconf_read   = params_test.get_string("gauge_config_type_input");
    const string        readfile         = params_test.get_string("config_filename_input");
    const string        str_gconf_write  = params_test.get_string("gauge_config_type_output");
    const string        writefile        = params_test.get_string("config_filename_output");
    const string        str_rand_type    = params_test.get_string("random_number_type");
    const unsigned long seed             = params_test.get_unsigned_long("seed_for_random_number");
    int                 i_conf           = params_test.get_int("trajectory_number");
    const int           Ntraj            = params_test.get_int("trajectory_number_step");
    const int           i_save_conf      = params_test.get_int("save_config_interval");
    const string        str_vlevel       = params_test.get_string("verbose_level");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

    const string           str_action_G_type = params_action_G.get_string("action_type");
    const int              Nlevels           = params_integrator.get_int("number_of_levels");
    const std::vector<int> level_action      = params_integrator.get_int_vector("level_of_actions");

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
    vout.general(vl, "  Ntraj        = %d\n", Ntraj);
    vout.general(vl, "  i_save_conf  = %d\n", i_save_conf);
    vout.general(vl, "  vlevel       = %s\n", str_vlevel.c_str());
    vout.general(vl, "\n");

    //- input parameter check
    int err = 0;
    err += ParameterCheck::non_NULL(str_gconf_status);
    err += ParameterCheck::non_negative(i_conf);
    err += ParameterCheck::non_negative(Ntraj);
    err += ParameterCheck::non_negative(i_save_conf);

    if (err) {
      vout.crucial(vl, "Error at %s: input parameters have not been set\n", test_name.c_str());
      exit(EXIT_FAILURE);
    }

    //- NB. Stout_SU3 is implicitly used in Action_F_Staggered_eo()
    if (CommonParameters::Nc() != 3) {
      vout.crucial(vl, "check skipped: Nc = 3 is needed, but Nc = %d.\n\n", CommonParameters::Nc());
      return EXIT_SKIP;
    }


    RandomNumberManager::initialize(str_rand_type, seed);


    // #####  object setup  #####
    Field_G U(Nvol, Ndim);

    if (str_gconf_status == "Continue") {
      GaugeConfig(str_gconf_read).read(U, readfile);
    } else if (str_gconf_status == "Cold_start") {
      GaugeConfig("Unit").read(U);
    } else if (str_gconf_status == "Hot_start") {
      GaugeConfig("Random").read(U);
    } else {
      vout.crucial(vl, "Error at %s: unsupported gconf status \"%s\"\n", test_name.c_str(), str_gconf_status.c_str());
      exit(EXIT_FAILURE);
    }

    GaugeConfig gconf_write(str_gconf_write);


    unique_ptr<Action> action_G(Action::New(str_action_G_type, params_action_G));

    unique_ptr<Fopr>  fopr_ks(new Fopr_Staggered_eo(params_staggered));
    unique_ptr<Force> force_fopr_ks(new Force_F_Staggered_eo(params_staggered));

    unique_ptr<Action> action_F2(new Action_F_Staggered_eo(fopr_ks.get(), force_fopr_ks.get(), params_action_F));

    ActionList actions(Nlevels);
    actions.append(level_action[0], action_F2.get());
    actions.append(level_action[1], action_G.get());

    unique_ptr<Builder_Integrator> builder(new Builder_Integrator(actions, params_integrator));
    Integrator *integrator = builder->build();

    unique_ptr<RandomNumbers> rand(new RandomNumbers_Mseries(i_conf));

    HMC_General hmc(actions, integrator, rand.get(), params_hmc);

    Timer timer(test_name);


    // ####  Execution main part  ####
    timer.start();

    vout.general(vl, "HMC: Ntraj = %d\n", Ntraj);

    double result = 0.0;
    for (int traj = 0; traj < Ntraj; ++traj) {
      vout.general(vl, "\n");
      vout.general(vl, "traj = %d\n", traj);

      result = hmc.update(U);

      if ((i_conf + traj + 1) % i_save_conf == 0) {
        std::string filename = FileUtils::generate_filename("%s-%06d", writefile.c_str(), (i_conf + traj + 1));
        gconf_write.write_file(U, filename);
      }
    }

    gconf_write.write_file(U, writefile);

    timer.report();

    RandomNumberManager::finalize();


    if (do_check) {
      return Test::verify(result, expected_result);
    } else {
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }
  }
} // namespace Test_HMC_Staggered
