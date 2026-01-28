/*!
        @file    test_HMC_Overlap_Nf2.cpp

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate: 2013-01-22 13:51:53 #$

        @version $LastChangedRevision: 2422 $
*/

#include "test.h"

#include "Action/Fermion/action_F_Overlap_Nf2.h"
#include "Action/Fermion/action_F_Standard_lex.h"

#include "Force/Fermion/force_F_Overlap_Nf2.h"

#include "HMC/hmc_Leapfrog.h"

#include "IO/gaugeConfig.h"

#include "Measurements/Fermion/fprop_Standard_lex.h"

#include "Tools/file_utils.h"
#include "Tools/randomNumberManager.h"
#include "Tools/randomNumbers_Mseries.h"

//====================================================================
//! Test of HMC update for overlap fermions.

/*!
    This class tests HMC update for dynamical overlap fermions.
    Smearing of gauge configuration for the fermion operator
    is incorporated.
                                         [12 Apr 2012 H.Matsufuru]
    (Coding history will be recovered from trac.)
    Implement YAML.                      [14 Nov 2012 Y.Namekawa]
    Implement Fprop and selectors.       [03 Mar 2013 Y.Namekawa]
    (Selectors are replaced with factories by Aoyama-san)
    Introduce unique_ptr to avoid memory leaks.
                                         [21 Mar 2015 Y.Namekawa]
 */

namespace Test_HMC_Overlap {
  const std::string test_name = "HMC.Overlap.Leapfrog_Nf2";

  //- test-private parameters
  namespace {
    const std::string filename_input = "test_HMC_Overlap_Nf2.yaml";
  }

  //- prototype declaration
  int leapfrog_Nf2(void);

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
#if defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
    static const bool is_registered = TestManager::RegisterTest(
      test_name,
      leapfrog_Nf2
      );
#endif
  }
#endif

  //====================================================================
  int leapfrog_Nf2(void)
  {
// multithread not supported

    // if (ThreadManager::get_num_threads_available() > 1) {
    //   vout.general("%s: multithread unsupported.\n\n", test_name.c_str());
    //   return EXIT_SKIP;
    // }


    // #####  parameter setup  #####
    const int Nc   = CommonParameters::Nc();
    const int Nvol = CommonParameters::Nvol();
    const int Ndim = CommonParameters::Ndim();
    const int NinG = 2 * Nc * Nc;

    const Parameters params_all = ParameterManager::read(filename_input);

    const Parameters params_test      = params_all.lookup("Test_HMC_Overlap");
    const Parameters params_action_G  = params_all.lookup("Action_G");
    const Parameters params_fopr      = params_all.lookup("Fopr");
    const Parameters params_overlap   = params_all.lookup("Fopr_Overlap");
    const Parameters params_solver_MD = params_all.lookup("Solver_MD");
    const Parameters params_solver_H  = params_all.lookup("Solver_H");
    const Parameters params_hmc       = params_all.lookup("HMC_Leapfrog");

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

    const string str_action_G_type  = params_action_G.get_string("action_type");
    const string str_fopr_type      = params_fopr.get_string("fermion_type");
    const string str_gmset_type     = params_fopr.get_string("gamma_matrix_type");
    const string str_solver_MD_type = params_solver_MD.get_string("solver_type");
    const string str_solver_H_type  = params_solver_H.get_string("solver_type");

    const Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

    //- print input parameters
    vout.general(vl, "  gconf_status = %s\n", str_gconf_status.c_str());
    vout.general(vl, "  gconf_read     = %s\n", str_gconf_read.c_str());
    vout.general(vl, "  readfile       = %s\n", readfile.c_str());
    vout.general(vl, "  gconf_write    = %s\n", str_gconf_write.c_str());
    vout.general(vl, "  writefile      = %s\n", writefile.c_str());
    vout.general(vl, "  rand_type      = %s\n", str_rand_type.c_str());
    vout.general(vl, "  seed           = %lu\n", seed);
    vout.general(vl, "  i_conf         = %d\n", i_conf);
    vout.general(vl, "  Ntraj          = %d\n", Ntraj);
    vout.general(vl, "  i_save_conf    = %d\n", i_save_conf);
    vout.general(vl, "  vlevel         = %s\n", str_vlevel.c_str());
    vout.general(vl, "  gmset_type     = %s\n", str_gmset_type.c_str());
    vout.general(vl, "  solver_MD_type = %s\n", str_solver_MD_type.c_str());
    vout.general(vl, "  solver_H_type  = %s\n", str_solver_H_type.c_str());
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

    //-- N_f=2 part
    unique_ptr<Fopr> fopr(Fopr::New(str_fopr_type, params_fopr));

    unique_ptr<Fopr> fopr_overlap(Fopr::New("Overlap", fopr.get(), params_overlap));

    unique_ptr<Force> force_fopr_overlap(new Force_F_Overlap_Nf2(params_overlap));

    unique_ptr<Solver> solver_MD(Solver::New(str_solver_MD_type, fopr_overlap.get(), params_solver_MD));
    unique_ptr<Fprop> fprop_MD(new Fprop_Standard_lex(solver_MD.get()));

    unique_ptr<Solver> solver_H(Solver::New(str_solver_H_type, fopr_overlap.get(), params_solver_H));
    unique_ptr<Fprop> fprop_H(new Fprop_Standard_lex(solver_H.get()));

    unique_ptr<Action> action_F(new Action_F_Overlap_Nf2(fopr_overlap.get(), force_fopr_overlap.get(), fprop_MD.get(), fprop_H.get()));

    ActionList actions(1);
    actions.append(0, action_F.get());
    actions.append(0, action_G.get());

    //- Random number is initialized with a parameter specified by i_conf
    unique_ptr<RandomNumbers> rand(new RandomNumbers_Mseries(i_conf));

    // define hmc_leapfrog (SA)
    HMC_Leapfrog hmc(actions, rand.get(), params_hmc);

    Timer timer(test_name);


    // ####  Execution main part  ####
    timer.start();

    vout.general(vl, "HMC: Ntraj = %d\n", Ntraj); // a number of trajectory (SA)

    double result = 0.0;
    for (int traj = 0; traj < Ntraj; ++traj) {
      vout.general(vl, "\n");
      vout.general(vl, "traj = %d\n", traj);

      result = hmc.update(U); // hmc update (SA)

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
} // namespace Test_HMC_Overlap
