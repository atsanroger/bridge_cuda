/*!
        @file    test_HMC_Wilson_Nf2.cpp

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate: 2013-01-22 13:51:53 #$

        @version $LastChangedRevision: 2520 $
*/

#include "test.h"

#include "lib/Action/Fermion/action_F_Standard_lex.h"

#include "lib/Fopr/fopr_Smeared.h"

#include "lib/Force/Fermion/force_F_Wilson_Nf2.h"
//#include "Force/Fermion/force_F_Smeared.h"
#include "lib_alt/Force/Fermion/aforce_F_Smeared_alt.h"

#include "lib/HMC/hmc_General.h"
#include "lib/HMC/builder_Integrator.h"

#include "lib/IO/gaugeConfig.h"

#include "lib/Measurements/Fermion/fprop_Standard_lex.h"

#include "lib/Tools/file_utils.h"
#include "lib/Tools/randomNumberManager.h"
#include "lib/Tools/randomNumbers_Mseries.h"

#ifdef USE_ALT_CODE
#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_lex.h"
#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_eo.h"
#include "lib_alt/Smear/director_alt_Smear.h"
#include "lib_alt/Force/Fermion/aforce_F_Wilson_Nf2.h"
#include "lib_alt/Force/Fermion/aforce_F_Wilson_eo_Nf2.h"
#include "lib_alt/Action/Fermion/action_F_alt_Standard_lex.h"
#include "lib_alt/Action/Fermion/action_F_alt_Standard_eo.h"
#include "lib_alt/Action/Fermion/action_F_alt_Ratio_eo.h"
#endif

#ifdef USE_ALT_QXS
#include "lib_alt_QXS/bridge_alt_qxs.h"
#endif


//====================================================================
//! Test of HMC update for Wilson fermions.

/*!
    This class tests HMC update for dynamical Wilson fermions.
    Smearing of gauge configuration for the fermion operator
    is incorporated.
                                          [12 Apr 2012 H.Matsufuru]
    (Coding history will be recovered from trac.)
    Implement YAML.                       [14 Nov 2012 Y.Namekawa]
    Implement Fprop and selectors.        [03 Mar 2013 Y.Namekawa]
    (Selectors are replaced with factories by Aoyama-san)
    Introduce unique_ptr to avoid memory leaks.
                                          [21 Mar 2015 Y.Namekawa]
    Add Nc check for USE_GROUP_SU_N.      [31 May 2021 Y.Namekawa]
 */

namespace Test_HMC_Wilson {
  const std::string test_name = "HMC.Wilson_eo.Nf2";

  //- test-private parameters
  namespace {
    const std::string filename_input = "test_HMC_Wilson_Nf2_alt_eo.yaml";
  }

  //- prototype declaration
  int update_Nf2_alt_eo(void);

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
#if defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
    static const bool is_registered = TestManager::RegisterTest(
      test_name,
      update_Nf2
      );
#endif
  }
#endif

  //====================================================================
  int update_Nf2_alt_eo(void)
  {
    //typedef  Field  AFIELD;
    typedef  AField<double,QXS>  AFIELD;


    // #####  parameter setup  #####
    const int Nc   = CommonParameters::Nc();
    const int Nvol = CommonParameters::Nvol();
    const int Ndim = CommonParameters::Ndim();

    const Parameters params_all = ParameterManager::read(filename_input);

    const Parameters params_test       = params_all.lookup("Test_HMC_Wilson");
    const Parameters params_action_G   = params_all.lookup("Action_G");
    const Parameters params_smear      = params_all.lookup("LinkSmearing");
    const Parameters params_fopr1      = params_all.lookup("Fopr1");
    // const Parameters params_fopr2      = params_all.lookup("Fopr2");
    const Parameters params_solver_MD  = params_all.lookup("Solver_MD");
    const Parameters params_solver_H   = params_all.lookup("Solver_H");
    const Parameters params_integrator = params_all.lookup("Builder_Integrator");
    const Parameters params_hmc        = params_all.lookup("HMC_General");

    const string     str_gconf_status = params_test.get_string("gauge_config_status");
    const string     str_gconf_read   = params_test.get_string("gauge_config_type_input");
    const string     readfile         = params_test.get_string("config_filename_input");
    const string     str_gconf_write  = params_test.get_string("gauge_config_type_output");
    const string     writefile        = params_test.get_string("config_filename_output");
    const string     str_rand_type    = params_test.get_string("random_number_type");
    const unsigned long seed             = params_test.get_unsigned_long("seed_for_random_number");
    int              i_conf           = params_test.get_int("trajectory_number");
    const int        Ntraj            = params_test.get_int("trajectory_number_step");
    const int        i_save_conf      = params_test.get_int("save_config_interval");
    const string     str_vlevel       = params_test.get_string("verbose_level");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

    const string str_action_G_type = params_action_G.get_string("action_type");
    const string str_gmset_type    = params_fopr1.get_string("gamma_matrix_type");
    const string           str_solver_MD_type = params_solver_MD.get_string("solver_type");
    const string           str_solver_H_type  = params_solver_H.get_string("solver_type");
    const int              Nlevels            = params_integrator.get_int("number_of_levels");
    const std::vector<int> level_action       = params_integrator.get_int_vector("level_of_actions");

    const Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

    //- print input parameters
    vout.general(vl, "  gconf_status   = %s\n", str_gconf_status.c_str());
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

    //- smearing director
    typedef Director_alt_Smear<AFIELD>  Director_alt;
    unique_ptr<Director_alt> dr_smear(new Director_alt(params_smear));

    //- N_f=2 preconditioner part

    string fopr1_type = params_fopr1.get_string("fermion_type");
    if(fopr1_type.substr(fopr1_type.size()-3 ,3) != "_eo"){
      vout.crucial(vl, "fop1_type = %s\n", fopr1_type.c_str());
      fopr1_type += "_eo";
      vout.crucial(vl, "fop1_type = %s\n", fopr1_type.c_str());
    }
    unique_ptr<AFopr<AFIELD> > fopr1(AFopr<AFIELD>::New(
                                           fopr1_type, params_fopr1));

    unique_ptr<AForce_F<AFIELD> > forceF1(
	   new AForce_F_Wilson_eo_Nf2<AFIELD>(params_fopr1));

    unique_ptr<AFopr<AFIELD> > fopr1_smr(AFopr<AFIELD>::New(
                      "Smeared_alt", fopr1.get(), dr_smear.get()));

    unique_ptr<AForce_F<AFIELD> > forceF1_smr(
      new AForce_F_Smeared_alt<AFIELD>(forceF1.get(), dr_smear.get()));

    unique_ptr<Fprop_alt<AFIELD> > fprop1_MD(
      new Fprop_alt_Standard_eo<AFIELD>(params_fopr1, params_solver_MD,
                                      dr_smear.get()));

    unique_ptr<Fprop_alt<AFIELD> > fprop1_H(
      new Fprop_alt_Standard_eo<AFIELD>(params_fopr1, params_solver_H,
                                        dr_smear.get()));

    unique_ptr<Action> action_F1(
      new Action_F_alt_Standard_eo<AFIELD>(
                            fopr1_smr.get(), forceF1_smr.get(),
                            fprop1_MD.get(), fprop1_H.get()));

    //- N_f=2 dynamical part
    /*
    string fopr2_type = params_fopr2.get_string("fermion_type");
    if(fopr2_type.substr(fopr2_type.size()-3 ,3) != "_eo"){
      vout.crucial(vl, "fop2_type = %s\n", fopr2_type.c_str());
      fopr2_type += "_eo";
      vout.crucial(vl, "fop2_type = %s\n", fopr2_type.c_str());
    }
    unique_ptr<AFopr<AFIELD> > fopr2(AFopr<AFIELD>::New(
                                           fopr2_type, params_fopr2));

    unique_ptr<AForce_F<AFIELD> > forceF2(
      new AForce_F_Wilson_eo_Nf2<AFIELD>(params_fopr2));

    unique_ptr<AFopr<AFIELD> > fopr2_smr(AFopr<AFIELD>::New(
                      "Smeared_alt", fopr2.get(), dr_smear.get()));

    unique_ptr<AForce_F<AFIELD> > forceF2_smr(
      new AForce_F_Smeared_alt<AFIELD>(forceF2.get(), dr_smear.get()));

    unique_ptr<Fprop_alt<AFIELD> > fprop2_MD(
      new Fprop_alt_Standard_eo<AFIELD>(params_fopr2, params_solver_MD,
                                      dr_smear.get()));

    unique_ptr<Fprop_alt<AFIELD> > fprop2_H(
      new Fprop_alt_Standard_eo<AFIELD>(params_fopr2, params_solver_H,
                                        dr_smear.get()));

    unique_ptr<Action> action_F2(
      new Action_F_alt_Ratio_eo<AFIELD>(
                            fopr1_smr.get(), forceF1_smr.get(),
                            fopr2_smr.get(), forceF2_smr.get(),
                            fprop1_H.get(),
                            fprop2_MD.get(), fprop2_H.get()));
    */

    //- set actions
    ActionList actions(Nlevels);
    actions.append(level_action[0], action_F1.get());
    actions.append(level_action[1], action_G.get());

    std::vector<Director *> directors(1);
    directors[0] = static_cast<Director *>(dr_smear.get());

    unique_ptr<Builder_Integrator> builder(
      new Builder_Integrator(actions, directors, params_integrator));
    Integrator *integrator = builder->build();

    //- Random number is initialized with a parameter specified by i_conf
    unique_ptr<RandomNumbers> rand(new RandomNumbers_Mseries(i_conf));


    HMC_General hmc(actions, directors, integrator, rand.get(), params_hmc);

    Timer timer(test_name);

    // ####  Execution main part  ####
    timer.start();

    vout.general(vl, "HMC: Ntraj = %d\n", Ntraj); // a number of trajectory (SA)

    double result = 0.0;
    for (int traj = 0; traj < Ntraj; ++traj) {
      vout.general(vl, "\n");
      vout.general(vl, "traj = %d\n", traj);

      result = hmc.update(U);

      if ((i_conf + traj + 1) % i_save_conf == 0) {
        std::string filename = FileUtils::generate_filename("%s-%06d",
                                writefile.c_str(), (i_conf + traj + 1));
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
} // namespace Test_HMC_Wilson
