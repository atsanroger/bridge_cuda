/*!
        @file    test_HMC_Domainwall_Nf2_PV.cpp

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate: 2013-04-27 12:28:50 #$

        @version $LastChangedRevision: 2520 $
*/

#include "test.h"

#include "lib/Action/Fermion/action_F_Standard_lex.h"
#include "lib/Action/Fermion/action_F_Ratio_lex.h"

#include "lib/Fopr/fopr_Smeared.h"
#include "lib/Fopr/fopr_Smeared_alt.h"

#include "lib/Force/Fermion/force_F_Domainwall.h"
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
#include "lib_alt/Force/Fermion/aforce_F_Domainwall.h"
#include "lib_alt/Force/Fermion/aforce_F_Domainwall_eo.h"
#include "lib_alt/Action/Fermion/action_F_alt_Ratio_lex.h"
#include "lib_alt/Action/Fermion/action_F_alt_Ratio_eo.h"
#endif

#ifdef USE_ALT_QXS
#include "lib_alt_QXS/bridge_alt_qxs.h"
#endif



//====================================================================
//! Test of HMC update for domain-wall fermion operator.

/*!
    This class tests HMC update for dynamical domain-wall fermions.
    Alt-code version.
                                          [12 Apr 2023 H.Matsufuru]
 */

namespace Test_HMC_Domainwall {
  const std::string test_name = "HMC.Domainwall_eo.Nf2_PV_nosmr";

  //- test-private parameters
  namespace {
    const std::string filename_input
                     = "test_HMC_Domainwall_Nf2_PV_alt_eo_nosmr.yaml";
  }

  //- prototype declaration
  int update_Nf2_PV_alt_eo_nosmr(void);

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
#if defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
    static const bool is_registered = TestManager::RegisterTest(
      test_name,
      update_Nf2_PV_alt_nosmr
      );
#endif
  }
#endif

  //====================================================================
  int update_Nf2_PV_alt_eo_nosmr(void)
  {
    vout.crucial("test update_Nf2_PV_alt_eo_nosmr start.\n");

    //typedef  Field  AFIELD;
    typedef  AField<double,QXS>  AFIELD;


    // #####  parameter setup  #####
    const int Nc   = CommonParameters::Nc();
    const int Nvol = CommonParameters::Nvol();
    const int Ndim = CommonParameters::Ndim();

    const Parameters params_all = ParameterManager::read(filename_input);

    const Parameters params_test       = params_all.lookup("Test_HMC_Domainwall");
    const Parameters params_action_G   = params_all.lookup("Action_G");
    const Parameters params_smear      = params_all.lookup("LinkSmearing");

    const Parameters params_dw         = params_all.lookup("Fopr_Domainwall");
    const Parameters params_dw_pv      = params_all.lookup("Fopr_Domainwall_PauliVillars");
    const Parameters params_solver_MD  = params_all.lookup("Solver_MD");
    const Parameters params_solver_H   = params_all.lookup("Solver_H");
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

    const string str_action_G_type = params_action_G.get_string("action_type");
    const string str_gmset_type    = params_dw.get_string("gamma_matrix_type");
    const string           str_solver_MD_type = params_solver_MD.get_string("solver_type");
    const string           str_solver_H_type  = params_solver_H.get_string("solver_type");
    const int              Nlevels            = params_integrator.get_int("number_of_levels");
    const std::vector<int> level_action       = params_integrator.get_int_vector("level_of_actions");

    const Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

    vout.crucial("parameter setup finished.\n");
    Communicator::sync();

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

    vout.crucial("parameter output finished.\n");
    Communicator::sync();

    //- input parameter check
    int err = 0;
    err += ParameterCheck::non_NULL(str_gconf_status);
    err += ParameterCheck::non_negative(i_conf);
    err += ParameterCheck::non_negative(Ntraj);
    err += ParameterCheck::non_negative(i_save_conf);

    if (err) {
      vout.crucial(vl, "Error at %s: input parameters have not been set\n",
                   test_name.c_str());
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
      vout.crucial(vl, "Error at %s: unsupported gconf status \"%s\"\n",
                   test_name.c_str(), str_gconf_status.c_str());
      exit(EXIT_FAILURE);
    }

    GaugeConfig gconf_write(str_gconf_write);

    vout.crucial("gauge config. finished\n");
    Communicator::sync();

    unique_ptr<Action> action_G(
                       Action::New(str_action_G_type, params_action_G));

    //typedef Director_alt_Smear<AFIELD>  Director_alt;
    //unique_ptr<Director_alt> dr_smear(new Director_alt(params_smear));

    unique_ptr<AFopr<AFIELD> > fopr_dw(
                           AFopr<AFIELD>::New("Domainwall_eo", params_dw));

    unique_ptr<AForce_F<AFIELD> > force_fopr_dw(
                          new AForce_F_Domainwall_eo<AFIELD>(params_dw));


    //  unique_ptr<AFopr<AFIELD> >  fopr_smear(
    //   AFopr<AFIELD>::New("Smeared_alt", fopr_dw.get(), dr_smear.get()));

    //unique_ptr<AForce_F<AFIELD>> force_fopr_smear(new
    //   AForce_F_Smeared_alt<AFIELD>(force_fopr_dw.get(), dr_smear.get()));


    unique_ptr<Fprop_alt<AFIELD> > fprop_MD(
       new Fprop_alt_Standard_eo<AFIELD>(params_dw, params_solver_MD));

    unique_ptr<Fprop_alt<AFIELD> > fprop_H(
       new Fprop_alt_Standard_eo<AFIELD>(params_dw, params_solver_H));

    unique_ptr<AFopr<AFIELD> > fopr_dw_pv(
                         AFopr<AFIELD>::New("Domainwall_eo", params_dw_pv));

    unique_ptr<AForce_F<AFIELD> > force_fopr_dw_pv(
                        new AForce_F_Domainwall_eo<AFIELD>(params_dw_pv));

    //unique_ptr<AFopr<AFIELD> >  fopr_smear_pv(
    //   AFopr<AFIELD>::New("Smeared_alt", fopr_dw_pv.get(), dr_smear.get()));

    //unique_ptr<AForce_F<AFIELD> > force_fopr_smear_pv(new
    //   AForce_F_Smeared_alt<AFIELD>(force_fopr_dw_pv.get(), dr_smear.get()));

    unique_ptr<Fprop_alt<AFIELD> > fprop_H_pv(
       new Fprop_alt_Standard_eo<AFIELD>(params_dw_pv, params_solver_H));

    unique_ptr<Action> action_F(new Action_F_alt_Ratio_eo<AFIELD>(
                        //fopr_smear_pv.get(), force_fopr_smear_pv.get(),
                        //fopr_smear.get(), force_fopr_smear.get(),
                        fopr_dw_pv.get(), force_fopr_dw_pv.get(),
                        fopr_dw.get(), force_fopr_dw.get(),
                        fprop_H_pv.get(), fprop_MD.get(), fprop_H.get()));

    ActionList actions(Nlevels);
    actions.append(level_action[0], action_F.get());
    actions.append(level_action[1], action_G.get());

    //  std::vector<Director *> directors(1);
    //  directors[0] = static_cast<Director *>(dr_smear.get());

    unique_ptr<Builder_Integrator> builder(
       //new Builder_Integrator(actions, directors, params_integrator));
       new Builder_Integrator(actions, params_integrator));
    Integrator *integrator = builder->build();

    unique_ptr<RandomNumbers> rand(new RandomNumbers_Mseries(i_conf));

    HMC_General hmc(actions, integrator, rand.get(), params_hmc);
    //HMC_General hmc(actions, directors, integrator, rand.get(), params_hmc);

    vout.crucial(vl, "Object construction finished.\n");
    Communicator::sync();

    Timer timer(test_name);


    // ####  Execution main part  ####
    timer.start();

    vout.general(vl, "HMC: Ntraj = %d\n", Ntraj);

    double result = 0.0;
    for (int traj = 0; traj < Ntraj; ++traj) {
      vout.general(vl, "\n");
      vout.general(vl, "---------------------------------------------------\n");
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
} // namespace Test_HMC_Domainwall
