/*!
    @file    test_HMC_Staggered_Leapfrog_alt.cpp
    @brief
    @author  Hideo Matsufuru <hideo.matsufuru@kek.jp> (matufuru)
             $LastChangedBy: matufuru $
    @date    $LastChangedDate: 2013-01-22 13:51:53 #$
    @version $LastChangedRevision: 2197 $
*/

// core library
#include "lib/IO/gaugeConfig.h"
#include "lib/Action/Fermion/action_F_Staggered_lex.h"
#include "lib/HMC/hmc_Leapfrog.h"
#include "lib/Fopr/fopr_Smeared.h"
#include "lib/Force/Fermion/force_F_Staggered.h"
#include "lib/Force/Fermion/force_F_Smeared.h"
//#include "lib/Smear/forceSmear.h"
#include "lib/Smear/forceSmear.h"
#include "lib/Solver/solver.h"
#include "lib/Smear/projection.h"
#include "lib/Smear/smear.h"
#include "lib/Smear/director_Smear.h"
#include "lib/Tools/randomNumbers_Mseries.h"
#include "lib/Tools/file_utils.h"
#include "lib/Measurements/Fermion/fprop_Standard_lex.h"
#include "lib/Measurements/Gauge/polyakovLoop.h"
#include "lib/Measurements/Gauge/staple_lex.h"

// alt-code
#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_lex.h"

#ifdef USE_ALT_ACCEL
#include "lib_alt_Accel/bridge_alt_accel.h"
#define IMPL ACCEL
#endif

#ifdef USE_ALT_QXS
#include "lib_alt_QXS/bridge_alt_qxs.h"
#define IMPL QXS
#endif

// local
#include "Tests/test.h"
#include "Tests/job_Utils.h"

//====================================================================
//! Test of HMC update for Staggered fermions.

/*!
    Staggeered HMC.
                                         [19 Aug 2019 H.Matsufuru]
 */

namespace Test_HMC_Staggered_Leapfrog
{



//====================================================================
int update_alt()
{
  string file_params = "test_HMC_Staggered_Leapfrog.yaml";
  string inputfile = "input";  // configuration numbers are stored

  string run_mode  = "test";
  //string run_mode  = "job";

  string code_mode = "org";
  //string code_mode = "alt";

  //  static const Impl IMPL = OPENACC;
  typedef AField<double, IMPL> AField_d;

    // #####  parameter setup  #####
    int Nc   = CommonParameters::Nc();
    int Nvol = CommonParameters::Nvol();
    int Ndim = CommonParameters::Ndim();

    Parameters params_all = ParameterManager::read(file_params);

    Parameters params_test      = params_all.lookup("Test_HMC_Staggered");
    Parameters params_action_G  = params_all.lookup("Action_G");
    Parameters params_action_F  = params_all.lookup("Action_F");
    Parameters params_fopr      = params_all.lookup("Fopr_Staggered");
    Parameters params_hmc       = params_all.lookup("HMC_Leapfrog");
    Parameters params_solver_MD = params_all.lookup("Solver_MD");
    Parameters params_solver_H  = params_all.lookup("Solver_H");

    string str_vlevel = params_test.get_string("verbose_level");
    Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

    // input config
    string config_type_input, config_file_input;
    if(run_mode == "test"){
      config_type_input  = params_test.get_string("config_type_test");
      config_file_input  = params_test.get_string("config_file_test");
    }else if(run_mode == "job"){
      config_type_input  = params_test.get_string("config_type_input");
      config_file_input  = params_test.get_string("config_file_input");
    }else{
      vout.crucial(vl, "%irrelevant run_mode =%s\n", run_mode.c_str());
      exit(EXIT_FAILURE);
    }

    // output config
    string config_type_output = params_test.get_string("config_type_output");
    string config_file_output = params_test.get_string("config_file_output");


    // for test mode
    bool   do_check         = false;
    double expected_result  = 0.0;
    if(run_mode == "test"){
      do_check         = params_test.is_set("expected_result");
      expected_result  =
             do_check ? params_test.get_double("expected_result") : 0.0;
    }


    //- print input parameters
    vout.general(vl, "  config_type_input  = %s\n", config_type_input.c_str());
    vout.general(vl, "  config_file_input  = %s\n", config_file_input.c_str());
    vout.general(vl, "  config_type_output = %s\n", config_type_output.c_str());
    vout.general(vl, "  config_file_output = %s\n", config_file_output.c_str());
    vout.general(vl, "  vlevel = %s\n", str_vlevel.c_str());
    vout.general(vl, "\n");


    //#####  setup of objects  #####

    unique_ptr<Field_G> U(new Field_G(Nvol, Ndim));

    unique_ptr<GaugeConfig> gconf_read(new GaugeConfig(config_type_input));
    unique_ptr<GaugeConfig> gconf_write(new GaugeConfig(config_type_output));

    // Gauge action
    string action_G_type = params_action_G.get_string("action_type");
    unique_ptr<Action> action_G(Action::New(action_G_type));
    action_G->set_parameters(params_action_G);


    // Fermion action: N_f=2 part
    unique_ptr<Fopr> fopr(new Fopr_Staggered());
    fopr->set_parameters(params_fopr);

    unique_ptr<Force> force_fopr(new Force_F_Staggered());
    force_fopr->set_parameters(params_fopr);


    // setup fermion propagators
    unique_ptr<Solver> solver_MD; // used only for code_mode = "org"
    unique_ptr<Solver> solver_H;  // used only for code_mode = "org"

    unique_ptr<Fprop> fprop_MD;
    unique_ptr<Fprop> fprop_H;

    if(code_mode == "org"){
      string solver_MD_type = params_solver_MD.get_string("solver_type");
      solver_MD.reset(Solver::New(solver_MD_type, fopr.get()));
      solver_MD->set_parameters(params_solver_MD);
      //solver_MD.reset(Solver::New(solver_MD_type, fopr, params_solver_MD);
      fprop_MD.reset(new Fprop_Standard_lex(solver_MD.get()));

      string solver_H_type = params_solver_H.get_string("solver_type");
      solver_H.reset(Solver::New(solver_H_type, fopr.get()));
      solver_H->set_parameters(params_solver_H);
      fprop_H.reset(new Fprop_Standard_lex(solver_H.get()));
    }else if(code_mode == "alt"){
      fprop_MD.reset(new Fprop_alt_Standard_lex<AField_d>(
                                       params_fopr, params_solver_MD));

      fprop_H.reset(new Fprop_alt_Standard_lex<AField_d>(
                                       params_fopr, params_solver_H));
    }else{
      vout.crucial(vl, "%irrelevant code_mode =%s\n", code_mode.c_str());
      exit(EXIT_FAILURE);
    }

    // fermion action
    unique_ptr<Action> action_F(new
          Action_F_Staggered_lex(fopr, force_fopr, fprop_MD, fprop_H));


    ActionList actions(1);        // one level
    actions.append(0, action_F.get());  // register actions at level0
    actions.append(0, action_G.get());

    //- Random number is initialized with a parameter specified by i_conf
    //  unique_ptr<RandomNumbers> rand(new RandomNumbers_Mseries(i_conf));

    unique_ptr<PolyakovLoop> ploop(new PolyakovLoop());

    unique_ptr<Staple> staple(new Staple_lex());

    unique_ptr<Timer> timer(new Timer);


    int iconf, nconf;

    if(run_mode == "test"){
      iconf = 1;
      nconf = 1;
      gconf_read->read_file(*U, config_file_input);
    }else{
      read_input(iconf, nconf, inputfile);
      vout.general(vl, "initial iconf   = %d\n", iconf);
      vout.general(vl, "number of confs = %d\n", nconf);
      if(iconf == 1){  // cold start
        U->set_unit();
      }else{           // continued
        gconf_read->read_file(*U, filename_config(config_file_input, iconf-1));
      }

    }

    double result = 0.0;

    // ####  Execution main part  ####
    for(int i = 0; i < nconf; ++i){

      vout.general(vl, "\n");
      vout.general(vl, "iconf = %d\n", iconf);

      unique_ptr<RandomNumbers> rand(new RandomNumbers_Mseries(iconf));

      //HMC_Leapfrog hmc(actions, rand);
      //hmc.set_parameters(params_hmc);
      HMC_Leapfrog hmc(actions, rand.get(), params_hmc);

      int Ntraj = params_hmc.get_int("number_of_trajectory");
      // this is not a genuine parameter of HMC class

      timer->start();

      vout.general(vl, "HMC: Ntraj = %d\n", Ntraj);

      for (int traj = 0; traj < Ntraj; ++traj) {

        vout.general(vl, "\n");
        vout.general(vl, "traj = %d\n", traj);

        hmc.update(*U);

        vout.general(vl, "\n");
        vout.general(vl, "Measurements: iconf = %d  traj = %d\n", iconf, traj);
        double plaq  = staple->plaquette(*U);
        vout.general(vl, "  Plaquette     = %16.12f\n", plaq);

        dcomplex plp = ploop->measure_ploop(*U);
        vout.general(vl, "  Polyakov loop = (%16.12f, %16.12f)\n",
                     real(plp), imag(plp));

        result = plaq;
      }

      if(run_mode == "job")
        gconf_write->write_file(*U,
                         filename_config(config_file_output, iconf));

      ++iconf;
      if(run_mode == "job") write_input(iconf, nconf, inputfile);

      timer->report();

    }  // i-loop (configuration)


    if(do_check){
      return Test::verify(result, expected_result);
    }else{
      if(run_mode == "test"){
        vout.general(vl, "check skipped: expected_result not set.\n\n");
      }
      return EXIT_SKIP;
    }

}

} // namespace Test_HMC_Wilson
