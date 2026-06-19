/*!
         @file    $Id:: main.cpp #$

         @brief

         @author  Hideo Matsufuru (matsufuru)
                  $LastChangedBy: matufuru $

         @date    $LastChangedDate:: 2023-07-02 11:40:28 #$

         @version $LastChangedRevision: 2528 $
*/

//#include "bridge.h"
//using Bridge::vout;

// corelib
#include "lib/Communicator/communicator.h"
#include "lib/ResourceManager/threadManager.h"
#include "lib/Parameters/commonParameters.h"
#include "lib/Parameters/parameters.h"
#include "lib/Parameters/parameterManager_YAML.h"
#include "lib/Tools/timer.h"
#include "lib/bridge_init_factory.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#ifdef USE_ALT_CODE
#include "lib_alt/bridge_alt_init.h"
#endif

// local
#include "Tests/testlist.h"

// device
//#include "lib_alt_OpenACC/ResourceManager/deviceManager_OpenACC.h"


const std::string filename_main_input = "main.yaml";
// const std::string filename_main_input = "stdin";

//- prototype declarations
int run_test();


//====================================================================
int main(int argc, char *argv[])
{
  // ###  initial setup  ###
  Bridge::VerboseLevel vl = Bridge::GENERAL;

  Communicator::init(&argc, &argv);

  // ####  banner  ####
  vout.general(vl, "Bridge++ %s\n\n", BRIDGE_VERSION);

  std::string filename_input = filename_main_input;
  if (filename_input == "stdin") {
    vout.general(vl, "input filename : ");
    std::cin >> filename_input;
    vout.general(vl, "%s\n", filename_input.c_str());
  } else {
    vout.general(vl, "input filename : %s\n", filename_input.c_str());
  }
  vout.general(vl, "\n");

  Parameters params_all  = ParameterManager::read(filename_input);
  Parameters params_main = params_all.lookup("Main");

  const vector<int> lattice_size = params_main.get_int_vector("lattice_size");
  const vector<int> grid_size    = params_main.get_int_vector("grid_size");
  const int number_of_thread     = params_main.get_int("number_of_thread");
  const int device_number_offset = params_main.get_int("device_number_offset");
  const int number_of_color      = params_main.get_int("number_of_color");
  const string str_logfile       = params_main.get_string("log_filename");
  const string str_ildg_logfile  = params_main.get_string("ildg_log_filename");
  const string str_vlevel        = params_main.get_string("verbose_level");

  //- initializations
  vl = vout.set_verbose_level(str_vlevel);
  CommonParameters::init_Vlevel(vl);

  if (str_logfile != "stdout") {
    vout.init(str_logfile);
  }

  if (str_ildg_logfile != "stdout") {
    vout.ildg_init(str_ildg_logfile);
  }


  // CommonParameters::init(lattice_size, grid_size);
  CommonParameters::init(lattice_size, grid_size, number_of_color);
  Communicator::setup();

  ThreadManager_OpenMP::init(number_of_thread);


  //- print input parameters
  vout.general(vl, "Main: input parameters\n");
  vout.general(vl, "  lattice_size     = %s\n", Parameters::to_string(lattice_size).c_str());
  vout.general(vl, "  grid_size        = %s\n", Parameters::to_string(grid_size).c_str());

  vout.general(vl, "  number of thread = %d\n", number_of_thread);
  vout.general(vl, "  device number offset = %d\n", device_number_offset);
  vout.general(vl, "  number of color  = %d\n", number_of_color);
  vout.general(vl, "  logfile          = %s\n", str_logfile.c_str());
  vout.general(vl, "  ildg_logfile     = %s\n", str_ildg_logfile.c_str());
  vout.general(vl, "  vlevel           = %s\n", str_vlevel.c_str());
  //vout.general(vl, "\n");
  vout.crucial(vl, "\n");

  //- input parameter check
  int err = 0;
  err += ParameterCheck::non_NULL(str_logfile);
  err += ParameterCheck::non_NULL(str_ildg_logfile);

  if (err) {
    vout.crucial(vl, "Error at main: input parameters have not been set.\n");
    exit(EXIT_FAILURE);
  }

#ifdef USE_FACTORY
#ifdef USE_FACTORY_AUTOREGISTER
#else
  bridge_init_factory();
  vout.crucial(vl, "bridge_init_factory finished.\n");
#endif


#ifdef DEBUG
  // bridge_report_factory();
#endif
#endif

  // initialization of alt-code
#ifdef USE_ALT_CODE
  bridge_alt_init(params_main);
  vout.crucial(vl, "bridge_alt_init finished.\n");
#endif
  
  //- timestamp (starting time)
  unique_ptr<Timer> timer(new Timer("Main"));
  timer->timestamp();
  timer->start();

  //####  here function is called explicitly  ####

  // this fails.
  //Test_HMC_Staggered_Leapfrog::update_alt();


  //###  Wilson  ###

  // this works
  //Test_HMC_Wilson::update_Nf2_nosmr();

  // this works
  //Test_HMC_Wilson::update_Nf2_alt_nosmr();

  // this works
  //Test_HMC_Wilson::update_Nf2();

  // this works
  //Test_HMC_Wilson::update_Nf2_alt_Field();

  // this works
  //Test_HMC_Wilson::update_Nf2_alt();

  // this works
  //Test_HMC_Wilson::update_Nf2_alt_eo_nosmr();

  // this works
  //Test_HMC_Wilson::update_Nf2_alt_eo();

  // this works
  //Test_HMC_Wilson::update_Nf2_alt_eo_ratio();

  // this works
  //Test_HMC_Wilson::RHMC_Nf2p1();

  // this works
  //Test_HMC_Wilson::RHMC_Nf2p1_alt();

  // this works
  //Test_HMC_Wilson::RHMC_Nf2p1_alt_eo();

  // this works
  //Test_HMC_Wilson::RHMC_Nf2p1_alt_eo_ratio();

  //###  Domainwall  ###

  // this is original test (with smearing)
  //Test_HMC_Domainwall::update_Nf2_PV();

  // alt_Field successfully works (with smearing)
  //Test_HMC_Domainwall::update_Nf2_PV_alt_Field();

  // this works
  //Test_HMC_Domainwall::update_Nf2_PV_alt();

  // this works
  //Test_HMC_Domainwall::update_Nf2_PV_nosmr();

  // this works
  //Test_HMC_Domainwall::update_Nf2_PV_alt_nosmr();

  // this works
  //Test_HMC_Domainwall::update_Nf2_PV_alt_eo_nosmr();

  // this works
  //Test_HMC_Domainwall::update_Nf2_PV_alt_eo();

  //
  Test_HMC_Domainwall::update_Nf2p1_PV_alt();

  // this works
  //Test_HMC_Domainwall::update_Nf2p1_PV_alt_eo();

  //- timestamp (end time)
  timer->report();
  timer->timestamp();

  ThreadManager_OpenMP::finalize();
  Communicator::finalize();

  return EXIT_SUCCESS;

}

//============================================================END=====
