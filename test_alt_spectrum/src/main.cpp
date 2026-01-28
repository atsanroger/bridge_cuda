/*!
         @file    main.cpp
         @brief
         @author  Hideo Matsufuru (matufuru)
                  $LastChangedBy: matufuru $
         @date    $LastChangedDate:: 2023-03-20 12:31:44 #$
         @version $LastChangedRevision: 2500 $
*/

#include <vector>
#include <string>
using namespace std;

//#include "bridge.h"

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

// alt-code
#ifdef USE_ALT_CODE
#include "lib_alt/bridge_alt_init.h"
#endif

// local
#include "Tests/testlist_alt.h"


const std::string filename_main_input = "main.yaml";
// const std::string filename_main_input = "stdin";

//- prototype declarations
//int run_test();

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

  const vector<int> lattice_size     = params_main.get_int_vector("lattice_size");
  const vector<int> grid_size        = params_main.get_int_vector("grid_size");
  const int         number_of_thread = params_main.get_int("number_of_thread");
  const int         number_of_color  = params_main.get_int("number_of_color");
  const string      str_logfile      = params_main.get_string("log_filename");
  const string      str_ildg_logfile = params_main.get_string("ildg_log_filename");
  const string      str_vlevel       = params_main.get_string("verbose_level");


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

  ThreadManager::init(number_of_thread);

  //- print input parameters
  vout.general(vl, "Main: input parameters\n");
  vout.general(vl, "  lattice_size     = %s\n", Parameters::to_string(lattice_size).c_str());
  vout.general(vl, "  grid_size        = %s\n", Parameters::to_string(grid_size).c_str());

  vout.general(vl, "  number of thread = %d\n", number_of_thread);
  vout.general(vl, "  number of color  = %d\n", number_of_color);
  vout.general(vl, "  logfile          = %s\n", str_logfile.c_str());
  vout.general(vl, "  ildg_logfile     = %s\n", str_ildg_logfile.c_str());
  vout.general(vl, "  vlevel           = %s\n", str_vlevel.c_str());
  vout.general(vl, "\n");

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
  // bridge_report_factory();
#endif
#endif


  // initialization of alt-code
#ifdef USE_ALT_CODE
  bridge_alt_init(params_main);
#endif

  //- timestamp (starting time)
  unique_ptr<Timer> timer(new Timer("Main"));
  timer->timestamp();
  timer->start();

  //####  here function is called explicitly  ####

  //  Test_alt_Corelib::test_all();

#ifdef USE_ALT_QXS
  Test_alt_QXS::test_all();
#endif

#ifdef USE_ALT_ACCEL
  Test_alt_Accel::test_all();
#endif

#ifdef USE_ALT_SIMD2
  Test_alt_SIMD2::test_all();
#endif

#ifdef USE_ALT_VECTOR
  Test_alt_Vector::test_all();
#endif

  //- timestamp (end time)
  timer->report();
  timer->timestamp();

  // finalization of alt-code
#ifdef USE_ALT_CODE
  bridge_alt_fin();
#endif


  ThreadManager::finalize();
  Communicator::finalize();

  return EXIT_SUCCESS;
}


//============================================================END=====
