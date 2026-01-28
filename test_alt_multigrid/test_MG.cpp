/*!
         @file    $Id:: test_MG.cpp #$

         @brief  test program for multigrid solver

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $

         @date    $LastChangedDate:: 2022-12-16 15:57:38 #$

         @version $LastChangedRevision: 2422 $
*/
//====================================================================

//#include "bridge.h"
#include "lib/bridge_setup.h"
#include "lib/Communicator/communicator.h"
#include "lib/ResourceManager/threadManager.h"

#include "lib/Parameters/commonParameters.h"
#include "lib/Parameters/parameters.h"
#include "lib/Parameters/parameterManager_YAML.h"

#include "lib/Tools/timer.h"
#ifdef USE_MPI
#include "mpi.h"
#endif

#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#ifdef USE_ALT_OPENACC
#include "lib_alt_OpenACC/init_alt_OpenACC.h"
#endif

#ifdef USE_ALT_QXS
#include "lib_alt_QXS/bridge_init_factory_alt_QXS.h"
#endif

#include "show_version.h"

const std::string filename_main_input = "main.yaml";
// const std::string filename_main_input = "stdin";

//- prototype declarations
int run_test_MG();
int run_test_afopr_coarse();

//====================================================================
int main(int argc, char *argv[])
{
  // ###  initial setup  ###
  Bridge::VerboseLevel vl = Bridge::GENERAL;
  bridge_initialize(&argc, &argv);
  Bridge::show_version();

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

  //  vout.general("hoge hoge\n");
  const std::vector<int> lattice_size     = params_main.get_int_vector("lattice_size");
  const std::vector<int> grid_size        = params_main.get_int_vector("grid_size");
  const int              number_of_thread = params_main.get_int("number_of_thread");
  const int              number_of_color  = params_main.get_int("number_of_color");
  const std::string      str_logfile      = params_main.get_string("log_filename");
  const std::string      str_ildg_logfile = params_main.get_string("ildg_log_filename");
  const std::string      str_vlevel       = params_main.get_string("verbose_level");


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
  vout.general("calling Communicator::setup\n");
  Communicator::setup();
  vout.general("calling Communicator::setup, done\n");

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

#ifdef USE_ALT_OPENACC
  init_alt_OpenACC(params_main);
#endif

#ifdef USE_ALT_QXS
  bridge_init_factory_alt_QXS();
#endif



  //- timestamp (starting time)
  unique_ptr<Timer> timer(new Timer("Main"));
  timer->timestamp();
  timer->start();

  //run_test_afopr_coarse();
  run_test_MG();

  //- timestamp (end time)
  timer->report();
  timer->timestamp();

  bridge_finalize();
  //ThreadManager::finalize();
  //  Communicator::finalize();

  return EXIT_SUCCESS;
}


//====================================================================
//============================================================END=====
