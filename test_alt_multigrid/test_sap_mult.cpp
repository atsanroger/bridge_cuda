/*!

        @file    test_sap_mult.cpp

        @brief   test for SAP operator

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: #$

        @version $LastChangedRevision: 2533 $
*/
//====================================================================
#include "test_sap_mult.h"


#include "lib/IO/gaugeConfig.h"
#include "lib/Tools/timer.h"
#include "lib/Field/field.h"
#include "lib/Field/field_F.h"
#include "lib/Field/field_G.h"


//#include "bridge.h"
#include "lib/Communicator/communicator.h"
#include "lib/ResourceManager/threadManager.h"

#include "lib/Parameters/commonParameters.h"
#include "lib/Parameters/parameters.h"
#include "lib/Parameters/parameterManager_YAML.h"

#include "lib/Tools/randomNumberManager.h"

#include "lib/Tools/timer.h"
#ifdef USE_MPI
#include "mpi.h"
#endif

#include "lib/IO/bridgeIO.h"
#include "show_version.h"

#ifdef USE_ALT_ACCEL
#include "lib_alt_Accel/bridge_alt_accel.h"
#endif

#ifdef USE_ALT_QXS
#include "lib_alt_QXS/bridge_alt_qxs.h"
#endif

#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_SAP.h"


// Registration as a test
namespace test_sap_mult {
  const string test_name = "Test_SapMult";

  int test(void)
  {
    Test_SapMult test;
    return test.test();
  }
}

// class name
const std::string Test_SapMult::class_name =
  "Test_SapMult";

//====================================================================
void Test_SapMult::init()
{
  // do nothing.
}


//====================================================================
int Test_SapMult::test()
{
  //    using namespace Test_Simd_Common;

  // ####  parameter setup  ####
  int Nc   = CommonParameters::Nc();
  int Nd   = CommonParameters::Nd();
  int Ndim = CommonParameters::Ndim();
  int Nvol = CommonParameters::Nvol();


  using real_t = float;
  //  using complex_t=typename AFIELD_f::complex_t;
  vout.general("\n*****************************************************\n");
  vout.general("*\n");
  vout.general("*\n");
  vout.general("*****************************************************\n\n");
  const std::string filename_input = "test_alt_SAP.yaml";
  const std::string test_name      = Test_SapMult::class_name;

  unique_ptr<Timer> timer(new Timer(test_name));
  timer->start();


  vout.general("\n");
  vout.general("test name: %s\n", test_name.c_str());


  Parameters params_all = ParameterManager::read(filename_input);

  Parameters params_gauge = params_all.lookup("Gauge");

  const string str_gconf_status = params_gauge.get_string("gauge_config_status");
  const string str_gconf_read   = params_gauge.get_string("gauge_config_type_input");
  const string readfile         = params_gauge.get_string("config_filename_input");

  Parameters params_mg_solver = params_all.lookup("MGSolver");
  //  const int    Nmult            = params_test.get_int("number_of_mult");
  const string str_vlevel = params_mg_solver.get_string("verbose_level");

  Parameters             params_coarse  = params_mg_solver.lookup("MultiGrid_Level1");
  const std::vector<int> sap_block_size = params_coarse.get_int_vector("sap_block");
  const int              num_vectors    = params_coarse.get_int("setup_number_of_vectors");
  //    const string solver_type = params_coarse.get_string("solver_type"); // MG

  Parameters   params_fopr = params_all.lookup("Fopr");
  const string fopr_type   = params_fopr.get_string("fermion_type");

  Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

  // setup random number manager
  RandomNumberManager::initialize("Mseries", 1234567UL);

  //- print input parameters
  vout.general(vl, "  gconf_status = %s\n", str_gconf_status.c_str());
  vout.general(vl, "  gconf_read   = %s\n", str_gconf_read.c_str());
  vout.general(vl, "  readfile     = %s\n", readfile.c_str());
  //  vout.general(vl, "  Nmult        = %d\n", Nmult);
  vout.general(vl, "  vlevel       = %s\n", str_vlevel.c_str());
  vout.general(vl, "  sap_block_size        = %s\n", Parameters::to_string(sap_block_size).c_str());
  vout.general(vl, "  num_vectors  = %d\n", num_vectors);
  //    vout.general(vl, "  solver_type  = %d\n", solver_type.c_str());
  vout.general(vl, "  Fopr         = %s\n", fopr_type.c_str());
  const double kappa = params_fopr.get_double("hopping_parameter");
  const double csw   = params_fopr.get_double("clover_coefficient");
  vout.general(vl, "     kappa  = %f\n", kappa);
  vout.general(vl, "     CSW    = %f\n", csw);

  //  const string fopr_type_dd = fopr_type + "_dd";
  if (fopr_type != "Clover") {
    vout.crucial("only for Clover fermion (with Fopr = Clover), sorry\n");
    exit(EXIT_FAILURE);
  }
  //  Parameters params_fopr_dd = params_fopr;
  //  params_fopr_dd.set_string("fopr_type", fopr_type_dd.c_str());
  //- input parameter check
  int err = 0;
  err += ParameterCheck::non_NULL(str_gconf_status);

  if (err) {
    vout.crucial(vl, "Error at %s: input parameters have not been set\n", test_name.c_str());
    exit(EXIT_FAILURE);
  }

  unique_ptr<Timer> timer_total(new Timer(test_name));
  timer_total->start();


  //RandomNumberManager::initialize("Mseries", 1234567UL);
  RandomNumbers *random = RandomNumberManager::getInstance();

  // ####  Setup gauge configuration  ####
  vout.general("Field, creating: U\n");
  U.reset(new Field_G(Nvol, Ndim));

  if (str_gconf_status == "Continue") {
    GaugeConfig(str_gconf_read).read(*U.get(), readfile);
  } else if (str_gconf_status == "Cold_start") {
    GaugeConfig("Unit").read(*U.get());
  } else if (str_gconf_status == "Hot_start") {
    GaugeConfig("Random").read(*U.get());
  } else {
    vout.crucial(m_vl, "Error at %s: unsupported gconf status \"%s\"\n",
                 test_name.c_str(), str_gconf_status.c_str());
    exit(EXIT_FAILURE);
  }

  // Source vector
  Field_F b;

  Field_F y;
  double  norm_b;
  //#pragma omp parallel
  {
    random->uniform_lex_global(b);
    norm_b = b.norm2();
  }
  vout.general(vl, "|source|   =%.8f\n", sqrt(norm_b));
  vout.general(vl, "|source|^2 =%.8f\n", norm_b);


  Parameters params_solver;  // for a moment this is dummy
  params_solver.set_int("maximum_number_of_iteration", 10);
  params_solver.set_double("convergence_criterion_squared", 0.0);
  // params_solver.set_string("verbose_level", "Detailed");

#ifdef USE_ALT_ACCEL
  const Impl IMPL = ACCEL;
#endif

#ifdef USE_ALT_QXS
  const Impl IMPL = QXS;
#endif

  typedef AField<float, IMPL> AFIELD_f;

  Fprop *fprop = new Fprop_alt_Standard_SAP<AFIELD_f>(
    params_fopr, params_solver);
  fprop->set_config(U.get());
  fprop->set_mode("D");

  int    nconv = -1;
  double diff  = 0.0;
  fprop->invert(y, b, nconv, diff);


  return EXIT_SUCCESS;
}


//============================================================END=====
