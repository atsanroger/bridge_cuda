/*!
        @file    test_Mult_Wilson.cpp

        @brief

        @author  Yusuke Namekawa  (namekawa)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2022-12-16 15:57:38 #$

        @version $LastChangedRevision: 2422 $
*/

#include "test.h"

#include "Fopr/fopr_Wilson.h"

#include "IO/gaugeConfig.h"

#include "Tools/randomNumberManager.h"

//- profiler
#ifdef FUJITSU_FX
#  ifdef __has_include
#    if __has_include(<fj_tool/fapp.h>)
#      include <fj_tool/fapp.h>
#    endif
#  endif
#endif

//====================================================================
//! Test of Mult with Wilson fermion, prepared for beginners.

/*!
                               [06 Jun 2014 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                               [21 Mar 2015 Y.Namekawa]
 */

namespace Test_Mult_Wilson {
  const std::string test_name = "Mult.Wilson";

  //- test-private parameters
  namespace {
    const std::string filename_input = "test_Mult_Wilson.yaml";
  }

  //- prototype declaration
  int mult(void);

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
    static const bool is_registered = TestManager::RegisterTest(
      test_name,
      mult
      );
  }
#endif

  //====================================================================
  int mult(void)
  {
    // ####  parameter setup  ####
    const int Nc   = CommonParameters::Nc();
    const int Nd   = CommonParameters::Nd();
    const int Ndim = CommonParameters::Ndim();
    const int Nvol = CommonParameters::Nvol();

    const int NPE     = CommonParameters::NPE();
    const int Nthread = ThreadManager::get_num_threads_available();

    const Parameters params_all = ParameterManager::read(filename_input);

    const Parameters params_test = params_all.lookup("Test_Mult");
    const Parameters params_fopr = params_all.lookup("Fopr");

    const string        str_gconf_status = params_test.get_string("gauge_config_status");
    const string        str_gconf_read   = params_test.get_string("gauge_config_type_input");
    const string        readfile         = params_test.get_string("config_filename_input");
    const string        str_rand_type    = params_test.get_string("random_number_type");
    const unsigned long seed             = params_test.get_unsigned_long("seed_for_random_number");
    const int           Nmult            = params_test.get_int("number_of_mult");
    const string        str_vlevel       = params_test.get_string("verbose_level");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;
    const double tolerance       = params_test.get_double("tolerance");

    const Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

    //- print input parameters
    vout.general(vl, "  gconf_status = %s\n", str_gconf_status.c_str());
    vout.general(vl, "  gconf_read   = %s\n", str_gconf_read.c_str());
    vout.general(vl, "  readfile     = %s\n", readfile.c_str());
    vout.general(vl, "  rand_type    = %s\n", str_rand_type.c_str());
    vout.general(vl, "  seed         = %lu\n", seed);
    vout.general(vl, "  Nmult        = %d\n", Nmult);
    vout.general(vl, "  vlevel       = %s\n", str_vlevel.c_str());

    //- input parameter check
    int err = 0;
    err += ParameterCheck::non_NULL(str_gconf_status);

    if (err) {
      vout.crucial(vl, "Error at %s: input parameters have not been set\n", test_name.c_str());
      exit(EXIT_FAILURE);
    }


    RandomNumberManager::initialize(str_rand_type, seed);


    // ####  Set up a gauge configuration  ####
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


    // ####  object setup  #####
    Fopr_Wilson fopr(params_fopr);
    fopr.set_config(&U);
    fopr.set_mode("D");

    Timer timer(test_name);

    vout.general(vl, "\n");


    // ####  Execution main part  ####
    Field_F b, y;
    b.set(1.0);

    double result = 0.0;

    timer.start();
#ifdef FUJITSU_FX
    //- profiler starts
    // fapp_start("Mult.Wilson",1,1);
#endif

#pragma omp parallel
    {
      for (int i = 0; i < Nmult; ++i) {
        fopr.mult(y, b);
      }
      result = y.norm();
    }

#ifdef FUJITSU_FX
    //- profiler ends
    // fapp_stop("Mult.Wilson",1,1);
#endif
    timer.stop();
    const double elapse_sec = timer.elapsed_sec();


    //- Flops counting in giga unit
    const double gflops_mult = fopr.flop_count() * Nmult / (elapse_sec * NPE * Nthread);

    vout.general(vl, "%s: %12.4e GFlops / core\n", test_name.c_str(), gflops_mult);
    vout.general(vl, "\n");


    RandomNumberManager::finalize();


    if (do_check) {
#if defined(USE_GROUP_SU2)
      // Nc=2 verificaction data is not available.
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
#else
      return Test::verify(result, expected_result, tolerance);
#endif
    } else {
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }
  }
} // namespace Test_Mult_Wilson
