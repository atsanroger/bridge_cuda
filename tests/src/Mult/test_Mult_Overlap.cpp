/*!
        @file    test_Mult_Overlap.cpp

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate: 2013-04-27 12:28:50 #$

        @version $LastChangedRevision: 2422 $
*/

#include "test.h"

#include "Eigen/eigensolver_IRLanczos.h"

#include "Fopr/fopr_Wilson.h"
#include "Fopr/fopr_Overlap.h"

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
//! Test of Mult with Overlap fermion.

/*!
                               [06 Jun 2014 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                               [21 Mar 2015 Y.Namekawa]
*/

namespace Test_Mult_Overlap {
  const std::string test_name = "Mult.Overlap";

  //- test-private parameters
  namespace {
    const std::string filename_input = "test_Mult_Overlap.yaml";
  }

  //- prototype declaration
  int mult(void);

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
#if defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
    static const bool is_registered = TestManager::RegisterTest(
      test_name,
      mult
      );
#endif
  }
#endif

  //====================================================================
  int mult(void)
  {
    // ####  parameter setup  ####
    const int Ndim = CommonParameters::Ndim();
    const int Nc   = CommonParameters::Nc();
    const int Nd   = CommonParameters::Nd();
    const int Nvol = CommonParameters::Nvol();

    const int NPE     = CommonParameters::NPE();
    const int Nthread = ThreadManager::get_num_threads_available();

    const Parameters params_all = ParameterManager::read(filename_input);

    const Parameters params_test      = params_all.lookup("Test_Mult_Overlap");
    const Parameters params_wilson    = params_all.lookup("Fopr_Wilson");
    const Parameters params_overlap   = params_all.lookup("Fopr_Overlap");
    const Parameters params_irlanczos = params_all.lookup("Eigensolver");

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


    // ####  object setup  ####
    unique_ptr<Fopr_Wilson> fopr_w(new Fopr_Wilson(params_wilson));
    fopr_w->set_config(&U);

    //-  with low-mode subtraction
    // int Nsbt = 0;
    // std::vector<double> TDa;
    // std::vector<Field>  vk;
    // calc_lowmodes(params_irlanczos, Nsbt, TDa, vk, fopr_w);

    unique_ptr<Fopr> fopr_overlap(new Fopr_Overlap(fopr_w.get(), params_overlap));
    fopr_overlap->set_config(&U);
    fopr_overlap->set_mode("D");
    // dynamic_cast<Fopr_Overlap*>(fopr_overlap)->set_lowmodes(Nsbt, &TDa, &vk);

    Timer timer(test_name);


    // ####  Execution main part  ####
    Field_F b, y;
    b.set(1.0);

    double result = 0.0;

    timer.start();
#ifdef FUJITSU_FX
    //- profiler starts
    // fapp_start("Mult.Overlap",1,1);
#endif

#pragma omp parallel
    {
      for (int i = 0; i < Nmult; ++i) {
        fopr_overlap->mult(y, b);
      }
      result = y.norm();
    }

#ifdef FUJITSU_FX
    //- profiler ends
    // fapp_stop("Mult.Overlap",1,1);
#endif
    timer.stop();
    const double elapse_sec = timer.elapsed_sec();


    //- Flops counting in giga unit
    vout.general(vl, "Warning at %s: flop_count has not been implemented.\n", test_name.c_str());
    vout.general(vl, "\n");


    RandomNumberManager::finalize();


    if (do_check) {
      return Test::verify(result, expected_result, tolerance);
    } else {
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }
  }
} // namespace Test_Mult_Overlap
