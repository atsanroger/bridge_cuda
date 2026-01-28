/*!
        @file    test_Spectrum_Overlap_checkSignFunction.cpp

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$

        @version $LastChangedRevision: 2492 $
*/

#include "test.h"

#include "Fopr/fopr_Wilson.h"
#include "Fopr/fopr_Sign.h"

#include "IO/gaugeConfig.h"

#include "Measurements/Fermion/corr2pt_4spinor.h"
#include "Measurements/Fermion/source.h"
#include "Measurements/Gauge/gaugeFixing.h"

#include "Tools/randomNumberManager.h"

//====================================================================
//! Test class of sign function for overlap fermion operators.

/*!
    This class tests sign function for overlap fermion operator.
                                          [13 Apr 2012 H.Matsufuru]
    (Coding history will be recovered from trac.)
    Implement Selectors.                  [03 Mar 2013 Y.Namekawa]
    (Selectors are replaced with factories by Aoyama-san)
    Introduce unique_ptr to avoid memory leaks.
                                          [21 Mar 2015 Y.Namekawa]
    Add Nc check for USE_GROUP_SU_N.      [31 May 2021 Y.Namekawa]
 */

namespace Test_Spectrum_Overlap {
  const std::string test_name = "Spectrum.Overlap.CheckSignFunction";

  //- test-private parameters
  namespace {
    const std::string filename_input = "test_Spectrum_Overlap_CheckSignFunction.yaml";
  }

  //- prototype declaration
  int check_sign(void);

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
#if defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
    static const bool is_registered = TestManager::RegisterTest(
      test_name,
      check_sign
      );
#endif
  }
#endif

  //====================================================================
  int check_sign(void)
  {
    // ####  parameter setup  ####
    const int Ndim = CommonParameters::Ndim();
    const int Nc   = CommonParameters::Nc();
    const int Nvol = CommonParameters::Nvol();

    const Parameters params_all = ParameterManager::read(filename_input);

    const Parameters params_test    = params_all.lookup("Test_Spectrum_Overlap");
    const Parameters params_gfix    = params_all.lookup("GaugeFixing");
    const Parameters params_wilson  = params_all.lookup("Fopr_Wilson");
    const Parameters params_overlap = params_all.lookup("Fopr_Overlap");
    const Parameters params_source  = params_all.lookup("Source");

    const string        str_gconf_status = params_test.get_string("gauge_config_status");
    const string        str_gconf_read   = params_test.get_string("gauge_config_type_input");
    const string        readfile         = params_test.get_string("config_filename_input");
    const string        str_rand_type    = params_test.get_string("random_number_type");
    const unsigned long seed             = params_test.get_unsigned_long("seed_for_random_number");
    const string        str_vlevel       = params_test.get_string("verbose_level");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

    const string str_gfix_type   = params_gfix.get_string("gauge_fixing_type");
    const string str_gmset_type  = params_wilson.get_string("gamma_matrix_type");
    const string str_source_type = params_source.get_string("source_type");

    const Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

    //- print input parameters
    vout.general(vl, "  gconf_status = %s\n", str_gconf_status.c_str());
    vout.general(vl, "  gconf_read   = %s\n", str_gconf_read.c_str());
    vout.general(vl, "  readfile     = %s\n", readfile.c_str());
    vout.general(vl, "  rand_type    = %s\n", str_rand_type.c_str());
    vout.general(vl, "  seed         = %lu\n", seed);
    vout.general(vl, "  vlevel       = %s\n", str_vlevel.c_str());
    vout.general(vl, "  gfix_type    = %s\n", str_gfix_type.c_str());
    vout.general(vl, "  gmset_type   = %s\n", str_gmset_type.c_str());
    vout.general(vl, "  source_type  = %s\n", str_source_type.c_str());

    //- input parameter check
    int err = 0;
    err += ParameterCheck::non_NULL(str_gconf_status);

    if (err) {
      vout.crucial(vl, "Error at %s: input parameters have not been set\n", test_name.c_str());
      exit(EXIT_FAILURE);
    }

    if ((str_gfix_type == "Coulomb") || (str_gfix_type == "Landau")) {
      if (CommonParameters::Nc() != 3) {
        vout.crucial(vl, "check skipped: Nc = 3 is needed, but Nc = %d.\n\n", CommonParameters::Nc());
        return EXIT_SKIP;
      }
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


    // ####  Gauge fixing  ####
    {
      Field_G                 Ufix(Nvol, Ndim);
      unique_ptr<GaugeFixing> gfix(GaugeFixing::New(str_gfix_type, params_gfix));

      gfix->fix(Ufix, U);
      copy(U, Ufix);
    }


    // ####  object setup  ####
    unique_ptr<Fopr> fopr_w(new Fopr_Wilson(params_wilson));
    fopr_w->set_config(&U);

    int                 Nsbt = 0;
    std::vector<double> TDa;
    std::vector<Field>  vk;

    //- with low-mode subtraction
    //double Vthreshold = 0.15;  // example value with low-mode subtraction
    //calc_lowmodes(Nsbt, TDa, vk, Vthreshold, fopr_w);

    unique_ptr<Fopr> fopr_sign(new Fopr_Sign(fopr_w.get(), params_overlap));
    fopr_sign->set_config(&U);
    // if (Nsbt > 0) fopr_sign->set_lowmodes(Nsbt, &TDa, &vk);

    unique_ptr<Source> source(Source::New(str_source_type, params_source));

    Timer timer(test_name);

    // ####  Execution main part  ####
    timer.start();

    int Nin = fopr_w->field_nin();
    int Nex = fopr_w->field_nex();
    // note that Nvol has been already set.

    Field b(Nin, Nvol, Nex);
    Field xq(Nin, Nvol, Nex);
    Field xq2(Nin, Nvol, Nex);

    double result = 1.0;

#pragma omp parallel
    {
      for (int ispin = 0; ispin < 1; ++ispin) {
        for (int icolor = 0; icolor < 1; ++icolor) {
          // for(int ispin = 0; ispin < Nd; ++ispin){
          //for(int icolor = 0; icolor < Nc; ++icolor){

          int icd = icolor + Nc * ispin;

          source->set(b, icd);

          fopr_sign->mult(xq, b);
          fopr_sign->mult(xq2, xq);

          axpy(xq2, -1.0, b);

          double xq_norm = xq2.norm2();
          if (icd == 0) {
            int ith = ThreadManager::get_thread_id();
            if (ith == 0) result = xq_norm;
          }

          vout.general(vl, " |sign^2(b) - b|^2 = %16.8e\n", xq_norm);
        }
      }
    }

    timer.report();

    RandomNumberManager::finalize();


    if (do_check) {
      return Test::verify(result, expected_result);
    } else {
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }
  }
} // namespace Test_Spectrum_Overlap
