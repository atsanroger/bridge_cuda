/*!
        @file    test_Spectrum_Nonrelativistic_2ptFunction.cpp

        @brief

        @author  Hideo Matsufuru  (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate: 2013-04-27 12:28:50 #$

        @version $LastChangedRevision: 2492 $
*/

#include "test.h"

#include "lib/Fopr/fopr_NonRelativistic.h"
#include "lib/IO/gaugeConfig.h"
#include "lib/Measurements/Gauge/gaugeFixing.h"
#include "lib/Tools/gammaMatrixSet.h"
#include "lib/Tools/randomNumberManager.h"
#include "lib/Measurements/Fermion/corr2pt_4spinor.h"
#include "lib/Measurements/Fermion/source.h"

//====================================================================
//! Test class of NonRelativistic fermion operators.

/*!
    This class tests NonRelativistic fermion.
    Now heavy-heavy correlator is implemented.
                                        [30 Dec 2022 H.Matsufuru]
*/

namespace Test_Spectrum_NonRelativistic {
  const std::string test_name
    = "Spectrum.NonRelativistic.heavy-heavy-2ptFunction";

  //- test-private parameters
  namespace {
    const std::string parameter_file
      = "test_Spectrum_NonRelativistic_heavy_heavy_2ptFunction.yaml";
  }

  //- prototype declaration
  int heavy_heavy_2ptFunction();

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
#if defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
    static const bool is_registered = TestManager::RegisterTest(
      test_name,
      heavy_heavy_2ptFunction
      );
#endif
  }
#endif

//====================================================================
  int heavy_heavy_2ptFunction()
  {
    using namespace std;

    // ####  parameter setup  ####
    const int Nc   = CommonParameters::Nc();
    const int Nd   = CommonParameters::Nd();
    const int Ndim = CommonParameters::Ndim();
    const int Nvol = CommonParameters::Nvol();

    Parameters params_all = ParameterManager::read(parameter_file);

    Parameters params_spectrum     = params_all.lookup("Spectrum");
    Parameters params_gfix         = params_all.lookup("GaugeFixing");
    Parameters params_fopr_heavy   = params_all.lookup("Fopr_heavy");
    Parameters params_source_heavy = params_all.lookup("Source_heavy");

    string gconf_status = params_spectrum.get_string("gauge_config_status");
    string gconf_type_read
      = params_spectrum.get_string("gauge_config_type_input");
    string readfile  = params_spectrum.get_string("config_filename_input");
    string rand_type = params_spectrum.get_string("random_number_type");

    unsigned long seed
      = params_spectrum.get_unsigned_long("seed_for_random_number");
    string vlevel = params_spectrum.get_string("verbose_level");

    const bool   do_check = params_spectrum.is_set("expected_result");
    const double expected_result
      = do_check ? params_spectrum.get_double("expected_result") : 0.0;

    string gfix_type         = params_gfix.get_string("gauge_fixing_type");
    string fopr_heavy_type   = params_fopr_heavy.get_string("fermion_type");
    string gmset_type        = params_fopr_heavy.get_string("gamma_matrix_type");
    string source_heavy_type = params_source_heavy.get_string("source_type");

    const Bridge::VerboseLevel vl = vout.set_verbose_level(vlevel);

    //- print input parameters
    vout.general(vl, "Measurement parameters\n");
    vout.general(vl, "  gconf_status    = %s\n", gconf_status.c_str());
    vout.general(vl, "  gconf_type_read = %s\n", gconf_type_read.c_str());
    vout.general(vl, "  readfile        = %s\n", readfile.c_str());
    vout.general(vl, "  rand_type    = %s\n", rand_type.c_str());
    vout.general(vl, "  seed         = %lu\n", seed);
    vout.general(vl, "  vlevel       = %s\n", vlevel.c_str());
    vout.general(vl, "  gfix_type    = %s\n", gfix_type.c_str());
    vout.general(vl, "  gmset_type   = %s\n", gmset_type.c_str());
    vout.general(vl, "  source_heavy_type  = %s\n", source_heavy_type.c_str());

    //- input parameter check
    int err = 0;
    err += ParameterCheck::non_NULL(gconf_status);

    if (err) {
      vout.crucial(vl, "Error at %s: input parameters have not been set\n",
                   test_name.c_str());
      exit(EXIT_FAILURE);
    }

    RandomNumberManager::initialize(rand_type, seed);


    // ####  Set up a gauge configuration  ####
    unique_ptr<Field_G> U(new Field_G(Nvol, Ndim));

    if (gconf_status == "Continue") {
      GaugeConfig(gconf_type_read).read(*U, readfile);
    } else if (gconf_status == "Cold_start") {
      GaugeConfig("Unit").read(*U);
    } else if (gconf_status == "Hot_start") {
      GaugeConfig("Random").read(*U);
    } else {
      vout.crucial(vl, "Error at %s: unsupported gconf status \"%s\"\n",
                   test_name.c_str(), gconf_status.c_str());
      exit(EXIT_FAILURE);
    }

    // ####  Gauge fixing  ####
    {
      unique_ptr<Field_G>           Ufix(new Field_G(Nvol, Ndim));
      const unique_ptr<GaugeFixing> gfix(GaugeFixing::New(gfix_type));
      gfix->set_parameters(params_gfix);
      // gfix->fix(*Ufix, *U);
      // copy(*U, *Ufix);
    }

    // object setup for light quark
    unique_ptr<GammaMatrixSet> gmset(GammaMatrixSet::New(gmset_type));

    // heavy quark
    unique_ptr<Fopr> fopr_heavy(Fopr::New(fopr_heavy_type, params_fopr_heavy));
    fopr_heavy->set_config(U.get());

    const unique_ptr<Source> source_heavy(Source::New(source_heavy_type));
    source_heavy->set_parameters(params_source_heavy);


    // meason correlator
    Corr2pt_4spinor corr(gmset.get());
    corr.set_parameters(params_all.lookup("Corr2pt_4spinor"));

    const unique_ptr<Timer> timer(new Timer(test_name));


    // ####  Execution main part  ####
    timer->start();

    // heavy quark propagator
    vout.general(vl, "Solving heavy quark propagator:\n");

    std::vector<Field_F> sq_heavy(Nc * Nd);
    for (int icd = 0; icd < Nc * Nd; ++icd) {
      sq_heavy[icd].set(0.0);
    }

    fopr_heavy->set_mode("Evolve");
    for (int id = 0; id < Nd / 2; ++id) {
      for (int ic = 0; ic < Nc; ++ic) {
        int icd = ic + Nc * id;

        Field_F b;  // b.set(0.0);
        source_heavy->set(b, icd);

        int Nconv;
#pragma omp parallel
        {
          fopr_heavy->mult(sq_heavy[icd], b);
        }
      }
    }

    //- heavy-heavy meson correlators
    vout.general(vl, "\n");
    vout.general(vl, "heavy-heavy meson correlator:\n");
    const double result = corr.meson_all(sq_heavy, sq_heavy);

    timer->report();

    RandomNumberManager::finalize();


    if (do_check) {
      return Test::verify(result, expected_result);
    } else {
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }
    return EXIT_SUCCESS;
  }
}  //namespace Test_Spectrum_Nonrelativistic

//============================================================END=====
