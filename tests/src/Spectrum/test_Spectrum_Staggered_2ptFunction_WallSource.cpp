/*!
        @file    test_Spectrum_Staggered_2ptFunction_WallSource.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$

        @version $LastChangedRevision: 2492 $
*/

#include "test.h"

#include "IO/gaugeConfig.h"

#include "Measurements/Fermion/corr2pt_Staggered.h"
#include "Fopr/fopr_Staggered.h"
#include "Measurements/Fermion/fprop_Standard_lex.h"
#include "Measurements/Fermion/source_Staggered_Wall.h"
#include "Measurements/Gauge/gaugeFixing.h"

#include "Tools/randomNumberManager.h"

//====================================================================
//! Spectrum test of staggered fermion.

/*!
    Staggered fermion spectroscopy with lexical version of Fopr.
    For ver.2.0.                            [21 Nov 2021 H.Matsufuru]
 */

namespace Test_Spectrum_Staggered {
  const std::string test_name =
    "Spectrum.Staggered.Hadron2ptFunction.WallSource";

  //- test-private parameters
  namespace {
    const std::string filename_input =
      "test_Spectrum_Staggered_Hadron2ptFunction_WallSource.yaml";
  }

  //- prototype declaration
  int hadron_2ptFunction_wallSource(void);

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
#if defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
    static const bool is_registered = TestManager::RegisterTest(
      test_name,
      hadron_2ptFunction_wallSource
      );
#endif
  }
#endif

//====================================================================
  int hadron_2ptFunction_wallSource(void)
  {
    // ####  parameter setup  ####
    const int Nc   = CommonParameters::Nc();
    const int Nvol = CommonParameters::Nvol();
    const int Ndim = CommonParameters::Ndim();

    const Parameters params_all = ParameterManager::read(filename_input);

    const Parameters params_test   = params_all.lookup("Test_Spectrum_Staggered");
    const Parameters params_gfix   = params_all.lookup("GaugeFixing");
    const Parameters params_fopr   = params_all.lookup("Fopr");
    const Parameters params_solver = params_all.lookup("Solver");
    const Parameters params_source = params_all.lookup("Source_Wall");

    const string gconf_status = params_test.get_string("gauge_config_status");
    const string gconf_read   = params_test.get_string("gauge_config_type_input");
    const string readfile     = params_test.get_string("config_filename_input");
    const string rand_type    = params_test.get_string("random_number_type");
    const string str_vlevel   = params_test.get_string("verbose_level");

    const unsigned long seed = params_test.get_unsigned_long("seed_for_random_number");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

    const string str_gfix_type = params_gfix.get_string("gauge_fixing_type");

    const Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

    //- print input parameters
    vout.general(vl, "  gconf_status = %s\n", gconf_status.c_str());
    vout.general(vl, "  gconf_read   = %s\n", gconf_read.c_str());
    vout.general(vl, "  readfile     = %s\n", readfile.c_str());
    vout.general(vl, "  rand_type    = %s\n", rand_type.c_str());
    vout.general(vl, "  seed         = %lu\n", seed);
    vout.general(vl, "  vlevel       = %s\n", str_vlevel.c_str());
    vout.general(vl, "  gfix_type    = %s\n", str_gfix_type.c_str());

    //- input parameter check
    int err = 0;
    err += ParameterCheck::non_NULL(gconf_status);

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

    RandomNumberManager::initialize(rand_type, seed);


    // ####  Set up a gauge configuration  ####
    Field_G U(Nvol, Ndim);

    if (gconf_status == "Continue") {
      GaugeConfig(gconf_read).read(U, readfile);
    } else if (gconf_status == "Cold_start") {
      GaugeConfig("Unit").read(U);
    } else if (gconf_status == "Hot_start") {
      GaugeConfig("Random").read(U);
    } else {
      vout.crucial(vl, "Error at %s: unsupported gconf status \"%s\"\n",
                   test_name.c_str(), gconf_status.c_str());
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

    std::string      fopr_type = params_fopr.get_string("fermion_type");
    unique_ptr<Fopr> fopr(Fopr::New(fopr_type, params_fopr));
    fopr->set_config(&U);

    std::string        solver_type = params_solver.get_string("solver_type");
    unique_ptr<Solver> solver(Solver::New(solver_type,
                                          fopr.get(), params_solver));

    unique_ptr<Fprop_Standard_lex> fprop(
      new Fprop_Standard_lex(solver.get()));

    unique_ptr<Source_Staggered_Wall> source(
      new Source_Staggered_Wall(params_source));

    Timer timer(test_name);

    // ####  Execution main part  ####
    timer.start();

    std::vector<Field_F_1spinor> sq(Nc * 2);
    for (int icd = 0; icd < 2 * Nc; ++icd) {
      sq[icd].set(0.0);
    }

    vout.general(vl, "    ic   isrc    Nconv       diff\n");
    for (int isrc = 0; isrc < 2; ++isrc) {
      for (int ic = 0; ic < Nc; ++ic) {
        Field_F_1spinor src;
        source->set(src, ic, isrc);

        int    icd = ic + Nc * isrc;
        int    nconv;
        double diff;
        fprop->invert_D(sq[icd], src, nconv, diff);
        vout.general(vl, "  %4d  %4d  %8d  %12.4e\n", ic, isrc, nconv, diff);
      }
    }

    //- meson correlators
    vout.general(vl, "\n");
    vout.general(vl, "meson correlator:\n");

    Corr2pt_Staggered   corr2pt;
    std::vector<double> meson_corr;
    corr2pt.meson(meson_corr, sq, sq);

    vout.general(vl, "PS <-- PS correlator:\n");
    for (int t = 0; t < meson_corr.size(); ++t) {
      vout.general(vl, "  %4d %20.12e\n", t, meson_corr[t]);
    }

    const double result = meson_corr[0];

    timer.report();

    RandomNumberManager::finalize();


    if (do_check) {
      return Test::verify(result, expected_result);
    } else {
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }
  }
} // namespace Test_Spectrum_Staggered
