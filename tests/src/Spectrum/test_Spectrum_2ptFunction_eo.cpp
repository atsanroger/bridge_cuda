/*!
        @file    test_Spectrum_2ptFunction_eo.cpp

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$

        @version $LastChangedRevision: 2492 $
*/

#include "test.h"

#include "Field/shiftField_lex.h"

#include "Fopr/fopr_Smeared.h"

#include "IO/gaugeConfig.h"

#include "Measurements/Fermion/corr2pt_4spinor.h"
#include "Measurements/Fermion/fprop_Standard_eo.h"
#include "Measurements/Fermion/source.h"
#include "Measurements/Gauge/gaugeFixing.h"

#include "Tools/filename.h"
#include "Tools/gammaMatrixSet.h"
#include "Tools/randomNumberManager.h"

//====================================================================
//! Test of spectroscopy with even-odd preconditioning.

/*!
    This class tests hadron spectroscopy on a given gauge configuration.
                                          [06 Jun 2014 Y.Namekawa]
    Introduce unique_ptr to avoid memory leaks.
                                          [21 Mar 2015 Y.Namekawa]
    Add shift_origin of gauge config.     [21 May 2017 Y.Namekawa]
    Modify to choose a solver type quark by quark.
                                          [23 Apr 2018 Y.Namekawa]
    Add Nc check for USE_GROUP_SU_N.      [31 May 2021 Y.Namekawa]
 */

namespace Test_Spectrum {
  const std::string test_name = "Spectrum.Hadron2ptFunction_eo";

  //- test-private parameters
  namespace {
    // const std::string filename_input  = "test_Spectrum_Clover_Hadron2ptFunction.yaml";
  }

  //- prototype declaration
  int hadron_2ptFunction_eo(const std::string&);

  //- hadron_2ptFunction for various fermions
  int hadron_2ptFunction_eo_Clover()
  {
    return hadron_2ptFunction_eo("test_Spectrum_Clover_Hadron2ptFunction.yaml");
  }


  int hadron_2ptFunction_eo_Wilson()
  {
    return hadron_2ptFunction_eo("test_Spectrum_Wilson_Hadron2ptFunction_2.yaml");
  }


#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
#if defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
    static const bool is_registered_Clover_eo = TestManager::RegisterTest(
      "Spectrum.Clover.Hadron2ptFunction_eo",
      hadron_2ptFunction_eo_Clover
      );

    //- NB. test_Spectrum_Wilson is implemented separately for beginners
    // static const bool is_registered_Wilson_eo = TestManager::RegisterTest(
    //   "Spectrum.Wilson.Hadron2ptFunction_eo",
    //   hadron_2ptFunction_eo_Wilson
    // );
#endif
  }
#endif

  //====================================================================
  int hadron_2ptFunction_eo(const std::string& filename_input)
  {
    // ####  parameter setup  ####
    const int Nc   = CommonParameters::Nc();
    const int Nd   = CommonParameters::Nd();
    const int Ndim = CommonParameters::Ndim();
    const int Nvol = CommonParameters::Nvol();

    const Parameters params_all = ParameterManager::read(filename_input);

    const Parameters params_test = params_all.lookup("Test_Spectrum");

    const int N_quark = params_test.get_int("number_of_valence_quarks");

    const std::string   str_gconf_status = params_test.get_string("gauge_config_status");
    const std::string   str_gconf_read   = params_test.get_string("gauge_config_type_input");
    const vector<int>   Nshift_origin    = params_test.get_int_vector("shift_origin");
    const std::string   readfile         = params_test.get_string("config_filename_input");
    const std::string   str_rand_type    = params_test.get_string("random_number_type");
    const unsigned long seed             = params_test.get_unsigned_long("seed_for_random_number");
    const std::string   str_vlevel       = params_test.get_string("verbose_level");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

    const std::string str_gfix_type  = params_all.lookup("GaugeFixing").get_string("gauge_fixing_type");
    const std::string str_proj_type  = params_all.lookup("Projection").get_string("projection_type");
    const std::string str_smear_type = params_all.lookup("Smear").get_string("smear_type");
    const int         Nsmear         = params_all.lookup("Director_Smear").get_int("number_of_smearing");

    const Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

    //- print input parameters
    vout.general(vl, "  gconf_status = %s\n", str_gconf_status.c_str());
    vout.general(vl, "  gconf_read   = %s\n", str_gconf_read.c_str());
    vout.general(vl, "  readfile     = %s\n", readfile.c_str());
    for (int mu = 0; mu < Ndim; ++mu) {
      vout.general(vl, "  shift_origin[%d] = %d\n", mu, Nshift_origin[mu]);
    }
    vout.general(vl, "  rand_type    = %s\n", str_rand_type.c_str());
    vout.general(vl, "  seed         = %lu\n", seed);
    vout.general(vl, "  vlevel       = %s\n", str_vlevel.c_str());
    vout.general(vl, "  gfix_type    = %s\n", str_gfix_type.c_str());
    vout.general(vl, "  proj_type    = %s\n", str_proj_type.c_str());
    vout.general(vl, "  smear_type   = %s\n", str_smear_type.c_str());

    vector<Parameters> params_quark;
    for (int iq = 0; iq < N_quark; ++iq) {
      std::string qlabel = Filename("Quark_{id}").format(iq + 1);
      params_quark.push_back(params_all.lookup(qlabel));
    }

    // NB. all str_gmset_type are supposed to be the same.
    const std::string str_gmset_type = params_quark[0].lookup("Fopr").get_string("gamma_matrix_type");
    vout.general(vl, "  gmset_type  = %s\n", str_gmset_type.c_str());

    for (int iq = 0; iq < N_quark; ++iq) {
      vout.general(vl, "  Quark_%d:\n", iq + 1);
      vout.general(vl, "    solver_type = %s\n",
                   params_quark[iq].lookup("Solver").get_string("solver_type").c_str());
      vout.general(vl, "    source_type = %s\n",
                   params_quark[iq].lookup("Source").get_string("source_type").c_str());

      std::string str_solver_type = params_quark[iq].lookup("Solver").get_string("solver_type");
      vout.general(vl, "    solver_type  = %s\n", str_solver_type.c_str());

      vector<int> source_position = params_quark[iq].lookup("Source").get_int_vector("source_position");
      for (int mu = 0; mu < Ndim; ++mu) {
        vout.general(vl, "    source_position[%d] = %d\n", mu, source_position[mu]);
      }
    }
    vout.general(vl, "\n");

    //- input parameter check
    int err = 0;
    err += ParameterCheck::non_NULL(str_gconf_status);

    if (err) {
      vout.crucial(vl, "Error at %s: input parameters have not been set\n", test_name.c_str());
      exit(EXIT_FAILURE);
    }

    for (int iq = 0; iq < N_quark; ++iq) {
      std::string str_solver_type = params_quark[iq].lookup("Solver").get_string("solver_type");

      if (str_solver_type == "CG") {
        vout.crucial(vl, "Error at %s: CG can not be adopted. Use CGNE,CGNR, instead.\n", test_name.c_str());
        exit(EXIT_FAILURE);
      }
    }

    if ((Nsmear > 0) && (str_proj_type == "Stout_SU3")) {
      if (CommonParameters::Nc() != 3) {
        vout.crucial(vl, "check skipped: Nc = 3 is needed, but Nc = %d.\n\n", CommonParameters::Nc());
        return EXIT_SKIP;
      }
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

    {
      ShiftField_lex shift;
      for (int i_dir = 0; i_dir < Ndim; ++i_dir) {
        for (int i_shift = 0; i_shift < Nshift_origin[i_dir]; ++i_shift) {
          Field_G Ushift(Nvol, Ndim);
          shift.backward(Ushift, U, i_dir);
          copy(U, Ushift);
        }
      }
    }

    if (Nsmear > 0) {
      unique_ptr<Projection>     proj(Projection::New(str_proj_type, params_all.lookup("Projection")));
      unique_ptr<Smear>          smear(Smear::New(str_smear_type, proj.get(), params_all.lookup("Smear")));
      unique_ptr<Director_Smear> dr_smear(new Director_Smear(smear.get(), params_all.lookup("Director_Smear")));
      dr_smear->set_config(&U);

      const Field_G *Usmear = (Field_G *)dr_smear->getptr_smearedConfig(Nsmear);
      copy(U, *Usmear);
    }

    // ####  Gauge fixing  ####
    {
      Field_G                 Ufix(Nvol, Ndim);
      unique_ptr<GaugeFixing> gfix(GaugeFixing::New(str_gfix_type, params_all.lookup("GaugeFixing")));

      gfix->fix(Ufix, U);
      copy(U, Ufix);
    }


    // ####  object setup  #####
    unique_ptr<GammaMatrixSet> gmset(GammaMatrixSet::New(str_gmset_type));

    struct QuarkType
    {
      unique_ptr<Fopr>   fopr;
      unique_ptr<Fopr>   fopr_eo;
      unique_ptr<Solver> solver;
      unique_ptr<Fprop>  fprop_eo;
      unique_ptr<Source> source;
    };
    // NB. Fopr is used only for check of diff2 below.

    std::vector<QuarkType> quark(N_quark);

    for (int iq = 0; iq < N_quark; ++iq) {
      const Parameters& params_fopr   = params_quark[iq].lookup("Fopr");
      const Parameters& params_solver = params_quark[iq].lookup("Solver");
      const Parameters& params_source = params_quark[iq].lookup("Source");

      quark[iq].fopr.reset(Fopr::New(params_fopr.get_string("fermion_type"),
                                     params_fopr));
      quark[iq].fopr->set_config(&U);

      quark[iq].fopr_eo.reset(Fopr::New(params_fopr.get_string("fermion_type") + "_eo",
                                        params_fopr));
      quark[iq].fopr_eo->set_config(&U);

      quark[iq].solver.reset(Solver::New(params_solver.get_string("solver_type"),
                                         quark[iq].fopr_eo.get(),
                                         params_solver));

      quark[iq].fprop_eo.reset(new Fprop_Standard_eo(quark[iq].solver.get()));

      quark[iq].source.reset(Source::New(params_source.get_string("source_type"),
                                         params_source));
    }

    Corr2pt_4spinor corr(gmset.get(), params_all.lookup("Corr2pt_4spinor"));

    Timer timer(test_name);


    // ####  Execution main part  ####
    timer.start();

    typedef std::vector<Field_F> PropagatorSet;
    std::vector<PropagatorSet> sq(N_quark);
    for (int iq = 0; iq < N_quark; ++iq) {
      sq[iq].resize(Nc * Nd);

      for (int i_cd = 0; i_cd < Nc * Nd; ++i_cd) {
        sq[iq][i_cd].set(0.0);
      }
    }

    for (int iq = 0; iq < N_quark; ++iq) {
      vout.general(vl, "Solving quark propagator, flavor = %d:\n", iq + 1);
      vout.general(vl, "  color spin   Nconv      diff           diff2\n");

      for (int ispin = 0; ispin < Nd; ++ispin) {
        for (int icolor = 0; icolor < Nc; ++icolor) {
          int i_cd = icolor + Nc * ispin;

          Field_F b;  // b.set(0.0);
          quark[iq].source->set(b, i_cd);

          int    Nconv;
          double diff;
          quark[iq].fprop_eo->invert_D(sq[iq][i_cd], b, Nconv, diff);

          Field_F y(b);
          quark[iq].fopr->set_mode("D");
          quark[iq].fopr->mult(y, sq[iq][i_cd]); // y  = fopr[iq]->mult(sq[iq][i_cd]);
          axpy(y, -1.0, b);                      // y -= b;
          double diff2 = y.norm2() / b.norm2();

          vout.general(vl, "   %2d   %2d   %6d   %12.4e   %12.4e\n",
                       icolor, ispin, Nconv, diff, diff2);
        }
      }

      vout.general(vl, "\n");
    }


    //- meson correlators
    std::vector<double> result(N_quark);

    vout.general(vl, "2-point correlator:\n");

    //- case(iq_1 == iq_2)
    for (int iq = 0; iq < N_quark; ++iq) {
      vout.general(vl, "Flavor combination = %d, %d\n", iq + 1, iq + 1);
      result[iq] = corr.meson_all(sq[iq], sq[iq]);
      vout.general(vl, "\n");
    }


    //- case(iq_1 < iq_2)
    for (int iq = 0; iq < N_quark; ++iq) {
      for (int jq = iq + 1; jq < N_quark; ++jq) {
        vout.general(vl, "Flavor combination = %d, %d\n", iq + 1, jq + 1);
        double result_2 = corr.meson_all(sq[iq], sq[jq]);
        vout.general(vl, "\n");
      }
    }


    timer.report();

    RandomNumberManager::finalize();


    if (do_check) {
      return Test::verify(result[0], expected_result);
    } else {
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }
  }
} // namespace Test_Spectrum
