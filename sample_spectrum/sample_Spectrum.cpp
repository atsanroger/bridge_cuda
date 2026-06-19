/*!
        @file    $Id: sample_Spectrum.cpp #$

        @brief   An example code to measure the spectrum

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-04-17 13:21:11 #$

        @version $LastChangedRevision: 2514 $
*/

#include "bridge.h"
using Bridge::vout;


namespace {
  const std::string test_name = "sample_Spectrum";
  const std::string parameter_file_default = "sample_Spectrum_Wilson_Hadron2ptFunction.yaml";
}

//====================================================================
//! An example code to use gauge field.

/*!
    This example is taken from test_Spectrum_Wilson_2ptFunction.cpp
                                      [25 Sep 2018 T.Aoyama]
    updated for version 2.0  [3 Apr 2023 I.Kanamori]
 */


//====================================================================
void usage(const char *program_name)
{
  vout.general("\n\nusage: %s parameter_file\n\n", program_name);
}


  //====================================================================
  int hadron_2ptFunction(const Parameters& params_all)
  {
    // ####  parameter setup  ####
    const int Nc   = CommonParameters::Nc();
    const int Nd   = CommonParameters::Nd();
    const int Ndim = CommonParameters::Ndim();
    const int Nvol = CommonParameters::Nvol();

    const Parameters params_spectrum = params_all.lookup("Spectrum");
    const Parameters params_gfix     = params_all.lookup("GaugeFixing");
    const Parameters params_fopr     = params_all.lookup("Fopr");
    const Parameters params_solver   = params_all.lookup("Solver");
    const Parameters params_source   = params_all.lookup("Source");

    const string        str_gconf_status = params_spectrum.get_string("gauge_config_status");
    const string        str_gconf_read   = params_spectrum.get_string("gauge_config_type_input");
    const string        readfile         = params_spectrum.get_string("config_filename_input");
    const string        str_rand_type    = params_spectrum.get_string("random_number_type");
    const unsigned long seed             = params_spectrum.get_unsigned_long("seed_for_random_number");
    const string        str_vlevel       = params_spectrum.get_string("verbose_level");

    //    const bool   do_check        = params_spectrum.is_set("expected_result");
    //    const double expected_result = do_check ? params_spectrum.get_double("expected_result") : 0.0;

    const string str_gfix_type   = params_gfix.get_string("gauge_fixing_type");
    const string str_fopr_type   = params_fopr.get_string("fermion_type");
    const string str_gmset_type  = params_fopr.get_string("gamma_matrix_type");
    const string str_solver_type = params_solver.get_string("solver_type");
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
    vout.general(vl, "  solver_type  = %s\n", str_solver_type.c_str());
    vout.general(vl, "  source_type  = %s\n", str_source_type.c_str());

    //- input parameter check
    int err = 0;
    err += ParameterCheck::non_NULL(str_gconf_status);

    if (err) {
      vout.crucial(vl, "Error at %s: input parameters have not been set\n", test_name.c_str());
      exit(EXIT_FAILURE);
    }

    if (str_solver_type == "CG") {
      vout.crucial(vl, "Error at %s: CG can not be adopted. Use CGNE,CGNR, instead.\n", test_name.c_str());
      exit(EXIT_FAILURE);
    }

    if ((str_gfix_type == "Coulomb") || (str_gfix_type == "Landau")) {
      if (CommonParameters::Nc() != 3) {
        vout.crucial(vl, "check skipped: Nc = 3 is needed, but Nc = %d.\n\n", CommonParameters::Nc());
        //return EXIT_SKIP;
        exit(EXIT_FAILURE);
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


    // ####  object setup  #####
    unique_ptr<GammaMatrixSet> gmset(GammaMatrixSet::New(str_gmset_type));

    unique_ptr<Fopr> fopr(Fopr::New(str_fopr_type, params_fopr));
    fopr->set_config(&U);

    unique_ptr<Solver> solver(Solver::New(str_solver_type, fopr.get(), params_solver));

    unique_ptr<Fprop> fprop_lex(new Fprop_Standard_lex(solver.get()));

    unique_ptr<Source> source(Source::New(str_source_type, params_source));

    Corr2pt_4spinor corr(gmset.get(), params_all.lookup("Corr2pt_4spinor"));

    Timer timer(test_name);


    // ####  Execution main part  ####
    timer.start();

    std::vector<Field_F> sq(Nc * Nd);
    for (int i_cd = 0; i_cd < Nc * Nd; ++i_cd) {
      sq[i_cd].set(0.0);
    }

    vout.general(vl, "\n");
    vout.general(vl, "Solving quark propagator:\n");
    vout.general(vl, "  color spin   Nconv      diff           diff2\n");

    for (int ispin = 0; ispin < Nd; ++ispin) {
      for (int icolor = 0; icolor < Nc; ++icolor) {
        int i_cd = icolor + Nc * ispin;

        Field_F b;  // b.set(0.0);
        source->set(b, i_cd);

        int    Nconv;
        double diff;
        fprop_lex->invert_D(sq[i_cd], b, Nconv, diff);

        Field_F y(b);
        fopr->set_mode("D");
        fopr->mult(y, sq[i_cd]); // y  = fopr->mult(sq[i_cd]);
        axpy(y, -1.0, b);        // y -= b;
        double diff2 = y.norm2() / b.norm2();

        vout.general(vl, "   %2d   %2d   %6d   %12.4e   %12.4e\n",
                     icolor, ispin, Nconv, diff, diff2);
      }
    }

    //- meson correlators
    vout.general(vl, "\n");
    vout.general(vl, "2-point correlator:\n");

    const double result = corr.meson_all(sq, sq);

    timer.report();

    RandomNumberManager::finalize();


    //    if (do_check) {
    //      return Test::verify(result, expected_result);
    //    } else {
    //      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
    //      return EXIT_SKIP;
    //    }
    return EXIT_SUCCESS;
  }


//====================================================================
int main(int argc, char *argv[])
{
  bridge_initialize(&argc, &argv);

  if (argc < 2) {
    usage(argv[0]);

    bridge_finalize();
    return EXIT_FAILURE;
  }

  std::string parameter_file = std::string(argv[1]);

  Parameters params = ParameterManager::read(parameter_file);
  bridge_setup(params.lookup("Main"));

  Timer timer("Main");
  timer.start();

  hadron_2ptFunction(params);

  timer.stop();
  timer.report();

  bridge_finalize();

  return EXIT_SUCCESS;
}
