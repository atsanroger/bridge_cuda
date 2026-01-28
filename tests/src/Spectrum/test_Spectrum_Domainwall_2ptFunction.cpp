/*!
        @file    test_Spectrum_Domainwall_2ptFunction.cpp

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$

        @version $LastChangedRevision: 2492 $
*/

#include "test.h"

#include "IO/gaugeConfig.h"

#include "Measurements/Fermion/corr2pt_4spinor.h"
#include "Measurements/Fermion/fprop_Standard_lex.h"
#include "Measurements/Fermion/source.h"
#include "Measurements/Gauge/gaugeFixing.h"
#include "Fopr/fopr_Domainwall.h"
#include "Tools/gammaMatrixSet.h"
#include "Tools/randomNumberManager.h"

//====================================================================
//! Test of spectroscopy for the domain-wall fermion operator.

/*!
    This class tests the domain-wall fermion propagator and
    construct hadron correlator.
                                          [12 Dec 2012 H.Matsufuru]
    (Coding history will be recovered from trac.)
    Implement YAML.                       [14 Nov 2012 Y.Namekawa]
    Implement Selectors.                  [03 Mar 2013 Y.Namekawa]
    (Selectors are replaced with factories by Aoyama-san)
    Introduce unique_ptr to avoid memory leaks.
                                          [21 Mar 2015 Y.Namekawa]
    Add Nc check for USE_GROUP_SU_N.      [31 May 2021 Y.Namekawa]

    New convention:
     -(P_ 0 ... 0 P_+) (D_5d)^{-1} D_- (P_+ 0 ... 0 P-)^T
     with  D_- = [ c_s Dwilson(-M0) - 1 ]
                                          [29 Dec 2022 I.Kanamori]

 */

namespace Test_Spectrum_Domainwall {
  const std::string test_name = "Spectrum.Domainwall.Hadron2ptFunction";

  //- test-private parameters
  namespace {
    const std::string filename_input = "test_Spectrum_Domainwall_Hadron2ptFunction.yaml";
  }

  //- prototype declaration
  int hadron_2ptFunction(void);

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
#if defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
    static const bool is_registered = TestManager::RegisterTest(
      test_name,
      hadron_2ptFunction
      );
#endif
  }
#endif

  //====================================================================
  int hadron_2ptFunction(void)
  {
    // ####  parameter setup  ####
    const int Nc   = CommonParameters::Nc();
    const int Nd   = CommonParameters::Nd();
    const int Nvol = CommonParameters::Nvol();
    const int Ndim = CommonParameters::Ndim();

    const Parameters params_all = ParameterManager::read(filename_input);

    const Parameters params_test   = params_all.lookup("Test_Spectrum_Domainwall");
    const Parameters params_gfix   = params_all.lookup("GaugeFixing");
    const Parameters params_dw     = params_all.lookup("Fopr_Domainwall");
    const Parameters params_solver = params_all.lookup("Solver");
    const Parameters params_source = params_all.lookup("Source");

    const string        str_gconf_status = params_test.get_string("gauge_config_status");
    const string        str_gconf_read   = params_test.get_string("gauge_config_type_input");
    const string        readfile         = params_test.get_string("config_filename_input");
    const string        str_rand_type    = params_test.get_string("random_number_type");
    const unsigned long seed             = params_test.get_unsigned_long("seed_for_random_number");
    const string        str_vlevel       = params_test.get_string("verbose_level");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

    const string str_gfix_type   = params_gfix.get_string("gauge_fixing_type");
    const string str_gmset_type  = params_dw.get_string("gamma_matrix_type");
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


    // ####  object setup  #####
    unique_ptr<GammaMatrixSet> gmset(GammaMatrixSet::New(str_gmset_type));

    // Fopr *fopr_w = Fopr::New("Wilson", str_gmset_type);
    // Fopr *fopr_dw = Fopr::New("Domainwall", fopr_w);
    // fopr_dw->set_parameters(params_dw);
    // fopr_dw->set_config(&U);

    // Fopr_Wilson *fopr_w = new Fopr_Wilson(str_gmset_type);
    // Fopr_Domainwall *fopr_dw = new Fopr_Domainwall(fopr_w, params_dw);
    // fopr_dw->set_config(&U);

    unique_ptr<Fopr_Domainwall> fopr_dw(new Fopr_Domainwall(params_dw));
    fopr_dw->set_mode("D");
    fopr_dw->set_config(&U);

    // kernel operator for 4d <--> 5d conversion
    Parameters params_dw_tmp;
    fopr_dw->get_parameters(params_dw_tmp);  // params_dw may not have coefficeint_c
    const string kernel_type = params_dw_tmp.get_string("kernel_type");
    const double M0          = params_dw_tmp.get_double("domain_wall_height");
    const double coeff_c     = params_dw_tmp.get_double("coefficient_c");
    const int    Ns          = params_dw_tmp.get_int("extent_of_5th_dimension");

    Parameters params_kernel   = params_dw;
    double     kappa           = 1.0 / (8.0 - 2.0 * M0);
    double     one_over_2kappa = 4.0 - M0;
    params_kernel.set_double("hopping_parameter", kappa);

    unique_ptr<Fopr> foprw(Fopr::New(kernel_type, params_kernel));
    foprw->set_mode("D");
    foprw->set_config(&U);

    unique_ptr<Solver> solver(Solver::New(str_solver_type, fopr_dw.get(), params_solver));

    //    unique_ptr<Fprop_Domainwall_4d> fprop(new Fprop_Domainwall_4d(fopr_dw.get(), solver.get()));
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

    Field_F b;
    Field_F vtmp_p, vtmp_m, vtmp;

    Field b5(b.nin(), b.nvol(), Ns);
    Field y5(b.nin(), b.nvol(), Ns);
    Field x5(b.nin(), b.nvol(), Ns);

    vout.general(vl, "\n");
    vout.general(vl, "Solving quark propagator:\n");
    vout.general(vl, "  color spin   Nconv      diff           diff2\n");

    for (int ispin = 0; ispin < Nd; ++ispin) {
      for (int icolor = 0; icolor < Nc; ++icolor) {
        int i_cd = icolor + Nc * ispin;

        source->set(b, i_cd);

        // set 5d source as
        //  - D_- ( P_+ 0 ....0 P_-)^T b4
#pragma omp parallel
        {
          // build 5dim vector
          b5.set(0.0);

          // s5=0
          fopr_dw->mult_chproj_4d(vtmp_p, b, +1);
          fopr_dw->mult_chproj_4d(vtmp_m, b, -1);

          foprw->mult(vtmp, vtmp_p);
          axpy(vtmp_p, -one_over_2kappa * coeff_c, vtmp);  // (-cD+1) (1+gm5)b/2

          foprw->mult(vtmp, vtmp_m);
          axpy(vtmp_m, -one_over_2kappa * coeff_c, vtmp);  // (-cD+1) (1-gm5)b/2

          copy(b5, 0, vtmp_p, 0);
          copy(b5, Ns - 1, vtmp_m, 0);
        } // omp parallel
        int    Nconv;
        double diff;
        fprop_lex->invert_D(x5, b5, Nconv, diff);
#pragma omp parallel
        {
          fopr_dw->set_mode("D");
          fopr_dw->mult(y5, x5);
          axpy(y5, -1.0, b5);
          double diff2 = y5.norm();

          vout.general(vl, "   %2d   %2d   %6d   %12.4e   %12.4e\n",
                       icolor, ispin, Nconv, diff, diff2);


          // convert to 4D propagator
          copy(vtmp, 0, x5, 0);
          fopr_dw->mult_chproj_4d(vtmp_m, vtmp, -1);
          copy(vtmp, 0, x5, Ns - 1);
          fopr_dw->mult_chproj_4d(vtmp_p, vtmp, +1);

          copy(sq[i_cd], vtmp_m);
          axpy(sq[i_cd], 1.0, vtmp_p);
        } // omp parallel
      }
    }

    //- hadron correlators
    vout.general(vl, "\n");
    vout.general(vl, "2-point correlator:\n");

    const double result = corr.meson_all(sq, sq);

    timer.report();

    RandomNumberManager::finalize();


    if (do_check) {
      return Test::verify(result, expected_result);
    } else {
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }
  }
} // namespace Test_Spectrum_Domainwall
