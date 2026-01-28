/*!
        @file    test_SF_fAfP_Boundary_Meson_2ptFunction.cpp

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$

        @version $LastChangedRevision: 2492 $
*/

#include "test.h"

#include "Fopr/fopr_Clover_SF.h"
#include "Fopr/fopr_Smeared.h"

#include "IO/gaugeConfig.h"

#include "Measurements/Fermion/corr2pt_Wilson_SF.h"
#include "Measurements/Fermion/fprop_Standard_lex.h"
#include "Measurements/Fermion/source_Wall_SF.h"

#include "Tools/gammaMatrixSet.h"
#include "Tools/randomNumberManager.h"

//====================================================================

/*!
    Calculate fA and fP for PCAC mass with SF BC.
    <ul>
      <li>The Dirac representaion of the gamma matrix is adopted.
      <li>Specified for
      <ul>
        <li>Number of smearing Nsmear=1 with parameter rho=0.1.
      </ul>
    </ul>
    (Coding history will be recovered from trac.)
    Implement YAML.                       [14 Nov 2012 Y.Namekawa]
    Implement Fprop and selectors.        [03 Mar 2013 Y.Namekawa]
    (Selectors are replaced with factories by Aoyama-san)
    Introduce unique_ptr to avoid memory leaks.
                                          [21 Mar 2015 Y.Namekawa]
    Add Nc check for USE_GROUP_SU_N.      [31 May 2021 Y.Namekawa]
*/

namespace Test_SF_fAfP {
  const std::string test_name = "SF_fAfP.Boundary_Meson_2ptFunction";

  //- test-private parameters
  namespace {
    const std::string filename_input = "test_SF_fAfP_Boundary_Meson_2ptFunction.yaml";
  }

  //- prototype declaration
  int boundary_meson_2ptFunction(void);

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
#if defined(USE_OPENMP) || defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
    static const bool is_registered = TestManager::RegisterTest(
      test_name,
      boundary_meson_2ptFunction
      );
#endif
  }
#endif

  //====================================================================
  int boundary_meson_2ptFunction(void)
  {
    // #####  parameter setup  #####
    const int Nc   = CommonParameters::Nc();
    const int Nd   = CommonParameters::Nd();
    const int Ndim = CommonParameters::Ndim();
    const int Nvol = CommonParameters::Nvol();

    const Parameters params_all = ParameterManager::read(filename_input);

    const Parameters params_test     = params_all.lookup("Test_SF_fAfP");
    const Parameters params_clover   = params_all.lookup("Fopr_Clover_SF");
    const Parameters params_proj     = params_all.lookup("Projection");
    const Parameters params_smear    = params_all.lookup("Smear_SF");
    const Parameters params_dr_smear = params_all.lookup("Director_Smear");
    const Parameters params_solver   = params_all.lookup("Solver");
    const Parameters params_source   = params_all.lookup("Source_Wall_SF");

    const string        str_gconf_status = params_test.get_string("gauge_config_status");
    const string        str_gconf_read   = params_test.get_string("gauge_config_type_input");
    const string        readfile         = params_test.get_string("config_filename_input");
    const string        str_rand_type    = params_test.get_string("random_number_type");
    const unsigned long seed             = params_test.get_unsigned_long("seed_for_random_number");
    const string        str_vlevel       = params_test.get_string("verbose_level");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

    const string str_gmset_type  = params_clover.get_string("gamma_matrix_type");
    const string str_proj_type   = params_proj.get_string("projection_type");
    const string str_smear_type  = params_smear.get_string("smear_type");
    const int    Nsmear          = params_dr_smear.get_int("number_of_smearing");
    const string str_solver_type = params_solver.get_string("solver_type");

    const Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

    //- print input parameters
    vout.general(vl, "  gconf_status = %s\n", str_gconf_status.c_str());
    vout.general(vl, "  gconf_read   = %s\n", str_gconf_read.c_str());
    vout.general(vl, "  readfile     = %s\n", readfile.c_str());
    vout.general(vl, "  rand_type    = %s\n", str_rand_type.c_str());
    vout.general(vl, "  seed         = %lu\n", seed);
    vout.general(vl, "  vlevel       = %s\n", str_vlevel.c_str());
    vout.general(vl, "  gmset_type   = %s\n", str_gmset_type.c_str());
    vout.general(vl, "  proj_type    = %s\n", str_proj_type.c_str());
    vout.general(vl, "  smear_type   = %s\n", str_smear_type.c_str());
    vout.general(vl, "  solver_type  = %s\n", str_solver_type.c_str());
    vout.general(vl, "\n");

    //- input parameter check
    int err = 0;
    err += ParameterCheck::non_NULL(str_gconf_status);

    if (err) {
      vout.crucial(vl, "Error at %s: input parameters have not been set\n", test_name.c_str());
      exit(EXIT_FAILURE);
    }

    if ((Nsmear > 0) && (str_proj_type == "Stout_SU3")) {
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

    if (Nsmear > 0) {
      unique_ptr<Projection>     proj(Projection::New(str_proj_type, params_all.lookup("Projection")));
      unique_ptr<Smear>          smear(Smear::New(str_smear_type, proj.get(), params_all.lookup("Smear_SF")));
      unique_ptr<Director_Smear> dr_smear(new Director_Smear(smear.get(), params_all.lookup("Director_Smear")));
      dr_smear->set_config(&U);

      const Field_G *Usmear = (Field_G *)dr_smear->getptr_smearedConfig(Nsmear);
      copy(U, *Usmear);
    }

    // ####  object setup  #####
    unique_ptr<GammaMatrixSet> gmset(GammaMatrixSet::New(str_gmset_type));

    //- NB. Chiral repr has not been implemented for SF, yet.
    unique_ptr<Fopr> fopr(new Fopr_Clover_SF(params_clover));
    fopr->set_config(&U);

    unique_ptr<Solver> solver(Solver::New(str_solver_type, fopr.get(), params_solver));

    unique_ptr<Fprop> fprop_lex(new Fprop_Standard_lex(solver.get()));

    unique_ptr<Source_Wall_SF> source(new Source_Wall_SF(params_source));
    source->set_config(&U);

    Timer timer(test_name);


    // ####  Execution main part  ####
    timer.start();

    std::vector<Field_F> H(Nc * Nd);
    std::vector<Field_F> Hpr(Nc * Nd);

    vout.general(vl, "\n");
    vout.general(vl, "Solving quark propagator:\n");
    vout.general(vl, "  color spin   Nconv      diff           diff2\n");

    for (int icolor = 0; icolor < Nc; ++icolor) {
      for (int ispin = 0; ispin < Nd / 2; ++ispin) {
        Field_F b;
        source->set_t0(b, icolor, ispin);

        int idx = icolor + Nc * ispin;

        Field_F xq;
        int     Nconv;
        double  diff;
        fprop_lex->invert_D(xq, b, Nconv, diff);

        Field_F y(b);
        fopr->set_mode("D");

        fopr->mult(y, xq);  // y  = fopr_w->mult(xq);
        axpy(y, -1.0, b);   // y -= b;
        double diff2 = y.norm2();

        vout.general(vl, "   %2d   %2d   %6d   %12.4e   %12.4e\n",
                     icolor, ispin, Nconv, diff, diff2);

        H[idx] = xq;
        xq.set(0.0);
        Hpr[idx] = xq;
      }

      for (int ispin = Nd / 2; ispin < Nd; ++ispin) {
        Field_F b;
        source->set_tT(b, icolor, ispin);

        int idx = icolor + Nc * ispin;

        Field_F xq;
        int     Nconv;
        double  diff;
        fprop_lex->invert_D(xq, b, Nconv, diff);

        Field_F y(b);
        fopr->set_mode("D");

        fopr->mult(y, xq);  // y  = fopr_w->mult(xq);
        axpy(y, -1.0, b);   // y -= b;
        double diff2 = y.norm2();

        vout.general(vl, "   %2d   %2d   %6d   %12.4e   %12.4e\n",
                     icolor, ispin, Nconv, diff, diff2);

        Hpr[idx] = xq;
        xq.set(0.0);
        H[idx] = xq;
      }
    }

#ifdef DEBUG
    Index_lex index;
    for (int t = 0; t < 2; ++t) {
      int site = index.site(0, 0, 0, t);
      for (int c1 = 0; c1 < Nc; ++c1) {
        for (int c0 = 0; c0 < Nc; ++c0) {
          for (int s1 = 0; s1 < Nd; ++s1) {
            for (int s0 = 0; s0 < Nd; ++s0) {
              vout.general(vl, "H[%d,%d,%d,%d,%d]=%lf %lf\n",
                           t, c1, c0, s1, s0,
                           H[c0 + Nc * s0].cmp_r(c1, s1, site),
                           H[c0 + Nc * s0].cmp_i(c1, s1, site));
            }
          }
        }
      }
    }
    for (int t = 0; t < 2; ++t) {
      int site = index.site(0, 0, 0, t);
      for (int c1 = 0; c1 < Nc; ++c1) {
        for (int c0 = 0; c0 < Nc; ++c0) {
          for (int s1 = 0; s1 < Nd; ++s1) {
            for (int s0 = 0; s0 < Nd; ++s0) {
              vout.general(vl, "Hpr[%d,%d,%d,%d,%d]=%lf %lf\n",
                           t, c1, c0, s1, s0,
                           Hpr[c0 + Nc * s0].cmp_r(c1, s1, site),
                           Hpr[c0 + Nc * s0].cmp_i(c1, s1, site));
            }
          }
        }
      }
    }
#endif

    //- meson correlators
    vout.general(vl, "\n");
    vout.general(vl, "boundary 2-point correlator with SF BC:\n");

    Corr2pt_Wilson_SF corr(gmset.get());
    const double      result = corr.fAfP(H, Hpr);

    timer.report();

    RandomNumberManager::finalize();


    if (do_check) {
      return Test::verify(result, expected_result);
    } else {
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }
  }
} // namespace Test_SF_fAfP
