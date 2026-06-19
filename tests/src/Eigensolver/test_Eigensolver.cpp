/*!
        @file    test_Eigensolver.cpp

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate: 2013-01-22 13:51:53 #$

        @version $LastChangedRevision: 2492 $
*/

#include "test.h"

#include "Eigen/eigensolver_IRLanczos.h"

#include "Field/field_F.h"
#include "Fopr/fopr_Smeared.h"

#include "IO/gaugeConfig.h"

#include "Tools/randomNumberManager.h"

//====================================================================
//! Test of eigenvalue solver.

/*!
                                          [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    Implement YAML.                       [14 Nov 2012 Y.Namekawa]
    Implement Selectors.                  [03 Mar 2013 Y.Namekawa]
    (Selectors are replaced with factories by Aoyama-san)
    Introduce unique_ptr to avoid memory leaks.
                                          [21 Mar 2015 Y.Namekawa]
    Add Nc check for USE_GROUP_SU_N.      [31 May 2021 Y.Namekawa]
 */

namespace Test_Eigensolver {
  const std::string test_name = "Eigensolver.Test";

  //- test-private parameters
  namespace {
    const std::string filename_input = "test_Eigensolver.yaml";
  }

  //- prototype declaration
  int solve(void);

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
#if defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
    static const bool is_registered = TestManager::RegisterTest(
      test_name,
      solve
      );
#endif
  }
#endif

  //====================================================================
  //- Check of eigenvalue solver
  int solve(void)
  {
    // #####  parameter setup  #####
    const int Ndim = CommonParameters::Ndim();
    const int Nvol = CommonParameters::Nvol();

    const Parameters params_all = ParameterManager::read(filename_input);

    const Parameters params_test      = params_all.lookup("Test_Eigensolver");
    const Parameters params_fopr      = params_all.lookup("Fopr");
    const Parameters params_proj      = params_all.lookup("Projection");
    const Parameters params_smear     = params_all.lookup("Smear");
    const Parameters params_dr_smear  = params_all.lookup("Director_Smear");
    const Parameters params_irlanczos = params_all.lookup("Eigensolver");

    const string        str_gconf_status = params_test.get_string("gauge_config_status");
    const string        str_gconf_read   = params_test.get_string("gauge_config_type_input");
    const string        readfile         = params_test.get_string("config_filename_input");
    const string        str_rand_type    = params_test.get_string("random_number_type");
    const unsigned long seed             = params_test.get_unsigned_long("seed_for_random_number");
    const string        str_vlevel       = params_test.get_string("verbose_level");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

    const string str_fopr_type      = params_fopr.get_string("fermion_type");
    const string str_gmset_type     = params_fopr.get_string("gamma_matrix_type");
    const string str_proj_type      = params_proj.get_string("projection_type");
    const string str_smear_type     = params_smear.get_string("smear_type");
    const int    Nsmear             = params_dr_smear.get_int("number_of_smearing");
    const string str_sortfield_type = params_irlanczos.get_string("eigensolver_mode");
    const int    Nk = params_irlanczos.get_int("number_of_wanted_eigenvectors");
    const int    Np = params_irlanczos.get_int("number_of_working_eigenvectors");

    const Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

    //- print input parameters
    vout.general(vl, "  gconf_status   = %s\n", str_gconf_status.c_str());
    vout.general(vl, "  gconf_read     = %s\n", str_gconf_read.c_str());
    vout.general(vl, "  readfile       = %s\n", readfile.c_str());
    vout.general(vl, "  rand_type      = %s\n", str_rand_type.c_str());
    vout.general(vl, "  seed           = %lu\n", seed);
    vout.general(vl, "  vlevel         = %s\n", str_vlevel.c_str());
    vout.general(vl, "  gmset_type     = %s\n", str_gmset_type.c_str());
    vout.general(vl, "  proj_type      = %s\n", str_proj_type.c_str());
    vout.general(vl, "  smear_type     = %s\n", str_smear_type.c_str());
    vout.general(vl, "  sortfield_type = %s\n", str_sortfield_type.c_str());
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
      unique_ptr<Smear>          smear(Smear::New(str_smear_type, proj.get(), params_all.lookup("Smear")));
      unique_ptr<Director_Smear> dr_smear(new Director_Smear(smear.get(), params_all.lookup("Director_Smear")));
      dr_smear->set_config(&U);

      const Field_G *Usmear = (Field_G *)dr_smear->getptr_smearedConfig(Nsmear);
      copy(U, *Usmear);
    }


    // ####  object setup  #####
    unique_ptr<Fopr> fopr(Fopr::New(str_fopr_type, params_fopr));
    fopr->set_config(&U);
    fopr->set_mode("H");

    unique_ptr<Eigensolver> eigen(new Eigensolver_IRLanczos(fopr.get(), params_irlanczos));

    Timer timer(test_name);


    // ####  Execution main part  ####
    timer.start();

    Field_F   b2;
    const int NFin  = b2.nin();
    const int NFvol = b2.nvol();
    const int NFex  = b2.nex();

    const int           Nm = Nk + Np;
    std::vector<double> TDa(Nm);

    std::vector<Field> vk(Nm);
    for (int k = 0; k < Nm; ++k) {
      vk[k].reset(NFin, NFvol, NFex);
    }

    int Nsbt  = -1;
    int Nconv = -100;
    eigen->solve(TDa, vk, Nsbt, Nconv, (Field)b2);

    for (int i = 0; i < Nsbt + 1; ++i) {
      Field v(NFin, NFvol, NFex);

      fopr->mult(v, vk[i]);
      axpy(v, -TDa[i], vk[i]);  // v -= TDa[i] * vk[i];
      double vv = v.norm2();    // vv = v * v;

      vout.general(vl, "Eigenvalues: %4d %20.14f %20.15e \n", i, TDa[i], vv);
    }

    const double result = TDa[0];

    timer.report();

    RandomNumberManager::finalize();


    if (do_check) {
      return Test::verify(result, expected_result);
    } else {
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }
  }
} // namespace Test_Eigensolver
