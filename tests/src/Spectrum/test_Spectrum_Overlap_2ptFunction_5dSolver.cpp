/*!
        @file    test_Spectrum_Overlap_2ptFunction_5dSolver.cpp

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate: 2013-04-27 12:28:50 #$

        @version $LastChangedRevision: 2359 $
*/

#include "test.h"

#include "Eigen/eigensolver_IRLanczos.h"

#include "IO/gaugeConfig.h"

#include "Measurements/Fermion/corr2pt_4spinor.h"
#include "Measurements/Fermion/fprop_Overlap_5d.h"
#include "Measurements/Fermion/source.h"
#include "Measurements/Gauge/gaugeFixing.h"

#include "Fopr/fopr_Wilson.h"

#include "Tools/gammaMatrixSet.h"
#include "Tools/randomNumberManager.h"

//====================================================================
//! Test class of overlap fermion operators with 5 dimensional solver.

/*!
    This class tests overlap fermion using 5d version of solver,
    with and without subtraction of low-lying modes of
    the base fermion operator. Whether subtraction is performed or not
    is controlled in the code explicitly by commenting out some lines.
    As the base fermion operator, Wilson fermion operator is used.
                                          [13 Apr 2012  H.Matsufuru]
    (Coding history will be recovered from trac.)
    Implement YAML.                       [14 Nov 2012 Y.Namekawa]
    Implement Selectors.                  [03 Mar 2013 Y.Namekawa]
    (Selectors are replaced with factories by Aoyama-san)
    Introduce unique_ptr to avoid memory leaks.
                                          [21 Mar 2015 Y.Namekawa]
    Add Nc check for USE_GROUP_SU_N.      [31 May 2021 Y.Namekawa]
*/

namespace Test_Spectrum_Overlap {
  const std::string test_name = "Spectrum.Overlap.Hadron2ptFunction_5dSolver";

  //- test-private parameters
  namespace {
    const std::string filename_input = "test_Spectrum_Overlap_Hadron2ptFunction_5dSolver.yaml";
  }

  //- prototype declaration
  int hadron_2ptFunction_5d_solver(void);
  void calc_lowmodes(Parameters& params_irlanczos,
                     int& Nsbt,
                     std::vector<double>& TDa,
                     std::vector<Field>& vk,
                     Fopr *fopr);

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
#if defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
    static const bool is_registered = TestManager::RegisterTest(
      test_name,
      hadron_2ptFunction_5d_solver
      );
#endif
  }
#endif

  //====================================================================
  int hadron_2ptFunction_5d_solver(void)
  {
    // ####  parameter setup  ####
    const int Nc   = CommonParameters::Nc();
    const int Nd   = CommonParameters::Nd();
    const int Nvol = CommonParameters::Nvol();
    const int Ndim = CommonParameters::Ndim();
    const int NinF = 2 * Nc * Nd;

    const Parameters params_all = ParameterManager::read(filename_input);

    const Parameters params_test    = params_all.lookup("Test_Spectrum_Overlap");
    const Parameters params_gfix    = params_all.lookup("GaugeFixing");
    const Parameters params_wilson  = params_all.lookup("Fopr_Wilson");
    const Parameters params_overlap = params_all.lookup("Fopr_Overlap_5d_solver");
    //- NB. params_irlanczos is set in the following calc_lowmodes
    Parameters       params_irlanczos = params_all.lookup("Eigensolver");
    const Parameters params_solver    = params_all.lookup("Solver");
    const Parameters params_source    = params_all.lookup("Source");

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

    if ( (str_gfix_type == "Coulomb") || (str_gfix_type == "Landau") ) {
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
      Field_G Ufix(Nvol, Ndim);
      unique_ptr<GaugeFixing> gfix(GaugeFixing::New(str_gfix_type, params_gfix));

      gfix->fix(Ufix, U);
      copy(U, Ufix);
    }


    // ####  object setup  ####
    unique_ptr<GammaMatrixSet> gmset(GammaMatrixSet::New(str_gmset_type));

    unique_ptr<Fopr_Wilson> fopr_w(new Fopr_Wilson(params_wilson));
    fopr_w->set_config(&U);

    int Nsbt = 0;

    std::vector<double> TDa;
    std::vector<Field>  vk;

    //- with low-mode subtraction
    calc_lowmodes(params_irlanczos, Nsbt, TDa, vk, fopr_w.get());

    unique_ptr<Fopr_Overlap> fopr_overlap(new Fopr_Overlap(fopr_w.get(), params_overlap));
    fopr_overlap->set_config(&U);
    fopr_overlap->set_lowmodes(Nsbt, &TDa, &vk);

    // TODO
    unique_ptr<Fprop_Overlap_5d> fprop_5d(new Fprop_Overlap_5d(fopr_overlap.get(), &U));
    fprop_5d->set_parameters(params_overlap, params_solver);
    fprop_5d->set_lowmodes(Nsbt, &TDa, &vk);

    unique_ptr<Source> source(Source::New(str_source_type, params_source));

    Corr2pt_4spinor corr(gmset.get(), params_all.lookup("Corr2pt_4spinor"));

    Timer timer(test_name);

    // ####  Execution main part  ####
    timer.start();

    std::vector<Field_F> sq(Nc * Nd);
    for (int icd = 0; icd < Nc * Nd; ++icd) {
      sq[icd].set(0.0);
    }

    Field b(NinF, Nvol, 1);

    vout.general(vl, "\n");
    vout.general(vl, "Solving quark propagator:\n");
    vout.general(vl, "  color spin   Nconv      diff           diff2\n");

    for (int ispin = 0; ispin < Nd; ++ispin) {
      for (int icolor = 0; icolor < Nc; ++icolor) {
        int icd = icolor + Nc * ispin;

        source->set(b, icd);

        fprop_5d->invert_D(sq[icd], b);
      }
    }

    //- meson correlators
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


  //====================================================================
  void calc_lowmodes(Parameters& params_irlanczos,
                     int& Nsbt,
                     std::vector<double>& TDa,
                     std::vector<Field>& vk,
                     Fopr *fopr)
  {
    // #####  parameter setup  #####
    const Bridge::VerboseLevel vl = CommonParameters::Vlevel();

    const string str_sortfield_type = params_irlanczos.get_string("eigensolver_mode");
    const int    Nk = params_irlanczos.get_int("number_of_wanted_eigenvectors");
    const int    Np = params_irlanczos.get_int("number_of_working_eigenvectors");

    Field_F   b;
    const int Nin  = b.nin();
    const int Nvol = b.nvol();
    const int Nex  = b.nex();

    const int Nmm = Nk + Np;

    TDa.resize(Nmm);
    vk.resize(Nmm);


    // ####  object setup  ####
    Eigensolver_IRLanczos eigen(fopr, params_irlanczos);


    // ####  Execution main part  ####
    for (int k = 0; k < Nmm; ++k) {
      vk[k].reset(Nin, Nvol, Nex);
    }
    int Nconv = -100;

    fopr->set_mode("H");
    eigen.solve(TDa, vk, Nsbt, Nconv, (Field)b);

    Field v;
    v.reset(Nin, Nvol, Nex);

    fopr->set_mode("H");
    for (int i = 0; i < Nsbt; ++i) {
      fopr->mult(v, vk[i]);
      axpy(v, -TDa[i], vk[i]);  // v -= TDa[i] * vk[i];
      double vv = v.norm2();    // double vv = v * v;

      vout.general(vl, "Eigenvalues: %4d %20.14f %20.15e\n", i, TDa[i], vv);
    }

  }
} // namespace Test_Spectrum_Overlap
