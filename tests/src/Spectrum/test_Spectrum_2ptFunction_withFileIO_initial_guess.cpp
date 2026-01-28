/*!
        @file    test_Spectrum_2ptFunction_withFileIO_initial_guess.cpp

        @brief

        @author  Yusuke Namekawa  (namekawa)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$

        @version $LastChangedRevision: 2492 $
*/

#include "test.h"

#include "Field/shiftField_lex.h"

#include "Fopr/fopr_Smeared.h"

#include "IO/gaugeConfig.h"

#include "Measurements/Fermion/corr2pt_4spinor.h"
#include "Measurements/Fermion/fprop_Standard_lex.h"
#include "Measurements/Fermion/source.h"
#include "Measurements/Gauge/gaugeFixing.h"

#include "Tools/filename.h"
#include "Tools/file_utils.h"
#include "Tools/randomNumberManager.h"
#include "Tools/gammaMatrixSet.h"

//====================================================================
//! Test of Spectroscopy with fileIO.

/*!
    This test class obtains the quark propagator for clover fermion
    and calculates typical hadron correlators.
    The quantum numbers of hadrons are specified with gamma matrices
    (GammaMatrix class instance) whose set is defined in a subclass
    of GammaMatrixSet class.
                                             [12 Apr 2012 H.Matsufuru]
    (Coding history will be recovered from trac.)
    Implement Selectors.                     [02 Feb 2013 Y.Namekawa]
    Implement Smearing.                      [03 Mar 2013 Y.Namekawa]
    (Selectors are replaced with factories by Aoyama-san)
    Implement Number_of_valence_quarks.      [15 May 2013 Y.Namekawa]
    Introduce unique_ptr to avoid memory leaks.
                                             [21 Mar 2015 Y.Namekawa]
    Add shift_origin of gauge config.        [21 May 2017 Y.Namekawa]
    Modify to choose a solver type quark by quark.
                                             [23 Apr 2018 Y.Namekawa]
    Employ initial guess for solver.         [ 1 May 2018 Y.Namekawa]
    Add Nc check for USE_GROUP_SU_N.         [31 May 2021 Y.Namekawa]
 */

namespace Test_Spectrum {
  const std::string test_name = "Spectrum.Hadron2ptFunction_withFileIO_initial_guess";

  //- test-private parameters
  namespace {
    // const std::string filename_input  = "test_Spectrum_Clover_Hadron2ptFunction.yaml";
  }

  //- prototype declaration
  int hadron_2ptFunction_withFileIO_initial_guess(const std::string&);

  int calculate_quark_propagator_withFileIO_initial_guess(const std::string&);
  int calculate_hadron_correlator_withFileIO_initial_guess(const std::string&);

  //- hadron_2ptFunction for various fermions
  int hadron_2ptFunction_withFileIO_Clover_initial_guess()
  {
    return hadron_2ptFunction_withFileIO_initial_guess("test_Spectrum_Clover_Hadron2ptFunction_initial_guess.yaml");
  }


  // int hadron_2ptFunction_withFileIO_CloverGeneral()
  // {
  //   return hadron_2ptFunction_withFileIO("test_Spectrum_CloverGeneral_Hadron2ptFunction.yaml");
  // }


  // int hadron_2ptFunction_withFileIO_Wilson_TwistedMass()
  // {
  //   return hadron_2ptFunction_withFileIO("test_Spectrum_Wilson_TwistedMass_Hadron2ptFunction.yaml");
  // }


  // int hadron_2ptFunction_withFileIO_Wilson()
  // {
  //   return hadron_2ptFunction_withFileIO("test_Spectrum_Wilson_Hadron2ptFunction.yaml");
  // }


  // int hadron_2ptFunction_withFileIO_WilsonGeneral()
  // {
  //   return hadron_2ptFunction_withFileIO("test_Spectrum_WilsonGeneral_Hadron2ptFunction.yaml");
  // }


#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
#if defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
    static const bool is_registered_Clover = TestManager::RegisterTest(
      "Spectrum.Clover.Hadron2ptFunction_withFileIO_initial_guess",
      hadron_2ptFunction_withFileIO_Clover_initial_guess
      );
#endif
  }
#endif

  //====================================================================
  int hadron_2ptFunction_withFileIO_initial_guess(const std::string& filename_input)
  {
    int result = 0;

    Timer timer(test_name);


    // ####  Execution main part  ####
    timer.start();

    result += calculate_quark_propagator_withFileIO_initial_guess(filename_input);
    result += calculate_hadron_correlator_withFileIO_initial_guess(filename_input);

    timer.report();

    return result;
  }


  //====================================================================
  int calculate_quark_propagator_withFileIO_initial_guess(const std::string& filename_input)
  {
    // ####  parameter setup  ####
    const int Nc   = CommonParameters::Nc();
    const int Nd   = CommonParameters::Nd();
    const int Ndim = CommonParameters::Ndim();
    const int Nvol = CommonParameters::Nvol();

    const Parameters params_all = ParameterManager::read(filename_input);

    const Parameters params_test = params_all.lookup("Test_Spectrum");

    const int N_quark = params_test.get_int("number_of_valence_quarks");

    const std::string str_gconf_status = params_test.get_string("gauge_config_status");
    const std::string str_gconf_read   = params_test.get_string("gauge_config_type_input");
    const std::string readfile         = params_test.get_string("config_filename_input");
    const vector<int> Nshift_origin    = params_test.get_int_vector("shift_origin");
    const std::string str_vlevel       = params_test.get_string("verbose_level");

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
    std::string str_gmset_type = params_quark[0].lookup("Fopr").get_string("gamma_matrix_type");
    vout.general(vl, "  gmset_type  = %s\n", str_gmset_type.c_str());

    std::vector<std::string> savefile_base(N_quark);
    for (int iq = 0; iq < N_quark; ++iq) {
      savefile_base[iq] = params_quark[iq].get_string("temporary_filename_base");
    }

    for (int iq = 0; iq < N_quark; ++iq) {
      std::string str_solver_type = params_quark[iq].lookup("Solver").get_string("solver_type");
      std::string str_source_type = params_quark[iq].lookup("Source").get_string("source_type");
      vector<int> source_position = params_quark[iq].lookup("Source").get_int_vector("source_position");

      vout.general(vl, "  Quark_%d:\n", iq + 1);
      vout.general(vl, "    savefile_base[iq] = %s[%d]\n", savefile_base[iq].c_str(), iq);
      vout.general(vl, "    solver_type       = %s\n", str_solver_type.c_str());
      vout.general(vl, "    source_type       = %s\n", str_source_type.c_str());
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
        //- NB. EXIT_SKIP is set in calculate_hadron_correlator_withFileIO()
        // return EXIT_SKIP;
        return EXIT_SUCCESS;
      }
    }

    if ((str_gfix_type == "Coulomb") || (str_gfix_type == "Landau")) {
      if (CommonParameters::Nc() != 3) {
        vout.crucial(vl, "check skipped: Nc = 3 is needed, but Nc = %d.\n\n", CommonParameters::Nc());
        return EXIT_SKIP;
      }
    }


    RandomNumberManager::initialize("Mseries", 1234567UL);


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

    struct QuarkType
    {
      unique_ptr<Fopr>   fopr;
      unique_ptr<Solver> solver;
      unique_ptr<Fprop>  fprop;
      unique_ptr<Source> source;
      bool               use_init_guess;
    };

    std::vector<QuarkType> quark(N_quark);

    for (int iq = 0; iq < N_quark; ++iq) {
      const Parameters& params_fopr   = params_quark[iq].lookup("Fopr");
      const Parameters& params_solver = params_quark[iq].lookup("Solver");
      const Parameters& params_source = params_quark[iq].lookup("Source");

      quark[iq].fopr.reset(Fopr::New(params_fopr.get_string("fermion_type"),
                                     params_fopr));
      quark[iq].fopr->set_config(&U);

      quark[iq].solver.reset(Solver::New(params_solver.get_string("solver_type"),
                                         quark[iq].fopr.get(),
                                         params_solver));
      quark[iq].use_init_guess = params_solver.get_bool("use_initial_guess");

      quark[iq].fprop.reset(new Fprop_Standard_lex(quark[iq].solver.get()));

      quark[iq].source.reset(Source::New(params_source.get_string("source_type"),
                                         params_source));
    }
    vout.general(vl, "\n");

    unique_ptr<FieldIO> field_io(new FieldIO_Binary(IO_Format::Trivial));


    // ####  Execution main part  ####
    for (int iq = 0; iq < N_quark; ++iq) {
      vout.general(vl, "Solving quark propagator, flavor = %d:\n", iq + 1);
      vout.general(vl, "  color spin   Nconv      diff           diff2\n");

      for (int ispin = 0; ispin < Nd; ++ispin) {
        for (int icolor = 0; icolor < Nc; ++icolor) {
          int i_cd = icolor + Nc * ispin;

          Field_F b;  // b.set(0.0);
          quark[iq].source->set(b, i_cd);

          Field_F xq;
          int     Nconv;
          double  diff;

          std::string filename = FileUtils::generate_filename(
            "%s_c%d_s%d.dat",
            savefile_base[iq].c_str(), (icolor + 1), (ispin + 1));

          if (quark[iq].use_init_guess) {
            // check if file exists, and read.
            std::ifstream ifs(filename);
            if (ifs.is_open()) {
              field_io->read_file(xq, filename);
            }
          }

          quark[iq].fprop->invert_D(xq, b, Nconv, diff);

          Field_F y(b);
          quark[iq].fopr->set_mode("D");
          quark[iq].fopr->mult(y, xq); // y  = fopr[iq]->mult(xq);
          axpy(y, -1.0, b);            // y -= b;
          double diff2 = y.norm2() / b.norm2();

          vout.general(vl, "   %2d   %2d   %6d   %12.4e   %12.4e\n",
                       icolor, ispin, Nconv, diff, diff2);

          field_io->write_file(xq, filename);
        }
      }

      vout.general(vl, "\n");
    }

    return EXIT_SUCCESS;
  }


  //====================================================================
  int calculate_hadron_correlator_withFileIO_initial_guess(const std::string& filename_input)
  {
    // ####  parameter setup  ####
    const int Nc   = CommonParameters::Nc();
    const int Nd   = CommonParameters::Nd();
    const int Ndim = CommonParameters::Ndim();
    const int Nvol = CommonParameters::Nvol();

    const Parameters params_all = ParameterManager::read(filename_input);

    const Parameters params_test = params_all.lookup("Test_Spectrum");

    const int N_quark = params_test.get_int("number_of_valence_quarks");

    const std::string str_vlevel = params_test.get_string("verbose_level");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

    const std::string str_proj_type = params_all.lookup("Projection").get_string("projection_type");
    const int         Nsmear        = params_all.lookup("Director_Smear").get_int("number_of_smearing");

    std::vector<std::string> savefile_base(N_quark);

    for (int iq = 0; iq < N_quark; ++iq) {
      savefile_base[iq] = params_all
                          .lookup(Filename("Quark_{id}").format(iq + 1))
                          .get_string("temporary_filename_base");
    }

    const Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

    //- print input parameters
    vout.general(vl, "  vlevel      = %s\n", str_vlevel.c_str());

    // NB. all str_gmset_type are supposed to be the same.
    const std::string str_gmset_type = params_all.lookup("Quark_1").lookup("Fopr").get_string("gamma_matrix_type");
    vout.general(vl, "  gmset_type  = %s\n", str_gmset_type.c_str());

    for (int iq = 0; iq < N_quark; ++iq) {
      vout.general(vl, "  savefile_base[iq] = %s[%d]\n", savefile_base[iq].c_str(), iq);
    }

    if ((Nsmear > 0) && (str_proj_type == "Stout_SU3")) {
      if (CommonParameters::Nc() != 3) {
        vout.crucial(vl, "check skipped: Nc = 3 is needed, but Nc = %d.\n\n", CommonParameters::Nc());
        return EXIT_SKIP;
      }
    }


    // ####  Set up objects  #####
    unique_ptr<GammaMatrixSet> gmset(GammaMatrixSet::New(str_gmset_type));

    unique_ptr<FieldIO> field_io(new FieldIO_Binary(IO_Format::Trivial));

    Corr2pt_4spinor corr(gmset.get(), params_all.lookup("Corr2pt_4spinor"));


    // ####  Execution main part  ####
    typedef std::vector<Field_F> PropagatorSet;
    std::vector<PropagatorSet> sq(N_quark);
    for (int iq = 0; iq < N_quark; ++iq) {
      sq[iq].resize(Nc * Nd);

      for (int i_cd = 0; i_cd < Nc * Nd; ++i_cd) {
        sq[iq][i_cd].set(0.0);
      }
    }


    for (int iq = 0; iq < N_quark; ++iq) {
      for (int ispin = 0; ispin < Nd; ++ispin) {
        for (int icolor = 0; icolor < Nc; ++icolor) {
          int i_cd = icolor + Nc * ispin;

          std::string filename = FileUtils::generate_filename("%s_c%d_s%d.dat", savefile_base[iq].c_str(), (icolor + 1), (ispin + 1));
          field_io->read_file(sq[iq][i_cd], filename);
        }
      }
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

    RandomNumberManager::finalize();


    if (do_check) {
      return Test::verify(result[0], expected_result);
    } else {
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }
  }
} // namespace Test_Spectrum
