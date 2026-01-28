/*!
        @file    test_Spectrum_CRSMatrix_Domainwall.cpp

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate: 2013-01-22 22:08:29 #$

        @version $LastChangedRevision: 2492 $
*/

#include "test.h"

#include "Fopr/fopr_CRS.h"
#include "Field/field_F.h"
#include "Fopr/fopr_Domainwall.h"

#include "IO/gaugeConfig.h"

#include "Measurements/Fermion/source.h"

#include "Solver/solver_CG.h"

//====================================================================
//! Test of CRS matrix format.

/*!
    This class generates and tests the CRS matrix format data.
    Matrix data in CRS format are generated for the following
    fermion operators together with sample source and solution
    vectors.

    The egenrated data can be tested using the function
    CRSsolver() whether really solved result satisfy the linear
    equation.
                                [24 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.        [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                [21 Mar 2015 Y.Namekawa]
 */

namespace Test_Spectrum_CRSMatrix {
  //- test-private parameters
  namespace {
    const std::string filename_input  = "test_Spectrum_CRSMatrix_Domainwall.yaml";
    const std::string filename_output = "stdout";
  }

  //- prototype declaration
  int domainwall(void);

  int CRSsolver(
    const string& solution,
    const string& matrix,
    const string& source,
    double& result       /* return value */
    );

  void write_text(const Field& f, const string& filename);
  void read_text(Field& f, const string& filename);

#ifdef USE_TESTMANAGER_AUTOREGISTER
#ifdef USE_MPI
  // this test runs only in single-node environment.
#else
  //- NB. CRS test is skipped, because it is time-consuming.
#if 0
  namespace {
    static const bool is_registered = TestManager::RegisterTest(
      "Spectrum.CRSMatrix.Domainwall",
      domainwall
      );
  }
#endif
#endif
#endif

  //====================================================================
  int domainwall(void)
  {
    // ####  parameter setup  ####
    const int Nc   = CommonParameters::Nc();
    const int Nd   = CommonParameters::Nd();
    const int Nvol = CommonParameters::Nvol();
    const int Ndim = CommonParameters::Ndim();

    const Parameters params_all = ParameterManager::read(filename_input);

    const Parameters params_test = params_all.lookup("Test_Spectrum_CRSMatrix");

    const string str_gconf_read = params_test.get_string("gauge_config_type_input");
    const string readfile       = params_test.get_string("config_filename_input");
    const string matrix_file    = params_test.get_string("matrix_output");
    const string source_file    = params_test.get_string("source_output");
    const string solution_file  = params_test.get_string("solution_output");
    const string str_vlevel     = params_test.get_string("verbose_level");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

    const Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

    //- print input parameters
    vout.general(vl, "  gconf_read      = %s\n", str_gconf_read.c_str());
    vout.general(vl, "  readfile        = %s\n", readfile.c_str());
    vout.general(vl, "  matrix_output   = %s\n", matrix_file.c_str());
    vout.general(vl, "  source_output   = %s\n", source_file.c_str());
    vout.general(vl, "  solution_output = %s\n", solution_file.c_str());
    vout.general(vl, "  vlevel          = %s\n", str_vlevel.c_str());
    vout.general(vl, "\n");

    //- input parameter check
    int err = 0;
    err += ParameterCheck::non_NULL(str_gconf_read);
    err += ParameterCheck::non_NULL(readfile);
    err += ParameterCheck::non_NULL(matrix_file);
    err += ParameterCheck::non_NULL(source_file);
    err += ParameterCheck::non_NULL(solution_file);

    if (err) {
      vout.crucial(vl, "Error: Test_Spectrum_CRSMatrix: input parameters have not been set\n");
      exit(EXIT_FAILURE);
    }


    // ####  Set up a gauge configuration  ####
    Field_G     U(Nvol, Ndim);
    GaugeConfig gconf_read(str_gconf_read);
    gconf_read.read_file(U, readfile);


    // ####  object setup  #####
    const double     mq = 0.1;
    const double     M0 = 1.6;
    const int        Ns = 8;
    std::vector<int> boundary(Ndim);
    boundary[0] = 1;
    boundary[1] = 1;
    boundary[2] = 1;
    boundary[3] = -1;


    Parameters params_dw;
    params_dw.set_double("quark_mass", mq);
    params_dw.set_double("domain_wall_height", M0);
    params_dw.set_int("extent_of_5th_dimension", Ns);
    params_dw.set_int_vector("boundary_condition", boundary);


    // Fopr_Wilson fopr_w("Dirac");
    // Fopr_Domainwall *fopr_dw = new Fopr_Domainwall(&fopr_w);
    // fopr_dw->set_parameters(params_dw);

    // Fopr_Wilson *fopr_w = new Fopr_Wilson("Dirac");
    // Fopr_Domainwall *fopr_dw = new Fopr_Domainwall(fopr_w, params_dw);

    unique_ptr<Fopr_Domainwall> fopr_dw(new Fopr_Domainwall(params_dw));

    fopr_dw->set_config(&U);
    fopr_dw->set_mode("D");

    unique_ptr<Fopr_CRS> fopr_crs(new Fopr_CRS(fopr_dw.get()));
    fopr_crs->write_matrix(matrix_file);

    // std::vector<int> source_position(Ndim);
    // source_position[0] = 0;
    // source_position[1] = 0;
    // source_position[2] = 0;
    // source_position[3] = 0;
    std::vector<int> source_position { 0, 0, 0, 0 };
    Parameters       params_source;
    params_source.set_int_vector("source_position", source_position);

    unique_ptr<Source> source(Source::New("Local", params_source));


    //- CGNE solver
    int    Niter          = 100;
    int    Nrestart       = 40;
    double Stop_cond      = 1.0e-28;
    bool   use_init_guess = false;

    Parameters params_solver;
    params_solver.set_int("maximum_number_of_iteration", Niter);
    params_solver.set_int("maximum_number_of_restart", Nrestart);
    params_solver.set_double("convergence_criterion_squared", Stop_cond);
    params_solver.set_double("use_initial_guess", use_init_guess);

    unique_ptr<Solver> solver(new Solver_CG(fopr_dw.get(), params_solver));


    // ####  Execution main part  ####
    const int Nin = 2 * Nc * Nd;
    const int Nex = Ns;

    for (int ispin = 0; ispin < Nd; ++ispin) {
      for (int icolor = 0; icolor < Nc; ++icolor) {
        for (int ch = -1; ch <= 1; ch = ch + 2) {
          int idx = icolor + Nc * ispin;

          Field_F b;
          source->set(b, idx);

          Field bt(Nin, Nvol, 1);
          fopr_dw->mult_chproj_4d(bt, b, ch);

          Field bt_5d(Nin, Nvol, Nex);
          bt_5d.set(0.0);
          if (ch == -1) {
            bt_5d.setpart_ex(Nex - 1, bt, 0);
          } else {
            bt_5d.setpart_ex(0, bt, 0);
          }

          //- save the source vector
          write_text(bt_5d, source_file);

          Field b_5d(Nin, Nvol, Nex);
          fopr_dw->set_mode("Ddag");
          fopr_dw->mult(b_5d, bt_5d);

          vout.general(vl, " b_5d norm2 = %14.6e\n", b_5d.norm2());

          Field  xq_5d(Nin, Nvol, Nex);
          int    Nconv;
          double diff;
          fopr_dw->set_mode("DdagD");
          solver->solve(xq_5d, b_5d, Nconv, diff);
          vout.general(vl, "  chirality: %2d\n", ch);
          vout.general(vl, "    Nconv = %d,  diff = %16.8e\n", Nconv, diff);

          //- save the solution vector
          write_text(xq_5d, solution_file);
        }
      }
    }

    double result = 0.0;
    CRSsolver(solution_file, matrix_file, source_file, result);


    if (do_check) {
      return Test::verify(result, expected_result);
    } else {
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }
  }
} // namespace Test_Spectrum_CRSMatrix
