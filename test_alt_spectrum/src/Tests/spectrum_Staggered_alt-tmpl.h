/*!
      @file    spectrum_Staggered_alt.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-10-02 13:50:31 #$
      @version $LastChangedRevision: 2543 $
*/

#include "spectrum_Staggered_alt.h"

#include "Tests/test.h"

//#include "lib/Field/field_F.h"
#include "lib/Field/field_F_1spinor.h"
#include "lib/Field/field_G.h"
#include "lib/Tools/gammaMatrixSet.h"
#include "lib/Tools/randomNumberManager.h"
#include "lib/Fopr/fopr_Staggered.h"
#include "lib/Solver/solver.h"
#include "lib/Measurements/Fermion/fprop_Standard_lex.h"
#include "lib/Measurements/Fermion/source.h"
//#include "lib/Measurements/Fermion/corr2pt_4spinor.h"
#include "lib/Measurements/Gauge/staple_lex.h"
#include "lib/Tools/filename.h"
#include "lib/IO/gaugeConfig.h"
#include "lib/Measurements/Fermion/corr2pt_Staggered_Local.h"
#include "lib/Measurements/Fermion/corr2pt_Staggered_Extended.h"
#include "lib/Measurements/Gauge/gaugeFixing.h"

// alt-code
#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_lex.h"
#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_eo.h"
#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_lex_Mixedprec.h"

#include "Tests/job_Utils.h"


// class name
template<Impl IMPL>
const std::string Spectrum_Staggered_alt<IMPL>::class_name =
  "Spectrum_Staggered_alt";

//====================================================================
template<Impl IMPL>
void Spectrum_Staggered_alt<IMPL>::init()
{
  // do nothing.
}


//====================================================================
template<Impl IMPL>
int Spectrum_Staggered_alt<IMPL>::hadron_2ptFunction_Evenodd(
  std::string file_params,
  std::string test_mode,
  std::string run_mode)
{
  std::string test_name = class_name + "hadron_2ptFunction";
  std::string inputfile = "input";  // configuration number is given

  vout.general("\n");
  vout.general("-------------------------------------------------"
               "-------------------\n");
  vout.general("test name: %s\n", test_name.c_str());
  vout.general("parameter file: %s\n", file_params.c_str());
  vout.general("test mode: %s\n", test_mode.c_str());
  vout.general("run mode:  %s\n", run_mode.c_str());

  int Nc   = CommonParameters::Nc();
  int Ndim = CommonParameters::Ndim();
  int Nvol = CommonParameters::Nvol();
  int Lt   = CommonParameters::Lt();

  params_all = ParameterManager::read(file_params);

  Parameters params_test   = params_all.lookup("Test_Spectrum");
  Parameters params_fopr   = params_all.lookup("Fopr");
  Parameters params_source = params_all.lookup("Source");

  string               str_vlevel = params_test.get_string("verbose_level");
  Bridge::VerboseLevel m_vl       = vout.set_verbose_level(str_vlevel);
  vout.general(m_vl, "  vlevel       = %s\n", str_vlevel.c_str());

  bool   do_check        = false;
  double expected_result = 0.0;
  if (run_mode == "test") {
    do_check = params_test.is_set("expected_result");
    if (do_check) expected_result = params_test.get_double("expected_result");
  }

  // setup random number manager
  //RandomNumberManager::initialize("Mseries", 1234567UL);

  // Setup gauge configuration
  unique_ptr<Field_G> U(new Field_G(Nvol, Ndim));
  // U.reset(new Field_G(Nvol, Ndim));
  // setup_config(U, params_test);

  string config_type_input
    = params_test.get_string("gauge_config_type_input");
  string config_file_input
    = params_test.get_string("config_filename_input");
  unique_ptr<GaugeConfig> gconf_read(new GaugeConfig(config_type_input));

  // test
  Staple_lex staple;
  //  double plaq = staple.plaquette((Field_G)*U);
  //  vout.general(m_vl, "Corelib:  plaq = %f\n", plaq);

  // Gauge fixing
  // Parameters params_gfix = params_all.lookup("GaugeFixing");
  // gauge_fixing(U, params_gfix);

  // setup Fprop: computer of fermion propagators
  unique_ptr<Fprop> fprop;

  // objects only used for org modes
  unique_ptr<Solver> solver;
  //  unique_ptr<Fopr>   fopr_eo;  // only for org_eo mode

  Parameters params_solver = params_all.lookup("Solver");

  // fermion operator (reference)
  // unique_ptr<Fopr> fopr(new Fopr_Staggered());
  // fopr->set_parameters(params_fopr);
  unique_ptr<Fopr> fopr(new Fopr_Staggered(params_fopr));

  if (test_mode == "double") {
    fprop.reset(new Fprop_alt_Standard_lex<AField<double, IMPL> >(
                  params_fopr, params_solver));
  } else if (test_mode == "float") {
    fprop.reset(new Fprop_alt_Standard_lex<AField<float, IMPL> >(
                  params_fopr, params_solver));
  } else if (test_mode == "double_eo") {
    fprop.reset(new Fprop_alt_Standard_eo<AField<double, IMPL> >(
                  params_fopr, params_solver));
  } else if (test_mode == "float_eo") {
    fprop.reset(new Fprop_alt_Standard_eo<AField<float, IMPL> >(
                  params_fopr, params_solver));
  } else if (test_mode == "org") {
    string solver_type = params_solver.get_string("solver_type");
    solver.reset(Solver::New(solver_type, fopr.get()));
    fprop.reset(new Fprop_Standard_lex(solver.get()));
  } else {
    vout.crucial(m_vl, "%s: irrelevant test_mode =%s\n",
                 class_name.c_str(), test_mode.c_str());
    exit(EXIT_FAILURE);
  }

  // source setup
  string source_type = params_source.get_string("source_type");
  vout.general(m_vl, "source_type  = %s\n", source_type.c_str());
  unique_ptr<Source> source(Source::New(source_type));
  source->set_parameters(params_source);


  unique_ptr<Timer> timer(new Timer(test_name));

  // Corr2pt_Staggered_Local corr2pt;
  Corr2pt_Staggered_Extended corr2pt;

  Filename filename("Data/meson_{iconf:6}_0.{imass:3}.dat");
  Filename filename2("Data/baryon_{iconf:6}_0.{imass:3}.dat");

  // ####  Execution main part  ####

  int iconf, nconf;
  iconf = 20000;
  nconf = 1;
  if (run_mode == "job") {
    read_input(iconf, nconf, inputfile);
    vout.general(m_vl, "initial iconf   = %d\n", iconf);
    vout.general(m_vl, "number of confs = %d\n", nconf);
  }

  double result = 0.0;

  // loop for configurations
  for (int i = 0; i < nconf; ++i) {
    vout.general(m_vl, "\n");
    vout.general(m_vl, "iconf = %d\n", iconf);

    timer->start();

    if (run_mode == "test") {
      gconf_read->read_file(*U, config_file_input);
    } else {
      gconf_read->read_file(*U, filename_config(config_file_input, iconf));
    }
    double plaq = staple.plaquette((Field_G) * U);
    vout.general(m_vl, "  plaq = %f\n", plaq);

    // gfix->fix(*Ufix, *U);
    // copy(*U, *Ufix);


    fopr->set_config(U.get());
    fprop->set_config(U.get());


    if (run_mode == "test") {
      vout.general(m_vl, "performance measurement start.\n");
      fprop->mult_performance("D", 10);
      fprop->mult_performance("D", 500);
      vout.general(m_vl, "performance measurement finished.\n");
    }

    fprop->set_mode("D");

    Field_F_1spinor              src;
    std::vector<Field_F_1spinor> sq(Nc * 2);

    int    nconv;
    double diff;


    vout.general(m_vl,
                 "   ic  isrc     nconv         diff          diff2\n"
                 " ---------------------------------------------------\n");
    for (int isrc = 0; isrc < 2; ++isrc) {
      for (int ic = 0; ic < Nc; ++ic) {
        source->set(src, ic, isrc);

        int idx = ic + Nc * isrc;
        fprop->invert(sq[idx], src, nconv, diff);

        // check
        Field_F_1spinor y;
        fopr->set_mode("D");
        fopr->mult(y, sq[idx]);
        axpy(y, -1.0, src);
        double diff2 = y.norm2() / src.norm2();

        vout.general(m_vl, "   %2d    %2d    %6d   %12.4e   %12.4e\n",
                     ic, isrc, nconv, diff, diff2);
      }
    }

    vout.general(m_vl,
                 " ---------------------------------------------------\n");

    double mass  = params_fopr.get_double("quark_mass");
    int    imass = int(mass * 1000);

    //- meson correlators
    vout.general(m_vl, "\n");
    vout.general(m_vl, "meson correlator:\n");

    if (run_mode == "job") {
      string file_meson = filename.format(iconf, imass);
      vout.general("    output file = %s\n", file_meson.c_str());
      corr2pt.set_output_file(file_meson);
    }

    std::vector<double> mcorr(Lt);
    corr2pt.meson_all(mcorr, sq, sq);
    //result = mcorr[0];

    int Lx = CommonParameters::Lx();
    int Ly = CommonParameters::Ly();
    int Lz = CommonParameters::Lz();
    result = mcorr[0] / double(Lx * Ly * Lz);

    vout.general(m_vl, "  PS <-- PS correlator:\n");
    for (int t = 0; t < mcorr.size(); ++t) {
      vout.general(m_vl, "  %4d %20.12e\n", t, mcorr[t]);
    }

    //- baryon correlators
    vout.general(m_vl, "\n");
    vout.general(m_vl, "baryon correlator:\n");

    if (run_mode == "job") {
      string file_baryon = filename2.format(iconf, imass);
      vout.general("    output file = %s\n", file_baryon.c_str());
      corr2pt.set_output_file(file_baryon);
    }

    std::vector<double> bcorr(Lt);
    corr2pt.baryon_all(bcorr, sq, sq);

    vout.general(m_vl, "  nucleon correlator:\n");
    for (int t = 0; t < bcorr.size(); ++t) {
      vout.general(m_vl, "  %4d %20.12e\n", t, bcorr[t]);
    }

    vout.general(m_vl, "    measurement done\n");

    timer->report();

    ++iconf;
    write_input(iconf, nconf, inputfile);
  }

  // RandomNumberManager::finalize();

  if (do_check) {
    return Test::verify(result, expected_result);
  } else {
    vout.detailed(m_vl, "check skipped: expected_result not set.\n\n");
    return EXIT_SKIP;
  }
}


//====================================================================
template<Impl IMPL>
int Spectrum_Staggered_alt<IMPL>::hadron_2ptFunction_Cube(
  std::string file_params,
  std::string test_mode,
  std::string run_mode)
{
  std::string test_name = class_name + "hadron_2ptFunction";
  std::string inputfile = "input";  // configuration number is given

  vout.general("\n");
  vout.general("-------------------------------------------------"
               "-------------------\n");
  vout.general("test name: %s\n", test_name.c_str());
  vout.general("parameter file: %s\n", file_params.c_str());
  vout.general("test mode: %s\n", test_mode.c_str());
  vout.general("run mode:  %s\n", run_mode.c_str());

  int Nc   = CommonParameters::Nc();
  int Ndim = CommonParameters::Ndim();
  int Nvol = CommonParameters::Nvol();
  int Lt   = CommonParameters::Lt();

  params_all = ParameterManager::read(file_params);

  Parameters params_test   = params_all.lookup("Test_Spectrum");
  Parameters params_fopr   = params_all.lookup("Fopr");
  Parameters params_source = params_all.lookup("Source");

  string               str_vlevel = params_test.get_string("verbose_level");
  Bridge::VerboseLevel m_vl       = vout.set_verbose_level(str_vlevel);
  vout.general(m_vl, "  vlevel       = %s\n", str_vlevel.c_str());

  bool   do_check        = false;
  double expected_result = 0.0;
  if (run_mode == "test") {
    do_check = params_test.is_set("expected_result");
    if (do_check) expected_result = params_test.get_double("expected_result");
  }

  // setup random number manager
  // RandomNumberManager::initialize("Mseries", 1234567UL);

  // Setup gauge configuration
  unique_ptr<Field_G> U(new Field_G(Nvol, Ndim));
  //  U.reset(new Field_G(Nvol, Ndim));
  // setup_config(U, params_test);

  string config_type_input
    = params_test.get_string("gauge_config_type_input");
  string config_file_input
    = params_test.get_string("config_filename_input");
  unique_ptr<GaugeConfig> gconf_read(new GaugeConfig(config_type_input));

  // test
  Staple_lex staple;
  //  double plaq = staple.plaquette((Field_G)*U);
  //  vout.general(m_vl, "Corelib:  plaq = %f\n", plaq);

  // Gauge fixing
  // Parameters params_gfix = params_all.lookup("GaugeFixing");
  // gauge_fixing(U, params_gfix);

  // setup Fprop: computer of fermion propagators
  unique_ptr<Fprop> fprop;

  // objects only used for org modes
  unique_ptr<Solver> solver;
  //  unique_ptr<Fopr>   fopr_eo;  // only for org_eo mode

  Parameters params_solver = params_all.lookup("Solver");

  // fermion operator (reference)
  unique_ptr<Fopr> fopr(new Fopr_Staggered());
  fopr->set_parameters(params_fopr);


  if (test_mode == "double") {
    fprop.reset(new Fprop_alt_Standard_lex<AField<double, IMPL> >(
                  params_fopr, params_solver));
  } else if (test_mode == "float") {
    fprop.reset(new Fprop_alt_Standard_lex<AField<float, IMPL> >(
                  params_fopr, params_solver));
  } else if (test_mode == "double_eo") {
    fprop.reset(new Fprop_alt_Standard_eo<AField<double, IMPL> >(
                  params_fopr, params_solver));
  } else if (test_mode == "float_eo") {
    fprop.reset(new Fprop_alt_Standard_eo<AField<float, IMPL> >(
                  params_fopr, params_solver));
  } else if (test_mode == "org") {
    string solver_type = params_solver.get_string("solver_type");
    solver.reset(Solver::New(solver_type, fopr.get()));
    solver->set_parameters(params_solver);
    fprop.reset(new Fprop_Standard_lex(solver.get()));
  } else {
    vout.crucial(m_vl, "%s: irrelevant test_mode =%s\n",
                 class_name.c_str(), test_mode.c_str());
    exit(EXIT_FAILURE);
  }

  // source setup
  string source_type = params_source.get_string("source_type");
  vout.general(m_vl, "source_type  = %s\n", source_type.c_str());
  unique_ptr<Source> source(Source::New(source_type));
  source->set_parameters(params_source);

  unique_ptr<Timer> timer(new Timer(test_name));

  Corr2pt_Staggered_Local corr2pt;

  Filename filename("Data/meson_{iconf:6}_0.{imass:3}.dat");
  Filename filename2("Data/baryon_{iconf:6}_0.{imass:3}.dat");

  // ####  Execution main part  ####

  int iconf, nconf;
  iconf = 20000;
  nconf = 1;
  if (run_mode == "job") {
    read_input(iconf, nconf, inputfile);
    vout.general(m_vl, "initial iconf   = %d\n", iconf);
    vout.general(m_vl, "number of confs = %d\n", nconf);
  }

  double result = 0.0;

  // loop for configurations
  for (int i = 0; i < nconf; ++i) {
    vout.general(m_vl, "\n");
    vout.general(m_vl, "iconf = %d\n", iconf);

    timer->start();

    if (run_mode == "test") {
      gconf_read->read_file(*U, config_file_input);
    } else {
      gconf_read->read_file(*U, filename_config(config_file_input, iconf));
    }
    double plaq = staple.plaquette((Field_G) * U);
    vout.general(m_vl, "  plaq = %f\n", plaq);

    // check of operator if necessary
    // check_operator(fopr.get(), params_fopr, U.get());
    //check_operator_eo(fopr.get(), params_fopr, U.get());
    // return 0;

    // gfix->fix(*Ufix, *U);
    // copy(*U, *Ufix);

    fopr->set_config(U.get());
    fprop->set_config(U.get());

    if (run_mode == "test") {
      vout.general(m_vl, "performance measurement start.\n");
      fprop->mult_performance("D", 10);
      fprop->mult_performance("D", 500);
      vout.general(m_vl, "performance measurement finished.\n");
    }

    //    return 0;

    fprop->set_mode("D");

    Field_F_1spinor              src;
    std::vector<Field_F_1spinor> sq(Nc * 2);

    int    nconv;
    double diff;

    vout.general(m_vl,
                 "   ic  isrc     nconv         diff          diff2\n"
                 " ---------------------------------------------------\n");
    // for (int isrc = 0; isrc < 2; ++isrc) {
    for (int isrc = 0; isrc < 1; ++isrc) {
      for (int ic = 0; ic < Nc; ++ic) {
        source->set(src, ic, isrc);

        int idx = ic + Nc * isrc;
        fprop->invert(sq[idx], src, nconv, diff);

        // check
        Field_F_1spinor y;
        fopr->set_mode("D");
        fopr->mult(y, sq[idx]);
        axpy(y, -1.0, src);
        double diff2 = y.norm2() / src.norm2();

        vout.general(m_vl, "   %2d    %2d    %6d   %12.4e   %12.4e\n",
                     ic, isrc, nconv, diff, diff2);
      }
    }

    vout.general(m_vl,
                 " ---------------------------------------------------\n");

    double mass  = params_fopr.get_double("quark_mass");
    int    imass = int(mass * 1000);

    //- meson correlators
    vout.general(m_vl, "\n");
    vout.general(m_vl, "meson correlator:\n");

    if (run_mode == "job") {
      string file_meson = filename.format(iconf, imass);
      vout.general("    output file = %s\n", file_meson.c_str());
      corr2pt.set_output_file(file_meson);
    }

    std::vector<double> mcorr(Lt);
    corr2pt.meson_all(mcorr, sq, sq);
    // result = mcorr[0];

    int Lx = CommonParameters::Lx();
    int Ly = CommonParameters::Ly();
    int Lz = CommonParameters::Lz();
    result = mcorr[0] / double(Lx * Ly * Lz);

    vout.general(m_vl, "  PS <-- PS correlator:\n");
    for (int t = 0; t < mcorr.size(); ++t) {
      vout.general(m_vl, "  %4d %20.12e\n", t, mcorr[t]);
    }

    //- baryon correlators
    vout.general(m_vl, "\n");
    vout.general(m_vl, "baryon correlator:\n");

    if (run_mode == "job") {
      string file_baryon = filename2.format(iconf, imass);
      vout.general("    output file = %s\n", file_baryon.c_str());
      corr2pt.set_output_file(file_baryon);
    }

    std::vector<double> bcorr(Lt);
    corr2pt.baryon_all(bcorr, sq, sq);

    vout.general(m_vl, "  nucleon correlator:\n");
    for (int t = 0; t < bcorr.size(); ++t) {
      vout.general(m_vl, "  %4d %20.12e\n", t, bcorr[t]);
    }

    vout.general(m_vl, "    measurement done\n");

    timer->report();

    ++iconf;
    write_input(iconf, nconf, inputfile);
  }

  // RandomNumberManager::finalize();

  if (do_check) {
    return Test::verify(result, expected_result);
  } else {
    vout.detailed(m_vl, "check skipped: expected_result not set.\n\n");
    return EXIT_SKIP;
  }
}


//====================================================================
template<Impl IMPL>
int Spectrum_Staggered_alt<IMPL>::check_operator(
  Fopr *fopr_ref,
  Parameters& params_fopr,
  Field_G *U)
{
  typedef AField<double, IMPL> AFIELD;

  vout.general(m_vl, "\n");
  vout.general(m_vl, "Check of fermion operator\n");

  int Nvol = CommonParameters::Nvol();
  int Nex  = fopr_ref->field_nex();


  //  fopr_ref->set_config(U.get());
  fopr_ref->set_config(U);
  fopr_ref->set_mode("D");

  Field_F_1spinor b(Nvol, Nex), bt(Nvol, Nex);
  b.set(0.0);
  b.set(0, 1.0);

  Field_F_1spinor y(Nvol, Nex), x(Nvol, Nex);

  // extend source vector
  for (int i = 0; i < 10; ++i) {
    fopr_ref->mult(y, b);
    fopr_ref->mult(b, y);
  }
  double bb = b.norm2();
  bb = 1.0 / sqrt(bb);
  scal(b, bb);
  vout.general("norm of source vector: %f\n", b.norm());

  //### from here, check main part ###

  // construction of alternative fermion object.
  string fopr_type = params_fopr.get_string("fermion_type");
  // fopr_type += "_eo";
  unique_ptr<AFopr<AFIELD> >
  fopr_alt(AFopr<AFIELD>::New(fopr_type, params_fopr));

  fopr_alt->set_config(U);

  // check
  int nin  = fopr_alt->field_nin();
  int nvol = fopr_alt->field_nvol();
  int nex  = fopr_alt->field_nex();
  vout.general(" nin = %d  nvol = %d  nex = %d\n", nin, nvol, nex);

  AFIELD abq(nin, Nvol, nex);
  AFIELD axq(nin, Nvol, nex);
  AFIELD ayq(nin, Nvol, nex);

  AIndex_lex<double, IMPL> index_alt;

  if (fopr_alt->needs_convert()) {
    vout.general("convert required.\n");
    fopr_alt->convert(abq, (Field&)b);
  } else {
    vout.general("convert not required.\n");
    convert(index_alt, abq, (Field&)b);
  }

  //std::string mode = "D";
  std::string mode = "DdagD";
  fopr_ref->set_mode(mode);
  fopr_alt->set_mode(mode);

  fopr_ref->mult(y, b);
  fopr_alt->mult(axq, abq);
  // copy(y, b);
  // copy(axq, abq);

  if (fopr_alt->needs_convert()) {
    fopr_alt->reverse((Field&)x, axq);
  } else {
    reverse(index_alt, (Field&)x, axq);
  }

  double yrnorm = y.norm2();
  double xxnorm = axq.norm2();
  double xrnorm = x.norm2();
  axpy(y, double(-1.0), x);

  double diff2 = y.norm2();
  //yy = yy/yynorm;

  vout.general(m_vl, "\n");
  vout.general(m_vl, "norm2(ref) = %f\n", yrnorm);
  vout.general(m_vl, "norm2(alt) = %f\n", xxnorm);
  vout.general(m_vl, "norm2(alt) = %f\n", xrnorm);
  vout.general(m_vl, "diff2 = %f\n", diff2);


  double epsilon = 1.e-16;
  if (diff2 < epsilon) return 0;

  return 1;
}


//====================================================================
template<Impl IMPL>
int Spectrum_Staggered_alt<IMPL>::check_operator_eo(
  Fopr *fopr_ref,
  Parameters& params_fopr,
  Field_G *U)
{
  typedef float                real_t;
  //typedef double real_t;
  typedef AField<real_t, IMPL> AFIELD;

  vout.general(m_vl, "\n");
  vout.general(m_vl, "Check of fermion operator\n");

  int Nvol = CommonParameters::Nvol();
  int Nex  = fopr_ref->field_nex();


  //  fopr_ref->set_config(U.get());
  fopr_ref->set_config(U);
  fopr_ref->set_mode("D");

  Field_F_1spinor b(Nvol, Nex), bt(Nvol, Nex);
  b.set(0.0);
  b.set(0, 1.0);

  Field_F_1spinor y(Nvol, Nex), x(Nvol, Nex);

  // extend source vector
  for (int i = 0; i < 10; ++i) {
    fopr_ref->mult(y, b);
    fopr_ref->mult(b, y);
  }
  double bb = b.norm2();
  bb = 1.0 / sqrt(bb);
  scal(b, bb);
  vout.general("norm of source vector: %f\n", b.norm());

  //### from here, check main part ###

  // construction of alternative fermion object.
  string fopr_type = params_fopr.get_string("fermion_type");
  fopr_type += "_eo";
  unique_ptr<AFopr<AFIELD> >
  fopr_alt(AFopr<AFIELD>::New(fopr_type, params_fopr));

  fopr_alt->set_config(U);

  // check
  int nin   = fopr_alt->field_nin();
  int nvol2 = fopr_alt->field_nvol();
  int nex   = fopr_alt->field_nex();
  vout.general(" nin = %d  nvol2 = %d  nex = %d\n", nin, nvol2, nex);

  AFIELD abq(nin, Nvol, nex);
  AFIELD axq(nin, Nvol, nex);
  AFIELD ayq(nin, Nvol, nex);

  AIndex_lex<real_t, IMPL> index_alt;

  if (fopr_alt->needs_convert()) {
    vout.general("convert required.\n");
    fopr_alt->convert(abq, (Field&)b);
  } else {
    vout.general("convert not required.\n");
    convert(index_alt, abq, (Field&)b);
  }

  // check of shift
  AFIELD be(nin, nvol2, nex), bo(nin, nvol2, nex);
  AFIELD xe(nin, nvol2, nex), xo(nin, nvol2, nex);
  AFIELD xt(nin, nvol2, nex);

  AIndex_eo<real_t, IMPL> index_eo;
  index_eo.split(be, bo, abq);

  fopr_alt->mult(xe, be, "Dee");
  fopr_alt->mult(xt, bo, "Deo");
  axpy(xe, 1.0, xt);
  fopr_alt->mult(xo, bo, "Doo");
  fopr_alt->mult(xt, be, "Doe");
  axpy(xo, 1.0, xt);

  //  return 0;

  index_eo.merge(axq, xe, xo);

  if (fopr_alt->needs_convert()) {
    fopr_alt->reverse((Field&)x, axq);
  } else {
    reverse(index_alt, (Field&)x, axq);
  }

  std::string mode = "D";
  fopr_ref->set_mode(mode);
  fopr_ref->mult(y, b);


  double yrnorm = y.norm2();
  double xxnorm = axq.norm2();
  double xrnorm = x.norm2();
  axpy(y, double(-1.0), x);

  double diff2 = y.norm2();
  //yy = yy/yynorm;

  vout.general(m_vl, "\n");
  vout.general(m_vl, "norm2(ref) = %f\n", yrnorm);
  vout.general(m_vl, "norm2(alt) = %f\n", xxnorm);
  vout.general(m_vl, "norm2(alt) = %f\n", xrnorm);
  vout.general(m_vl, "diff2 = %f\n", diff2);


  double epsilon = 1.e-16;
  if (diff2 < epsilon) return 0;

  return 1;
}


//============================================================END=====
