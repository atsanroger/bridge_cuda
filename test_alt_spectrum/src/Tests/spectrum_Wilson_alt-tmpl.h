/*!
        @file    spectrum_Wilson_alt-tmpl.h

        @brief

        @author  Hideo Matsufuru
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2024-01-18 21:39:51 #$

        @version $LastChangedRevision: 2570 $
*/

#include "spectrum_Wilson_alt.h"
//#include "spectrum_alt.inc"
#include "test.h"

#include "lib/Field/field_G.h"
#include "lib/Field/field_F.h"
#include "lib/Fopr/afopr.h"
#include "lib/Fopr/fopr.h"
#include "lib/Tools/gammaMatrixSet.h"
#include "lib/Tools/randomNumberManager.h"
#include "lib/Solver/solver.h"
#include "lib/Measurements/Fermion/fprop_Standard_lex.h"
#include "lib/Measurements/Fermion/fprop_Standard_eo.h"
#include "lib/Measurements/Fermion/source.h"
#include "lib/Measurements/Fermion/corr2pt_4spinor.h"

#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_lex.h"
#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_lex_Mixedprec.h"
#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_eo.h"
#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_eo_Mixedprec.h"
#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_eo_Richardson.h"


template<Impl IMPL>
const std::string Spectrum_Wilson_alt<IMPL>::class_name =
  "Spectrum_Wilson_alt";

//====================================================================
template<Impl IMPL>
void Spectrum_Wilson_alt<IMPL>::init()
{
  // do nothing.
}


//====================================================================
template<Impl IMPL>
int Spectrum_Wilson_alt<IMPL>::hadron_2ptFunction(
                           std::string file_params, std::string mode)
{
  std::string test_name = class_name + "hadron_2ptFunction";

  vout.general("\n");
  vout.general("-------------------------------------------------"
               "-------------------\n");
  vout.general("test name: %s\n", test_name.c_str());
  vout.general("test file: %s\n", file_params.c_str());
  vout.general("test mode: %s\n", mode.c_str());

  int Nc   = CommonParameters::Nc();
  int Nd   = CommonParameters::Nd();
  int Ndim = CommonParameters::Ndim();
  int Nvol = CommonParameters::Nvol();

  params_all = ParameterManager::read(file_params);

  Parameters params_test   = params_all.lookup("Test_Spectrum");
  Parameters params_fopr   = params_all.lookup("Fopr");
  Parameters params_source = params_all.lookup("Source");

  string               str_vlevel = params_test.get_string("verbose_level");
  Bridge::VerboseLevel m_vl       = vout.set_verbose_level(str_vlevel);
  vout.general(m_vl, "  vlevel       = %s\n", str_vlevel.c_str());

  bool   do_check        = params_test.is_set("expected_result");
  double expected_result = 0.0;
  if (do_check) expected_result = params_test.get_double("expected_result");

  // setup random number manager
  RandomNumberManager::initialize("Mseries", 1234567UL);

  // Setup gauge configuration
  U.reset(new Field_G(Nvol, Ndim));
  setup_config(U, params_test);

  // Gauge fixing
  Parameters params_gfix = params_all.lookup("GaugeFixing");
  gauge_fixing(U, params_gfix);
  vout.general("\n");

  // gamma matrix set
  string gmset_type = params_fopr.get_string("gamma_matrix_type");
  vout.general(m_vl, "gmset_type   = %s\n", gmset_type.c_str());
  unique_ptr<GammaMatrixSet> gmset(GammaMatrixSet::New(gmset_type));

  // fermion operator (reference)
  string fopr_type = params_fopr.get_string("fermion_type");
  if (fopr_type == "Clover_dd") fopr_type = "Clover";
  //  unique_ptr<Fopr> fopr(Fopr::New(fopr_type, gmset_type));
  unique_ptr<Fopr> fopr(Fopr::New(fopr_type, params_fopr));
  //fopr->set_parameters(params_fopr);
  //unique_ptr<Fopr> fopr(Fopr::New(fopr_type, params_fopr));
  fopr->set_config(U.get());

  // To confirm the alt operator works same as corelib.
  if(mode == "double"){
    check_operator<AField<double, IMPL> >(fopr, params_fopr);
  }else if (mode == "double_eo"){
    check_operator_eo<AField<double, IMPL> >(fopr, params_fopr);
  }else if (mode == "float"){
    check_operator<AField<float, IMPL> >(fopr, params_fopr);
  }else if (mode == "float_eo"){
    check_operator_eo<AField<float, IMPL> >(fopr, params_fopr);
#ifdef USE_FP16
  }else if (mode == "half"){
    check_operator<AField<half, IMPL> >(fopr, params_fopr);
  }else if (mode == "half_eo"){
    check_operator_eo<AField<half, IMPL> >(fopr, params_fopr);
#endif
  }

  //return 0;

  // setup Fprop: computer of fermion propagators
  unique_ptr<Fprop> fprop;

  // objects only used for org modes
  unique_ptr<Solver> solver;
  unique_ptr<Fopr>   fopr_eo;  // only for org_eo mode

  Parameters params_solver = params_all.lookup("Solver");

  if (mode == "double") {
    fprop.reset(new Fprop_alt_Standard_lex<AField<double, IMPL> >(
                  params_fopr, params_solver));
  } else if (mode == "double_eo") {
    fprop.reset(new Fprop_alt_Standard_eo<AField<double, IMPL> >(
                  params_fopr, params_solver));
  } else if (mode == "float") {
    params_solver.set_double("convergence_criterion_squared", 1.0e-14);
    fprop.reset(new Fprop_alt_Standard_lex<AField<float, IMPL> >(
                  params_fopr, params_solver));
  } else if (mode == "float_eo") {
    params_solver.set_double("convergence_criterion_squared", 1.0e-14);
    fprop.reset(new Fprop_alt_Standard_eo<AField<float, IMPL> >(
                  params_fopr, params_solver));
#ifdef USE_FP16
  } else if (mode == "half") {
    params_solver.set_double("convergence_criterion_squared", 1.0e-7);
    fprop.reset(new Fprop_alt_Standard_lex<AField<half, IMPL> >(
                  params_fopr, params_solver));
  } else if (mode == "half_eo") {
    params_solver.set_double("convergence_criterion_squared", 1.0e-7);
    fprop.reset(new Fprop_alt_Standard_eo<AField<half, IMPL> >(
                  params_fopr, params_solver));
#endif
  } else if (mode == "mixed") {
    fprop.reset(new Fprop_alt_Standard_lex_Mixedprec<AField<double, IMPL>,
                                                     AField<float, IMPL> >(
                  params_fopr, params_solver));
  } else if (mode == "mixed_eo") {
    fprop.reset(new Fprop_alt_Standard_eo_Mixedprec<AField<double, IMPL>,
                                                    AField<float, IMPL> >(
                  params_fopr, params_solver));
  } else if (mode == "richardson_eo") {
    fprop.reset(new Fprop_alt_Standard_eo_Richardson<AField<double, IMPL>,
                                                     AField<float, IMPL> >(
                  params_fopr, params_solver));
  } else if (mode == "org") {
    string solver_type = params_solver.get_string("solver_type");
    if (solver_type == "BiCGStab") solver_type += "_Cmplx";
    solver.reset(Solver::New(solver_type, fopr.get()));
    solver->set_parameters(params_solver);
    fprop.reset(new Fprop_Standard_lex(solver.get()));
  } else if (mode == "org_eo") {
    string fopr_type_eo = fopr_type + "_eo";
    fopr_eo.reset(Fopr::New(fopr_type_eo, gmset_type));
    fopr_eo->set_parameters(params_fopr);
    fopr_eo->set_config(U.get());
    string solver_type = params_solver.get_string("solver_type");
    if (solver_type == "BiCGStab") solver_type += "_Cmplx";
    solver.reset(Solver::New(solver_type, fopr_eo.get()));
    solver->set_parameters(params_solver);
    fprop.reset(new Fprop_Standard_eo(solver.get()));
  } else {
    vout.crucial(m_vl, "%s: irrelevant mode =%s\n",
                 class_name.c_str(), mode.c_str());
    exit(EXIT_FAILURE);
  }

  fprop->set_config(U.get());

  fprop->mult_performance("D", 10);

  fprop->mult_performance("D", 500);

  fprop->set_mode("D");

  // source setup
  string str_source_type = params_source.get_string("source_type");
  vout.general(m_vl, "source_type  = %s\n", str_source_type.c_str());
  unique_ptr<Source> source(Source::New(str_source_type));
  source->set_parameters(params_source);

  //  check<double>(fopr.get(), source.get(), params_fopr);

  // getting fermion propagator
  vout.general(m_vl, "\n");
  vout.general(m_vl, "Solving quark propagator:\n");
  vout.general(m_vl, "  color spin   Nconv      diff           diff2\n");

  std::vector<Field_F> sq(Nc * Nd);

  Field_F b, y;
  b.set(0.0);

  int    nconv;
  double diff;

  for (int id = 0; id < Nd; ++id) {
    for (int ic = 0; ic < Nc; ++ic) {
      int idx = ic + Nc * id;
      sq[idx].set(0.0);
      source->set(b, idx);

      // fprop->invert_D(sq[idx], b, nconv, diff);
      fprop->invert(sq[idx], b, nconv, diff);

      fopr->set_mode("D");
      fopr->mult(y, sq[idx]);
      axpy(y, -1.0, b);
      double diff2 = y.norm2() / b.norm2();

      vout.general(m_vl, "   %2d   %2d   %6d   %12.4e   %12.4e\n",
                   ic, id, nconv, diff, diff2);
    }
  }

  fprop->report_performance();


  unique_ptr<Timer> timer(new Timer(test_name));
  timer->start();

  // meson correlators
  vout.general(m_vl, "\n");
  vout.general(m_vl, "2-point correlator:\n");

  Corr2pt_4spinor corr(gmset.get());
  corr.set_parameters(params_all.lookup("Corr2pt_4spinor"));

  double result = corr.meson_all(sq, sq);

  timer->report();

  RandomNumberManager::finalize();

  if (do_check) {
    return Test::verify(result, expected_result);
  } else {
    vout.detailed(m_vl, "check skipped: expected_result not set.\n\n");
    return EXIT_SKIP;
  }
}


//====================================================================
template<Impl IMPL>
template<typename AFIELD>
int Spectrum_Wilson_alt<IMPL>::check_operator(
                                        unique_ptr<Fopr>& fopr_ref,
                                        Parameters& params_fopr)
{
  typedef typename AFIELD::real_t real_t;

  vout.general(m_vl, "\n");
  vout.general(m_vl, "Check of fermion operator\n");
  vout.general(m_vl, "  sizeof(real_t) = %d\n", sizeof(real_t));

  int Nvol = CommonParameters::Nvol();
  int Nex  = fopr_ref->field_nex();

  fopr_ref->set_config(U.get());
  fopr_ref->set_mode("D");

  Field_F b(Nvol, Nex), bt(Nvol, Nex);
  b.set(0.0);
  b.set(0, 1.0);

  Field_F y(Nvol, Nex), x(Nvol, Nex);

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
  unique_ptr<AFopr<AFIELD> >
     fopr_alt(AFopr<AFIELD>::New(fopr_type, params_fopr));

  fopr_alt->set_config(U.get());

  // check
  int nin  = fopr_alt->field_nin();
  int nvol = fopr_alt->field_nvol();
  int nex  = fopr_alt->field_nex();

  vout.general(" nin = %d  nvol = %d  nex = %d\n", nin, nvol, nex);

  AFIELD abq(nin, Nvol, nex);
  AFIELD axq(nin, Nvol, nex);
  AFIELD ayq(nin, Nvol, nex);

  AIndex_lex<real_t, IMPL> index_alt;

  if (fopr_alt->needs_convert()) {
    vout.general("convert required.\n");
    fopr_alt->convert(abq, (Field&)b);
  } else {
    vout.general("convert not required.\n");
    convert_spinor(index_alt, abq, (Field&)b);
  }

  double bn2 = b.norm2();
  vout.crucial("b norm2 = %e\n", bn2);

  double abqn2 = abq.norm2();
  vout.crucial("abq norm2 = %e\n", abqn2);

  std::string fopr_mode = "DdagD";
  fopr_ref->set_mode(fopr_mode);
  fopr_alt->set_mode(fopr_mode);

  y.set(0.0);
  fopr_ref->mult(y, b);

  axq.set(0.0);
  fopr_alt->mult(axq, abq);

  /*
  for(int i = 0; i < 1; ++i){
    fopr_alt->mult(ayq, axq);
    fopr_alt->mult(axq, ayq);
    fopr_ref->mult(x, y);
    fopr_ref->mult(y, x);
  }
  */

  if (fopr_alt->needs_convert()) {
    fopr_alt->reverse((Field&)x, axq);
  } else {
    reverse_spinor(index_alt, (Field&)x, axq);
  }

  double xxnorm = axq.norm2();

  double yrnorm = y.norm2();
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
template<typename AFIELD>
int Spectrum_Wilson_alt<IMPL>::check_operator_eo(
                                         unique_ptr<Fopr>& fopr_ref,
                                         Parameters& params_fopr)
{
  typedef typename AFIELD::real_t real_t;

  vout.general(m_vl, "\n");
  vout.general(m_vl, "Check of fermion operator\n");
  vout.general(m_vl, "  sizeof(real_t) = %d\n", sizeof(real_t));

  int Nvol = CommonParameters::Nvol();
  int Nex  = fopr_ref->field_nex();

  fopr_ref->set_config(U.get());
  fopr_ref->set_mode("D");

  Field_F b(Nvol, Nex), bt(Nvol, Nex);
  b.set(0.0);
  b.set(0, 1.0);

  Field_F y(Nvol, Nex), x(Nvol, Nex);

  // extend source vector
  for (int i = 0; i < 10; ++i) {
    //for (int i = 0; i < 2; ++i) {  // for half multiply less times
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

  fopr_alt->set_config(U.get());

  // check
  int nin  = fopr_alt->field_nin();
  int nvol = fopr_alt->field_nvol();
  int nex  = fopr_alt->field_nex();

  vout.general(" nin = %d  nvol = %d  nex = %d\n", nin, nvol, nex);

  AFIELD abq(nin, Nvol, nex);
  AFIELD axq(nin, Nvol, nex);
  AFIELD ayq(nin, Nvol, nex);

  AIndex_lex<real_t, IMPL> index_alt;

  if (fopr_alt->needs_convert()) {
    vout.general("convert required.\n");
    fopr_alt->convert(abq, (Field&)b);
  } else {
    vout.general("convert not required.\n");
    convert_spinor(index_alt, abq, (Field&)b);
  }

  double bn2 = b.norm2();
  vout.crucial("b norm2 = %e\n", bn2);

  double abqn2 = abq.norm2();
  vout.crucial("abq norm2 = %e\n", abqn2);

  int nvol2 = Nvol/2;
  AFIELD be(nin, nvol2, nex), bo(nin, nvol2, nex);
  AFIELD xe(nin, nvol2, nex), xo(nin, nvol2, nex);
  AFIELD xt(nin, nvol2, nex);

  AIndex_eo<real_t, IMPL> index_eo;
  index_eo.split(be, bo, abq);

  fopr_alt->mult(xe, be, "Dee");
  fopr_alt->mult(xt, bo, "Deo");
  axpy(xe, real_t(1.0), xt);
  fopr_alt->mult(xo, bo, "Doo");
  fopr_alt->mult(xt, be, "Doe");
  axpy(xo, real_t(1.0), xt);

  index_eo.merge(axq, xe, xo);

  if (fopr_alt->needs_convert()) {
    fopr_alt->reverse((Field&)x, axq);
  } else {
    reverse_spinor(index_alt, (Field&)x, axq);
  }

  std::string fopr_mode = "D";
  fopr_ref->set_mode(fopr_mode);
  y.set(0.0);
  fopr_ref->mult(y, b);

  //double xxenorm = xe.norm2();
  //double xxonorm = xo.norm2();
  //vout.general(m_vl, "norm2_e(alt) = %f\n", xxenorm);
  //vout.general(m_vl, "norm2_o(alt) = %f\n", xxonorm);

  double xxnorm = axq.norm2();
  double yrnorm = y.norm2();
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
