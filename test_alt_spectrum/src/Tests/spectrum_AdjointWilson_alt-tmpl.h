/*!
        @file    $Id:: spectrum_AdjointWilson_alt-tmpl#$

        @brief

        @author  Hideo Matsufuru
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2022-12-16 15:57:38 #$

        @version $LastChangedRevision: 2422 $
*/

#include "spectrum_AdjointWilson_alt.h"

//#include "spectrum_alt.inc"

#include "Tests/test.h"

#include "lib/Field/field_G.h"
#include "lib/Tools/gammaMatrixSet.h"
#include "lib/Solver/solver.h"
#include "lib/Tools/randomNumberManager.h"
#include "lib/Measurements/Fermion/fprop_Standard_lex.h"
#include "lib/Measurements/Fermion/source.h"
#include "lib/Eigen/eigensolver_IRLanczos.h"
#include "lib/Measurements/Fermion/corr2pt_4spinor.h"

//#include "set_impl.h"

#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_lex.h"
#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_lex_Mixedprec.h"

/*
// Registration as a test
namespace Test_Spectrum_AdjointWilson_alt {
  const string test_name = "Spectrum.AdjointWilson.Hadron2ptFunction_alt";

  int hadron_2ptFunction_alt(void){
    Spectrum_AdjointWilson_alt test_wilson;
    int result = 0;
    // result += test_wilson.hadron_2ptFunction("double");
    //  test_wilson.hadron_2ptFunction("float");  // precision not enough
    // result += test_wilson.hadron_2ptFunction("mixed");
    result += test_wilson.hadron_2ptFunction("org");
    return result;
  }

  int eigenspectrum_alt(void){
    Spectrum_AdjointWilson_alt test_wilson;
    int result = 0;
    // result += test_wilson.hadron_2ptFunction("double");
    //  test_wilson.hadron_2ptFunction("float");  // precision not enough
    // result += test_wilson.hadron_2ptFunction("mixed");
    result += test_wilson.eigenspectrum("org");
    return result;
  }

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
    static const bool is_registered = TestManager::RegisterTest(
      test_name,
      hadron_2ptFunction_alt);
  }
#endif
};
*/

// class name
template<Impl IMPL>
const std::string Spectrum_AdjointWilson_alt<IMPL>::class_name =
                                        "Spectrum_AdjointWilson_alt";
//====================================================================
template<Impl IMPL>
void Spectrum_AdjointWilson_alt<IMPL>::init()
{
  // do nothing.
}

//====================================================================
template<Impl IMPL>
int Spectrum_AdjointWilson_alt<IMPL>::hadron_2ptFunction(
                            std::string test_file, std::string mode)
{
  const std::string test_name = class_name + "hadron_2ptFunction";

  vout.general("\n");
  vout.general("test name: %s\n", test_name.c_str());
  vout.general("test file: %s\n", test_file.c_str());

  // ####  parameter setup  ####
  int Nc   = CommonParameters::Nc();
  int Nd   = CommonParameters::Nd();
  int Ndim = CommonParameters::Ndim();
  int Nvol = CommonParameters::Nvol();

  params_all = ParameterManager::read(test_file);

  Parameters params_test = params_all.lookup("Test_Spectrum");
  Parameters params_fopr = params_all.lookup("Fopr");
  Parameters params_source = params_all.lookup("Source");

  string str_vlevel = params_test.get_string("verbose_level");
  Bridge::VerboseLevel m_vl = vout.set_verbose_level(str_vlevel);
  vout.general(m_vl, "  vlevel       = %s\n", str_vlevel.c_str());

  const bool   do_check        = params_test.is_set("expected_result");
  const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

  // setup random number manager
  RandomNumberManager::initialize("Mseries", 1234567UL);

  // Setup gauge configuration 
  U.reset(new Field_G(Nvol, Ndim));
  setup_config(U, params_test);

  // Gauge fixing
  Parameters params_gfix = params_all.lookup("GaugeFixing");
  gauge_fixing(U, params_gfix);

  // gamma matrix set
  string str_gmset_type  = params_fopr.get_string("gamma_matrix_type");
  vout.general(m_vl, "gmset_type   = %s\n", str_gmset_type.c_str());
  unique_ptr<GammaMatrixSet> gmset(GammaMatrixSet::New(str_gmset_type));

  // fermion operator (reference)
  string str_fopr_type   = params_fopr.get_string("fermion_type");
  unique_ptr<Fopr> fopr(Fopr::New(str_fopr_type, str_gmset_type));
  fopr->set_parameters(params_fopr);
  fopr->set_config(U.get());

  // setup Fprop: computer of fermion propagators
  unique_ptr<Solver> solver;
  unique_ptr<Fprop>  fprop;

  Parameters params_solver = params_all.lookup("Solver");

  if(mode == "double"){
    vout.crucial(m_vl, "not supported yet.\n");
    // fprop.reset(new
    //    Fprop_alt_Standard_lex<AField<double> >(params_fopr,
    //                                           params_solver) );
  }else if(mode == "float"){
    vout.crucial(m_vl, "not supported yet.\n");
  }else if(mode == "mixed"){
    vout.crucial(m_vl, "not supported yet.\n");
  }else if(mode == "org"){
    string solver_type = params_solver.get_string("solver_type");
    //    solver_type += "_Cmplx";
    solver.reset(Solver::New(solver_type, fopr.get()));
    solver->set_parameters(params_solver);
    fprop.reset(new Fprop_Standard_lex(solver.get()));
  }else{
    vout.crucial(m_vl, "%s: irrelevant mode =%s\n",
		 class_name.c_str(), mode.c_str());
    exit(EXIT_FAILURE);
  }

  fprop->set_config(U.get());

  fprop->mult_performance("D", 100);

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

  //  Field_F b, y;
  int ninF = fopr->field_nin();
  int nvol = fopr->field_nvol();
  int nex  = fopr->field_nex();

  //  std::vector<Field_F> sq(Nc * Nd);
  std::vector<Field> sq(Nc * Nd);
  for(int i = 0; i < Nc*Nd; ++i){
    sq[i].reset(ninF, nvol, nex);
  }

  Field b(ninF, nvol, nex);
  Field y(ninF, nvol, nex);
  b.set(0.0);

  int    nconv;
  double diff;

  for (int ispin = 0; ispin < Nd; ++ispin) {
    for (int icolor = 0; icolor < Nc; ++icolor) {

      int idx = icolor + Nc * ispin;
      sq[idx].set(0.0);
      source->set(b, idx);

      // fprop->invert_D(sq[idx], b, nconv, diff);
      fprop->invert(sq[idx], b, nconv, diff);

      fopr->set_mode("D");
      fopr->mult(y, sq[idx]);
      axpy(y, -1.0, b);
      double diff2 = y.norm2() / b.norm2();

      vout.general(m_vl, "   %2d   %2d   %6d   %12.4e   %12.4e\n",
                   icolor, ispin, nconv, diff, diff2);

    }
  }

  fprop->report_performance();


  unique_ptr<Timer> timer(new Timer(test_name));
  timer->start();

  // meson correlators
  vout.general(m_vl, "\n");
  vout.general(m_vl, "2-point correlator:\n");

  double result = 0.0;
  /*
  Corr2pt_4spinor   corr(gmset);
  corr.set_parameters(params_all.lookup("Corr2pt_4spinor"));

  result = corr.meson_all(sq, sq);
  */

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
int Spectrum_AdjointWilson_alt<IMPL>::eigenspectrum(
                            std::string test_file, std::string mode)
{
  const std::string test_name = class_name + ".eigenspectrum";

  vout.general("\n");
  vout.general("test name: %s\n", test_name.c_str());
  vout.general("test file: %s\n", test_file.c_str());

  // ####  parameter setup  ####
  int Nc   = CommonParameters::Nc();
  int Nd   = CommonParameters::Nd();
  int Ndim = CommonParameters::Ndim();
  int Nvol = CommonParameters::Nvol();

  params_all = ParameterManager::read(test_file);

  Parameters params_test = params_all.lookup("Test_Spectrum");
  Parameters params_fopr = params_all.lookup("Fopr");

  string str_vlevel = params_test.get_string("verbose_level");
  Bridge::VerboseLevel m_vl = vout.set_verbose_level(str_vlevel);
  vout.general(m_vl, "  vlevel       = %s\n", str_vlevel.c_str());

  const bool   do_check        = params_test.is_set("expected_result");
  const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

  // setup random number manager
  RandomNumberManager::initialize("Mseries", 1234567UL);

  // Setup gauge configuration 
  U.reset(new Field_G(Nvol, Ndim));
  setup_config(U, params_test);

  // Gauge fixing
  Parameters params_gfix = params_all.lookup("GaugeFixing");
  gauge_fixing(U, params_gfix);

  // gamma matrix set
  string str_gmset_type  = params_fopr.get_string("gamma_matrix_type");
  vout.general(m_vl, "gmset_type   = %s\n", str_gmset_type.c_str());
  unique_ptr<GammaMatrixSet> gmset(GammaMatrixSet::New(str_gmset_type));

  // fermion operator (reference)
  string str_fopr_type   = params_fopr.get_string("fermion_type");
  unique_ptr<Fopr> fopr(Fopr::New(str_fopr_type, str_gmset_type));
  fopr->set_parameters(params_fopr);
  fopr->set_config(U.get());
  fopr->set_mode("H");

  Parameters params_eigen = params_all.lookup("Eigensolver");

  unique_ptr<Eigensolver> eigen(new Eigensolver_IRLanczos(fopr.get()));
  eigen->set_parameters(params_eigen);



  if(mode == "double"){
    vout.crucial(m_vl, "not supported yet.\n");
    // fprop.reset(new
    //    Fprop_Standard_lex_alt<AField<double> >(params_fopr,
    //                                           params_solver) );
  }else if(mode == "float"){
    vout.crucial(m_vl, "not supported yet.\n");
     // fprop.reset(new
     //    Fprop_Standard_lex_alt<AField<float> >(params_fopr,
    //                                          params_solver) );
  }else if(mode == "mixed"){
    vout.crucial(m_vl, "not supported yet.\n");
     // fprop.reset(new
     //    Fprop_Standard_lex_alt_Mixedprec(params_fopr, params_solver) );
  }else if(mode == "org"){
    // do nothing.
  }else{
    vout.crucial(m_vl, "%s: irrelevant mode =%s\n",
		 class_name.c_str(), mode.c_str());
    exit(EXIT_FAILURE);
  }


  unique_ptr<Timer> timer(new Timer(test_name));


  // ####  Execution main part  ####                                          
  timer->start();


  const int Nk = params_eigen.get_int("number_of_wanted_eigenvectors");
  const int Np = params_eigen.get_int("number_of_working_eigenvectors");

  int Nm = Nk + Np;
  std::vector<double> TDa(Nm);
  std::vector<Field>  vk(Nm);

  int     NFin  = fopr->field_nin();
  int     NFvol = fopr->field_nvol();
  int     NFex  = fopr->field_nex();
  for (int k = 0; k < Nm; ++k) {
    vk[k].reset(NFin, NFvol, NFex);
  }
  Field b2(NFin, NFvol, NFex);
  b2.set(0.0);
  b2.set(0, 0, 0, 1.0);

  int Nsbt  = -1;
  int Nconv = -100;
  eigen->solve(TDa, vk, Nsbt, Nconv, b2);

  Field v(NFin, NFvol, NFex);
  double vv = 0.0;  // superficial initialization                             

  for (int i = 0; i < Nsbt + 1; ++i) {
    fopr->mult(v, vk[i]);
    axpy(v, -TDa[i], vk[i]);  // v -= TDa[i] * vk[i];                         
    vv = v.norm2();           // vv = v * v;                                  

    vout.general(m_vl, "Eigenvalues: %4d %20.14f %20.15e \n",
                 i, TDa[i], vv);
  }

  double result = TDa[0];

  timer->report();

  RandomNumberManager::finalize();


  if (do_check) {
    return Test::verify(result, expected_result);
  } else {
    vout.detailed(m_vl, "check skipped: expected_result not set.\n\n");
    return EXIT_SKIP;
  }

}

//============================================================END=====
