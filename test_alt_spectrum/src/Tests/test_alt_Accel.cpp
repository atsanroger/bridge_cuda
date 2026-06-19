/*!
      @File    test_alt_Accel.cpp
      @brief
      @author  Hideo Matsufuru
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2024-10-27 21:04:29 #$
      @version $LastChangedRevision: 2610 $
*/

#ifdef USE_ALT_ACCEL

#include "lib_alt_Accel/bridge_alt_accel.h"
#define IMPL ACCEL


#include "spectrum_Wilson_alt-tmpl.h"
#include "spectrum_Domainwall_alt-tmpl.h"
#include "spectrum_Staggered_alt-tmpl.h"
#include "spectrum_LinearAlgebra_alt-tmpl.h"
#include "eigenvalue_alt-tmpl.h"
#include "spectrum_Domainwall_PVprec_alt-tmpl.h"
#include "spectrum_AdjointWilson_alt-tmpl.h"

namespace Test_alt_Accel{
  int test_all();
  int test_Spectrum_Wilson();
  int test_Spectrum_Clover();
  int test_Spectrum_Domainwall();
  int test_Spectrum_Domainwall_5din();
  int test_Spectrum_OptimalDomainwall();
  int test_Spectrum_Staggered();
  int test_Spectrum_LinearAlgebra();
  int test_Eigenvalue_Wilson_Lanczos();
  int test_Eigenvalue_Clover_Lanczos();
  int test_Eigenvalue_Wilson_Arnoldi();
  int test_Eigenvalue_Domainwall_PVprec_Arnoldi();
  int test_Solve_Domainwall_PVprec_BiCGStab();
  int test_Solve_Domainwall_PVprec_Propagator();
  int test_Solve_Domainwall_PVprec_2pt();
  int test_Eigenvalue_Staggered_Lanczos();
  int test_Eigenvalue_Staggered_Arnoldi();
  int test_Spectrum_AdjointWilson_Lanczos();
  int test_Eigenvalue_AdjointWilson_Lanczos();
}

//====================================================================
int Test_alt_Accel::test_all()
{
  int result = 0;
  // result += test_Spectrum_Wilson();
  // result += test_Spectrum_Clover();
  // result += test_Spectrum_Staggered();
  // result += test_Eigenvalue_Domainwall_PVprec_Arnoldi();
  // result += test_Solve_Domainwall_PVprec_BiCGStab();
  // result += test_Solve_Domainwall_PVprec_Propagator(); // FLOAT GMRES-on-A: convergence (nconv) + true residual vs alpha
   result += test_Spectrum_Domainwall_5din();           // CGNE 2pt + m_res (reference; non-eo vs eo alpha cross-check)
   result += test_Solve_Domainwall_PVprec_2pt();        // A-based GMRES 2pt + m_res: confirm vs preset 4.58e-01 / mres 1.87e-3
  // result += test_Spectrum_Domainwall();
  // result += test_Spectrum_OptimalDomainwall();
  // result += test_Eigenvalue_Wilson_Lanczos();
  // result += test_Eigenvalue_Clover_Lanczos();
  // result += test_Eigenvalue_Wilson_Arnoldi();
  // result += test_Eigenvalue_Staggered_Arnoldi();
  // result += test_Spectrum_LinearAlgebra();
  // result += test_Spectrum_AdjointWilson_Lanczos();
  // result += test_Eigenvalue_AdjointWilson_Lanczos();

  vout.general("Test result for Accel = %d\n", result);
  return result;
}

//====================================================================
int Test_alt_Accel::test_Spectrum_Wilson()
{
  Spectrum_Wilson_alt<IMPL> test_wilson;
  string file_params = "test_alt_Spectrum_Wilson_Hadron2ptFunction.yaml";

  int result = 0;
  result += test_wilson.hadron_2ptFunction(file_params, "double");
  //test_wilson.hadron_2ptFunction(file_params, "float");
  //result += test_wilson.hadron_2ptFunction(file_params, "double_eo");
  //test_wilson.hadron_2ptFunction(file_params, "float_eo");

  //CLE test
  //string file_params = "test_alt_Spectrum_Wilson_CLE_Hadron2ptFunction.yaml";
  //result += test_wilson.hadron_2ptFunction(file_params, "double");


  /*
  result += test_wilson.hadron_2ptFunction(file_params, "mixed");
  result += test_wilson.hadron_2ptFunction(file_params, "mixed_eo");
  result += test_wilson.hadron_2ptFunction(file_params, "org");
  result += test_wilson.hadron_2ptFunction(file_params, "org_eo");
  */

  return result;
}

//====================================================================
int Test_alt_Accel::test_Spectrum_Clover()
{
  Spectrum_Wilson_alt<IMPL> test_wilson;
  string file_params = "test_alt_Spectrum_Clover_Hadron2ptFunction.yaml";

  int result = 0;
  //result += test_wilson.hadron_2ptFunction(file_params, "double");
  //test_wilson.hadron_2ptFunction(file_params, "float");
  result += test_wilson.hadron_2ptFunction(file_params, "double_eo");
  //test_wilson.hadron_2ptFunction(file_params, "float_eo");
  /*
  result += test_wilson.hadron_2ptFunction(file_params, "mixed");
  result += test_wilson.hadron_2ptFunction(file_params, "mixed_eo");
  result += test_wilson.hadron_2ptFunction(file_params, "org");
  result += test_wilson.hadron_2ptFunction(file_params, "org_eo");
  */
  return result;
}

//====================================================================
int Test_alt_Accel::test_Spectrum_Staggered()
{
  Spectrum_Staggered_alt<IMPL> test_staggered;

  int result = 0;

  //string run_mode = "job";
  string run_mode = "test";

  string file_gfix = "test_alt_Spectrum_Staggered_GaugeFixing.yaml";
  //result += test_staggered.gauge_fixing(file_gfix, run_mode);

  string file_eo =
     "test_alt_Spectrum_Staggered_Hadron2ptFunction_Wall_Evenodd.yaml";
  // "job_alt_Spectrum_Staggered_Hadron2ptFunction_Wall_Evenodd.yaml";
  // result += test_staggered.hadron_2ptFunction_Evenodd(file_eo,
  //                                               "double", run_mode);

  string file_cube =
    "test_alt_Spectrum_Staggered_Hadron2ptFunction_Wall_Cube.yaml";
    // "job_alt_Spectrum_Staggered_Hadron2ptFunction_Wall_Cube.yaml";
  result += test_staggered.hadron_2ptFunction_Cube(file_cube,
                                               "double", run_mode);

  result += test_staggered.hadron_2ptFunction_Cube(file_cube,
                                               "float", run_mode);

  result += test_staggered.hadron_2ptFunction_Cube(file_cube,
                                               "double_eo", run_mode);

  result += test_staggered.hadron_2ptFunction_Cube(file_cube,
                                               "float_eo", run_mode);


  // following tests have not confirmed.
  //test_staggered.hadron_2ptFunction(test_file, "float");
  //result += test_staggered.hadron_2ptFunction(test_file, "double_eo");
  // test_staggered.hadron_2ptFunction("float");  // precision not enough
  // result += test_staggered.hadron_2ptFunction("mixed");
  //result += test_staggered.hadron_2ptFunction(test_file, "org");

  return result;
}

//====================================================================
int Test_alt_Accel::test_Spectrum_Domainwall()
{
  Spectrum_Domainwall_alt<IMPL> test_dw;
  string test_file
     = "test_alt_Spectrum_Domainwall_Hadron2ptFunction.yaml";

  int result = 0;
  result += test_dw.hadron_2ptFunction(test_file, "double"); //ok
  test_dw.hadron_2ptFunction(test_file, "float");
  // result += test_dw.hadron_2ptFunction(test_file, "org"); //ok
  // result += test_dw.hadron_2ptFunction(test_file, "double_prec"); //ok

  //failed
  result += test_dw.hadron_2ptFunction(test_file, "double_eo");
  result += test_dw.hadron_2ptFunction(test_file, "float_eo");

  // unsupported
  // test_dw.hadron_2ptFunction(test_file, "mixed");
  return result;
}

//====================================================================
int Test_alt_Accel::test_Spectrum_Domainwall_5din()
{
  Spectrum_Domainwall_alt<IMPL> test_dw;
  string test_file
     = "test_alt_Spectrum_Domainwall_5din_Hadron2ptFunction.yaml";

  Parameters params_all  = ParameterManager::read(test_file);
  Parameters params_test = params_all.lookup("Test_Spectrum");
  string field_type = "double_eo";
  params_test.fetch_string("field_type", field_type);

  int result = 0;
  result += test_dw.hadron_2ptFunction(test_file, field_type);
  return result;
}

//====================================================================
int Test_alt_Accel::test_Spectrum_OptimalDomainwall()
{
  Spectrum_Domainwall_alt<IMPL> test_dw;
  string test_file
     = "test_alt_Spectrum_OptimalDomainwall_Hadron2ptFunction.yaml";

  int result = 0;
  result += test_dw.hadron_2ptFunction(test_file, "double");

  // unchecked:
  // result += test_optdw.hadron_2ptFunction("org");
  // test_optdw.hadron_2ptFunction("float");
  // test_optdw.hadron_2ptFunction("mixed");

  return result;
}

//====================================================================
int Test_alt_Accel::test_Spectrum_LinearAlgebra()
{
  Spectrum_LinearAlgebra_alt<IMPL> test_la;
  string test_file = "test_alt_Spectrum_Wilson_Hadron2ptFunction.yaml";

  int result = 0;
  result += test_la.hadron_2ptFunction(test_file, "double");

  return result;
}

//====================================================================
int Test_alt_Accel::test_Eigenvalue_Wilson_Lanczos()
{
  Eigenvalue_alt<AField<double, IMPL> > test_wilson;
  // Eigenvalue_alt<IMPL> test_wilson;
  string test_file = "test_alt_Eigenvalue_Lanczos_Wilson.yaml";
  // string test_file = "test_alt_Eigenvalue_Lanczos_Wilson_CLE.yaml";

  int result = 0;
  result += test_wilson.eigenvalue_Lanczos(test_file);

  return result;
}

//====================================================================
int Test_alt_Accel::test_Eigenvalue_Clover_Lanczos()
{
  Eigenvalue_alt<AField<double, IMPL> > test_wilson;
  // Eigenvalue_alt<IMPL> test_wilson;
  string test_file = "test_alt_Eigenvalue_Lanczos_Clover.yaml";

  int result = 0;
  result += test_wilson.eigenvalue_Lanczos(test_file);

  return result;
}

//====================================================================
int Test_alt_Accel::test_Eigenvalue_Wilson_Arnoldi()
{
  Eigenvalue_alt<AField<double, IMPL> > test_wilson;
  // Eigenvalue_alt<IMPL> test_wilson;
  // string test_file = "test_alt_Eigenvalue_Arnoldi_Wilson.yaml";
  string test_file = "test_alt_Eigenvalue_Arnoldi_Wilson_CLE.yaml";

  int result = 0;
  result += test_wilson.eigenvalue_Arnoldi(test_file);

  return result;
}

//====================================================================
int Test_alt_Accel::test_Eigenvalue_Domainwall_PVprec_Arnoldi()
{
  // FP (single precision) spectrum of the Pauli-Villars + site-diagonal
  // preconditioned Moebius DW operator  (D_PV C_PV^{-1})^dagger C^{-1} D
  // [Kanamori, Chen, Matsufuru, PoS(LATTICE2024)414].
  Eigenvalue_alt<AField<double, IMPL> > test_pvprec;
  string test_file = "test_alt_Eigenvalue_Arnoldi_Domainwall_PVprec.yaml";

  int result = 0;
  result += test_pvprec.eigenvalue_Arnoldi(test_file);

  return result;
}

//====================================================================
// Read the working precision ("double" or "float") from the Test_Eigensolver
// section of the PVprec parameter file. Defaults to "float" (the RTX 3080
// research path). Lets the yaml pick FP64 vs FP32 without a recompile.
static std::string pvprec_precision(const std::string& test_file)
{
  Parameters params_all = ParameterManager::read(test_file);
  std::string prec = "float";
  // 2pt / propagator yamls use the Test_Spectrum section; the eigenvalue
  // diagnostic yaml uses Test_Eigensolver.
  if (params_all.is_set("Test_Spectrum"))
    params_all.lookup("Test_Spectrum").fetch_string("precision", prec);
  else if (params_all.is_set("Test_Eigensolver"))
    params_all.lookup("Test_Eigensolver").fetch_string("precision", prec);
  return prec;
}

//====================================================================
int Test_alt_Accel::test_Solve_Domainwall_PVprec_BiCGStab()
{
  // Solve A x = b on the positive MG target operator
  // A = (D_PV C_PV^{-1})^dagger C^{-1} D  [Kanamori, Chen, Matsufuru].
  // (GMRES(m); BiCGStab breaks down on the non-normal A.)
  string test_file = "test_alt_Eigenvalue_Arnoldi_Domainwall_PVprec.yaml";
  int result = 0;
  if (pvprec_precision(test_file) == "double") {
    Spectrum_Domainwall_PVprec_alt<AField<double, IMPL> > t;
    result += t.solve_A(test_file);
  } else {
    Spectrum_Domainwall_PVprec_alt<AField<float, IMPL> > t;
    result += t.solve_A(test_file);
  }
  return result;
}

//====================================================================
int Test_alt_Accel::test_Solve_Domainwall_PVprec_Propagator()
{
  // Solve the physical propagator D psi = eta via the preconditioned operator
  // A = (D_PV C_PV^{-1})^dag C^{-1} D, and certify ||D psi - eta||/||eta||.
  // psi == D^{-1} eta == CGNE solution -> identical 2pt functions.
  // Precision (double/float) is selected by the yaml "precision" field.
  string test_file = "test_alt_Spectrum_Domainwall_PVprec_Hadron2ptFunction.yaml";
  int result = 0;
  if (pvprec_precision(test_file) == "double") {
    Spectrum_Domainwall_PVprec_alt<AField<double, IMPL> > t;
    result += t.solve_propagator(test_file);
  } else {
    Spectrum_Domainwall_PVprec_alt<AField<float, IMPL> > t;
    result += t.solve_propagator(test_file);
  }
  return result;
}

//====================================================================
int Test_alt_Accel::test_Solve_Domainwall_PVprec_2pt()
{
  // Full hadron 2pt function with every quark propagator solved through the
  // preconditioned operator A (psi = D^{-1} eta exactly), i.e. the SAME 2pt a
  // direct CGNE measurement gives. Precision selected by the yaml "precision".
  string test_file = "test_alt_Spectrum_Domainwall_PVprec_Hadron2ptFunction.yaml";
  int result = 0;
  if (pvprec_precision(test_file) == "double") {
    Spectrum_Domainwall_PVprec_alt<AField<double, IMPL> > t;
    result += t.solve_2pt(test_file);
  } else {
    Spectrum_Domainwall_PVprec_alt<AField<float, IMPL> > t;
    result += t.solve_2pt(test_file);
  }
  return result;
}

//====================================================================
int Test_alt_Accel::test_Eigenvalue_Staggered_Lanczos()
{
  Eigenvalue_alt<AField<double, IMPL> > test_staggered;
  // Eigenvalue_alt<IMPL> test_staggered;
  string test_file = "test_alt_Eigenvalue_Lanczos_Staggered.yaml";

  int result = 0;
  result += test_staggered.eigenvalue_Lanczos(test_file);

  return result;
}

//====================================================================
int Test_alt_Accel::test_Eigenvalue_Staggered_Arnoldi()
{
  Eigenvalue_alt<AField<double, IMPL> > test_staggered;
  // Eigenvalue_alt<IMPL> test_staggered;
  // string test_file = "test_alt_Eigenvalue_Arnoldi_Staggered.yaml";
  string test_file = "test_alt_Eigenvalue_Arnoldi_Staggered_CLE.yaml";

  int result = 0;
  result += test_staggered.measure(test_file, "test");
  // result += test_staggered.measure(test_file, "job");

  return result;
}

//====================================================================
int Test_alt_Accel::test_Spectrum_AdjointWilson_Lanczos()
{
  Spectrum_AdjointWilson_alt<IMPL> test_wilson;

  string test_file =
          "test_alt_Spectrum_AdjointWilson_Hadron2ptFunction.yaml";

  int result = 0;
  result += test_wilson.hadron_2ptFunction(test_file, "org");

  return result;
}

//====================================================================
int Test_alt_Accel::test_Eigenvalue_AdjointWilson_Lanczos()
{
  Spectrum_AdjointWilson_alt<IMPL> test_wilson;

  string test_file =
              "test_alt_Spectrum_AdjointWilson_Eigenspectrum.yaml";

  int result = 0;
  result += test_wilson.eigenspectrum(test_file, "org");

  return result;
}


#endif
//============================================================END=====
