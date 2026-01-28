/*!
      @File    spectrum_Domainwall_alt.cpp
      @brief
      @author  Hideo Matsufuru
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-03-22 13:04:01 #$
      @version $LastChangedRevision: 2503 $
*/

#include "spectrum_Domainwall_alt.h"
//#include "spectrum_alt.inc"
#include "test.h"

#include "lib/Field/field_F.h"
#include "lib/Field/field_G.h"
#include "lib/Fopr/afopr.h"
#include "lib/Fopr/fopr_Domainwall.h"
#include "lib/Measurements/Fermion/fprop_Standard_lex.h"
#include "lib/Tools/gammaMatrixSet.h"
#include "lib/Tools/randomNumberManager.h"
#include "lib/Solver/solver.h"
#include "lib/Measurements/Fermion/source.h"
#include "lib/Measurements/Fermion/fprop_Standard_Precond.h"
#include "lib/Measurements/Fermion/corr2pt_4spinor.h"

#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_lex.h"
#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_eo.h"
#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_lex_Mixedprec.h"

// class name
template<Impl IMPL>
const std::string Spectrum_Domainwall_alt<IMPL>::class_name =
                                    "Spectrum_Domainwall_alt";

//====================================================================
template<Impl IMPL>
void Spectrum_Domainwall_alt<IMPL>::init()
{
  // do nothing.
}

//====================================================================
template<Impl IMPL>
int Spectrum_Domainwall_alt<IMPL>::hadron_2ptFunction(
                            std::string test_file, std::string mode)
{
  const std::string test_name = class_name + ".hadron_2ptFunction";

  vout.general("\n");
  vout.general("-------------------------------------------------"
               "-------------------\n");
  vout.general("test name: %s\n", test_name.c_str());
  vout.general("test file: %s\n", test_file.c_str());
  vout.general("test mode: %s\n", mode.c_str());

  // parameter setup.
  int Nc   = CommonParameters::Nc();
  int Nd   = CommonParameters::Nd();
  int Ndim = CommonParameters::Ndim();
  int Nvol = CommonParameters::Nvol();

  params_all = ParameterManager::read(test_file);

  Parameters params_test = params_all.lookup("Test_Spectrum");
  Parameters params_fopr = params_all.lookup("Fopr");
  Parameters params_fopr_ref = params_all.lookup("Fopr_ref");
  Parameters params_source = params_all.lookup("Source");

  const string str_vlevel       = params_test.get_string("verbose_level");
  Bridge::VerboseLevel m_vl = vout.set_verbose_level(str_vlevel);

  const bool   do_check        = params_test.is_set("expected_result");
  const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

  const string str_fopr_type   = params_fopr.get_string("fermion_type");
  const string str_source_type = params_source.get_string("source_type");

  vout.general(m_vl, "  vlevel       = %s\n", str_vlevel.c_str());
  vout.general(m_vl, "  source_type  = %s\n", str_source_type.c_str());

  RandomNumberManager::initialize("Mseries", 1234567UL);

  // configuration setup.
  U.reset(new Field_G(Nvol, Ndim));
  setup_config(U, params_test);

  Parameters params_gfix = params_all.lookup("GaugeFixing");
  gauge_fixing(U, params_gfix);


  // Reference fermion operator.
  string gmset_type  = params_fopr.get_string("gamma_matrix_type");
  vout.general(m_vl, "  gmset_type   = %s\n", gmset_type.c_str());
  unique_ptr<GammaMatrixSet> gmset(GammaMatrixSet::New(gmset_type));

  // domain-wall operator
  //unique_ptr<Fopr> fopr(new Fopr_Domainwall(gmset_type));
  string fopr_ref_type = params_fopr_ref.get_string("fermion_type");
  //unique_ptr<Fopr> fopr(Fopr::New(fopr_ref_type, gmset_type));
  unique_ptr<Fopr> fopr(Fopr::New(fopr_ref_type, params_fopr_ref));

  fopr->set_parameters(params_fopr);
  fopr->set_config(U.get());

  // kernel operator for 4d <--> 5d conversion
  std::string kernel_type;
  double M0;
  double coeff_c;
  int err = params_fopr.fetch_string("kernel_type", kernel_type);

  if(err > 0){
    vout.crucial(m_vl, "%s: Error: kernel_type is not specified.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
  err += params_fopr.fetch_double("domain_wall_height", M0);
  if(err > 0){
    vout.crucial(m_vl, "Error at %s: domain_wall_height is not specified.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
  err += params_fopr.fetch_double("coefficient_c", coeff_c);
  if(err > 0){
    vout.crucial(m_vl, "Error at %s: coefficient_c is not specified.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  Parameters params_kernel = params_fopr;
  double kappa = 1.0/(8.0 - 2.0 * M0);
  double one_over_2kappa = 4.0 - M0;
  params_kernel.set_double("hopping_parameter", kappa);

  unique_ptr<Fopr> foprw(Fopr::New(kernel_type, params_kernel));
  foprw->set_mode("D");
  foprw->set_config(U.get());



  unique_ptr<Source> source(Source::New(str_source_type));
  source->set_parameters(params_source);

  unique_ptr<Timer> timer(new Timer(test_name));

  // setup fermion propgator computer.
  Parameters params_solver = params_all.lookup("Solver");

  //  unique_ptr<Solver> solver; // only used for original Fprop.

  // check(fopr, params_fopr);
  // vout.general(m_vl, "  check finished\n\n");
  //return 0;

  unique_ptr<Fprop> fprop;

  if(mode == "org"){
    fprop.reset(new Fprop_Standard_Precond(fopr.get()));
    fprop->set_mode("D");

  }else if(mode == "double"){
    fprop.reset(new Fprop_alt_Standard_lex<AField<double,IMPL> >(
                                     params_fopr, params_solver) );
    fprop->set_mode("D");

  }else if(mode == "double_prec"){
    fprop.reset(new Fprop_alt_Standard_lex<AField<double,IMPL> >(
                                     params_fopr, params_solver) );
    fprop->set_mode("D_prec");

  }else if(mode == "double_eo"){
    fprop.reset(new Fprop_alt_Standard_eo<AField<double,IMPL> >(
                                     params_fopr, params_solver) );
    fprop->set_mode("D");

  }else if(mode == "float"){
    fprop.reset(new Fprop_alt_Standard_lex<AField<float,IMPL> >(
                                     params_fopr, params_solver) );
    fprop->set_mode("D_prec");

  }else if(mode == "float_eo"){
    fprop.reset(new Fprop_alt_Standard_eo<AField<float,IMPL> >(
                                     params_fopr, params_solver) );
    fprop->set_mode("D");

    /*
  }else if(mode == "mixed"){
    fprop.reset(new Fprop_Standard_lex_alt_Mixedprec(
                                     params_fopr, params_solver) );
    */
  }else{
    vout.crucial(m_vl, "%s: irrelevant mode =%s\n",
                 class_name.c_str(), mode.c_str());
    exit(EXIT_FAILURE);
  }

  fprop->set_config(U.get());

  //  fprop->mult_performance("DdagD_prec", 100);
  fprop->mult_performance("DdagD", 100);

  // setup source and progator.
  int Ns = params_fopr.get_int("extent_of_5th_dimension");
  params_fopr.fetch_int("extent_of_5th_dimension", Ns);

  std::vector<Field_F> sq(Nc * Nd);
  for (int i = 0; i < Nc * Nd; ++i) {
    sq[i].set(0.0);
  }

  Field_F b, bt;
  Field_F vt1, vt2;
  Field_F b5(Nvol, Ns), x5(Nvol, Ns), y5(Nvol, Ns);
  b.set(0.0);
  source->set(b, 0);

  // Solver main part.
  timer->start();

  int nconv;
  double diff;

  vout.general(m_vl, "\n");
  vout.general(m_vl, "Solving quark propagator:\n");
  vout.general(m_vl, "  color spin   Nconv      diff           diff2\n");


  for (int ispin = 0; ispin < Nd; ++ispin) {
    for (int icolor = 0; icolor < Nc; ++icolor) {
      int idx = icolor + Nc * ispin;
      source->set(b, idx);

      // set 5d source as
      //   D_- ( P_+ 0 ....0 P_-)^T b4
#pragma omp parallel
      {
        // build 5dim vector
        b5.set(0.0);

        // s5=0
        foprw->mult_gm5(vt1, b);
        axpy(vt1, +1.0, b);       // (1+gm5)ab
        foprw->mult(vt2, vt1);
        aypx(-one_over_2kappa * coeff_c, vt2, vt1); // (-cD+1) (1+gm5)ab
        scal(vt2, -0.5);          // 0.5*(cD-1) (1+gm5)ab
        copy(b5, 0, vt2, 0);

        foprw->mult_gm5(vt1, b);
        axpy(vt1, -1.0, b);      // (-1+gm5)ab  [ = -(1-gm5)ab ]
        foprw->mult(vt2, vt1);
        aypx(-one_over_2kappa * coeff_c, vt2, vt1); // (cD-1) (1-gm5)ab
        scal(vt2, 0.5);          // 0.5*(cD-1) (1-gm5)ab
        copy(b5, Ns-1, vt2, 0);
      }

      fprop->invert(x5, b5, nconv, diff);

#pragma omp parallel
      {
        {// check the solution
          copy(y5, b5);
          fopr->set_mode("D");
          fopr->mult(y5, x5);
          axpy(y5, -1.0, b5);
          double diff2 = y5.norm2() / b5.norm2();
          vout.general(m_vl, "   %2d   %2d   %6d   %12.4e   %12.4e\n",
                       icolor, ispin, nconv, diff, diff2);

        }
      }
#pragma omp parallel
      {
        // convert to 4D propagator
        copy(vt1, 0, x5, 0);
        copy(vt2, 0, x5, Ns-1);
        copy(sq[idx], vt1);
#pragma omp barrier

        axpy(sq[idx], 1.0, vt2);    // sq[idx] = (x5[0]+x5[Ls-1])
        axpy(vt1, -1.0, vt2);       //     vt1 = (x5[0]-x5[Ls-1])
#pragma omp barrier
        foprw->mult_gm5(vt2, vt1);  //     vt2 = gm5*(x5[0]-x5[Ls-1])
        axpy(sq[idx], -1.0, vt2);   // sq[idx] = (1-gm5) x5[0] + (1+gm5) x5[Ls-1]
        scal(sq[idx], 0.5);
      }
    }
  }

  fprop->report_performance();


  //- meson correlators
  Corr2pt_4spinor corr(gmset.get());
  corr.set_parameters(params_all.lookup("Corr2pt_4spinor"));

  vout.general(m_vl, "\n");
  vout.general(m_vl, "2-point correlator:\n");

  double result = corr.meson_all(sq, sq);

  timer->report();

  RandomNumberManager::finalize();

  vout.general(m_vl, " alt code.\n");


  if (do_check) {
    return Test::verify(result, expected_result);
  } else {
    vout.detailed(m_vl, "check skipped: expected_result not set.\n\n");
    return EXIT_SKIP;
  }

}

//====================================================================
template<Impl IMPL>
void Spectrum_Domainwall_alt<IMPL>::check(
                                    unique_ptr<Fopr>& fopr_ref,
                                    Parameters& params_fopr)
{
  //  static const Impl IMPL = OPENACC;

  typedef AField<double,IMPL> AFIELD;

  vout.general(m_vl, "\n");
  vout.general(m_vl, "Check of fermion operator\n");

  int Nvol = CommonParameters::Nvol();
  int Ns = params_fopr.get_int("extent_of_5th_dimension");

  fopr_ref->set_config(U.get());
  fopr_ref->set_mode("D");

  Field_F b(Nvol, 1), bt(Nvol, 1);
  b.set(0.0);
  b.set(0, 1.0);

  Field_F b5(Nvol, Ns), y5(Nvol, Ns);

  // set 5D source vector                                                   
  b5.set(0.0);
  fopr_ref->mult_gm5(bt, b);
  axpy(b5, Ns-1,  0.5, b,  0);
  axpy(b5, Ns-1, -0.5, bt, 0);
  axpy(b5, 0,     0.5, b,  0);
  axpy(b5, 0,     0.5, bt, 0);

  // extend source vector
  for(int i = 0; i < 10; ++i){
    fopr_ref->mult(y5, b5);
    fopr_ref->mult(b5, y5);
  }

  double bb = b5.norm2();
  bb = 1.0/sqrt(bb);
  scal(b5, bb);
  vout.general("norm of source vector: %f\n", b5.norm());

  //### from here, check main part ###

  // construction of alternative fermion object.
  string fopr_type = params_fopr.get_string("fermion_type");
  //  fopr_type += "_5din";
  //  fopr_type += "_eo";
  unique_ptr<AFopr<AFIELD> >
               fopr_alt(AFopr<AFIELD>::New(fopr_type, params_fopr));

  fopr_alt->set_config(U.get());

  // check
  Field_F x5(Nvol, Ns), z5(Nvol, Ns);

  int nin  = fopr_alt->field_nin();
  int nex  = fopr_alt->field_nex();
  AFIELD abq(nin, Nvol, nex);
  AFIELD axq(nin, Nvol, nex);
  AFIELD ayq(nin, Nvol, nex);

  AIndex_lex<double,IMPL> index_alt;

  if(fopr_alt->needs_convert()){
    vout.general(m_vl, "convert required.\n");
    fopr_alt->convert(abq, (Field&)b5);
    vout.general(m_vl, "convert performed.\n");
  }else{
    convert_spinor(index_alt, abq, (Field&)b5);
    vout.general(m_vl, "convert does not performed.\n");
  }

  // check of lexical operator from here ...
  std::string mode = "D";
  //  mode = "DdagD_prec";
  fopr_ref->set_mode(mode);
  fopr_alt->set_mode(mode);

  fopr_alt->mult(axq, abq);
  fopr_ref->mult(y5, b5);
  /*
  for(int i = 0; i < 1; ++i){
    fopr_alt->mult(ayq, axq);
    fopr_alt->mult(axq, ayq);
    fopr_ref->mult(x5, y5);
    fopr_ref->mult(y5, x5);
  }
  */
  // ... to here

  // check of even-odd operator from here ... 
  /*
  int nvol2 = Nvol/2;
  AFIELD be(nin, nvol2, nex), bo(nin, nvol2, nex);
  AFIELD xe(nin, nvol2, nex), xo(nin, nvol2, nex);
  AFIELD xt(nin, nvol2, nex);

  AIndex_eo<double, IMPL> index_eo;
  index_eo.split(be, bo, abq);

  fopr_alt->mult(xe, be, "Dee");
  fopr_alt->mult(xt, bo, "Deo");
  axpy(xe, 1.0, xt);
  fopr_alt->mult(xo, bo, "Doo");
  fopr_alt->mult(xt, be, "Doe");
  axpy(xo, 1.0, xt);

  index_eo.merge(axq, xe, xo);

  fopr_ref->mult(y5, b5);
  */
  // ... to here

  if(fopr_alt->needs_convert()){
    fopr_alt->convert(ayq, (Field&)y5);
  }else{
    convert_spinor(index_alt, ayq, (Field&)y5);
  }

  double y5norm = y5.norm2();
  double yynorm = ayq.norm2();
  double xxnorm = axq.norm2();
  axpy(ayq, double(-1.0), axq);

  double yy = ayq.norm2();
  //yy = yy/yynorm;

  vout.general(m_vl, "\n");
  vout.general(m_vl, "norm2(y5) = %f\n", y5norm);
  vout.general(m_vl, "norm2(ayq) = %f\n", yynorm);
  vout.general(m_vl, "norm2(axq) = %f\n", xxnorm);
  vout.general(m_vl, "diff2 = %f\n", yy);

  /*
  reverse_spinor(index_alt, (Field&)x5, axq);
  axpy(y5, -1.0, x5);
  double yy2 = y5.norm2();

  vout.general(m_vl, "diff2 = %f\n", yy);
  */

}

//============================================================END=====
