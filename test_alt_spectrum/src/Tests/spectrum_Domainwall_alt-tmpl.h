/*!
      @File    spectrum_Domainwall_alt.cpp
      @brief
      @author  Hideo Matsufuru
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2025-01-15 22:38:16 #$
      @version $LastChangedRevision: 2630 $
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
#include "lib/Tools/randomNumbers_Mseries.h"
#include "lib/Solver/solver.h"
#include "lib/Measurements/Fermion/source.h"
#include "lib/Measurements/Fermion/fprop_Standard_Precond.h"
#include "lib/Measurements/Fermion/corr2pt_4spinor.h"
#include "lib/Smear/director_Smear.h"
#include "lib/Measurements/Gauge/staple_lex.h"

#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_lex.h"
#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_eo.h"
#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_lex_Mixedprec.h"


double calc_plaquette(const Field_G *u){
  // parepare parameter
  Parameters params_staple;
  params_staple.set_string("verbose_level", "general");
  Timer timer("calc_plqauette");
  timer.start();

  // calculate plaquette
  const unique_ptr<Staple> staple(Staple::New("Lexical"));
  staple->set_parameters(params_staple);
  const double result = staple->plaquette(*u);

  // done!
  timer.stop();
  timer.report();
  return result;
}

// template<typename AFIELD>
void apply_smearing(Field_G *u, const Parameters& params){
  const string str_smear_type = params.get_string("smear_type");
  const string str_proj_type  = params.get_string("projection_type");
  //    const double rho_uniform    = params.get_double("rho_uniform");
  //    const int number_of_smearing= params.get_int("number_of_smearing");
  //    const string str_vlevel     = params.get_string("verbose_level");

  const Parameters &params_proj=params;
  const Parameters &params_smear=params;
  const Parameters &params_dr_smear=params;
  // projection
  unique_ptr<Projection> proj(Projection::New(str_proj_type));
  proj->set_parameters(params_proj);

  // smear
  unique_ptr<Smear> smear(Smear::New(str_smear_type, proj.get()));
  smear->set_parameters(params_smear);

  // director
  const unique_ptr<Director_Smear> dr_smear(new Director_Smear(smear.get()));
  dr_smear->set_parameters(params_dr_smear);
  dr_smear->set_config(u);

  // run smearing
  const int     Nsmear  = dr_smear->get_Nsmear();
  const Field_G *Usmear = (Field_G *)dr_smear->getptr_smearedConfig(Nsmear);
  copy(*u, *Usmear);

  // done!
}

// TWDIAG helper: apply one named eo sub-operator to the odd half of a source,
// reverse the result to a double Field_F. Works for any AField type / mw mode.
template<typename AF, Impl AIMPL>
static Field_F twdiag_subop(AFopr<AF>* fopr, const Field& src,
                            int Nvol, int Ns, int nex, int mwi, const char* op)
{
  typedef typename AF::real_t real_t;
  AIndex_eo<real_t, AIMPL> index_eo;
  const int nin2  = fopr->field_nin();
  const int nvol2 = Nvol / 2;
  AF aq(nin2, Nvol, nex);
  fopr->convert(aq, src);
  AF be(nin2, nvol2, nex), bo(nin2, nvol2, nex);
  AF out(nin2, nvol2, nex), zero(nin2, nvol2, nex);
  index_eo.split(be, bo, aq, mwi);
  zero.set(0.0);
  fopr->mult(out, bo, op);
  AF axq(nin2, Nvol, nex);
  index_eo.merge(axq, out, zero, mwi);
  Field_F z(Nvol, Ns);
  fopr->reverse((Field&)z, axq);
  return z;
}

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

  Parameters params_test     = params_all.lookup("Test_Spectrum");
  Parameters params_fopr     = params_all.lookup("Fopr");
  Parameters params_fopr_ref = params_all.lookup("Fopr_ref");
  Parameters params_source   = params_all.lookup("Source");

  const string         str_vlevel = params_test.get_string("verbose_level");
  Bridge::VerboseLevel m_vl       = vout.set_verbose_level(str_vlevel);

  const bool   do_check        = params_test.is_set("expected_result");
  const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

  const string str_fopr_type   = params_fopr.get_string("fermion_type");
  const string str_source_type = params_source.get_string("source_type");

  vout.general(m_vl, "  vlevel       = %s\n", str_vlevel.c_str());
  vout.general(m_vl, "  source_type  = %s\n", str_source_type.c_str());
  // for smearing
  const Parameters params_smearing   = params_all.lookup("Smearing");
  const string str_proj_type         = params_smearing.get_string("projection_type");
  const string str_smear_type        = params_smearing.get_string("smear_type");
  const double rho_uniform           = params_smearing.get_double("rho_uniform");
  const int number_of_smearing       = params_smearing.get_int("number_of_smearing");

  //RandomNumberManager::initialize("Mseries", 1234567UL);
  vout.general(m_vl, "Smearing parameters (assumes APE-Stout):\n");
  vout.general(m_vl, "   proj_type:  %s\n", str_proj_type.c_str());
  vout.general(m_vl, "   smear_type: %s\n", str_smear_type.c_str());
  vout.general(m_vl, "   rho_uniform: %f\n", rho_uniform);
  vout.general(m_vl, "   number_of_smearing: %d\n", number_of_smearing);
  vout.general(m_vl, "  vlevel       = %s\n", str_vlevel.c_str());

  if(str_smear_type == "None"){
    vout.general(m_vl, "skip smearing\n");
    return;
  }
  // configuration setup.
  U.reset(new Field_G(Nvol, Ndim));
  setup_config(U, params_test);

  Parameters params_gfix = params_all.lookup("GaugeFixing");
  gauge_fixing(U, params_gfix);

  double plaq=calc_plaquette(U.get());
  vout.general(m_vl, "plaquette (before smearing): %23.15f\n", plaq);
  vout.general(m_vl, "applying smearing\n");
  Timer timer_smearing("smearing");
  timer_smearing.reset();
  timer_smearing.start();

  apply_smearing(U.get(), params_smearing);
  double plaq_smear=calc_plaquette(U.get());

  timer_smearing.stop();
  timer_smearing.report();

  vout.general(m_vl, "plaquette (after smearing):  %23.15f\n", plaq_smear);

  // Reference fermion operator.
  string gmset_type = params_fopr.get_string("gamma_matrix_type");
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
  double      M0;
  double      coeff_c;
  int         err = params_fopr.fetch_string("kernel_type", kernel_type);

  if (err > 0) {
    vout.crucial(m_vl, "%s: Error: kernel_type is not specified.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
  err += params_fopr.fetch_double("domain_wall_height", M0);
  if (err > 0) {
    vout.crucial(m_vl, "Error at %s: domain_wall_height is not specified.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
  err += params_fopr.fetch_double("coefficient_c", coeff_c);
  if (err > 0) {
    vout.crucial(m_vl, "Error at %s: coefficient_c is not specified.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  Parameters params_kernel   = params_fopr;
  double     kappa           = 1.0 / (8.0 - 2.0 * M0);
  double     one_over_2kappa = 4.0 - M0;
  params_kernel.set_double("hopping_parameter", kappa);

  unique_ptr<Fopr> foprw(Fopr::New(kernel_type, params_kernel));
  foprw->set_mode("D");
  foprw->set_config(U.get());



  unique_ptr<Source> source(Source::New(str_source_type));
  source->set_parameters(params_source);

  unique_ptr<Timer> timer(new Timer(test_name));

  // setup fermion propgator computer.
  Parameters params_solver = params_all.lookup("Solver");

  // To confirm the alt operator works same as corelib.
  if(mode == "double"){
    check_operator<AField<double, IMPL> >(fopr, params_fopr);
  }else if (mode == "double_eo"){
    check_operator_eo<AField<double, IMPL> >(fopr, params_fopr);
  }else if (mode == "float"){
    check_operator<AField<float, IMPL> >(fopr, params_fopr);
  }else if (mode == "float_eo"){
    check_operator_eo<AField<float, IMPL> >(fopr, params_fopr);
  }else if (mode == "double_prec"){
    check_operator<AField<double, IMPL> >(fopr, params_fopr);
#ifdef USE_FP16
  }else if (mode == "half"){
    check_operator<AField<half, IMPL> >(fopr, params_fopr);
  }else if (mode == "half_eo"){
    check_operator_eo<AField<half, IMPL> >(fopr, params_fopr);
#endif
  }

  //  return 0;
  //  unique_ptr<Solver> solver; // only used for original Fprop.

  unique_ptr<Fprop> fprop;

  if (mode == "org") {
    fprop.reset(new Fprop_Standard_Precond(fopr.get()));
    fprop->set_mode("D");
  } else if (mode == "double") {
    fprop.reset(new Fprop_alt_Standard_lex<AField<double, IMPL> >(
                  params_fopr, params_solver));
    fprop->set_mode("D");
  } else if (mode == "double_prec") {
    fprop.reset(new Fprop_alt_Standard_lex<AField<double, IMPL> >(
                  params_fopr, params_solver));
    fprop->set_mode("D_prec");
  } else if (mode == "double_eo") {
    fprop.reset(new Fprop_alt_Standard_eo<AField<double, IMPL> >(
                  params_fopr, params_solver));
    fprop->set_mode("D");
  } else if (mode == "float") {
    params_solver.set_double("convergence_criterion_squared", 1.0e-14);
    fprop.reset(new Fprop_alt_Standard_lex<AField<float, IMPL> >(
                  params_fopr, params_solver));
    fprop->set_mode("D");
  } else if (mode == "float_eo") {
    // TW operator now runs GENUINE triple-word. Use the SAME convergence
    // criterion as the yaml (and as double_eo) so the float-triple vs FP64
    // comparison is fair: previously this was hardcoded to 1e-20 (rel.resid
    // ~1e-10), which capped the QTW solve 6 orders above the FP64 floor.
    fprop.reset(new Fprop_alt_Standard_eo<AField<float, IMPL> >(
                  params_fopr, params_solver));
    fprop->set_mode("D");
#ifdef USE_FP16
  } else if (mode == "half") {
    params_solver.set_double("convergence_criterion_squared", 1.0e-7);
    fprop.reset(new Fprop_alt_Standard_lex<AField<half, IMPL> >(
                  params_fopr, params_solver));
    fprop->set_mode("D");
  } else if (mode == "half_eo") {
    params_solver.set_double("convergence_criterion_squared", 1.0e-7);
    fprop.reset(new Fprop_alt_Standard_eo<AField<half, IMPL> >(
                  params_fopr, params_solver));
    fprop->set_mode("D");
#endif
    /*
  }else if(mode == "mixed"){
    fprop.reset(new Fprop_Standard_lex_alt_Mixedprec(
                                     params_fopr, params_solver) );
    */
  } else {
    vout.crucial(m_vl, "%s: irrelevant mode =%s\n",
                 class_name.c_str(), mode.c_str());
    exit(EXIT_FAILURE);
  }

  fprop->set_config(U.get());
  int measure_times = 200;

  if (mode == "double_eo" || mode == "float_eo") {
    fprop->mult_performance("Dee_inv", measure_times);
    fprop->mult_performance("Doe",     measure_times);
    fprop->mult_performance("D",       measure_times);
    fprop->mult_performance("Ddag",    measure_times);
    fprop->mult_performance("DdagD",   measure_times);
  }else{
    fprop->mult_performance("D",       measure_times);
    fprop->mult_performance("Ddag",    measure_times);
    fprop->mult_performance("DdagD",   measure_times);
  }

  // setup source and progator.
  int Ns = params_fopr.get_int("extent_of_5th_dimension");
  params_fopr.fetch_int("extent_of_5th_dimension", Ns);

  // boundary-projected propagator sq (physical surface field) and the
  // MIDPOINT-projected sq_mid (5th-dim midpoint field) for the residual mass.
  std::vector<Field_F> sq(Nc * Nd), sq_mid(Nc * Nd);
  for (int i = 0; i < Nc * Nd; ++i) {
    sq[i].set(0.0);
    sq_mid[i].set(0.0);
  }
  const int s_mid_lo = Ns / 2 - 1;   // lower-half top slice
  const int s_mid_hi = Ns / 2;       // upper-half bottom slice

  Field_F b, bt;
  Field_F vt1, vt2;
  Field_F b5(Nvol, Ns), x5(Nvol, Ns), y5(Nvol, Ns);
  b.set(0.0);
  source->set(b, 0);

  // Solver main part.
  timer->start();

  int    nconv;
  double diff;

  vout.general(m_vl, "\n");
  vout.general(m_vl, "Solving quark propagator:\n");
  vout.general(m_vl, "  color spin   Nconv      diff           diff2\n");


    for (int ispin = 0; ispin < Nd; ++ispin) {
      for (int icolor = 0; icolor < Nc; ++icolor) {
        int idx = icolor + Nc * ispin;
        source->set(b, idx);

        // set 5d source as
        // D_- ( P_+ 0 ....0 P_-)^T b4
  #pragma omp parallel
        {
          // build 5dim vector
          b5.set(0.0);

          // s5=0
          foprw->mult_gm5(vt1, b);
          axpy(vt1, +1.0, b);                         // (1+gm5)ab
          foprw->mult(vt2, vt1);
          aypx(-one_over_2kappa * coeff_c, vt2, vt1); // (-cD+1) (1+gm5)ab
          scal(vt2, -0.5);                            // 0.5*(cD-1) (1+gm5)ab
          copy(b5, 0, vt2, 0);

          foprw->mult_gm5(vt1, b);
          axpy(vt1, -1.0, b);                         // (-1+gm5)ab  [ = -(1-gm5)ab ]
          foprw->mult(vt2, vt1);
          aypx(-one_over_2kappa * coeff_c, vt2, vt1); // (cD-1) (1-gm5)ab
          scal(vt2, 0.5);                             // 0.5*(cD-1) (1-gm5)ab
          copy(b5, Ns - 1, vt2, 0);

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
          copy(vt2, 0, x5, Ns - 1);
          copy(sq[idx], vt1);
  #pragma omp barrier
          axpy(sq[idx], 1.0, vt2);    // sq[idx] = (x5[0]+x5[Ls-1])
          axpy(vt1, -1.0, vt2);       //     vt1 = (x5[0]-x5[Ls-1])
  #pragma omp barrier
          foprw->mult_gm5(vt2, vt1);  //     vt2 = gm5*(x5[0]-x5[Ls-1])
          axpy(sq[idx], -1.0, vt2);   // sq[idx] = (1-gm5) x5[0] + (1+gm5) x5[Ls-1]
          scal(sq[idx], 0.5);
        }

        // midpoint-projected propagator (for the residual mass): SAME
        // reconstruction at the s = Ns/2-1 / Ns/2 cut.
  #pragma omp parallel
        {
          copy(vt1, 0, x5, s_mid_hi);
          copy(vt2, 0, x5, s_mid_lo);
          copy(sq_mid[idx], vt1);
  #pragma omp barrier
          axpy(sq_mid[idx], 1.0, vt2);   // (x5[Ns/2]+x5[Ns/2-1])
          axpy(vt1, -1.0, vt2);          // (x5[Ns/2]-x5[Ns/2-1])
  #pragma omp barrier
          foprw->mult_gm5(vt2, vt1);     // gm5 (x5[Ns/2]-x5[Ns/2-1])
          axpy(sq_mid[idx], -1.0, vt2);  // (1-gm5)x5[Ns/2] + (1+gm5)x5[Ns/2-1]
          scal(sq_mid[idx], 0.5);
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

  vout.general(m_vl, "RESULT_HIPREC: %.17e\n", result);

  // ---- residual mass --------------------------------------------------------
  // m_res(t) = <J5q(t) P(0)> / <P(t) P(0)> = C_PP(sq_mid,sq_mid)/C_PP(sq,sq),
  // both pseudoscalar (g5 x g5) correlators from the same boundary point source.
  // Plateau (large-t) value is the residual mass. CGNE/reference path -> serves
  // as the trusted cross-check for the A-based solve_2pt m_res on the same setup.
  {
    GammaMatrix gm5 = gmset->get_GM(gmset->GAMMA5);
    const int Lt = CommonParameters::Lt();
    std::vector<dcomplex> c_PP(Lt), c_J5q(Lt);   // meson_correlator needs them pre-sized
    corr.meson_correlator(c_PP,  gm5, gm5, sq,     sq);
    corr.meson_correlator(c_J5q, gm5, gm5, sq_mid, sq_mid);
    vout.general(m_vl, "\n");
    vout.general(m_vl, "Residual mass (J5q/PP ratio):\n");
    vout.general(m_vl, "    t      C_PP(t)         C_J5q(t)        m_res(t)\n");
    double mres_sum = 0.0; int mres_n = 0;
    for (int t = 0; t < Lt; ++t) {
      double cpp  = real(c_PP[t]);
      double cj5q = real(c_J5q[t]);
      double r    = (cpp != 0.0) ? cj5q / cpp : 0.0;
      vout.general(m_vl, " %4d  %15.7e  %15.7e  %15.7e\n", t, cpp, cj5q, r);
      if (t >= Lt / 4 && t < 3 * Lt / 4) { mres_sum += r; ++mres_n; }
    }
    double m_res = (mres_n > 0) ? mres_sum / mres_n : 0.0;
    vout.general(m_vl, "RESULT_HIPREC_mres (plateau avg over t in [%d,%d)): %.10e\n",
                 Lt / 4, 3 * Lt / 4, m_res);
  }
  // ---------------------------------------------------------------------------

  timer->report();
  //RandomNumberManager::finalize();
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
template<typename AFIELD>
int Spectrum_Domainwall_alt<IMPL>::check_operator(
                                        unique_ptr<Fopr>& fopr_ref,
                                        Parameters& params_fopr)
{
  typedef typename AFIELD::real_t real_t;

  vout.general(m_vl, "\n");
  vout.general(m_vl, "Check of fermion operator (lexical)\n");

  int Nvol = CommonParameters::Nvol();
  int Ns   = params_fopr.get_int("extent_of_5th_dimension");

  fopr_ref->set_config(U.get());
  fopr_ref->set_mode("D");

  Field_F b(Nvol, 1), bt(Nvol, 1);
  b.set(0.0);
  b.set(0, 1.0);

  Field_F b5(Nvol, Ns), y5(Nvol, Ns);

  // set 5D source vector
  b5.set(0.0);
  fopr_ref->mult_gm5(bt, b);
  axpy(b5, Ns - 1,  0.5, b,  0);
  axpy(b5, Ns - 1, -0.5, bt, 0);
  axpy(b5, 0, 0.5, b,  0);
  axpy(b5, 0, 0.5, bt, 0);

  // extend source vector
  for (int i = 0; i < 10; ++i) {
    fopr_ref->mult(y5, b5);
    fopr_ref->mult(b5, y5);

    double bb = b5.norm2();
    bb = 1.0 / sqrt(bb);
    scal(b5, bb);
  }

  double bb = b5.norm2();
  vout.general("norm of source vector: %f\n", bb);

  bb = 1.0 / sqrt(bb);
  scal(b5, bb);
  bb = b5.norm2();
  vout.general("normalized: %f\n", bb);

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

  int    nin = fopr_alt->field_nin();
  int    nex = fopr_alt->field_nex();
  AFIELD abq(nin, Nvol, nex);
  AFIELD axq(nin, Nvol, nex);
  AFIELD ayq(nin, Nvol, nex);

  AIndex_lex<real_t, IMPL> index_alt;

  if (fopr_alt->needs_convert()) {
    vout.general(m_vl, "convert required.\n");
    fopr_alt->convert(abq, (Field&)b5);
    vout.general(m_vl, "convert performed.\n");
  } else {
    convert_spinor(index_alt, abq, (Field&)b5);
    vout.general(m_vl, "convert does not performed.\n");
  }

  double abqnorm = abq.norm2();
  
  //std::string mode = "D";
  //std::string mode = "Ddag";
  std::string mode = "DdagD";
  //std::string mode = "D_prec";
  //std::string mode = "Ddag_prec";
  //std::string mode = "DdagD_prec";
  fopr_ref->set_mode(mode);
  fopr_alt->set_mode(mode);
  vout.general(m_vl, "check mode = %s\n", mode.c_str());

  fopr_alt->mult(axq, abq);
  fopr_ref->mult(y5, b5);

  //for(int i = 0; i < 5; ++i){
  for(int i = 0; i < 0; ++i){
    double aa = axq.norm2();
    scal(axq, 1.0/sqrt(aa));
    fopr_alt->mult(ayq, axq);
    fopr_alt->mult(axq, ayq);

    double yy = y5.norm2();
    scal(y5, 1.0/sqrt(yy));
    fopr_ref->mult(x5, y5);
    fopr_ref->mult(y5, x5);
  }

  if (fopr_alt->needs_convert()) {
    fopr_alt->reverse((Field&)z5, axq);
  } else {
    reverse_spinor(index_alt, (Field&)z5, axq);
  }

  double y5norm = y5.norm2();
  double z5norm = z5.norm2();

  double axqnorm = axq.norm2();

  axpy(y5, double(-1.0), z5);
  double diff2 = y5.norm2();

  vout.general(m_vl, "\n");
  vout.general(m_vl, "norm2(abq) = %f\n", abqnorm);
  vout.general(m_vl, "norm2(y5)  = %f\n", y5norm);
  vout.general(m_vl, "norm2(z5)  = %f\n", z5norm);
  vout.general(m_vl, "norm2(axq) = %f\n", axqnorm);
  vout.general(m_vl, "diff2 = %f\n", diff2);

  double epsilon = 1.e-16;
  if (diff2 < epsilon) return 0;

  return 1;

}

//====================================================================
template<Impl IMPL>
template<typename AFIELD>
int Spectrum_Domainwall_alt<IMPL>::check_operator_eo(
    unique_ptr<Fopr>& fopr_ref,
    Parameters& params_fopr)
{
  typedef typename AFIELD::real_t real_t;

  vout.general(m_vl, "\n");
  vout.general(m_vl, "Check of fermion operator (even-odd)\n");

  int Nvol = CommonParameters::Nvol();
  int Ns   = params_fopr.get_int("extent_of_5th_dimension");

  fopr_ref->set_config(U.get());
  fopr_ref->set_mode("D");

  Field_F b(Nvol, 1), bt(Nvol, 1);
  b.set(0.0);
  b.set(0, 1.0);

  Field_F b5(Nvol, Ns), y5(Nvol, Ns), z5(Nvol, Ns);
  b5.set(0.0);
  fopr_ref->mult_gm5(bt, b);
  axpy(b5, Ns - 1, 0.5, b, 0);
  axpy(b5, Ns - 1, -0.5, bt, 0);
  axpy(b5, 0, 0.5, b, 0);
  axpy(b5, 0, 0.5, bt, 0);

  for (int i = 0; i < 10; ++i) {
    fopr_ref->mult(y5, b5);
    fopr_ref->mult(b5, y5);
  }

  double bb = b5.norm2();
  bb = 1.0 / sqrt(bb);
  scal(b5, bb);
  vout.general("norm of source vector: %f\n", b5.norm());

  string fopr_type = params_fopr.get_string("fermion_type") + "_eo";
  unique_ptr<AFopr<AFIELD>> fopr_alt(AFopr<AFIELD>::New(fopr_type, params_fopr));
  fopr_alt->set_config(U.get());

  Field_F x5(Nvol, Ns);

  int nin = fopr_alt->field_nin();
  int nex = fopr_alt->field_nex();
  AFIELD abq(nin, Nvol, nex), axq(nin, Nvol, nex);

  AIndex_lex<real_t, IMPL> index_alt;
  if (fopr_alt->needs_convert()) {
    vout.general(m_vl, "convert required.\n");
    fopr_alt->convert(abq, (Field&)b5);
    vout.general(m_vl, "convert performed.\n");
  } else {
    convert_spinor(index_alt, abq, (Field&)b5);
    vout.general(m_vl, "convert not performed.\n");
  }

  double abqnorm = (double)abq.norm2();

  // --- TWDIAG: localize the structural FP-vs-TW operator gap ---
  // The FP alt operator is validated against the reference (the D check below
  // passes < 1e-16). So comparing each TW sub-operator against its FP version
  // (both reversed to double) isolates which TW kernel deviates. Run with
  // TWDIAG=1 in the environment; it prints per-sub-op rel diffs and returns.
  if (getenv("TWDIAG")) {
    vout.general(m_vl, "\n=== TWDIAG: TW sub-ops vs EXACT double sub-ops (CG operator pieces) ===\n");

    // TWDIAG_RAND: replace the (smooth, operator-evolved) source with a random
    // Gaussian vector so the sub-op comparison sees ROUGH/full-spectrum modes,
    // which is what the CG Krylov space explores. If the TW operator deviates
    // from QDW much more on a random vector than on the smooth source, the
    // operator (not the BLAS) is what stalls the CG true residual at ~1e-8.
    if (getenv("TWDIAG_RAND")) {
      RandomNumbers_Mseries rng(1234567);
      rng.gauss_lex_global((Field&)b5);
      double bn = b5.norm2(); scal(b5, 1.0 / sqrt(bn));
      vout.general(m_vl, "  [TWDIAG_RAND] using RANDOM Gaussian source vector\n");
    }

    // Two references: plain double (FP) and double-double (QDW, DW mode).
    // QDW (~1e-30) is effectively exact; double (FP) carries its own rounding.
    string fopr_type2 = params_fopr.get_string("fermion_type") + "_eo";
    unique_ptr<AFopr<AField<double, IMPL>>>
        fopr_dbl(AFopr<AField<double, IMPL>>::New(fopr_type2, params_fopr));
    fopr_dbl->set_config(U.get());           // FP double

    unique_ptr<AFopr<AField<double, IMPL>>>
        fopr_qdw(AFopr<AField<double, IMPL>>::New(fopr_type2, params_fopr));
    fopr_qdw->set_mw_mode(AFopr<AField<double, IMPL>>::MWMode::DW);
    fopr_qdw->set_config(U.get());           // double-double reference

    // TW float operator (the one the CG uses): set mode, then (re)convert gauge.
    fopr_alt->set_mw_mode(AFopr<AFIELD>::MWMode::TW);
    fopr_alt->set_config(U.get());

    // The CG operator D() uses only these: D_eo (Deo/Doe) and LU_inv (Dee_inv/Doo_inv).
    // Compare both double(FP) and TW(float-triple) against the QDW(dd) reference.
    // If |TW-QDW| << |double-QDW|, TW genuinely beats double (triple-word OK).
    const char* ops[] = {"Deo", "Doe", "Dee_inv", "Doo_inv", "Dee", "Doo"};
    for (const char* op : ops) {
      Field_F zq = twdiag_subop<AField<double, IMPL>, IMPL>(fopr_qdw.get(), (Field&)b5, Nvol, Ns, nex, 1, op);
      Field_F zd = twdiag_subop<AField<double, IMPL>, IMPL>(fopr_dbl.get(), (Field&)b5, Nvol, Ns, nex, 0, op);
      Field_F zt = twdiag_subop<AFIELD, IMPL>(fopr_alt.get(),               (Field&)b5, Nvol, Ns, nex, 2, op);
      double nq = zq.norm2();
      Field_F ed = zq; axpy(ed, -1.0, zd);   // QDW - double(FP)
      Field_F et = zq; axpy(et, -1.0, zt);   // QDW - TW
      vout.general(m_vl, "  %-8s : |QDW-double|/|QDW|=%.6e   |QDW-TW|/|QDW|=%.6e\n",
                   op, sqrt(ed.norm2()/nq), sqrt(et.norm2()/nq));
    }
    fopr_alt->set_mw_mode(AFopr<AFIELD>::MWMode::FP);
    fopr_alt->set_config(U.get());
    return 0;
  }

  // --- D check ---
  {
    int nvol2 = Nvol / 2;
    AFIELD be(nin, nvol2, nex), bo(nin, nvol2, nex);
    AFIELD xe(nin, nvol2, nex), xo(nin, nvol2, nex), xt(nin, nvol2, nex);

    AIndex_eo<real_t, IMPL> index_eo;
    index_eo.split(be, bo, abq);

    fopr_alt->mult(xe, be, "Dee");
    fopr_alt->mult(xt, xe, "Dee_inv");
    fopr_alt->mult(xe, xt, "Dee");
    fopr_alt->mult(xt, bo, "Deo");
    axpy(xe, 1.0, xt);
    fopr_alt->mult(xo, bo, "Doo");
    fopr_alt->mult(xt, xo, "Doo_inv");
    fopr_alt->mult(xo, xt, "Doo");
    fopr_alt->mult(xt, be, "Doe");
    axpy(xo, 1.0, xt);

    index_eo.merge(axq, xe, xo);

    fopr_ref->set_mode("D");
    fopr_ref->mult(y5, b5);
    if (fopr_alt->needs_convert()) {
      fopr_alt->reverse((Field&)z5, axq);
    } else {
      reverse_spinor(index_alt, (Field&)z5, axq);
    }

    double y5norm_D = y5.norm2();
    double z5norm_D = z5.norm2();
    double axqnorm_D = axq.norm2();

    axpy(y5, -1.0, z5);
    double diff2_D = y5.norm2();

    vout.general(m_vl, "\n[D check]\n");
    vout.general(m_vl, "norm2(abq) = %f\n", abqnorm);
    vout.general(m_vl, "norm2(y5)  = %f\n", y5norm_D);
    vout.general(m_vl, "norm2(z5)  = %f\n", z5norm_D);
    vout.general(m_vl, "norm2(axq) = %f\n", axqnorm_D);
    vout.general(m_vl, "diff2 = %f\n", diff2_D);

    if (diff2_D > 1e-16) return 1;
  }

  // --- Ddag check ---
  {
    int nvol2 = Nvol / 2;
    AFIELD be(nin, nvol2, nex), bo(nin, nvol2, nex);
    AFIELD xe(nin, nvol2, nex), xo(nin, nvol2, nex), xt(nin, nvol2, nex);

    AIndex_eo<real_t, IMPL> index_eo;
    index_eo.split(be, bo, abq);

    fopr_alt->mult_dag(xe, be, "Dee");
    fopr_alt->mult_dag(xt, xe, "Dee_inv");
    fopr_alt->mult_dag(xe, xt, "Dee");
    fopr_alt->mult_dag(xt, bo, "Doe");
    axpy(xe, 1.0, xt);
    fopr_alt->mult_dag(xo, bo, "Doo");
    fopr_alt->mult_dag(xt, xo, "Doo_inv");
    fopr_alt->mult_dag(xo, xt, "Doo");
    fopr_alt->mult_dag(xt, be, "Deo");
    axpy(xo, 1.0, xt);

    index_eo.merge(axq, xe, xo);

    fopr_ref->set_mode("Ddag");
    fopr_ref->mult(y5, b5);
    if (fopr_alt->needs_convert()) {
      fopr_alt->reverse((Field&)z5, axq);
    } else {
      reverse_spinor(index_alt, (Field&)z5, axq);
    }

    double y5norm_Ddag = y5.norm2();
    double z5norm_Ddag = z5.norm2();
    double axqnorm_Ddag = axq.norm2();

    axpy(y5, -1.0, z5);
    double diff2_Ddag = y5.norm2();

    vout.general(m_vl, "\n[Ddag check]\n");
    vout.general(m_vl, "norm2(abq) = %f\n", abqnorm);
    vout.general(m_vl, "norm2(y5)  = %f\n", y5norm_Ddag);
    vout.general(m_vl, "norm2(z5)  = %f\n", z5norm_Ddag);
    vout.general(m_vl, "norm2(axq) = %f\n", axqnorm_Ddag);
    vout.general(m_vl, "diff2 = %f\n", diff2_Ddag);

    if (diff2_Ddag > 1e-16) return 1;
  }

  return 0;
}

//============================================================END=====
