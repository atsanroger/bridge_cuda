/*!
      @file    fprop_alt_Standard_SAP.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
      @version $LastChangedRevision: 2499 $
*/

#include "Measurements/Fermion/fprop_alt_Standard_SAP.h"

#include "lib/ResourceManager/threadManager.h"


template<typename AFIELD>
const std::string Fprop_alt_Standard_SAP<AFIELD>::class_name
  = "Fprop_alt_Standard_SAP";

//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_SAP<AFIELD>::init(
  const Parameters& params_fopr,
  const Parameters& params_solver)
{
  vout.general(m_vl, "\n");
  vout.general(m_vl, "%s: being setup (without link smearing).\n",
               class_name.c_str());

  m_dr_smear = 0;

  string fopr_type = params_fopr.get_string("fermion_type");
  vout.general("%s: fopr_type: %s\n", class_name.c_str(), fopr_type.c_str());
  string fopr_type_dd = fopr_type + "_dd";
  m_kernel = 0;

  m_afopr = dynamic_cast<AFopr_dd<AFIELD> *>(
    AFopr<AFIELD>::New(fopr_type_dd, params_fopr));

  if (m_afopr == nullptr) {
    vout.crucial(m_vl, "%s: failed to create AFopr_dd with fermion type %s.\n", class_name.c_str(), fopr_type.c_str());
    exit(EXIT_FAILURE);
  }
  init_common(params_fopr, params_solver);

  vout.general(m_vl, "%s: setup finished.\n", class_name.c_str());
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_SAP<AFIELD>::init(
  const Parameters& params_fopr,
  const Parameters& params_solver,
  Director_Smear *dr_smear)
{
  vout.general(m_vl, "\n");
  vout.general(m_vl, "%s: being setup (with link smearing).\n",
               class_name.c_str());

  vout.crucial(m_vl, "sorry, not available.\n");
  exit(EXIT_FAILURE);

  /*
  m_dr_smear = dr_smear;

  string fopr_type = params_fopr.get_string("fermion_type");

  m_kernel_d = AFopr<AFIELD_d>::New(fopr_type, params_fopr);
  m_fopr_d = new AFopr_Smeared<AFIELD_d>(m_kernel_d, m_dr_smear);

  m_kernel_f = AFopr<AFIELD_f>::New(fopr_type, params_fopr);
  m_fopr_f = new AFopr_Smeared<AFIELD_f>(m_kernel_f, m_dr_smear);

  init_common(params_fopr, params_solver);

  vout.general(m_vl, "%s: setup finished.\n", class_name.c_str());
  */
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_SAP<AFIELD>::init_common(
  const Parameters& params_fopr,
  const Parameters& params_solver)
{
  int Ndim = CommonParameters::Ndim();
  m_fine_lattice.resize(Ndim);
  m_coarse_lattice.resize(Ndim);

  m_fine_lattice[0] = CommonParameters::Nx();
  m_fine_lattice[1] = CommonParameters::Ny();
  m_fine_lattice[2] = CommonParameters::Nz();
  m_fine_lattice[3] = CommonParameters::Nt();

  std::vector<int> block_size;
  int              err = 0;
  err += params_fopr.fetch_int_vector("block_size", block_size);
  if (err) {
    vout.crucial(m_vl, "Error at %s: block_size not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  for (int mu = 0; mu < Ndim; ++mu) {
    m_coarse_lattice[mu] = m_fine_lattice[mu] / block_size[mu];
  }

  for (int mu = 0; mu < Ndim; ++mu) {
    vout.general(m_vl, "  block_size[%d] = %3d"
                       "   coarse_lattice[%d] = %3d\n",
                 mu, block_size[mu], mu, m_coarse_lattice[mu]);
  }

  m_block_index = new AIndex_block_lex<real_t, AFIELD::IMPL>(
    m_coarse_lattice,
    m_fine_lattice);

  m_asolver = new ASolver_SAP<AFIELD>(m_afopr, m_block_index);
  vout.general(m_vl, "SAP is initiated, done\n");

  int    niter;
  double stop_cond;
  err += params_solver.fetch_int("maximum_number_of_iteration",
                                 niter);
  err += params_solver.fetch_double("convergence_criterion_squared",
                                    stop_cond);
  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  m_asolver->set_parameters(niter, real_t(stop_cond));

  reset_performance();
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_SAP<AFIELD>::tidyup()
{
  delete m_asolver;
  delete m_afopr;
  delete m_block_index;
  if (m_kernel != 0) delete m_kernel;
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_SAP<AFIELD>::set_config(Field *U)
{
  m_afopr->set_config(U);
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_SAP<AFIELD>::invert(
  Field& xq, const Field& b,
  int& nconv, double& diff)
{
  vout.paranoiac(m_vl, "%s: invert is called.\n", class_name.c_str());
  vout.paranoiac(m_vl, "mode = %s.\n", m_mode.c_str());

  m_timer.reset();
  m_timer.start();

  if (m_mode == "D") {
    invert_D(xq, b, nconv, diff);
  } else if (m_mode == "DdagD") {
    invert_DdagD(xq, b, nconv, diff);
  } else {
    vout.crucial(m_vl, "%s: unsupported mode: %s\n",
                 class_name.c_str(), m_mode.c_str());
    exit(EXIT_FAILURE);
  }

  m_timer.stop();
  m_elapsed_time += m_timer.elapsed_sec();
  // m_flop_count   += m_solver->flop_count();
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_SAP<AFIELD>::invert_D(
  Field& xq, const Field& b,
  int& nconv, double& diff)
{
  int Nin  = m_afopr->field_nin();
  int Nvol = m_afopr->field_nvol();
  int Nex  = m_afopr->field_nex();

  // Source vector
  AFIELD ab(Nin, Nvol, Nex);

  // Solution vector
  AFIELD ax(Nin, Nvol, Nex);

  // conversion
  AIndex_lex<float, AFIELD::IMPL> index_alt;
  vout.paranoiac(m_vl, "index is ready\n");

#pragma omp parallel
  {
    if (m_afopr->needs_convert()) {
      vout.detailed(m_vl, "convert required.\n");
      m_afopr->convert(ab, b);
    } else {
      vout.detailed(m_vl, "convert not required.\n");
      convert(index_alt, ab, b);
    }
  }
  vout.general(m_vl, "converting source, done\n");

  real_t nr = ab.norm2();
  nr = sqrt(1.0 / nr);
  ab.scal(nr);

  invert_D(ax, ab, nconv, diff);

#pragma omp parallel
  {
    if (m_afopr->needs_convert()) {
      vout.detailed(m_vl, "convert required.\n");
      m_afopr->reverse(xq, ax);
    } else {
      vout.detailed(m_vl, "convert not required.\n");
      reverse(index_alt, xq, ax);
    }
  }
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_SAP<AFIELD>::invert_DdagD(
  Field& xq, const Field& b,
  int& nconv, double& diff)
{
  vout.crucial(m_vl, "%s: invert_DdagD is not available.\n",
               class_name.c_str());
  exit(EXIT_FAILURE);
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_SAP<AFIELD>::invert(
  AFIELD& xq, const AFIELD& b,
  int& nconv, double& diff)
{
  vout.paranoiac(m_vl, "%s: invert is called.\n", class_name.c_str());
  vout.paranoiac(m_vl, "mode = %s.\n", m_mode.c_str());

  if (m_mode == "D") {
    invert_D(xq, b, nconv, diff);
  } else if (m_mode == "DdagD") {
    invert_DdagD(xq, b, nconv, diff);
  } else {
    vout.crucial(m_vl, "%s: unsupported mode: %s\n",
                 class_name.c_str(), m_mode.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_SAP<AFIELD>::invert_D(
  AFIELD& ax, const AFIELD& ab,
  int& nconv, double& diff)
{
  int Nin  = m_afopr->field_nin();
  int Nvol = m_afopr->field_nvol();
  int Nex  = m_afopr->field_nex();

  m_afopr->set_mode(m_mode);

  real_t diff1, diff2;

#pragma omp parallel
  {
    m_asolver->solve(ax, ab, nconv, diff1);
  }

  AFIELD ay(Nin, Nvol, Nex);

#pragma omp parallel
  {
    m_afopr->mult(ay, ax);
#pragma omp barrier

    axpy(ay, real_t(-1.0), ab);
#pragma omp barrier

    diff2 = ay.norm2();
  }

  diff = double(diff2);

  real_t norm2 = ab.norm2();

  vout.general("AField solver done:\n");
  vout.general("  solver residue squared = %22.15e\n", diff1);
  vout.general("  real residue squared   = %22.15e\n", diff2);
  vout.general("  relative residue       = %22.15e\n", diff2 / norm2);
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_SAP<AFIELD>::invert_DdagD(
  AFIELD& xq, const AFIELD& b,
  int& nconv, double& diff)
{
  vout.crucial(m_vl, "%s: invert_DdagD is not available.\n",
               class_name.c_str());
  exit(EXIT_FAILURE);
}


//====================================================================
template<typename AFIELD>
double Fprop_alt_Standard_SAP<AFIELD>::flop_count()
{
  return m_asolver->flop_count();
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_SAP<AFIELD>::reset_performance()
{
  m_flop_count   = 0.0;
  m_elapsed_time = 0.0;
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_SAP<AFIELD>::report_performance()
{
  double flops  = m_flop_count / m_elapsed_time;
  double gflops = flops * 1.0e-9;

  vout.general(m_vl, "\n");
  vout.general(m_vl, "%s: solver performance:\n", class_name.c_str());
  vout.general(m_vl, "  NOT YET implelemted.\n");

  /*
  vout.general(m_vl, "  Elapsed time = %14.6f sec\n", m_elapsed_time);
  vout.general(m_vl, "  Flop(total)  = %18.0f\n",m_flop_count);
  vout.general(m_vl, "  Performance  = %11.3f GFlops\n", gflops);
  */
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_SAP<AFIELD>::mult_performance(
  const std::string mode,
  const int Nrepeat)
{
  /*
  int nin  = m_afopr->field_nin();
  int nvol = m_afopr->field_nvol();
  int nex  = m_afopr->field_nex();

  unique_ptr<Timer> timer(new Timer);

  std::string mode_prev_d = m_fopr_d->get_mode();
  std::string mode_prev_f = m_fopr_f->get_mode();

  m_fopr_d->set_mode(mode);
  m_fopr_f->set_mode(mode);

 {
  AFIELD axq(nin, nvol, nex), abq(nin, nvol, nex);
  abq.set(0.0);
  abq.set(0, 1.0);

  timer->start();

#pragma omp parallel
  {
    for(int i = 0; i < Nrepeat; ++i){
      m_fopr_d->mult(axq, abq);
      m_fopr_d->mult(abq, axq);
    }
  }

  timer->stop();
 }

  double flop_fopr = m_fopr_d->flop_count();
  double flop_total = flop_fopr * double(2 * Nrepeat);

  double elapsed_time = timer->elapsed_sec();
  double flops  = flop_total/elapsed_time;
  double gflops = flops * 1.0e-9;

  vout.general(m_vl, "\n");
  vout.general(m_vl, "%s: mult performance:\n", class_name.c_str());
  vout.general(m_vl, "  Double precision:\n");
  vout.general(m_vl, "   Elapsed time = %14.6f sec\n", elapsed_time);
  vout.general(m_vl, "   Flop(Fopr)   = %18.0f\n",flop_fopr);
  vout.general(m_vl, "   Flop(total)  = %18.0f\n",flop_total);
  vout.general(m_vl, "   Performance  = %11.3f GFlops\n", gflops);

 {
  AFIELD axq(nin, nvol, nex), abq(nin, nvol, nex);
  abq.set(0.0);
  abq.set(0, 1.0);

  timer->reset();
  timer->start();

#pragma omp parallel
  {
    for(int i = 0; i < Nrepeat; ++i){
      m_fopr_f->mult(axq, abq);
      m_fopr_f->mult(abq, axq);
    }
  }

  timer->stop();
 }

  flop_fopr = m_fopr_f->flop_count();
  flop_total = flop_fopr * double(2 * Nrepeat);

  elapsed_time = timer->elapsed_sec();
  flops  = flop_total/elapsed_time;
  gflops = flops * 1.0e-9;

  vout.general(m_vl, "  Single precision:\n");
  vout.general(m_vl, "   Elapsed time = %14.6f sec\n", elapsed_time);
  vout.general(m_vl, "   Flop(Fopr)   = %18.0f\n",flop_fopr);
  vout.general(m_vl, "   Flop(total)  = %18.0f\n",flop_total);
  vout.general(m_vl, "   Performance  = %11.3f GFlops\n", gflops);

  m_fopr_d->set_mode(mode_prev_d);
  m_fopr_f->set_mode(mode_prev_f);
  */
}


//============================================================END=====
