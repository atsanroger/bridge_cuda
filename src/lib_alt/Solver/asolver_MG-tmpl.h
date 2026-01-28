/*!

        @file    asolver_MG_tmpl.h
        @brief   multigrid solver (template version)
        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate::  $
        @version $LastChangedRevision: 2499 $

 */

#include "lib_alt/Solver/asolver_MG.h"

#include "lib/ResourceManager/threadManager.h"
#include "lib/Tools/randomNumberManager.h"

//====================================================================
template<typename AFIELD>
const std::string ASolver_MG<AFIELD>::class_name = "ASolver_MG";

/*  for book keeping
  // inner prod:   num_vector*(num_vector-1)/2 times
  // complex axpy: num_vector*(num_vector-1)/2 times
  // norm:         num_vector times
  // real scal:    num_vector times (GS)
  // convert rep: 2*num_vectors times
  // real scal:    num_vector_times (normalization for covert rep.)
  //   if one counts 0+a as one 1 flop, block inner prod has the same flop
  //   counting as the full inner prod
  const int flop_gs_site
    = 8*Nc*Nd*num_vectors*(num_vectors-1)/2   // inner prod
    +  8*Nc*Nd*num_vectors*(num_vectors-1)/2  // axpy
    +  4*Nc*Nd*num_vectors                    // norm
    +  2*Nc*Nd*num_vectors                    // scal (GS)
    +  2*Nc*Nd*2*num_vectors                  // convert_rep
    +  2*Nc*Nd*num_vectors;                   // scal (normalization)
  size_t Lvol=CommonParameters::Lvol();
  const double flop_gs=static_cast<double>(flop_gs_site)*static_cast<double>(Lvol);
*/

//====================================================================
template<typename AFIELD>
void ASolver_MG<AFIELD>::init()
{
  ThreadManager::assert_single_thread(class_name);
  constexpr int size = sizeof(typename AFIELD::real_t);
  if (size != sizeof(double)) {
    vout.crucial("ASolver_MG must be instanced with double prec. field\n");
    abort();
  }
}


//====================================================================
template<typename AFIELD>
void ASolver_MG<AFIELD>::tidyup()
{
}


//====================================================================
template<typename AFIELD>
void ASolver_MG<AFIELD>::set_parameters(const Parameters& params)
{
  // example of the parameter
  //   each of outer solver (in level 0) and coarse solver (in level 1)
  //   may futher has additional parameters like Omega_tolerance.
  //
  //  # general parameters, passed to the outer solver
  //  maximum_number_of_iteration : 100
  //  maximum_number_of_restart:     10
  //  convergence_criterion_squared: 1e-28
  //  verbose_level : Detailed
  //  # spedfice to MG solver
  //  MultiGrid_Level1:
  //    sap_block : [4,2,2,2]
  //    number_of_vectors: 4
  //    setup_number_of_setp  : 4
  //    maximum_number_of_iteration : 1000
  //    maximum_number_of_restart :   1
  //    convergence_criterion_squared: 2.0e-5
  //    smoother_number_of_iteration : 3
  //    smoother_convergence_criterion_squared: 1.0e-3
  //    verbose_level : General

  //  Parameters params_outer_solver=params.lookup("Level0");
  Parameters params_coarse = params.lookup("MultiGrid_Level1");

  // level0: outer solver
  set_parameters_level0(params);

  // level1: coarse grid solver and smoother
  set_parameters_level1(params_coarse);
}


//====================================================================
template<typename AFIELD>
void ASolver_MG<AFIELD>::set_parameters_level0(const Parameters& params)
{
  m_params_asolver_outer = params;
}


//====================================================================
template<typename AFIELD>
void ASolver_MG<AFIELD>::set_parameters_level1(const Parameters& params)
{
  // for setup
  m_sap_block_size = params.get_int_vector("sap_block");
  m_nvec           = params.get_int("setup_number_of_vectors");
  m_nsetup         = params.get_int("setup_number_of_step");

  // for smoother
  m_smoother_niter     = params.get_int("smoother_number_of_iteration");
  m_smoother_stop_cond = params.get_double("smoother_convergence_criterion_squared");

  // for coarse solver
  //   I assume that the above parameters harmless for the coarse grid solver
  m_params_asolver_coarse = params;

  set_lattice(m_sap_block_size);
}


//====================================================================
template<typename AFIELD>
void ASolver_MG<AFIELD>::set_parameters(
  const int Niter,
  const real_t Stop_cond,
  const std::string& outer_vlevel,
  const std::vector<int>& sap_block_size,
  const int nvec,
  const int nsetup,
  const int coarse_niter,
  const real_t coarse_stop_cond,
  const std::string& coarse_vlevel,
  const int smoother_niter,
  const real_t smoother_stop_cond)
{
  ThreadManager::assert_single_thread(class_name);

  // for outer solver
  Parameters param_level0;
  param_level0.set_int("maximum_number_of_iteration", Niter);
  param_level0.set_int("maximum_number_of_restart", 1);
  param_level0.set_double("convergence_criterion_squared", Stop_cond);
  param_level0.set_string("verbose_level", outer_vlevel);
  set_parameters_level0(param_level0);

  // coarse grid
  Parameters param_level1;

  // for building coarse grid
  param_level1.set_int_vector("sap_block", sap_block_size);
  param_level1.set_int("number_of_vectors", nvec);
  // coarse grid solver
  param_level1.set_int("maximum_number_of_iteration", coarse_niter);
  param_level1.set_int("maximum_number_of_restart", 1);
  param_level1.set_double("convergence_criterion_squared",
                          coarse_stop_cond);
  param_level1.set_string("verbose_level", coarse_vlevel);
  // smoother
  param_level1.set_int("smoother_number_of_iteration", smoother_niter);
  param_level1.set_int("smoother_convergence_criterion_squared",
                       smoother_stop_cond);

  set_parameters_level1(param_level1);
}


//====================================================================
template<typename AFIELD>
bool ASolver_MG<AFIELD>::use_fopr_for_smoother() const
{
#ifdef USE_FOPR_FOR_SMOOTHER
  return true;
#else
  return false;
#endif
}


//====================================================================
template<typename AFIELD>
AFopr_dd<typename ASolver_MG<AFIELD>::AFIELD_f> *ASolver_MG<AFIELD>::new_fopr_smoother(const Parameters& param) const
{
#ifdef USE_FOPR_FOR_SMOOTHER
  return new FoprSmoother_t(param);
#else
  return nullptr;
#endif
}


//====================================================================
template<typename AFIELD>
void ASolver_MG<AFIELD>::set_foprD(AFopr<AFIELD> *foprD)
{
  int th = ThreadManager::get_thread_id();
#pragma omp barrier
  if (th == 0) {
    m_afopr_fineD = dynamic_cast<FoprD_t *>(foprD);
  }
#pragma omp barrier
  if (m_afopr_fineD == nullptr) {
    vout.crucial("%s: bad afopr: only AFopr_Clover is vaild for FoprD]\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void ASolver_MG<AFIELD>::set_foprF(AFopr_dd<AFIELD_f> *foprF)
{
  m_afopr_fineF = dynamic_cast<FoprF_t *>(foprF);
  if (m_afopr_fineF == nullptr) {
    vout.crucial("%s: bad afopr: only AFopr_Clover_dd is vaild for FoprF]\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
#ifndef USE_FOPR_FOR_SMOOTHER
  m_afopr_smoother = m_afopr_fineF;
#endif
}


//====================================================================
template<typename AFIELD>
void ASolver_MG<AFIELD>::set_fopr_smoother(AFopr_dd<AFIELD_f> *foprF)
{
#ifdef USE_FOPR_FOR_SMOOTHER
  m_afopr_smoother = dynamic_cast<FoprSmoother_t *>(foprF);
  if (m_afopr_fineF == nullptr) {
    vout.crucial("%s: bad afopr: only AFopr_Clover_dd is vaild for Fopr smoother]\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
#else
  vout.general("%s: set_foprF_smoother is called but ignored.  The same one given with set_foprF will be used.\n", class_name.c_str());
#endif
}


//====================================================================
template<typename AFIELD>
void ASolver_MG<AFIELD>::init_solver(std::string mode)
{
  ThreadManager::assert_single_thread(class_name);
  m_mode = mode;

  // check fopr is registered
  if (!m_afopr_fineF) {
    vout.crucial(m_vl, "%s: init_solver, single prec. afopr is not yet set.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
  if (!m_afopr_fineD) {
    vout.crucial(m_vl, "%s: init_solver, double prec. afopr is not yet set.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  // fine grid operator
  m_afopr_fineD->set_mode(m_mode);
  m_afopr_fineF->set_mode(m_mode);

  // coarse grid solvers/smoothers
  init_coarse_grid();

  // outer solver
  m_outer_solver.reset(new OuterSolver_t(m_afopr_fineD, m_prec_mg.get()));
  m_outer_solver->set_parameters(m_params_asolver_outer);
}


//====================================================================
#ifndef SKIP_INIT_COARSE
template<typename AFIELD>
void ASolver_MG<AFIELD>::init_coarse_grid()
{
  vout.detailed(m_vl, "%s: init_coarse_grid() [template version] is called\n",
                class_name.c_str());

  // multigrid converter
  std::vector<int> fine_lattice = { CommonParameters::Nx(),
                                    CommonParameters::Ny(),
                                    CommonParameters::Nz(),
                                    CommonParameters::Nt() };
  const int        Nc = CommonParameters::Nc();
  const int        Nd = CommonParameters::Nd();

  MultiGrid_t *multigrid = new MultiGrid_t();
  multigrid->init(m_coarse_lattice, fine_lattice, 2 * Nc * Nd, m_nvec);
  multigrid->set_afopr_fine(m_afopr_fineF);
  m_multigrid.reset(multigrid);

  // afopr: downcasted objects
  //   downcast safty is guaranteed in the setter
  FoprF_t *afopr_fineF = static_cast<FoprF_t *>(m_afopr_fineF);
  FoprD_t *afopr_fineD = static_cast<FoprD_t *>(m_afopr_fineD);

  // initialize coase grid operator
  FoprCoarse_t *afopr_coarse = new FoprCoarse_t();
  afopr_coarse->set_parameters(m_nvec, m_coarse_lattice);
  afopr_coarse->set_mode(m_mode);
  m_afopr_coarse.reset(afopr_coarse);

  vout.general(m_vl, "afopr_coarse version is created\n");

  // initialize coarse grid solver
  CoarseSolver_t *asolver_coarse = new CoarseSolver_t(m_afopr_coarse.get());
  asolver_coarse->set_parameters(m_params_asolver_coarse);
  asolver_coarse->set_init_mode(ASolver<AFIELD_f>::InitialGuess::ZERO);
  m_asolver_coarse.reset(asolver_coarse);

  // initialize smoother
#ifdef USE_SAP_FOR_SMOOTHER
  Smoother_t *asolver_smoother =
    new Smoother_t(m_afopr_smoother, m_multigrid->get_block_index());
#else
  Smoother_t *asolver_smoother = new Smoother_t(m_afopr_smoother);
#endif
  asolver_smoother->set_parameters(m_smoother_niter, m_smoother_stop_cond);
  m_asolver_smoother.reset(asolver_smoother);

  // combine everything and initialize multgrid preconditionor
  m_prec_mg.reset(new APrecond_MG<AFIELD, AFIELD_f>());
  m_prec_mg->set_coarse_solver(m_asolver_coarse.get());
  m_prec_mg->set_smoother(m_asolver_smoother.get());
  m_prec_mg->set_multigrid(m_multigrid.get());
  m_prec_mg->set_fopr(m_afopr_fineD, m_afopr_fineF);
}


#endif

//====================================================================
template<typename AFIELD>
void ASolver_MG<AFIELD>::set_lattice(const vector<int>& sap_block_size)
{
  ThreadManager::assert_single_thread(class_name);

  // set coarse lattice
  assert(CommonParameters::Ndim() == 4);
  m_coarse_lattice.resize(4);
  std::vector<int> fine_lattice = { CommonParameters::Nx(),
                                    CommonParameters::Ny(),
                                    CommonParameters::Nz(),
                                    CommonParameters::Nt() };

  // sanity checks
  if (sap_block_size.size() != 4) {
    vout.crucial(m_vl, "%s: bad sap_block_size:  Must be 4-dim, but the given dimension is %.\n",
                 class_name.c_str(), sap_block_size.size());
    abort();
  }

  for (int i = 0; i < 4; ++i) {
    m_coarse_lattice[i] = fine_lattice[i] / sap_block_size[i];
    if (m_coarse_lattice[i] * sap_block_size[i] != fine_lattice[i]) {
      vout.crucial(m_vl, "bad sap_block_size: i=%d, sap_block_size=%d, fine_lattice=%d, coarse_lattice=%d\n",
                   i, sap_block_size[i], fine_lattice[i], m_coarse_lattice[i]);
      exit(EXIT_FAILURE);
    }
  }

  vout.general(m_vl, "  fine_lattice        = %s\n",
               Parameters::to_string(fine_lattice).c_str());
  vout.general(m_vl, "  coarse_lattice      = %s\n",
               Parameters::to_string(m_coarse_lattice).c_str());
  vout.general(m_vl, "  sap_block_size      = %s\n",
               Parameters::to_string(sap_block_size).c_str());
}


//====================================================================
template<typename AFIELD>
void ASolver_MG<AFIELD>::run_setup(std::vector<AFIELD_f>& testvec_work)
{
  assert(m_nvec == testvec_work.size());

  Timer timer_setup("setup total (work vectors given)");
#pragma omp master
  {
    timer_setup.start();
    m_timer_gramschmidt.reset(new Timer("Gramschmidt in the setup"));
    m_timer_generate_coarse_op.reset(new Timer("generate coarse op"));
  }

  run_setup_initial(testvec_work);
  run_setup_iterative(m_nsetup, testvec_work);

#pragma omp master
  {
    timer_setup.stop();

    m_timer_gramschmidt->report();
    m_timer_generate_coarse_op->report();
    m_prec_mg->report_timer();
    m_prec_mg->reset_flop_count();
    timer_setup.report();
  }
  vout.general("setup is done!\n");
}


//====================================================================
template<typename AFIELD>
void ASolver_MG<AFIELD>::run_setup()
{
  ThreadManager::assert_single_thread(class_name);

  Timer timer_setup("setup total");
  timer_setup.start();

  const int num_vectors = m_nvec;
  const int Nin         = m_afopr_fineD->field_nin();
  const int Nvol        = m_afopr_fineD->field_nvol();
  const int Nex         = m_afopr_fineD->field_nex();

  //  std::vector<double> flop_setup(m_nsetup+1);

  std::vector<AFIELD_f> testvec_work(num_vectors);
  for (int i = 0; i < num_vectors; ++i) {
    testvec_work[i].reset(Nin, Nvol, Nex);
  }
#pragma omp parallel
  {
    run_setup(testvec_work);
  }
  timer_setup.stop();
  timer_setup.report();
}


//====================================================================
template<typename AFIELD>
void ASolver_MG<AFIELD>::run_setup_initial(
  std::vector<AFIELD_f>& testvec_work)
{
  // Todo: add flops count

  //  ThreadManager::assert_single_thread(class_name);
  //#pragma omp parallel
  {
    vout.detailed("run_setup: using single precision Gramschmidt\n");


    assert(testvec_work.size() == m_nvec);
    unique_ptr<Timer> timer_initial_setup;
#pragma omp master
    {
      timer_initial_setup.reset(new Timer("initial setup"));
      timer_initial_setup->start();
    }

    ASolver<AFIELD_f> *asolver_setup = m_asolver_smoother.get();


    // random vectors are not thread paraleized....
    m_multigrid->set_testvectors(); // generate initial random vectors


    // initial setup
    for (int i = 0; i < m_nvec; ++i) {
      int   nconv = -1;
      float diff  = -1.0;
      asolver_setup->solve(testvec_work[i],
                           (*m_multigrid->get_testvectors())[i],
                           nconv, diff);
    } // i
#pragma omp master
    {
      m_timer_gramschmidt->start();
    }
    m_multigrid->gramschmidt(testvec_work);
#pragma omp master
    {
      m_timer_gramschmidt->stop();
    }
    m_multigrid->set_testvectors(testvec_work); // update the testvectors
    vout.general("initial testvectors are ready\n");
    FoprCoarse_t *afopr_coarse =
      static_cast<FoprCoarse_t *>(m_afopr_coarse.get());
#pragma omp master
    {
      m_timer_generate_coarse_op->start();
    }
    afopr_coarse->generate_coarse_op(m_afopr_fineF,
                                     *m_multigrid->get_testvectors());
#pragma omp master
    {
      m_timer_generate_coarse_op->stop();
    }
    vout.general(m_vl, "afopr_coarse version is ready\n");

#pragma omp master
    {
      timer_initial_setup->stop();
      timer_initial_setup->report();
    }
  }
}


//====================================================================
template<typename AFIELD>
void ASolver_MG<AFIELD>::run_setup_iterative(int niter,
                                             std::vector<AFIELD_f>& testvec_work)
{
  // Todo: add flops count

  // N.B.
  // working vectors [= testvec_work] are needed because
  // the "solver" uses the multigrid precondinor

  //  ThreadManager::assert_single_thread(class_name);
  //#pragma omp parallel
  {
    assert(testvec_work.size() == m_nvec);
    unique_ptr<Timer> timer_setup;
#pragma omp master
    {
      timer_setup.reset(new Timer("each setup"));
    }
    for (int n = 0; n < niter; ++n) {
#pragma omp master
      {
        timer_setup->reset();
        timer_setup->start();
      }
      for (int i = 0; i < m_nvec; ++i) {
        //#pragma omp barrier
        m_prec_mg->mult_as_setup(testvec_work[i],
                                 (*m_multigrid->get_testvectors())[i]);
      } // i

#pragma omp barrier

#pragma omp master
      {
        m_timer_gramschmidt->start();
      }
      m_multigrid->gramschmidt(testvec_work);
#pragma omp master
      {
        m_timer_gramschmidt->stop();
      }
      m_multigrid->set_testvectors(testvec_work); // update the testvectors
      FoprCoarse_t *afopr_coarse =
        static_cast<FoprCoarse_t *>(m_afopr_coarse.get());
#pragma omp master
      {
        m_timer_generate_coarse_op->start();
      }
      afopr_coarse->generate_coarse_op(m_afopr_fineF,
                                       *m_multigrid->get_testvectors());
#pragma omp master
      {
        m_timer_generate_coarse_op->stop();
      }
      vout.general(m_vl, "renewed afopr_coarse: n=%d / %d\n", n, niter);
#pragma omp master
      {
        timer_setup->stop();
        timer_setup->report();
      }
    } //n
  }   // omp parallel
}


//====================================================================
template<typename AFIELD>
void ASolver_MG<AFIELD>::solve(AFIELD& x, const AFIELD& b,
                               int& nconv, real_t& diff)
{
  m_outer_solver->solve(x, b, nconv, diff);
}


//====================================================================
template<typename AFIELD>
void ASolver_MG<AFIELD>::reset_flop_count()
{
  m_nconv = 0;
  m_prec_mg->reset_flop_count();
}


//====================================================================
template<typename AFIELD>
double ASolver_MG<AFIELD>::flop_count()
{
  double flop_solve    = m_outer_solver->flop_count();
  double flop_outer    = flop_solve - m_prec_mg->flop_count();
  double flop_coarse   = m_prec_mg->flop_count_coarse();
  double flop_smoother = m_prec_mg->flop_count_smoother();
  double flop_other    = m_prec_mg->flop_count_other();
  double flop_double   = m_prec_mg->flop_count_double();
  m_prec_mg->report_timer();

  vout.general(m_vl, "flop count (MG solver) [GFlop]:\n");
  vout.general(m_vl, "   solve total (double+float): %e\n", flop_solve * 1.0e-9);
  vout.general(m_vl, "   solve coarse (float):   %e\n", flop_coarse * 1.0e-9);
  vout.general(m_vl, "   solve smoother (float): %e\n", flop_smoother * 1.0e-9);
  vout.general(m_vl, "   solve other (float):    %e\n", flop_other * 1.0e-9);
  vout.general(m_vl, "   solve other (double):   %e\n", flop_double * 1.0e-9);

  return flop_solve;
}


//============================================================END=====
