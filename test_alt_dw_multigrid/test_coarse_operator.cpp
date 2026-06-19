/*!
        @file    test_coarse_operator.cpp

        @brief   test for coarse grid operator

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: #$

        @version $LastChangedRevision: 2598 $
*/

#include "test_coarse_operator.h"

//#include "test.h"

#include "lib/IO/gaugeConfig.h"
#include "lib/Tools/randomNumberManager.h"
#include "lib/Tools/randomNumbers_Mseries.h"
#include "lib/Tools/timer.h"
#include "lib/Fopr/fopr.h"
#include "lib/Field/field_F.h"
#include "lib/Field/field_G.h"

#include "lib/Field/shiftField_lex.h"
#include "lib/Parameters/parameterManager.h"


#ifdef USE_ALT_SIMD
#include "lib_alt_SIMD/bridge_alt_simd.h"
#include "lib_alt_SIMD/Fopr/afopr_Domainwall_5din.h"
#include "lib_alt_SIMD/Fopr/afopr_Domainwall_5din_dd.h"
#include "lib_alt_SIMD/Fopr/afopr_Domainwall_coarse.h"
const Impl IMPL = SIMD;
#endif

#ifdef USE_ALT_QXS
#include "lib_alt_QXS/bridge_alt_qxs.h"
#include "lib_alt_QXS/Field/aindex_block_lex.h"
#include "lib_alt_QXS/Field/aindex_coarse_lex.h"
#include "lib_alt_QXS/Fopr/afopr_Domainwall_5din.h"
#include "lib_alt_QXS/Fopr/afopr_Domainwall_5din_dd.h"
#include "lib_alt_QXS/Fopr/afopr_Domainwall_coarse.h"

const Impl IMPL = QXS;
#endif

#ifdef USE_ALT_ACCEL
#include "lib_alt_Accel/bridge_alt_accel.h"
#include "lib_alt_Accel/Field/aindex_block_lex.h"
#include "lib_alt_Accel/Field/aindex_coarse_lex.h"
#include "lib_alt_Accel/Fopr/afopr_Domainwall_5din.h"
#include "lib_alt_Accel/Fopr/afopr_Domainwall_5din_dd.h"
#include "lib_alt_Accel/Fopr/afopr_Domainwall_coarse.h"
const Impl IMPL = ACCEL;
#endif

typedef float real_t;


//#include "lib_alt/Solver/MultiGrid_Domainwall.h"
#include "MultiGrid_Domainwall.h"

#include "lib_alt/Solver/asolver.h"
#include "lib_alt/Solver/asolver_BiCGStab.h"
#include "lib_alt/Solver/asolver_BiCGStab_Cmplx.h"
#include "lib_alt/Solver/asolver_GMRES_m_Cmplx.h"
#include "lib_alt/Solver/asolver_CG.h"

#include "lib_alt/Eigen/aeigensolver_IRArnoldi_block.h"

#include "asolver_SAP_dw.h"
#include "asolver_SAP_MINRES_dw.h"


// Registration as a test
namespace test_coarse {
  const string test_name = "Test_Coarse";

  int test(void)
  {
    Test_Coarse test;
    return test.test();
  }
}

// class name
const std::string Test_Coarse::class_name =
  "Test_Coarse";

//====================================================================
void Test_Coarse::init()
{
  // do nothing.
}


//====================================================================
namespace {

//====================================================================
  void set_coarse_lattice(vector<int>& coarse_lattice, const vector<int>& sap_block_size, const Bridge::VerboseLevel vl)
  {
    assert(CommonParameters::Ndim() == 4);
    coarse_lattice.resize(4);
    vector<int> fine_lattice(4);
    fine_lattice[0] = CommonParameters::Nx();
    fine_lattice[1] = CommonParameters::Ny();
    fine_lattice[2] = CommonParameters::Nz();
    fine_lattice[3] = CommonParameters::Nt();
    for (int i = 0; i < 4; i++) {
      coarse_lattice[i] = fine_lattice[i] / sap_block_size[i];
      if (coarse_lattice[i] * sap_block_size[i] != fine_lattice[i]) {
        vout.crucial("bad sap_block_size: i=%d, sap_block_size=%d, fine_lattice=%d, coarse_lattice=%d\n",
                     i, sap_block_size[i], fine_lattice[i], coarse_lattice[i]);
        abort();
      }

      /*
      if(coarse_lattice[i] % CommonParameters::NPEsize(i) != 0){
        vout.crucial("bad sap_block_size: i=%d, grid_size=%d, coarse_lattice=%d\n",
                     i, sap_block_size[i], CommonParameters::NPEsize(i), coarse_lattice[i]);
        abort();
      }
      */
    }

    vout.general(vl, "  fine_lattice        = %s\n", Parameters::to_string(fine_lattice).c_str());
    vout.general(vl, "  coarse_lattice      = %s\n", Parameters::to_string(coarse_lattice).c_str());
  }


  typedef AField<float, IMPL> AFIELD_f;
  typedef AField<double, IMPL> AFIELD_d;

  using complexF_t = typename AFIELD_f::complex_t;

  void make_random_coarse_vector(AFIELD_f& vec, std::vector<int> coarse_lattice, int nvec)
  {
    const int             seed = 314159;
    RandomNumbers_Mseries random(seed);

    const int Nx = coarse_lattice[0];
    const int Ny = coarse_lattice[1];
    const int Nz = coarse_lattice[2];
    const int Nt = coarse_lattice[3];

    const int Lx = Nx * Communicator::npe(0);
    const int Ly = Ny * Communicator::npe(1);
    const int Lz = Nz * Communicator::npe(2);
    const int Lt = Nt * Communicator::npe(3);

    const int ipe_x = Communicator::ipe(0);
    const int ipe_y = Communicator::ipe(1);
    const int ipe_z = Communicator::ipe(2);
    const int ipe_t = Communicator::ipe(3);

    const int Nex = vec.nex();
    const int Nin = vec.nin();

    AIndex_coarse_lex<AFIELD_f::real_t, AFIELD_f::IMPL> index(
      Nx, Ny, Nz, Nt, nvec, 2);

    for (int j = 0; j < Nex; ++j) {
      bool in_j = true;

      for (int t = 0; t < Lt; ++t) {
        bool in_t = in_j && (t >= ipe_t * Nt) && (t < (ipe_t + 1) * Nt);

        for (int z = 0; z < Lz; ++z) {
          bool in_z = in_t && (z >= ipe_z * Nz) && (z < (ipe_z + 1) * Nz);

          for (int y = 0; y < Ly; ++y) {
            bool in_y = in_z && (y >= ipe_y * Ny) && (y < (ipe_y + 1) * Ny);

            for (int x = 0; x < Lx; ++x) {
              bool in_x = in_y && (x >= ipe_x * Nx) && (x < (ipe_x + 1) * Nx);

              //              int isite = index.site(x%Nx,y%Ny,z%Nz,t%Nt);
              //              int isite = x + Nx * (y + Ny * (z + Nz * t));
              int isite = (x % Nx) + Nx * ((y % Ny) + Ny * ((z % Nz) + Nz * (t % Nt)));
              for (int ic = 0; ic < nvec; ++ic) {
                for (int id = 0; id < 2; ++id) {
                  double re = random.get();
                  double im = random.get();
                  if (in_x) {
                    int idx_re = index.idx_SPr(ic, id, isite, j);
                    int idx_im = index.idx_SPi(ic, id, isite, j);
                    //                    vout.general("hoge: j, isite, ic, id: %d %d %d %d : idx = %d %d\n",
                    //                                j,isite,ic,id,idx_re, idx_im);
                    vec.set(idx_re, re);
                    vec.set(idx_im, im);
                  }
                }
              }    // ic, id
            }
          }
        }
      }          // x,y,z,t
    } // j
  }


  void make_random_fine_vector(AFIELD_f& vec)
  {
    const int             seed = 54321;
    RandomNumbers_Mseries random(seed);

    const int Nx = CommonParameters::Nx();
    const int Ny = CommonParameters::Ny();
    const int Nz = CommonParameters::Nz();
    const int Nt = CommonParameters::Nt();

    const int Lx = Nx * Communicator::npe(0);
    const int Ly = Ny * Communicator::npe(1);
    const int Lz = Nz * Communicator::npe(2);
    const int Lt = Nt * Communicator::npe(3);

    const int ipe_x = Communicator::ipe(0);
    const int ipe_y = Communicator::ipe(1);
    const int ipe_z = Communicator::ipe(2);
    const int ipe_t = Communicator::ipe(3);

    const int Nex = vec.nex();
    const int Nin = vec.nin();

    int Nc = 3;
    AIndex_lex<AFIELD_f::real_t, AFIELD_f::IMPL> index(
      Nx, Ny, Nz, Nt);

    for (int j = 0; j < Nex; ++j) {
      bool in_j = true;

      for (int t = 0; t < Lt; ++t) {
        bool in_t = in_j && (t >= ipe_t * Nt) && (t < (ipe_t + 1) * Nt);

        for (int z = 0; z < Lz; ++z) {
          bool in_z = in_t && (z >= ipe_z * Nz) && (z < (ipe_z + 1) * Nz);

          for (int y = 0; y < Ly; ++y) {
            bool in_y = in_z && (y >= ipe_y * Ny) && (y < (ipe_y + 1) * Ny);

            for (int x = 0; x < Lx; ++x) {
              bool in_x = in_y && (x >= ipe_x * Nx) && (x < (ipe_x + 1) * Nx);

              int isite = index.site(x % Nx, y % Ny, z % Nz, t % Nt);
              for (int ic = 0; ic < Nc; ++ic) {
                for (int id = 0; id < 4; ++id) {
                  double re = random.get();
                  double im = random.get();
                  if (in_x) {
                    int idx_re = index.idx_SPr(ic, id, isite, j);
                    int idx_im = index.idx_SPi(ic, id, isite, j);
                    vec.set(idx_re, re);
                    vec.set(idx_im, im);
                  }
                }
              }    // ic, id
            }
          }
        }
      }          // x,y,z,t
    } // j
  }
} // anonymous namespace

//====================================================================
int Test_Coarse::test()
{
  const std::string filename_input = "test_alt_Multigrid.yaml";

  const std::string test_name = class_name;

  vout.general("\n");
  vout.general("test name: %s\n", test_name.c_str());

  // ####  parameter setup  ####
  const int Nc   = CommonParameters::Nc();
  const int Nd   = CommonParameters::Nd();
  const int Ndim = CommonParameters::Ndim();
  const int Nvol = CommonParameters::Nvol();

  params_all = ParameterManager::read(filename_input);

  Parameters params_gauge = params_all.lookup("Gauge");

  const string str_gconf_status = params_gauge.get_string("gauge_config_status");
  const string str_gconf_read   = params_gauge.get_string("gauge_config_type_input");
  const string readfile         = params_gauge.get_string("config_filename_input");

  Parameters   params_mg_solver = params_all.lookup("MGSolver");
  const string str_vlevel       = params_mg_solver.get_string("verbose_level");

  Parameters             params_coarse  = params_mg_solver.lookup("MultiGrid_Level1");
  const std::vector<int> sap_block_size = params_coarse.get_int_vector("sap_block");
  const int              num_vectors    = params_coarse.get_int("setup_number_of_vectors");
  //    const string solver_type = params_coarse.get_string("solver_type"); // MG


  Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

  // setup random number manager
  RandomNumberManager::initialize("Mseries", 1234567UL);

  //- print input parameters
  vout.general(vl, "  gconf_status = %s\n", str_gconf_status.c_str());
  vout.general(vl, "  gconf_read   = %s\n", str_gconf_read.c_str());
  vout.general(vl, "  readfile     = %s\n", readfile.c_str());
  //  vout.general(vl, "  Nmult        = %d\n", Nmult);
  vout.general(vl, "  vlevel       = %s\n", str_vlevel.c_str());
  vout.general(vl, "  sap_block_size        = %s\n", Parameters::to_string(sap_block_size).c_str());
  vout.general(vl, "  num_vectors  = %d\n", num_vectors);

  Parameters params_fopr = params_all.lookup("Fopr");

  //- input parameter check
  int err = 0;
  err += ParameterCheck::non_NULL(str_gconf_status);

  if (err) {
    vout.crucial(vl, "Error at %s: input parameters have not been set\n", test_name.c_str());
    exit(EXIT_FAILURE);
  }

  unique_ptr<Timer> timer_total(new Timer(test_name));
  timer_total->start();

  //RandomNumberManager::initialize("Mseries", 1234567UL);
  RandomNumbers *random = RandomNumberManager::getInstance();

  // ####  Setup gauge configuration  ####
  vout.general("Field, creating: U\n");
  U.reset(new Field_G(Nvol, Ndim));


  if (str_gconf_status == "Continue") {
    GaugeConfig(str_gconf_read).read(*U.get(), readfile);
  } else if (str_gconf_status == "Cold_start") {
    GaugeConfig("Unit").read(*U.get());
  } else if (str_gconf_status == "Hot_start") {
    GaugeConfig("Random").read(*U.get());
  } else {
    vout.crucial(m_vl, "Error at %s: unsupported gconf status \"%s\"\n",
                 test_name.c_str(), str_gconf_status.c_str());
    exit(EXIT_FAILURE);
  }

  bool run_fopr       = false;
  bool run_afopr      = false;
  bool run_comparison = false;

  // prepareing fine grid operator w/ block op. extension
  Parameters params_fopr_dd = params_fopr;
  params_fopr_dd.set_int_vector("block_size", sap_block_size);
  AFopr_Domainwall_5din_dd<AFIELD_f> *afopr_fineF = new AFopr_Domainwall_5din_dd<AFIELD_f>(params_fopr_dd);
  afopr_fineF->set_config(U.get());

  int nin  = afopr_fineF->field_nin();
  int nvol = afopr_fineF->field_nvol();
  int nex  = afopr_fineF->field_nex();
  int Ls=nin/(2*Nc*Nd);

  vout.general("hoge: nin=%d, nvol=%d, nex=%d, Ls=%d\n", nin, nvol, nex, Ls);

  // Source vector
  Field_F b(nvol, Ls);

  Field_F y(nvol, Ls);
  double  norm_b;
  //#pragma omp parallel
  {
    b.set(1.0);
    //random->uniform_lex_global(b);
    //b.set(0.0);
    //b.set(0,1.0);
    //if(Communicator::nodeid()==0){
    //  b.set(12,1.0);  // fot HOP_TP
    //  b.set(14+24*(1),1.0);  // fot HOP_TP
    //}
    norm_b = b.norm();
  }
  vout.general(vl, "|source|   =%.8f\n", norm_b);
  vout.general(vl, "|source|^2 =%.8f\n", norm_b * norm_b);

  // Source vector
  AFIELD_f ab;
  vout.general(vl, "resizing soruce\n");
  ab.reset(nin, nvol, nex);

  // Solution vector and residual vector
  AFIELD_f ax, az, ar;
  AFIELD_f ay;
  vout.general(vl, "resizing solution\n");
  ax.reset(nin, nvol, nex);
  az.reset(nin, nvol, nex);
  ar.reset(nin, nvol, nex);
  ay.reset(nin, nvol, nex);

  // information of the coarse lattice/blcok on the fine lattice
  vector<int> coarse_lattice(Ndim);
  set_coarse_lattice(coarse_lattice, sap_block_size, vl);
  size_t coarse_nvol = 1;
  for (int i = 0; i < 4; i++) {
    coarse_nvol *= (size_t)coarse_lattice[i];
  }
  std::vector<int> fine_lattice
    = { CommonParameters::Nx(),
        CommonParameters::Ny(),
        CommonParameters::Nz(),
        CommonParameters::Nt() };


  // conversion
  AIndex_lex<float, AFIELD_f::IMPL> index_alt;
  vout.paranoiac(vl, "index is ready\n");
  vout.general(vl, "converting soruce\n");
  //  #pragma omp parallel
  {
    if (afopr_fineF->needs_convert()) {
      vout.detailed(m_vl, "convert in rep. required.\n");
      afopr_fineF->convert(ab, b);
    } else {
      vout.detailed(m_vl, "convert in rep. not required.\n");
      convert(index_alt, ab, b);
    }
  }
  vout.general(vl, "converting source, done\n");

  vout.general(m_vl, "norm check: %f  %f\n", ab.norm2(), b.norm2());


  vout.general("generating mitligrid: num_vectors=%d\n", num_vectors);
  MultiGrid_Domainwall<AFIELD_f, AFIELD_f> multigrid(coarse_lattice, fine_lattice, nin, num_vectors);

  vout.general("hoge: calling set_afopr_fine of MultiGrid_Domainwall\n");
  multigrid.set_afopr_fine(afopr_fineF);

  // test vectors
  std::vector<AFIELD_f>& atestvec = *(multigrid.get_testvectors());
  multigrid.set_testvectors();


  vout.general(m_vl, "Check of testvectors\n");
  vout.general(m_vl, " size of atestvec = %d\n", atestvec.size());
  for (int i = 0; i < num_vectors; ++i) {
    vout.general(m_vl, " norm2 of atestvec[%d] = %f\n",
                 i, atestvec[i].norm2());
  }

  multigrid.gramschmidt();

  vout.general(m_vl, "Check of testvectors\n");
  vout.general(m_vl, " size of atestvec = %d\n", atestvec.size());
  for (int i = 0; i < num_vectors; ++i) {
    vout.general(m_vl, " norm2 of atestvec[%d] = %f\n",
                 i, atestvec[i].norm2());
  }

  //

  // HM added
  params_fopr.set_int_vector("block_size", sap_block_size);
  // up to here

  AFopr_Domainwall_5din_dd<AFIELD_f> *afopr_fineF_dd =
    new AFopr_Domainwall_5din_dd<AFIELD_f>(params_fopr);
  vout.general("hoge: calling set_config for afopr_fineF_dd\n");
  afopr_fineF_dd->set_config(U.get());
  vout.general("hoge: calling set_mode for afopr_fineF_dd\n");
  afopr_fineF_dd->set_mode("D");
  int fine_nin  = afopr_fineF_dd->field_nin();
  int fine_nvol = afopr_fineF_dd->field_nvol();
  int fine_nex  = afopr_fineF_dd->field_nex();

  vout.general(vl, "fine grid operator (float) is ready\n");

  // preparing coase grid operator
  vout.general(vl, "afopr_coarse version\n");
  Parameters params_afopr_coarse;
  params_afopr_coarse.set_int("number_of_testvectors", num_vectors);
  params_afopr_coarse.set_int_vector("coarse_lattice_size", coarse_lattice);
  params_afopr_coarse.set_int("fine_operator_nin", fine_nin);
  params_afopr_coarse.set_int("fine_operator_nvol", fine_nvol);
  params_afopr_coarse.set_int("fine_operator_nex", fine_nex);

  AFopr_Domainwall_coarse<AFIELD_f> *afopr_coarse = new AFopr_Domainwall_coarse<AFIELD_f>();

  afopr_coarse->set_parameters(params_afopr_coarse);
  //afopr_coarse->set_parameters(num_vectors, coarse_lattice,
  //                             fine_nin, fine_nvol, fine_nex);
  AIndex_block_lex<real_t, AFIELD_f::IMPL> block_index(coarse_lattice, fine_lattice);

  Parameters params_fopr_PV = params_fopr;
  double m_PV = params_fopr.get_double("quark_mass_PauliVillars");
  params_fopr_PV.set_double("quark_mass", m_PV);
  AFopr_Domainwall_5din_dd<AFIELD_f> *afopr_fineF_PV =
    new AFopr_Domainwall_5din_dd<AFIELD_f>(params_fopr_PV);
  afopr_fineF_PV->set_mode("D");

  ASolver_SAP_MINRES_dw<AFIELD_f> solver_in_sap_PV(afopr_fineF_PV, &block_index);
  Parameters params_solver_in_sap;
  params_solver_in_sap.set_int("maximum_number_of_iteration", 10);
  params_solver_in_sap.set_int("maximum_number_of_restart", 1);
  params_solver_in_sap.set_double("convergence_criterion_squared", 1.0e-30);
  params_solver_in_sap.set_string("verbose_level", "Detailed");
  solver_in_sap_PV.set_parameters(params_solver_in_sap);

  ASolver_SAP_dw<AFIELD_f> solver_PV(afopr_fineF_PV, &block_index, &solver_in_sap_PV);
  Parameters params_solver_PV;
  params_solver_PV.set_string("verbose_level", "Detailed");
  params_solver_PV.set_int("maximum_number_of_iteration", 1);
  params_solver_PV.set_int("maximum_number_of_restart", 1);
  params_solver_PV.set_double("convergence_criterion_squared", 0.0);
  solver_PV.set_parameters(params_solver_PV);

  vout.general("hoge: calling afopr_coarse->generate_coarse_op\n");
#pragma omp parallel
  {
    //afopr_coarse->generate_coarse_op(afopr_fineF_dd, &solver_PV, atestvec);
    afopr_coarse->generate_coarse_op(afopr_fineF_dd, &solver_in_sap_PV, atestvec);
  }
  vout.general("hoge: calling afopr_coarse->generate_coarse_op, done\n");
  afopr_coarse->set_mode("D");
  vout.general(vl, "afopr_coarse version is ready\n");
  //  afopr_coarse->dump();

  // coarse grid solver
  ASolver_BiCGStab_Cmplx<AFIELD_f> *asolver_coarse =
    new ASolver_BiCGStab_Cmplx<AFIELD_f>(afopr_coarse);
  //  ASolver_GMRES_m_Cmplx<AFIELD_f> *asolver_coarse
  //    = new ASolver_GMRES_m_Cmplx<AFIELD_f>(afopr_coarse);
  afopr_coarse->set_mode("D");
  int   coarse_niter     = 50;
  float coarse_stop_cond = 1e-5;
  int   coarse_nrestart = 1;
  int   coarse_gmres_m = 3;

  //  asolver_coarse->set_parameters(coarse_niter, coarse_stop_cond);
  Parameters params_gmres;
  params_gmres.set_int("maximum_number_of_iteration", coarse_niter);
  params_gmres.set_int("maximum_number_of_restart", coarse_nrestart);
  params_gmres.set_double("convergence_criterion_squared", coarse_stop_cond);
  params_gmres.set_int("number_of_orthonormal_vectors", coarse_gmres_m);
  params_gmres.set_string("verbose_level", "Detailed");
  asolver_coarse->set_parameters(params_gmres);
  int   nconv = -1;
  float diff  = -1.0;


  // coarse source vector
  AFIELD_f bcoarse;
  bcoarse.reset(2 * 2 * num_vectors, coarse_nvol, 1);

  AFIELD_f bcoarse_precond;
  bcoarse_precond.reset(2 * 2 * num_vectors, coarse_nvol, 1);

  // coarse solution vector
  AFIELD_f xcoarse;
  xcoarse.reset(2 * 2 * num_vectors, coarse_nvol, 1);

  AFIELD_f ycoarse;
  ycoarse.reset(2 * 2 * num_vectors, coarse_nvol, 1);


  AFopr_Domainwall_5din<AFIELD_f> *afopr_fineF_tmp = new AFopr_Domainwall_5din<AFIELD_f>(params_fopr_dd);
  afopr_fineF_tmp->set_config(U.get());
  afopr_fineF_tmp->set_mode("D");

  /*
  //  make_random_fine_vector(ax);
  for (int ii = 0; ii < ax.ntot(); ++ii) {
    ax.set(0.0);
    ax.set(ii, 1.0);
    {
      //      const real_t *ptr=ax.ptr(0);
      //      for(int i=0; i<ax.ntot(); ++i){
      //        vout.general("__ax: %d %d %18.10e\n",ii, i, ptr[i]);
      //      }
    }

#pragma omp parallel
    {
      double x2 = ax.norm2();
      vout.general("  (fine) x2:                  %23.15e\n", x2);
      afopr_fineF_dd->mult(ay, ax);
      double y2 = ay.norm2();
      vout.general("  (fine) |full op (dd)  x|^2: %23.15e\n", y2);
      afopr_fineF_tmp->mult(ar, ax);
      double r2 = ar.norm2();
      vout.general("  (fine) |full op (qxs) x|^2: %23.15e\n", r2);
      axpy(ar, real_t(-1.0), ay);
      double diff2 = ar.norm2();
      vout.general("  (file) diff2 = %23.15e\n", diff2);
    }
    {
      const real_t *ptr = ay.ptr(0);
      for (int i = 0; i < ay.ntot(); ++i) {
        if (fabs(ptr[i]) > 1.0e-15) {
          vout.general("__ay: %d %d %18.10e\n", ii, i, ptr[i]);
        }
      }
    }
  }
  */
  // solve on the coarse grid
#pragma omp parallel
  {
    double ab2 = ab.norm2();
    vout.general("calling make_coarse_vector:        input norm2=%22.15e\n", ab2);
    multigrid.make_coarse_vector(bcoarse, ab);
    double bcoarse2 = bcoarse.norm2();
    vout.general("calling make_coarse_vector, done: output norm2=%22.15e\n", bcoarse2);
    afopr_coarse->mult(xcoarse, bcoarse);
    double xcoarse2 = xcoarse.norm2();
    vout.general("  x2: %23.15e\n", xcoarse2);
  }
  // xcoarse.set(0.0);
  // afopr_coarse->mult_dn(0, xcoarse, bcoarse);

  int num_mult = 1000;
  vout.general("applying afopr_coarse %d times\n", num_mult);
#pragma omp parallel
  {
    double xcoarse2 = xcoarse.norm2();

    xcoarse.scal(1.0 / sqrt(xcoarse2));
    Timer timer_mult("mult of coarse operator");
    timer_mult.reset();
    timer_mult.start();
    for (int i = 0; i < num_mult / 2; i++) {
      afopr_coarse->mult(ycoarse, xcoarse);
      afopr_coarse->mult(xcoarse, ycoarse);
    }
    timer_mult.stop();
    xcoarse2 = xcoarse.norm2();
    //vout.general("  x2: %23.15e\n", xcoarse2);
    vout.general("after  %d mult of the coarse op. x2: %23.15e\n", num_mult, xcoarse2);
    timer_mult.report();
    double flop_count = afopr_coarse->flop_count();
    double gflops     = num_mult * flop_count / timer_mult.elapsed_sec() * 1.0e-9;
    vout.general("  coaree op flop_count/mult:  %f\n", flop_count);
    vout.general("  coarse op Flops (float):  %15.7f GFlop/sec\n", gflops);
  }

  vout.general("applying afopr_dd (full op.) %d times\n", num_mult);
#pragma omp parallel
  {
    copy(ax, ab);
    double x2 = ax.norm2();
    ax.scal(1.0 / sqrt(x2));
    Timer timer_mult("mult of full operator");
#pragma omp barrier
    afopr_fineF_dd->mult(ar, ax);
#pragma omp barrier
    double norm2 = ar.norm2();
    vout.general("mult of full op: norm2=%23.15e\n", norm2);
    timer_mult.reset();
    timer_mult.start();
    for (int i = 0; i < num_mult / 2; i++) {
      afopr_fineF_dd->mult(ar, ax);
      afopr_fineF_dd->mult(ax, ar);
    }
    timer_mult.stop();
    //    xcoarse2 = xcoarse.norm2();
    //    vout.general("  x2: %23.15e\n", xcoarse2);
    //    vout.general("after  %d mult of the coarse op. x2: %23.15e\n", num_mult, xcoarse2);
    timer_mult.report();
    double flop_count = afopr_fineF_dd->flop_count();
    double gflops     = num_mult * flop_count / timer_mult.elapsed_sec() * 1.0e-9;
    vout.general("  full op flop_count/mult:  %f\n", flop_count);
    vout.general("  full op Flops (float):  %15.7f GFlop/sec\n", gflops);
  }

  vout.general("applying afopr_dd (sap op.) %d times\n", num_mult);

  //#pragma omp parallel
  {
    copy(ax, ab);
    double x2 = ax.norm2();
    ax.scal(1.0 / sqrt(x2));
    {
      ar.set(0.0);
      afopr_fineF_dd->mult_sap(ar, ax, 0);
      double norm2 = ar.norm2();
      vout.general("mult of SAP op[0]: norm2=%23.15e\n", norm2);
    }
    {
      ar.set(0.0);
      afopr_fineF_dd->mult_sap(ar, ax, 1);
      double norm2 = ar.norm2();
      vout.general("mult of SAP op[1]: norm2=%23.15e\n", norm2);
    }
    {
      ar.set(0.0);
      afopr_fineF_dd->mult_dd(ar, ax);
      double norm2 = ar.norm2();
      vout.general("mult of dd op: norm2=%23.15e\n", norm2);
    }
    {
      ar.set(0.0);
      afopr_fineF_dd->mult(ar, ax);
      afopr_fineF_tmp->mult(ay, ax);
      double norm2     = ar.norm2();
      double norm2_tmp = ay.norm2();
      axpy(ar, -1.0, ay);
      double diff2 = ar.norm2();
      vout.general("mult of full op (dd) : norm2=%23.15e\n", norm2);
      vout.general("mult of full op (qxs): norm2=%23.15e\n", norm2_tmp);
      vout.general("                     : diff2=%23.15e\n", diff2);
    }


    Timer timer_mult("mult of SAP operator");
    timer_mult.reset();
    timer_mult.start();
    for (int i = 0; i < num_mult / 2; i++) {
      afopr_fineF_dd->mult_sap(ar, ax, 0);
      afopr_fineF_dd->mult_sap(ax, ar, 0);
    }
    timer_mult.stop();
    //    xcoarse2 = xcoarse.norm2();
    //    vout.general("  x2: %23.15e\n", xcoarse2);
    //    vout.general("after  %d mult of the coarse op. x2: %23.15e\n", num_mult, xcoarse2);
    timer_mult.report();
    double flop_count = afopr_fineF_dd->flop_count_sap();
    double gflops     = num_mult * flop_count / timer_mult.elapsed_sec() * 1.0e-9;
    vout.general("  SAP op flop_count/mult:  %f\n", flop_count);
    vout.general("  SAP op Flops (float):  %15.7f GFlop/sec\n", gflops);
  }

#pragma omp parallel
  {
    bcoarse_precond.set(0.0);
    afopr_coarse->mult(bcoarse_precond, bcoarse, "prec");
    {
      afopr_coarse->mult(xcoarse, bcoarse_precond, "clov");
      axpy(xcoarse, -1.0, bcoarse);
      double b2=bcoarse.norm2();
      double diff2 = xcoarse.norm2();
      vout.general("checking prec/clov of the coarse operator\n");
      vout.general("      b2 = %23.15e\n", b2);
      vout.general("   diff2 = %23.15e\n", diff2);
    }
    asolver_coarse->solve(xcoarse, bcoarse_precond, nconv, diff);
    double xc2 = xcoarse.norm2();
#pragma omp barrier
    multigrid.make_fine_vector(ax, xcoarse);
#pragma omp barrier
    double ax2 = ax.norm2();

    //    afopr_fineF_dd->set_mode("D");
    //    afopr_fineF_dd->mult(ar,ax);
    afopr_fineF_tmp->set_mode("D");
    ar.set(0.0);
    int nconv_pv=-1;
    real_t diff=-1.0;
    //    solver_PV.solve(az, ax, nconv_pv, diff);
    solver_in_sap_PV.solve(az, ax, nconv_pv, diff);
    afopr_fineF_tmp->mult(ar, az);

    vout.general(" after PV solver: nconv=%d, diff=%23.15e\n", nconv_pv, diff);
    double ab2 = ab.norm2();
    double az2 = az.norm2();
    double ar2 = ar.norm2();
    vout.general("  |b|^2 :               %23.15e\n", ab2);
    vout.general("  |xcoarse|^2 :         %23.15e\n", xc2);
    vout.general("  |ax(fine)|^2:         %23.15e\n", ax2);
    vout.general("  |  PVinv ax(fine)|^2: %23.15e\n", az2);
    vout.general("  |D PVinv ax(fine)|^2: %23.15e\n", ar2);
    axpy(ar, (float)(-1.0), ab);
    double r2 = ar.norm2();
    vout.general("after coarse solver: fine r2 = %23.15e  (relative: %23.15e)\n",
                 r2, r2/ab2);

    using complex_t = AFIELD_f::complex_t;
    const std::vector<AFIELD_f>& testvectors = *(multigrid.get_testvectors());
    for (int i = 0; i < num_vectors; i++) {
      double    norm2 = testvectors[i].norm2();
      complex_t tmp   = dotc(testvectors[i], ar);
      double    re    = real(tmp);
      double    im    = imag(tmp);
      vout.general("  i=%d, <i|i> = %e,  <i|r>=%e %e  (must be almost zero)\n", i, norm2, re, im);
    }
    vout.general("hoge!\n");
  }


  delete afopr_coarse;

  delete asolver_coarse;
  timer_total->report();

  RandomNumberManager::finalize();

  return EXIT_SUCCESS;
}


//============================================================END=====
