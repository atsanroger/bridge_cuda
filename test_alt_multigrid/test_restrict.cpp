/*!
        @file    test_restrict.cpp

        @brief   test for restrition

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: #$

        @version $LastChangedRevision: 2533 $
*/

#include "test_restrict.h"

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
#include "lib_alt_SIMD/Field/aindex_block_lex.h"
#include "lib_alt_SIMD/Field/aindex_coarse_lex.h"
#include "lib_alt_SIMD/Fopr/afopr_Clover.h"
#include "lib_alt_SIMD/Fopr/afopr_Clover_dd.h"
#include "lib_alt_SIMD/Fopr/afopr_Clover_coarse.h"
const Impl IMPL = SIMD;
#endif

#ifdef USE_ALT_QXS
#include "lib_alt_QXS/bridge_alt_qxs.h"
#include "lib_alt_QXS/Field/aindex_block_lex.h"
#include "lib_alt_QXS/Field/aindex_coarse_lex.h"
#include "lib_alt_QXS/Fopr/afopr_Clover.h"
#include "lib_alt_QXS/Fopr/afopr_Clover_dd.h"
#include "lib_alt_QXS/Fopr/afopr_Clover_coarse.h"
const Impl IMPL = QXS;
#endif

#ifdef USE_ALT_ACCEL
#include "lib_alt_Accel/bridge_alt_accel.h"
#include "lib_alt_Accel/Field/aindex_block_lex.h"
#include "lib_alt_Accel/Field/aindex_coarse_lex.h"
#include "lib_alt_Accel/Fopr/afopr_Clover.h"
#include "lib_alt_Accel/Fopr/afopr_Clover_dd.h"
#include "lib_alt_Accel/Fopr/afopr_Clover_coarse.h"
const Impl IMPL = ACCEL;
#endif

typedef float real_t;


#include "lib_alt/Solver/MultiGrid_Clover.h"

//#include "lib_alt/Solver/asolver.h"
//#include "lib_alt/Solver/asolver_BiCGStab.h"
//#include "lib_alt/Solver/asolver_BiCGStab_Cmplx.h"
//#include "lib_alt/Solver/asolver_CG.h"



// Registration as a test
namespace test_restrict {
  const string test_name = "Test_restrict";

  int test(void)
  {
    Test_Restrict test;
    return test.test("restrict");
  }
}

// Registration as a test
namespace test_prolong {
  const string test_name = "Test_prolong";

  int test(void)
  {
    Test_Restrict test;
    return test.test("prolong");
  }
}

// class name
const std::string Test_Restrict::class_name =
  "Test_Restrict";

//====================================================================
void Test_Restrict::init()
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
        fflush(stdout);
        fflush(stderr);
        exit(EXIT_FAILURE);
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


  void make_random_fine_vector_keep(AFIELD_f& vec, AFopr<AFIELD_f> *afopr_fineF)
  {
    Field_F tmp_vec;
    //    RandomNumberManager::initialize("Mseries", 1234567UL);

    const int             seed = 54321;
    RandomNumbers_Mseries random(seed);
    random.uniform_lex_global(tmp_vec);

    // conversion
    AIndex_lex<float, AFIELD_f::IMPL> index_alt;
    //    vout.paranoiac(vl, "index is ready\n");
    //    vout.general(vl, "converting soruce\n");
#pragma omp parallel
    {
      if (afopr_fineF->needs_convert()) {
        vout.detailed("convert in rep. required.\n");
        afopr_fineF->convert(vec, tmp_vec);
      } else {
        vout.detailed("convert in pep. not required.\n");
        convert(index_alt, vec, tmp_vec);
      }

      double norm2_tmp = tmp_vec.norm2();
      double norm2_vec = vec.norm2();
      vout.general("norm check: Field, AFeild : %f,  %f\n", norm2_tmp, norm2_vec);
    } // omp parallel
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

              //              int isite = index.site(x%Nx,y%Ny,z%Nz,t%Nt);
              int isite = (x % Nx) + Nx * ((y % Ny) + Ny * ((z % Nz) + Nz * (t % Nt)));
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
} // anonymous namespace

//====================================================================
int Test_Restrict::test(const std::string mode)
{
  const std::string filename_input = "test_alt_Restrict.yaml";

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

  Parameters params_mg_solver = params_all.lookup("MGSolver");
  //  const int    Nmult            = params_test.get_int("number_of_mult");
  const string str_vlevel = params_mg_solver.get_string("verbose_level");

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

  /*
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

  //  bool run_fopr=false;
  //  bool run_afopr=false;
  //  bool run_comparison=false;
  */

  // fine vector
  AFIELD_f fine_vec;
  fine_vec.reset(2 * Nc * Nd, Nvol, 1);
  vout.general("fine vec: nin=%d, nvol=%d, nex=%d\n", fine_vec.nin(), fine_vec.nvol(), fine_vec.nex());

  // Solution vector
  //  AFIELD_f ax;
  //  vout.general(vl, "resizing solution\n");
  //  ax.reset(2*Nc*Nd, Nvol,1);

  std::vector<int> fine_lattice
    = { CommonParameters::Nx(),
        CommonParameters::Ny(),
        CommonParameters::Nz(),
        CommonParameters::Nt() };

  // information of the coarse lattice/blcok on the fine lattice
  vector<int> coarse_lattice(Ndim);
  set_coarse_lattice(coarse_lattice, sap_block_size, vl);
  size_t coarse_nvol = 1;
  for (int i = 0; i < 4; i++) {
    coarse_nvol *= (size_t)coarse_lattice[i];
  }

  const int coarse_nin = 4 * num_vectors; // 4 is for complex x chirality
  vout.general("hoge: nin=%d, nvol=d, nex=%d\n",
               coarse_nin, coarse_nvol, 1);
  fflush(0);
// coarse vector
  AFIELD_f coarse_vec;
  coarse_vec.reset(coarse_nin, coarse_nvol, 1);
  vout.general("coarse vec: nin=%d, nvol=d, nex=%d\n",
               coarse_vec.nin(), coarse_vec.nvol(), coarse_vec.nex());
  fflush(0);
  AIndex_block_lex<typename AFIELD_f::real_t, AFIELD_f::IMPL>
  block_index(coarse_lattice, fine_lattice);


  // fopr is only to convert fine spinor
  string fopr_type    = params_fopr.get_string("fermion_type");
  string fopr_type_dd = fopr_type + "_dd";
  params_fopr.set_int_vector("block_size", sap_block_size);

  // no factor for _dd operator
  unique_ptr<AFopr_Clover_dd<AFIELD_f> > afopr_fineF(new AFopr_Clover_dd<AFIELD_f>(params_fopr));

  //  unique_ptr< AFopr_dd<AFIELD_f> > afopr_fineF(AFopr_dd<AFIELD_f>::New(fopr_type, params_fopr));
  //  AFopr_dd<AFIELD_f> *afopr_dd = dynamic_cast<AFopr_dd<AFIELD_f>* >(afopr_fineF.get());
  //  AFopr_dd<AFIELD_f>  *afopr_dd = dynamic_cast<AFopr_dd<AFIELD_f>* >(afopr_fineF.get());
  //  if(afopr_dd == nullptr){
  //    vout.crucial(m_vl, "failed to generated a _dd type operator\n");
  //  }

  vout.general("generating mitligrid: num_vectors=%d\n", num_vectors);
  MultiGrid_Clover<AFIELD_f, AFIELD_f> multigrid(coarse_lattice, fine_lattice, 2 * Nc * Nd, num_vectors);

  vout.general("calling set_afopr_fine of MultiGrid_Clover\n");
  multigrid.set_afopr_fine(afopr_fineF.get());

  // test vectors
  multigrid.set_testvectors();
  std::vector<AFIELD_f>& atestvec = *(multigrid.get_testvectors());

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
    float *tmp = atestvec[i].ptr(0);
    for (int j = 0; j < 600; ++j) {
      vout.general("  %d:  %13.8e\n", j, tmp[j]);
    }
  }
  //

  double result = -1.0;
  if (mode == "restrict") {
    make_random_fine_vector(fine_vec);
    //    make_random_fine_vector(fine_vec, afopr_fineF.get());
    //    fine_vec.set(0,1.0);
#pragma omp parallel
    {
      multigrid.make_coarse_vector(coarse_vec, fine_vec);
      double fine2   = fine_vec.norm2();
      double coarse2 = coarse_vec.norm2();
      vout.general("calling make_coarse_vector, done\n");
      vout.general("  |fine|^2:    %24.15e\n", fine2);
      vout.general("  |coarse|^2:  %24.15e\n", coarse2);
#pragma omp master
      {
        result = coarse2;
//        float *tmp=coarse_vec.ptr(0);
//        int nodeid=Communicator::nodeid();
//        int num_nodes=Communicator::size();
//        for(int node=0; node<num_nodes; ++node){
//          fflush(stdout);
//          Communicator::sync();
//          if(nodeid==node){
//            for(int j=0; j<coarse_vec.size(); ++j){
//              int id=nodeid*coarse_vec.size()+j;
//              printf(" hoge_ %d  %15.10e : %d %d\n", id,  tmp[j],node, j);
//            }
//          }
//          fflush(stdout);
//          Communicator::sync();
//       }
      }
    } // omp parallel
  } else if (mode == "prolong") {
    make_random_coarse_vector(coarse_vec, coarse_lattice, num_vectors);
#pragma omp parallel
    {
      multigrid.make_fine_vector(fine_vec, coarse_vec);
      double fine2   = fine_vec.norm2();
      double coarse2 = coarse_vec.norm2();
      vout.general("calling make_fine_vector, done\n");
      vout.general("  |fine|^2:    %24.15e\n", fine2);
      vout.general("  |coarse|^2:  %24.15e\n", coarse2);
#pragma omp master
      {
        result = fine2;
      }
    }
    float *tmp      = coarse_vec.ptr(0);
    int   nodeid    = Communicator::nodeid();
    int   num_nodes = Communicator::size();
    for (int node = 0; node < 1; ++node) {
      fflush(stdout);
      Communicator::sync();
      if (nodeid == node) {
        printf(" coarse\n");
        for (int j = 0; j < coarse_vec.size(); ++j) {
          int id = nodeid * coarse_vec.size() + j;
          printf(" hoge_ %d  %15.10e : %d %d\n", id, tmp[j], node, j);
        }
      }
      fflush(stdout);
      Communicator::sync();
    }

    for (int node = 0; node < 1; ++node) {
      float *tmp = fine_vec.ptr(0);
      fflush(stdout);
      Communicator::sync();
      if (nodeid == node) {
        printf(" fine\n");
        for (int j = 0; j < fine_vec.size(); ++j) {
          int id = nodeid * fine_vec.size() + j;
          printf(" fuga_ %d  %15.10e : %d %d\n", id, tmp[j], node, j);
        }
      }
      fflush(stdout);
      Communicator::sync();
    }
  }
  timer_total->report();

  RandomNumberManager::finalize();

  return EXIT_SUCCESS;
}


//============================================================END=====
