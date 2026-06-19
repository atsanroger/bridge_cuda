/*!
      @File    eigenvalue_alt-tmpl.h
      @brief
      @author  Hideo Matsufuru
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-03-22 13:04:01 #$
      @version $LastChangedRevision: 2503 $
*/

#include "eigenvalue_alt.h"
#include "test.h"

#include "lib/Field/field_G.h"
#include "lib/Fopr/afopr.h"
#include "lib/Fopr/fopr.h"

#include "lib/Eigen/aeigensolver.h"
#include "lib/Eigen/aeigensolver_IRLanczos.h"
#include "lib/Eigen/aeigensolver_IRArnoldi.h"

// class name

template<typename AFIELD>
const std::string Eigenvalue_alt<AFIELD>::class_name = "Eigenvalue_alt";

//====================================================================
template<typename AFIELD>
void Eigenvalue_alt<AFIELD>::init()
{
 // do nothing
}

//====================================================================
template<typename AFIELD>
int Eigenvalue_alt<AFIELD>::measure(std::string file_params,
                                  std::string run_mode)
{
  typedef AFopr<AFIELD> AFOPR;
  typedef AEigensolver<AFIELD,AFOPR> AEIGENSOLVER;

  const std::string test_name = class_name + ".eigenspectrum";

  unique_ptr<Timer> timer(new Timer());
  timer->start();

  vout.general("\n");
  vout.general("-------------------------------------------------"
               "-------------------\n");
  vout.general("test name: %s\n", test_name.c_str());
  vout.general("parameter file: %s\n", file_params.c_str());
  vout.general("run_mode: %s\n", run_mode.c_str());
  vout.general("\n");

  string inputfile = "input";  // configuration numbers are stored
  if(run_mode == "job")
    vout.general("input file: %s\n", inputfile.c_str());

  int Nc   = CommonParameters::Nc();
  int Nd   = CommonParameters::Nd();
  int Ndim = CommonParameters::Ndim();
  int Nvol = CommonParameters::Nvol();

  // ####  parameter setup  ####
  params_all = ParameterManager::read(file_params);

  Parameters params_test = params_all.lookup("Test_Eigensolver");

  string str_vlevel = params_test.get_string("verbose_level");
  Bridge::VerboseLevel m_vl = vout.set_verbose_level(str_vlevel);
  vout.general(m_vl, "  vlevel       = %s\n", str_vlevel.c_str());

  bool do_check = false;
  double expected_result = 0.0;
  if(run_mode != "job"){
    do_check = params_test.is_set("expected_result");
    if(do_check) expected_result = params_test.get_double("expected_result");
  }

  // Field for gauge configuration
  U.reset(new Field_G(Nvol, Ndim));

  // ### Fopr setup ###
  Parameters params_fopr = params_all.lookup("Fopr");
  string fopr_type = params_fopr.get_string("fermion_type");
  vout.general(m_vl, "\n");
  vout.general(m_vl, "fermion_type = %s\n", fopr_type.c_str());
  unique_ptr<AFopr<AFIELD> > fopr(
                            AFopr<AFIELD>::New(fopr_type, params_fopr));

  // ### Eigensolver setup ###
  Parameters params_eigen = params_all.lookup("Eigensolver");
  string eigensolver_type = params_eigen.get_string("eigensolver_type");
  string fermion_mult_mode = params_eigen.get_string("fermion_mult_mode");
  // string eigensolver_type = "IRArnoldi";
  // string mult_mode = "D";
  vout.general(m_vl, "\n");
  vout.general(m_vl, "eigensolver_type = %s\n", eigensolver_type.c_str());
  vout.general(m_vl, "fermion_mult_mode = %s\n", fermion_mult_mode.c_str());
  //unique_ptr<AEIGENSOLVER> eigen(
  //                            AEIGENSOLVER::New(eigensolver_type, fopr));
  // eigen->set_parameters(params_eigen);
  unique_ptr<AEIGENSOLVER> eigen(
            AEIGENSOLVER::New(eigensolver_type, fopr.get(), params_eigen));

  // ### eigenmodes preparation ###
  int Nk = params_eigen.get_int("number_of_wanted_eigenvectors");
  int Np = params_eigen.get_int("number_of_working_eigenvectors");
  int Nm = Nk + Np;
  std::vector<dcomplex> TDa(Nm);
  std::vector<AFIELD>  vk(Nm);

  int NFin  = fopr->field_nin();
  int NFvol = fopr->field_nvol();
  int NFex  = fopr->field_nex();
  for (int k = 0; k < Nm; ++k) {
    vk[k].reset(NFin, NFvol, NFex);
  }

  AFIELD b2(NFin, NFvol, NFex);
  // b2.set(0.0);
  b2.set(0, 1.0);

  timer->stop();
  vout.general(m_vl, "\n");
  vout.general(m_vl, "object setup finished\n");
  timer->report();
  vout.general(m_vl, "\n");

  timer->reset();

  // input file
  int iconf = 1;
  int nconf = 1;
  if(run_mode == "job"){
    read_input(iconf, nconf, inputfile);
    vout.general(m_vl, "initial iconf   = %d\n", iconf);
    vout.general(m_vl, "number of confs = %d\n", nconf);
    vout.general(m_vl, "\n");
  }

  string filename_base = params_test.get_string("config_filename_input");
  string data_filename_base;
  if(run_mode == "job"){
     data_filename_base = params_test.get_string("data_filename_output");
  }
  double result = 0.0;

  // #### iconf loop: measurement main part ####
  for(int iconf2 = 0; iconf2 < nconf; ++iconf2){

    timer->start();
    vout.general(m_vl, "--------------------------------------"
                       "-----------\n");
    vout.general(m_vl, "iconf = %d\n", iconf);

    if(run_mode == "job"){
      string filename = filename_config(filename_base, iconf);
      params_test.set_string("config_filename_input", filename);
    }
    setup_config(U, params_test); // read config

    for (int k = 0; k < Nm; ++k) {
      vk[k].set(0.0);
      TDa[k] = cmplx(0.0, 0.0);
    }

    fopr->set_config(U.get());
    // fopr->set_mode("D");
    fopr->set_mode(fermion_mult_mode);

    int Nsbt  = -1;
    int Nconv = -100;
    vout.general(m_vl, "\n");
    eigen->solve(TDa, vk, Nsbt, Nconv, b2);

    AFIELD v(NFin, NFvol, NFex);
    double vv = 0.0;  // superficial initialization

    vout.general(m_vl, "Eigenvalues:\n");
    vout.general(m_vl, " index                 lambda    "
                       "                abs(lambda)        diff2\n");
    for (int i = 0; i < Nsbt + 1; ++i) {
      fopr->mult(v, vk[i]);
      axpy(v, -TDa[i], vk[i]);  // v -= TDa[i] * vk[i];
      vv = v.norm2();           // vv = v * v;
      vout.general(m_vl, " %4d  (%17.14f, %17.14f)  %16.14f  %10.4e\n",
                   i, real(TDa[i]), imag(TDa[i]), abs(TDa[i]), vv);
    }

    // output to file
    if(run_mode == "job"){
      string filename_data = filename_config(data_filename_base, iconf);
      vout.general("  output file = %s\n", filename_data.c_str());
      vout.init(filename_data);

      vout.general(m_vl, "  Nsbt = %4d\n", Nsbt);
      for(int i = 0; i < Nsbt + 1; ++i){
        vout.general(m_vl, " %4d  (%17.14f, %17.14f)  %16.14f\n",
                     i, real(TDa[i]), imag(TDa[i]), abs(TDa[i]));
      }

      vout.unset();
      vout.init("stdout");
      vout.general(m_vl, "  output to file done.\n");
    }

    vout.general(m_vl, "\n");

    result = abs(TDa[0]);

    ++iconf;
    if(run_mode == "job") write_input(iconf, nconf, inputfile);

    timer->stop();

  }

  timer->report();

  if(run_mode == "job"){
    return EXIT_SUCCESS;
  }else{
   if(do_check){
      return Test::verify(result, expected_result);
    }else{
      vout.detailed(m_vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }
  }

}

//====================================================================
template<typename AFIELD>
int Eigenvalue_alt<AFIELD>::eigenvalue_Arnoldi(std::string file_params)
{
  return measure(file_params, "test");
}

//====================================================================
template<typename AFIELD>
int Eigenvalue_alt<AFIELD>::eigenvalue_Lanczos(std::string file_params)
{
  return measure(file_params, "test");
}
//============================================================END=====
