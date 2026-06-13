/*!
        @file    testvector_eigen-tmpl.h
        @brief   test for coarse grid operator
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: #$
        @version $LastChangedRevision: 2605 $
*/

#ifndef TESTVECTOR_EIGEN_INCLUDED
#define TESTVECTOR_EIGEN_INCLUDED

#include "lib/Tools/timer.h"
#include "lib/Fopr/fopr.h"
#include "lib/Field/field_F.h"

#include "lib/Eigen/aeigensolver.h"
#include "lib/Eigen/aeigensolver_IRArnoldi.h"
#include "lib/Eigen/aeigensolver_IRLanczos.h"
#include "lib_alt/Eigen/aeigensolver_IRArnoldi_block.h"

#include "lib/Fopr/fopr_Domainwall.h"


#define COMPLEX_t typename ComplexTraits<typename AFIELD::real_t>::complex_t
// to be undefed at the end of this file.


// blocked eigensolver implemented as a functional template
//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
void set_testvectors_eigen(AFopr_dd<AFIELD_d>* fopr_d,
                           std::vector<AFIELD_f>& testvectors,
                           const Parameters& params_eigen,
                           std::string filename)
{
  Bridge::VerboseLevel vl = vout.set_verbose_level("Detailed");

  vout.crucial(vl, "set_testvectors_eigen start.\n");
  vout.increase_indent();

  typedef typename AFIELD_d::real_t real_t;
  typedef typename ComplexTraits<real_t>::complex_t complex_t;
  typedef AFopr<AFIELD_d> AFOPR_d;
  typedef AEigensolver<AFIELD_d, AFOPR_d> AEIGEN;

  int Ndim = CommonParameters::Ndim();

  int num_vectors = testvectors.size();
  //vout.general(vl, "Number of testvectors = %d\n", num_vectors);
  vout.crucial(vl, "Number of testvectors = %d\n", num_vectors);
  
  std::vector<int> fine_lattice(4);
  fine_lattice[0] = CommonParameters::Nx();
  fine_lattice[1] = CommonParameters::Ny();
  fine_lattice[2] = CommonParameters::Nz();
  fine_lattice[3] = CommonParameters::Nt();


  int Nk = params_eigen.get_int("number_of_wanted_eigenvectors");
  int Np = params_eigen.get_int("number_of_working_eigenvectors");

  int Nin  = fopr_d->field_nin();
  int Nvol = fopr_d->field_nvol();
  int Nex  = fopr_d->field_nex();
  std::vector<AFIELD_d> testvectors_d(Nk + Np);
  for(int ivec = 0; ivec < Nk + Np; ++ivec){
    testvectors_d[ivec].reset(Nin, Nvol, Nex);
  }

  std::string fopr_mode = params_eigen.get_string("fopr_mode");
  vout.general("fopr_mode = %s\n", fopr_mode.c_str());
  fopr_d->set_mode(fopr_mode);
  
  std::vector<complex_t> TDa(Nk + Np);

  // setup eigensolver
  AEIGEN* eigen;
  std::string eigensolver_type =
                         params_eigen.get_string("eigensolver_type");
  if (eigensolver_type == "IRArnoldi"){
    eigen = new AEigensolver_IRArnoldi<AFIELD_d, AFOPR_d>(
                                               fopr_d, params_eigen);
  }else if (eigensolver_type == "IRLanczos"){
    eigen = new AEigensolver_IRLanczos<AFIELD_d, AFOPR_d>(
                                               fopr_d, params_eigen);
  }else{
    vout.crucial(vl, "unsupported eigensolver_type = %s\n",
                 eigensolver_type.c_str());
    exit(EXIT_FAILURE);
  }

  int Nsbt  =  0;
  int Nconv = -1;
  AFIELD_d x(Nin, Nvol, Nex);
  x.set(0.0);

  eigen->solve(TDa, testvectors_d, Nsbt, Nconv, x);

  delete eigen;

  // check
  for(int ik = 0; ik < Nk; ++ik){
    fopr_d->mult(x, testvectors_d[ik]);
    complex_t fac = cmplx(-real(TDa[ik]), -imag(TDa[ik]));
    axpy(x, fac, testvectors_d[ik]);
    real_t anorm2 = x.norm2();
    vout.general(vl, " %4d  (%16.8e, %16.8e) %16.8e %12.4e\n",
                 ik, real(TDa[ik]), imag(TDa[ik]), abs(TDa[ik]), anorm2);
  }


  // check of eigenvectors

  std::string file_base = filename;
  write_eigenvectors(fopr_d, TDa, testvectors_d, Nk, file_base);

  // set testvectors
  int num_testvectors = testvectors.size();
  if(num_testvectors > Nk){
    vout.crucial(vl, "too large num_testvectors\n");
    exit(EXIT_FAILURE);
  }

  //for(int ik = 0; ik < num_testvectors; ++ik){
  //  copy(testvectors[ik], testvectors_d[ik]);
  //}

  vout.decrease_indent();
  vout.crucial(vl, "set_testvectors_eigen finished.\n");
  Communicator::sync();

}

//====================================================================
template<typename AFIELD>
void write_eigenvectors(AFopr<AFIELD>* fopr,
                        std::vector<COMPLEX_t>& TDa,
                        std::vector<AFIELD>& testvectors,
                        int Nk,
                        std::string filename)
{
  Bridge::VerboseLevel vl = vout.set_verbose_level("Detailed");

  int Ndim = CommonParameters::Ndim();

  int Nin  = testvectors[0].nin();
  int Nvol = testvectors[0].nvol();
  int Nex  = testvectors[0].nex();

  if(Nk > testvectors.size()){
    vout.crucial("size of testvectors mismatch.\n");
    exit(EXIT_FAILURE);
  }

  // output eigenvectors

  vout.detailed(vl, "eigenvectors are output.\n");

  unique_ptr<FieldIO> field_io(new FieldIO_Binary(IO_Format::Trivial));
  int NinF = 2 * CommonParameters::Nc() * CommonParameters::Nd();
  int NexF = Nin / NinF;
  Field xout(NinF, Nvol, NexF);

  std::string file_base = filename;

  if (fopr->needs_convert()) {

    vout.detailed(vl, "convert in rep. required.\n");
    for(int ik = 0; ik < Nk; ++ik){
      vout.detailed(vl, "norm of testvec = %e\n", testvectors[ik].norm2());
      fopr->reverse(xout, testvectors[ik]);
      vout.detailed(vl, "convert finished: %d\n", ik);

      std::string filename = FileUtils::generate_filename(
                               "%s_%04d.dat", file_base.c_str(), ik);
      field_io->write_file(xout, filename);
      vout.detailed(vl, "output finished: %d\n", ik);
    }

  } else {

    vout.detailed(vl, "convert in rep. not required.\n");
    for(int ik = 0; ik < Nk; ++ik){
      AIndex_lex<real_t, AFIELD::IMPL> index_alt;
      reverse(index_alt, xout, testvectors[ik]);

      std::string filename = FileUtils::generate_filename(
                               "%s_%04d.dat", file_base.c_str(), ik);
      field_io->write_file(xout, filename);
    }
  }

  // output eigenvalues

  for(int ik = 0; ik < Nk; ++ik){

    if(Communicator::nodeid() == 0){
      if(ThreadManager::get_thread_id() == 0){

        std::string filename2 = FileUtils::generate_filename(
                            "%s_%04d.eigen", file_base.c_str(), ik);

        FILE* fp = fopen(filename2.c_str(), "w");
        fprintf(fp, " %20.14e  %20.14e\n",
                real(TDa[ik]), imag(TDa[ik]));
        fclose(fp);

	/*
	std::fstream write_file(filename2.c_str(), std::ios::out);
        if (!write_file.is_open()) {
          vout.crucial(vl, "Error at opening file %s\n",
		       filename2.c_str());
          exit(EXIT_FAILURE);
        }
	write_file << real(TDa[ik]) << "  " << imag(TDa[ik]) << std::endl;
        write_file.close();
	*/
      }
    }

  }

  vout.crucial(vl, "write_testvectors finished.\n");
  Communicator::sync();

}

//====================================================================
template<typename AFIELD>
void read_eigenvectors(AFopr<AFIELD>* fopr,
                       std::vector<COMPLEX_t>& TDa,
                       std::vector<AFIELD>& testvectors,
                       std::string filename)
{
  typedef typename AFIELD::real_t real_t;
  typedef typename ComplexTraits<real_t>::complex_t complex_t;

  Bridge::VerboseLevel vl = vout.set_verbose_level("Detailed");

  vout.general(vl, "read_testvectors start.\n");

  int Ndim = CommonParameters::Ndim();

  int Nin  = testvectors[0].nin();
  int Nvol = testvectors[0].nvol();
  int Nex  = testvectors[0].nex();
  int Nk   = testvectors.size();

  unique_ptr<FieldIO> field_io(new FieldIO_Binary(IO_Format::Trivial));
  int NinF = 2 * CommonParameters::Nc() * CommonParameters::Nd();
  int NexF = Nin / NinF;

  Field xin(NinF, Nvol, NexF);
  AFIELD x(Nin, Nvol, Nex);

  //std::string file_base = "testvector";
  std::string file_base = filename;

  vout.detailed(vl, "  Nk = %d\n", Nk);
  vout.detailed(vl, "  Nin = %d  Nvol = %d  Nex = %d\n", Nin, Nvol, Nex);
  vout.detailed(vl, "  NinF = %d  Nvol = %d  NexF = %d\n", NinF, Nvol, NexF);

  // read eigenvectors

  if (fopr->needs_convert()) {

    vout.detailed(vl, "convert in rep. required.\n");
    for(int ik = 0; ik < Nk; ++ik){
      vout.detailed(vl, "ik = %d\n", ik);
      Communicator::sync();
      std::string filename = FileUtils::generate_filename(
                               "%s_%04d.dat", file_base.c_str(), ik);
      //vout.detailed(vl, "filename = %s\n", filename.c_str());
      field_io->read_file(xin, filename);
#pragma omp parallel
      {
       fopr->convert(testvectors[ik], xin);
      }
    }

  } else {

    vout.detailed(vl, "convert in rep. not required.\n");
    for(int ik = 0; ik < Nk; ++ik){
      std::string filename = FileUtils::generate_filename(
                               "%s_%04d.dat", file_base.c_str(), ik);
      //vout.detailed(vl, "filename = %s\n", filename.c_str());
      field_io->read_file(xin, filename);
      AIndex_lex<real_t, AFIELD::IMPL> index_alt;
#pragma omp parallel
      {
       convert(index_alt, testvectors[ik], xin);
      }
    }
  }

  // read eigenvalues
#pragma omp parallel
  if(ThreadManager::get_thread_id() == 0){

    for(int ik = 0; ik < Nk; ++ik){
      double TDa1[2];

      if(Communicator::nodeid() == 0){
        std::string filename2 = FileUtils::generate_filename(
                              "%s_%04d.eigen", file_base.c_str(), ik);

        std::fstream read_file(filename2.c_str(), std::ios::in);
        if (!read_file.is_open()) {
          vout.crucial(vl, "Error at opening file %s\n",
                       filename2.c_str());
          Communicator::abort();
        }
	read_file >> TDa1[0] >> TDa1[1];
        read_file.close();
      }
      Communicator::broadcast(2, TDa1, 0);

      TDa[ik] = cmplx(real_t(TDa1[0]), real_t(TDa1[1]));
    }

  }

  vout.detailed(vl, "read_testvectors finished.\n");

}

//====================================================================
template<typename AFIELD>
void check_eigenvectors(AFopr<AFIELD>* fopr,
                        std::vector<COMPLEX_t>& TDa,
                        std::vector<AFIELD>& testvectors)
{
  Bridge::VerboseLevel vl = vout.set_verbose_level("Detailed");

  vout.detailed(vl, "check_eigen_testvectors start\n");

  typedef typename AFIELD::real_t real_t;
  typedef typename ComplexTraits<real_t>::complex_t complex_t;

  int Ndim = CommonParameters::Ndim();

  int Nin  = testvectors[0].nin();
  int Nvol = testvectors[0].nvol();
  int Nex  = testvectors[0].nex();
  int Nk   = testvectors.size();

  AFIELD x(Nin, Nvol, Nex);

  vout.general(vl, "   ik                   TDa                  "
                   "      abs           diff2\n");
#pragma omp parallel
  for(int ik = 0; ik < Nk; ++ik){

    fopr->mult(x, testvectors[ik]);

    complex_t fac = cmplx(-real(TDa[ik]), -imag(TDa[ik]));
    axpy(x, fac, testvectors[ik]);

    real_t anorm2 = x.norm2();

    vout.general(vl, " %4d  (%16.8e, %16.8e) %16.8e %12.4e\n",
                 ik, real(TDa[ik]), imag(TDa[ik]), abs(TDa[ik]), anorm2);
  }

  vout.detailed(vl, "check_eigen_testvectors finished\n");

}

#undef COMPLEX_t

#endif
//============================================================END=====
