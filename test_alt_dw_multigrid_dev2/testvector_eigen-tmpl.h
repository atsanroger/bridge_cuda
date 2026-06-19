/*!
        @file    testvector_eigen-tmpl.h
        @brief   test for coarse grid operator
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: #$
        @version $LastChangedRevision: 2599 $
*/

#ifndef TESTVECTOR_EIGEN_INCLUDED
#define TESTVECTOR_EIGEN_INCLUDED

#include "lib/Tools/timer.h"
#include "lib/Fopr/fopr.h"
#include "lib/Field/field_F.h"

#include "lib_alt/Eigen/aeigensolver_IRArnoldi_block.h"

// blocked eigensolver implemented as a functional template
//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
void set_testvectors_eigen(AFopr_dd<AFIELD_d>* fopr_d,
                           std::vector<AFIELD_f>& testvectors,
                           std::vector<int>& coarse_lattice,
			   const Parameters& params_eigen)
{
  Bridge::VerboseLevel vl = vout.set_verbose_level("Detailed");

  vout.crucial(vl, "Setting testvectors with eigensolver\n");

  typedef typename AFIELD_d::real_t real_t;
  typedef typename ComplexTraits<real_t>::complex_t complex_t;
  typedef AFopr<AFIELD_d> AFOPR_d;
  typedef AEigensolver_IRArnoldi_block<AFIELD_d, AFOPR_d> AEIGEN;

  int Ndim = CommonParameters::Ndim();

  int num_vectors = testvectors.size();
  //vout.general(vl, "  Number of testvectors = %d\n", num_vectors);
  vout.crucial(vl, "  Number of testvectors = %d\n", num_vectors);
  
  std::vector<int> fine_lattice(4);
  fine_lattice[0] = CommonParameters::Nx();
  fine_lattice[1] = CommonParameters::Ny();
  fine_lattice[2] = CommonParameters::Nz();
  fine_lattice[3] = CommonParameters::Nt();

  int num_blocks = 1;
  for(int mu = 0; mu < 4; ++mu){
    vout.general(vl, "  fine_lattice[%d] = %d  coarse_lattice[%d] = %d\n",
		 mu, fine_lattice[mu], mu, coarse_lattice[mu]);
    num_blocks *= coarse_lattice[mu];
  }
  //vout.general(vl, "  Number of blocks = %d\n", num_vectors);
  vout.crucial(vl, "  Number of blocks = %d\n", num_vectors);

  // int    Nk = num_vectors + 4;
  // int    Np = 16;
  // int    Niter_eigen = 100;
  // double Enorm_eigen = 1.e-24;
  // double Vthreshold  = 1.0;
  /*
  Parameters params_eigen;
  params_eigen.set_string("eigensolver_mode", "abs_ascending");
  params_eigen.set_int("number_of_wanted_eigenvectors", Nk);
  params_eigen.set_int("number_of_working_eigenvectors", Np);
  params_eigen.set_int("maximum_number_of_iteration", Niter_eigen);
  params_eigen.set_double("convergence_criterion_squared", Enorm_eigen);
  params_eigen.set_double("threshold_value", Vthreshold);
  params_eigen.set_int_vector("fine_lattice", fine_lattice);
  params_eigen.set_int_vector("coarse_lattice", coarse_lattice);
  params_eigen.set_string("verbose_level", "Detailed");
  */
  int Nk = params_eigen.get_int("number_of_wanted_eigenvectors");
  int Np = params_eigen.get_int("number_of_working_eigenvectors");

  int Nin  = fopr_d->field_nin();
  int Nvol = fopr_d->field_nvol();
  int Nex  = fopr_d->field_nex();
  std::vector<AFIELD_d> testvectors_d(Nk + Np);
  for(int ivec = 0; ivec < Nk + Np; ++ivec){
    testvectors_d[ivec].reset(Nin, Nvol, Nex);
  }
  
  std::vector<complex_t> TDa((Nk + Np) * num_blocks);

  AEIGEN* eigen = new AEIGEN(fopr_d, params_eigen);

  int Nsbt  =  0;
  int Nconv = -1;
  AFIELD_d x(Nin, Nvol, Nex);
  x.set(0.0);

  eigen->solve(TDa, testvectors_d, Nsbt, Nconv, x);

  delete eigen;

  // check of eigenvectors
/*
  AIndex_block_lex<real_t, AFIELD_d::IMPL>
                     index_block(coarse_lattice, fine_lattice);

  std::vector<complex_t> TDa1(num_blocks);
  std::vector<real_t> anorm2(num_blocks);

  vout.general(vl, "    ik    diff2(ave)  "
                   "  diff2(max)          TDa[0]\n");
  for(int ik = 0; ik < Nk; ++ik){

    fopr_d->mult_dd(x, testvectors_d[ik]);

    for(int block = 0; block < num_blocks; ++block){
      TDa1[block] = TDa[ik + (Nk+Np)*block];
    }
    block_axpy(x, &TDa1[0], testvectors_d[ik], real_t(-1.0), index_block);

    block_norm2(&anorm2[0], x, index_block);

    real_t anorm2_ave = 0.0;
    real_t anorm2_max = 0.0;
    for(int block = 0; block < num_blocks; ++block){
      anorm2_ave += anorm2[block];
      anorm2_max = std::max(anorm2[block], anorm2_max);
    }
    anorm2_ave = Communicator::reduce_sum(anorm2_ave);
    anorm2_max = Communicator::reduce_max(anorm2_max);
    anorm2_ave = anorm2_ave/real_t(num_blocks * CommonParameters::NPE());
    vout.general(vl, "  %4d   %14.6e  %14.6e   (%14.6e, %14.6e)\n",
                 ik, anorm2_ave, anorm2_max, real(TDa1[0]), imag(TDa1[0]));
  }
*/

  // field output
  {
    vout.crucial(vl, "eigenvectors are output.\n");

    unique_ptr<FieldIO> field_io(new FieldIO_Binary(IO_Format::Trivial));
    int NinF = 2 * CommonParameters::Nc() * CommonParameters::Nd();
    int NexF = Nin / NinF;
    Field xout(NinF, Nvol, NexF);

    std::string file_base = "testvector";

    if (fopr_d->needs_convert()) {

      //vout.detailed(vl, "convert in rep. required.\n");
      vout.crucial(vl, "convert in rep. required.\n");
      for(int ik = 0; ik < Nk; ++ik){
	//for(int ik = Nk-1; ik >= 0; --ik){
        vout.crucial(vl, "norm of testvec = %e\n", testvectors_d[ik].norm2());
        fopr_d->reverse(xout, testvectors_d[ik]);
        vout.crucial(vl, "convert finished: %d\n", ik);
        std::string filename = FileUtils::generate_filename(
                                 "%s_%04d.dat", file_base.c_str(), ik);
        vout.crucial(vl, "file = %s\n", filename.c_str());
        vout.crucial(vl, "norm of xout = %e\n", xout.norm2());
        field_io->write_file(xout, filename);
        vout.crucial(vl, "output finished: %d\n", ik);
      }

    } else {

      vout.detailed(vl, "convert in rep. not required.\n");
      for(int ik = 0; ik < Nk; ++ik){
        AIndex_lex<real_t, AFIELD_d::IMPL> index_alt;
        reverse(index_alt, xout, testvectors_d[ik]);
        std::string filename = FileUtils::generate_filename(
                                 "%s_%4d.dat", file_base.c_str(), ik);
        field_io->write_file(xout, filename);
      }
    }

  } // field output

  Communicator::sync();

  // set testvectors in float
  int num_testvectors = testvectors.size();
  for(int ik = 0; ik < Nk; ++ik){
        copy(testvectors[ik], testvectors_d[ik]);
  }

  Communicator::sync();

  // field input
  /*
  {
    vout.crucial(vl, "eigenvectors are output.\n");

    unique_ptr<FieldIO> field_io(new FieldIO_Binary(IO_Format::Trivial));
    int NinF = 2 * CommonParameters::Nc() * CommonParameters::Nd();
    int NexF = Nin / NinF;
    Field xin(NinF, Nvol, NexF);

    std::string file_base = "testvector";

    if (fopr_d->needs_convert()) {

      vout.detailed(vl, "convert in rep. required.\n");
      for(int ik = 0; ik < Nk; ++ik){
        std::string filename = FileUtils::generate_filename(
                                 "%s_%04d.dat", file_base.c_str(), ik);
        field_io->read_file(xin, filename);
        fopr_d->convert(testvectors_d[ik], xin);
      }

    } else {

      vout.detailed(vl, "convert in rep. not required.\n");
      for(int ik = 0; ik < Nk; ++ik){
        std::string filename = FileUtils::generate_filename(
                                 "%s_%4d.dat", file_base.c_str(), ik);
        field_io->read_file(xin, filename);
        AIndex_lex<real_t, AFIELD_d::IMPL> index_alt;
        convert(index_alt, testvectors_d[ik], xin);
      }
    }
  } // field input
  */

  // check of eigenvectors
  /*
  {
    vout.crucial(vl, "check of read data.\n");

    AIndex_block_lex<real_t, AFIELD_d::IMPL>
                       index_block(coarse_lattice, fine_lattice);

    std::vector<complex_t> TDa1(num_blocks);
    std::vector<real_t> anorm2(num_blocks);

    vout.general(vl, "    ik    diff2(ave)  "
                     "  diff2(max)          TDa[0]\n");
    for(int ik = 0; ik < Nk; ++ik){

      fopr_d->mult_dd(x, testvectors_d[ik]);

      for(int block = 0; block < num_blocks; ++block){
        TDa1[block] = TDa[ik + (Nk+Np)*block];
      }
      block_axpy(x, &TDa1[0], testvectors_d[ik], real_t(-1.0), index_block);

      block_norm2(&anorm2[0], x, index_block);

      real_t anorm2_ave = 0.0;
      real_t anorm2_max = 0.0;
      for(int block = 0; block < num_blocks; ++block){
        anorm2_ave += anorm2[block];
        anorm2_max = std::max(anorm2[block], anorm2_max);

        if(block == 1 || block == num_blocks-1){
          vout.general(vl, "    %4d  %4d   (%14.6e, %14.6e)   %14.6e\n",
                       ik, block, real(TDa1[block]), imag(TDa1[block]),
                       anorm2[block]);
        }
      }
      anorm2_ave = Communicator::reduce_sum(anorm2_ave);
      anorm2_max = Communicator::reduce_max(anorm2_max);
      anorm2_ave = anorm2_ave/real_t(num_blocks * CommonParameters::NPE());
      vout.general(vl, "  %4d   %14.6e  %14.6e   (%14.6e, %14.6e)\n",
                   ik, anorm2_ave, anorm2_max, real(TDa1[0]), imag(TDa1[0]));
    }
  }
  */

}

//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
void read_testvectors_eigen(AFopr_dd<AFIELD_d>* fopr_d,
                            std::vector<AFIELD_f>& testvectors,
                            std::vector<int>& coarse_lattice)
{
  Bridge::VerboseLevel vl = vout.set_verbose_level("Detailed");

  vout.crucial(vl, "Setting testvectors with eigensolver\n");

  typedef typename AFIELD_d::real_t real_t;
  typedef typename ComplexTraits<real_t>::complex_t complex_t;
  typedef AFopr<AFIELD_d> AFOPR_d;
  typedef AEigensolver_IRArnoldi_block<AFIELD_d, AFOPR_d> AEIGEN;

  int Ndim = CommonParameters::Ndim();

  int num_vectors = testvectors.size();
  vout.general(vl, "  Number of testvectors = %d\n", num_vectors);

  
  std::vector<int> fine_lattice(4);
  fine_lattice[0] = CommonParameters::Nx();
  fine_lattice[1] = CommonParameters::Ny();
  fine_lattice[2] = CommonParameters::Nz();
  fine_lattice[3] = CommonParameters::Nt();

  int Nin  = fopr_d->field_nin();
  int Nvol = fopr_d->field_nvol();
  int Nex  = fopr_d->field_nex();

  int Nk = testvectors.size();  //this may differ from Nk of eigensolver

  // field input
  {
    vout.crucial(vl, "eigenvectors are output.\n");

    unique_ptr<FieldIO> field_io(new FieldIO_Binary(IO_Format::Trivial));
    int NinF = 2 * CommonParameters::Nc() * CommonParameters::Nd();
    int NexF = Nin / NinF;
    Field xin(NinF, Nvol, NexF);
    AFIELD_d x(Nin, Nvol, Nex);

    std::string file_base = "testvector";

    if (fopr_d->needs_convert()) {

      vout.detailed(vl, "convert in rep. required.\n");
      for(int ik = 0; ik < Nk; ++ik){
        std::string filename = FileUtils::generate_filename(
                                 "%s_%04d.dat", file_base.c_str(), ik);
        field_io->read_file(xin, filename);
        fopr_d->convert(x, xin);
        copy(testvectors[ik], x);
      }

    } else {

      vout.detailed(vl, "convert in rep. not required.\n");
      for(int ik = 0; ik < Nk; ++ik){
        std::string filename = FileUtils::generate_filename(
                                 "%s_%4d.dat", file_base.c_str(), ik);
        field_io->read_file(xin, filename);
        AIndex_lex<real_t, AFIELD_d::IMPL> index_alt;
        convert(index_alt, x, xin);
        copy(testvectors[ik], x);
      }
    }
  } // field input

}

#endif
//============================================================END=====
