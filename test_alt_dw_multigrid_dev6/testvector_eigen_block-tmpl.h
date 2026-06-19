/*!
        @file    testvector_eigen-tmpl.h
        @brief   test for coarse grid operator
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: #$
        @version $LastChangedRevision: 2604 $
*/

#ifndef TESTVECTOR_EIGEN_BLOCK_INCLUDED
#define TESTVECTOR_EIGEN_BLOCK_INCLUDED

#include "lib/Tools/timer.h"
#include "lib/Fopr/fopr.h"
#include "lib/Field/field_F.h"

#include "lib/Eigen/aeigensolver.h"
#include "lib/Eigen/aeigensolver_IRArnoldi.h"
#include "lib_alt/Eigen/aeigensolver_IRArnoldi_block.h"


#define COMPLEX_t typename ComplexTraits<typename AFIELD::real_t>::complex_t
// to be undefed at the end of this file.

// blocked eigensolver implemented as a functional template
//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
void set_testvectors_eigen_block(AFopr_dd<AFIELD_d>* fopr_d,
                                 std::vector<AFIELD_f>& testvectors,
                                 std::vector<int>& coarse_lattice,
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

  int num_blocks = 1;
  for(int mu = 0; mu < 4; ++mu){
    vout.general(vl, "fine_lattice[%2d] = %4d  coarse_lattice[%2d] = %4d\n",
		 mu, fine_lattice[mu], mu, coarse_lattice[mu]);
    num_blocks *= coarse_lattice[mu];
  }
  //vout.general(vl, "  Number of blocks = %d\n", num_vectors);
  vout.crucial(vl, "Number of blocks = %d\n", num_vectors);

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
  
  std::vector<complex_t> TDa(num_blocks * (Nk + Np));

  // setup eigensolver
  AEIGEN* eigen;
  std::string eigensolver_type =
                         params_eigen.get_string("eigensolver_type");
  if(eigensolver_type == "IRArnoldi_block"){
    eigen = new AEigensolver_IRArnoldi_block<AFIELD_d, AFOPR_d>(
                                               fopr_d, params_eigen);
  }else if (eigensolver_type == "IRArnoldi"){
    eigen = new AEigensolver_IRArnoldi<AFIELD_d, AFOPR_d>(
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
  
  // check of eigenvectors

  //std::string file_base = "testvector";
  std::string file_base = filename;
  write_eigenvectors_block(fopr_d, TDa, testvectors_d, Nk,
                           coarse_lattice, file_base);

  // set testvectors in float
  int num_testvectors = testvectors.size();
  if(num_testvectors > Nk){
    vout.crucial(vl, "too large num_testvectors\n");
    exit(EXIT_FAILURE);
  }
  for(int ik = 0; ik < num_testvectors; ++ik){
        copy(testvectors[ik], testvectors_d[ik]);
  }

  vout.decrease_indent();
  vout.crucial(vl, "set_testvectors_eigen finished.\n");
  Communicator::sync();

}

//====================================================================
template<typename AFIELD>
void check_eigenvectors_block(AFopr_dd<AFIELD>* fopr,
                              std::vector<COMPLEX_t>& TDa,
                              std::vector<AFIELD>& testvectors,
                              std::vector<int>& coarse_lattice)
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

  std::vector<int> fine_lattice(Ndim);
  fine_lattice[0] = CommonParameters::Nx();
  fine_lattice[1] = CommonParameters::Ny();
  fine_lattice[2] = CommonParameters::Nz();
  fine_lattice[3] = CommonParameters::Nt();

  int Nblock = 1;
  for(int mu = 0; mu < Ndim; ++ mu){
    Nblock *= coarse_lattice[mu];
  }

  AIndex_block_lex<real_t, AFIELD::IMPL>
                       index_block(coarse_lattice, fine_lattice);

  std::vector<complex_t> TDa1(Nblock);
  std::vector<real_t> anorm2(Nblock);

  AFIELD x(Nin, Nvol, Nex);

  vout.general(vl, "    ik    diff2(ave)  "
                   "  diff2(max)          TDa[0]\n");
  for(int ik = 0; ik < Nk; ++ik){

    fopr->mult_dd(x, testvectors[ik]);

    for(int block = 0; block < Nblock; ++block){
      TDa1[block] = TDa[block + Nblock * ik];
    }
    block_axpy(x, &TDa1[0], testvectors[ik], real_t(-1.0), index_block);

    block_norm2(&anorm2[0], x, index_block);

    real_t anorm2_ave = 0.0;
    real_t anorm2_max = 0.0;
    for(int block = 0; block < Nblock; ++block){
      anorm2_ave += anorm2[block];
      anorm2_max = std::max(anorm2[block], anorm2_max);

      if(block == 1 || block == Nblock-1){
        vout.paranoiac(vl, "    %4d  %4d   (%14.6e, %14.6e)   %14.6e\n",
                       ik, block, real(TDa1[block]), imag(TDa1[block]),
                       anorm2[block]);
      }
    }
    anorm2_ave = Communicator::reduce_sum(anorm2_ave);
    anorm2_max = Communicator::reduce_max(anorm2_max);
    anorm2_ave = anorm2_ave/real_t(Nblock * CommonParameters::NPE());
    vout.general(vl, "  %4d   %14.6e  %14.6e   (%14.6e, %14.6e)\n",
                 ik, anorm2_ave, anorm2_max, real(TDa1[0]), imag(TDa1[0]));
  }

  vout.detailed(vl, "check_eigen_testvectors finished\n");

}


//====================================================================
template<typename AFIELD>
void write_eigenvectors_block(AFopr<AFIELD>* fopr,
                              std::vector<COMPLEX_t>& TDa,
                              std::vector<AFIELD>& testvectors,
                              int Nk,
                              std::vector<int>& coarse_lattice,
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

  // field output
  vout.crucial(vl, "eigenvectors are output.\n");

  unique_ptr<FieldIO> field_io(new FieldIO_Binary(IO_Format::Trivial));
  int NinF = 2 * CommonParameters::Nc() * CommonParameters::Nd();
  int NexF = Nin / NinF;
  Field xout(NinF, Nvol, NexF);

  //  std::string file_base = "testvector";
  std::string file_base = filename;

  if (fopr->needs_convert()) {

    //vout.detailed(vl, "convert in rep. required.\n");
    vout.crucial(vl, "convert in rep. required.\n");
    for(int ik = 0; ik < Nk; ++ik){
      vout.crucial(vl, "norm of testvec = %e\n", testvectors[ik].norm2());
      fopr->reverse(xout, testvectors[ik]);
      vout.crucial(vl, "convert finished: %d\n", ik);

      std::string filename = FileUtils::generate_filename(
                               "%s_%04d.dat", file_base.c_str(), ik);
      //vout.crucial(vl, "file = %s\n", filename.c_str());
      field_io->write_file(xout, filename);
      vout.crucial(vl, "output finished: %d\n", ik);
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

  std::vector<int> fine_lattice(4);
  fine_lattice[0] = CommonParameters::Nx();
  fine_lattice[1] = CommonParameters::Ny();
  fine_lattice[2] = CommonParameters::Nz();
  fine_lattice[3] = CommonParameters::Nt();

  int Nx = CommonParameters::Nx();
  int Ny = CommonParameters::Ny();
  int Nz = CommonParameters::Nz();
  int Nt = CommonParameters::Nt();

  std::vector<int> block_size(Ndim);
  int Nblock = 1;
  for(int mu = 0; mu < Ndim; ++ mu){
    block_size[mu] = fine_lattice[mu]/coarse_lattice[mu];
    Nblock *= coarse_lattice[mu];
  }
  Field eigenvalues(2, Nvol, 1);
  eigenvalues.set(0.0);

  int Nbx = coarse_lattice[0];
  int Nby = coarse_lattice[1];
  int Nbz = coarse_lattice[2];
  int Nbt = coarse_lattice[3];

  for(int ik = 0; ik < Nk; ++ik){

    for(int ibt = 0; ibt < Nbt; ++ibt){
      for(int ibz = 0; ibz < Nbz; ++ibz){
        for(int iby = 0; iby < Nby; ++iby){
          for(int ibx = 0; ibx < Nbx; ++ibx){
            int block = ibx + Nbx * (iby + Nby *(ibz + Nbz * ibt));
            int ix = ibx * block_size[0];
            int iy = iby * block_size[1];
            int iz = ibz * block_size[2];
            int it = ibt * block_size[3];
            int site = ix + Nx * (iy + Ny * (iz + Nz * it));
            double eigr = double(real(TDa[block + Nblock * ik]));
            double eigi = double(imag(TDa[block + Nblock * ik]));
            eigenvalues.set(0, site, 0, eigr);
            eigenvalues.set(1, site, 0, eigi);
	  }
        }
      }
    }

    std::string filename2 = FileUtils::generate_filename(
                              "%s_%04d.eigen", file_base.c_str(), ik);
    field_io->write_file(eigenvalues, filename2);
  }

  vout.crucial(vl, "write_testvectors finished.\n");
  Communicator::sync();

}

//====================================================================
template<typename AFIELD>
void read_eigenvectors_block(AFopr<AFIELD>* fopr,
                             std::vector<COMPLEX_t>& TDa,
                             std::vector<AFIELD>& testvectors,
                             std::vector<int>& coarse_lattice,
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

  if (fopr->needs_convert()) {

    vout.detailed(vl, "convert in rep. required.\n");
    for(int ik = 0; ik < Nk; ++ik){
      vout.detailed(vl, "ik = %d\n", ik);
      Communicator::sync();
      std::string filename = FileUtils::generate_filename(
                               "%s_%04d.dat", file_base.c_str(), ik);
      // vout.detailed(vl, "filename = %s\n", filename.c_str());
      field_io->read_file(xin, filename);
      fopr->convert(testvectors[ik], xin);
    }

  } else {

    vout.detailed(vl, "convert in rep. not required.\n");
    for(int ik = 0; ik < Nk; ++ik){
      std::string filename = FileUtils::generate_filename(
                               "%s_%04d.dat", file_base.c_str(), ik);
      field_io->read_file(xin, filename);
      AIndex_lex<real_t, AFIELD::IMPL> index_alt;
      convert(index_alt, testvectors[ik], xin);
    }
  }

  // read eigenvalues

  std::vector<int> fine_lattice(4);
  fine_lattice[0] = CommonParameters::Nx();
  fine_lattice[1] = CommonParameters::Ny();
  fine_lattice[2] = CommonParameters::Nz();
  fine_lattice[3] = CommonParameters::Nt();

  int Nx = CommonParameters::Nx();
  int Ny = CommonParameters::Ny();
  int Nz = CommonParameters::Nz();
  int Nt = CommonParameters::Nt();

  std::vector<int> block_size(Ndim);
  int Nblock = 1;
  for(int mu = 0; mu < Ndim; ++ mu){
    block_size[mu] = fine_lattice[mu]/coarse_lattice[mu];
    Nblock *= coarse_lattice[mu];
  }
  Field eigenvalues(2, Nvol, 1);
  //  eigenvalues.set(0.0);

  int Nbx = coarse_lattice[0];
  int Nby = coarse_lattice[1];
  int Nbz = coarse_lattice[2];
  int Nbt = coarse_lattice[3];

  for(int ik = 0; ik < Nk; ++ik){

    std::string filename2 = FileUtils::generate_filename(
                              "%s_%04d.eigen", file_base.c_str(), ik);
    field_io->read_file(eigenvalues, filename2);

    for(int ibt = 0; ibt < Nbt; ++ibt){
      for(int ibz = 0; ibz < Nbz; ++ibz){
        for(int iby = 0; iby < Nby; ++iby){
          for(int ibx = 0; ibx < Nbx; ++ibx){
            int block = ibx + Nbx * (iby + Nby * (ibz + Nbz * ibt));
            int ix = ibx * block_size[0];
            int iy = iby * block_size[1];
            int iz = ibz * block_size[2];
            int it = ibt * block_size[3];
            int site = ix + Nx * (iy + Ny * (iz + Nz * it));

            double eigr = eigenvalues.cmp(0, site, 0);
            double eigi = eigenvalues.cmp(1, site, 0);
            TDa[block + Nblock * ik] = cmplx(real_t(eigr), real_t(eigi));
	  }
        }
      }
    }

  }

  vout.detailed(vl, "read_testvectors finished.\n");

}

#undef COMPLEX_t

#endif
//============================================================END=====
