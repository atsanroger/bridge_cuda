/*!
        @file    aeigensolver_IRArnoldi_block-tmpl.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2024-06-01 20:47:37 #$
        @version $LastChangedRevision: 2604 $
*/

#include "lib_alt/Eigen/aeigensolver_IRArnoldi_block.h"

#ifdef __PGI
#pragma global opt=1
// This is becaouse for PGI compiler on some environment,
// compilation with optimization level higher than O2 fails.
//                                  [12 Feb 2021 H.Matsufuru]
#endif

using std::string;

#include "lib/ResourceManager/threadManager.h"

template<typename FIELD, typename FOPR>
const string AEigensolver_IRArnoldi_block<FIELD, FOPR>::
class_name = "AEigensolver_IRArnoldi_block";

//====================================================================
namespace{
  template<typename real_t>
  real_t pythag(real_t a, real_t b){
    real_t absa, absb, v;
    absa = abs(a);
    absb = abs(b);
    if(absa > absb){
      real_t rat = absb/absa;
      v = absa * sqrt(1.0 + rat*rat);
    }else{
      if(absb == 0.0){
        v = 0.0;
      }else{
        real_t rat = absa/absb;
        v = absb * sqrt(1.0 + rat*rat);
      }
    }
    return v;
  }

}

//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRArnoldi_block<FIELD, FOPR>::init(
                                             const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  } else {
    m_vl = CommonParameters::Vlevel();
  }

  vout.general(m_vl, "%s: construction\n", class_name.c_str());
  vout.increase_indent();

  m_Ndim = CommonParameters::Ndim();

  m_fine_lattice.resize(m_Ndim);
  m_coarse_lattice.resize(m_Ndim);

  m_ecrit = 1.0e-14;  //!< criteria of too small number

  set_parameters(params);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}

//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRArnoldi_block<FIELD, FOPR>::tidyup()
{
  if(m_index_block) delete m_index_block;
  if (m_sorter) delete m_sorter;
}


//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRArnoldi_block<FIELD, FOPR>::set_parameters(
  const Parameters& params)
{
  string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  string sort_type;
  int    Nk, Np;
  int    Niter_eigen;
  double Enorm_eigen, Vthreshold;
  std::vector<int> fine_lattice(m_Ndim);
  std::vector<int> coarse_lattice(m_Ndim);

  int err = 0;
  err += params.fetch_string("eigensolver_mode", sort_type);
  err += params.fetch_int("number_of_wanted_eigenvectors", Nk);
  err += params.fetch_int("number_of_working_eigenvectors", Np);
  err += params.fetch_int("maximum_number_of_iteration", Niter_eigen);
  err += params.fetch_double("convergence_criterion_squared", Enorm_eigen);
  err += params.fetch_double("threshold_value", Vthreshold);
  err += params.fetch_int_vector("fine_lattice", fine_lattice);
  err += params.fetch_int_vector("coarse_lattice", coarse_lattice);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(sort_type, Nk, Np, Niter_eigen, Enorm_eigen,
                 Vthreshold, fine_lattice, coarse_lattice);
}


//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRArnoldi_block<FIELD, FOPR>::get_parameters(
                                          Parameters& params) const
{
  params.set_string("eigensolver_mode", m_sort_type);
  params.set_int("number_of_wanted_eigenvectors", m_Nk);
  params.set_int("number_of_working_eigenvectors", m_Np);
  params.set_int("maximum_number_of_iteration", m_Niter_eigen);
  params.set_double("convergence_criterion_squared", m_Enorm_eigen);
  params.set_double("threshold_value", m_Vthreshold);
  params.set_int_vector("fine_lattice", m_fine_lattice);
  params.set_int_vector("coarse_lattice", m_coarse_lattice);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRArnoldi_block<FIELD, FOPR>::set_parameters(
                                     const string& sort_type,
                                     int Nk,
                                     int Np,
                                     int Niter_eigen,
                                     double Enorm_eigen,
                                     double Vthreshold,
                                     std::vector<int> fine_lattice,
                                     std::vector<int> coarse_lattice)
{
  if (m_sorter) delete m_sorter;

  m_sort_type = sort_type;
  vout.general("  sort_type = %s\n", m_sort_type.c_str());
  m_sorter    = new Sorter<complex_t>(sort_type);

  set_parameters(Nk, Np, Niter_eigen, Enorm_eigen, Vthreshold,
                 fine_lattice, coarse_lattice);
}


//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRArnoldi_block<FIELD, FOPR>::set_parameters(
                                     int Nk,
                                     int Np,
                                     int Niter_eigen,
                                     double Enorm_eigen,
                                     double Vthreshold,
                                     std::vector<int> fine_lattice,
                                     std::vector<int> coarse_lattice)
{
  //- range check
  int err = 0;
  err += ParameterCheck::non_negative(Nk);
  err += ParameterCheck::non_negative(Np);
  err += ParameterCheck::non_negative(Niter_eigen);
  err += ParameterCheck::square_non_zero(Enorm_eigen);
  // NB. Vthreshold == 0 is allowed.

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  //- store values
  m_Nk          = Nk;
  m_Np          = Np;
  m_Nm          = m_Nk + m_Np;
  m_Niter_eigen = Niter_eigen;
  m_Enorm_eigen = real_t(Enorm_eigen);
  m_Vthreshold  = real_t(Vthreshold);
  m_fine_lattice = fine_lattice;
  m_coarse_lattice = coarse_lattice;

  if(m_index_block != 0) delete m_index_block;
  m_index_block = new AIndex_block(m_coarse_lattice, m_fine_lattice);
  m_num_blocks  = m_index_block->coarse_nvol();

  //- print input parameters
  vout.general(m_vl, "%s: set paramters\n", class_name.c_str());
  vout.general(m_vl, "  Nk          = %d\n", m_Nk);
  vout.general(m_vl, "  Np          = %d\n", m_Np);
  vout.general(m_vl, "  Niter_eigen = %d\n", m_Niter_eigen);
  vout.general(m_vl, "  Enorm_eigen = %16.8e\n", m_Enorm_eigen);
  vout.general(m_vl, "  Vthreshold  = %16.8e\n", m_Vthreshold);
  for(int mu = 0; mu < m_Ndim; ++mu){
    vout.general(m_vl, "  fine_lattice[%1d] = %4d"
		 "   coarse_lattice[%1d] = %d\n",
                 mu, m_fine_lattice[mu], mu, m_coarse_lattice[mu]);
  }
  vout.general(m_vl, "  num_blocks  = %d\n", m_num_blocks);

  m_TDa2.resize(m_num_blocks);
  m_HtB.resize(m_num_blocks);
  m_Ht2B.resize(m_num_blocks);
  m_YtB.resize(m_num_blocks);
  m_Iconv.resize(m_num_blocks);
  m_QtB.resize(m_num_blocks);
  for(int block = 0; block < m_num_blocks; ++block){
    m_TDa2[block].resize(m_Nm);
    m_HtB[block].resize((m_Nm + 1) * m_Nm);  // (Nm+1) * Nm matrix
    m_Ht2B[block].resize((m_Nm + 1) * m_Nm);
    m_YtB[block].resize(m_Nm * m_Nm);
    m_Iconv[block].resize(m_Nm);
    m_QtB[block].resize(m_Nm * m_Nm);
  }

  m_Kdis.resize(m_num_blocks);
  m_Kthreshold.resize(m_num_blocks);
  m_Kalmost.resize(m_num_blocks);

  m_reorder.resize(m_num_blocks);

  int Nin  = m_fopr->field_nin();
  int Nvol = m_fopr->field_nvol();
  int Nex  = m_fopr->field_nex();

  m_B.resize(m_Nm);
  for (int k = 0; k < m_Nm; ++k) {
    m_B[k].reset(Nin, Nvol, Nex);
  }

  m_f.reset(Nin, Nvol, Nex);
  m_v.reset(Nin, Nvol, Nex);

  m_bvec1.resize(m_num_blocks);

  m_rvec1.resize(m_num_blocks);
  m_rvec2.resize(m_num_blocks);

}

//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRArnoldi_block<FIELD, FOPR>::solve(
                                       std::vector<complex_t>& TDk,
                                       std::vector<FIELD>& vk,
                                       int& Nsbt,
                                       int& Nconv,
                                       const FIELD& b)
{
  ThreadManager::assert_single_thread(class_name);

  int blockF = m_num_blocks-1;  // monitored block
  //int blockF = 0;  // monitored block

  int    Nk          = m_Nk;
  int    Np          = m_Np;
  int    Niter_eigen = m_Niter_eigen;
  real_t Enorm_eigen = m_Enorm_eigen;
  real_t Vthreshold  = m_Vthreshold;

  real_t Enorm_almost = 1.e-12; // this is an experimental value

  int Nm = Nk + Np;

  if (Nm * m_num_blocks > TDk.size()) {
    vout.crucial(m_vl, "Error at %s: std::vector TDa is too small.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  } else if (Nk + Np > vk.size()) {
    vout.crucial(m_vl, "Error at %s: std::vector vk is too small.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  std::vector<complex_t> TDa(Nm * m_num_blocks);

  vout.general(m_vl, "Implicitly Restarted Arnoldi_block algorithm start\n");
  vout.general(m_vl, "  Nk = %d  Np = %d\n", Nk, Np);
  vout.general(m_vl, "  Nm = %d\n", Nm);
  vout.general(m_vl, "  size of TDa  = %d\n", TDk.size());
  vout.general(m_vl, "  size of vk   = %d\n", vk.size());
  vout.general(m_vl, "  Enorm_eigen  = %e\n", Enorm_eigen);
  vout.general(m_vl, "  Enorm_almost = %e\n", Enorm_almost);

  Nconv = -1;
  Nsbt  = 0;

  int k1 = 0;
  int k2 = Nk;

  for(int block = 0; block < m_num_blocks; ++block){
    m_Kdis[block]       = 0;
    m_Kthreshold[block] = 0;
    m_Kalmost[block]    = 0;
  }
  int Kdis       = 0;
  int Kthreshold = 0;
  int Kalmost    = 0;

  int nconv;

  std::vector<complex_t> beta(m_num_blocks);
  for(int block = 0; block < m_num_blocks; ++block){
    beta[block] = cmplx(real_t(0.0), real_t(0.0));
  }

#pragma omp parallel
  {
    int ith, nth, iblock, nblock;
    set_threadtask(ith, nth, iblock, nblock, m_num_blocks);

    //- set initial vector
    real_t bnorm2 = b.norm2();
    if (bnorm2 > 1.e-12) { // b is regarded as a nonzero vector
      copy(vk[0], b);
    } else {
      vk[0].set(1.0);  // start with a uniform vector
    }
#pragma omp barrier

    block_norm2(&m_rvec1[0], vk[0], *m_index_block);

#pragma omp barrier

    if(ith == 0){
      for(int block = 0; block < m_num_blocks; ++block){
        m_rvec1[block] = 1.0/sqrt(m_rvec1[block]);
      }
    }
#pragma omp barrier

    block_scal(vk[0], &m_rvec1[0], *m_index_block);

#pragma omp barrier

    //if(ith == 0){
      for (int block = iblock; block < nblock; ++block) {
        for (int j = 0; j < Nm; ++j) {
          for (int k = 0; k < Nm + 1; ++k) {
            m_HtB[block][ index(k, j)] = cmplx(real_t(0.0), real_t(0.0));
            m_Ht2B[block][index(k, j)] = cmplx(real_t(0.0), real_t(0.0));
          }
        }
      }
      //}
#pragma omp barrier

    //- initial Nk steps
    for (int k = 0; k < k2; ++k) {
      step(Nm, k, vk, m_f);
    }

#pragma omp barrier

    //- restarting loop begins
    for (int iter = 0; iter <= Niter_eigen; ++iter) {
      vout.detailed(m_vl, "\n iter = %d\n", iter);

      for(int block = iblock; block < nblock; ++block){

        // here size Nk eigenproblem is solved.

        beta[block] = m_HtB[block][index(k2, k2 - 1)];
        m_HtB[block][index(k2, k2 - 1)] = cmplx(real_t(0.0), real_t(0.0));

        m_Ht2B[block] = m_HtB[block];
        setUnit_Qt(Nm, m_QtB[block]);

        tqri(m_Ht2B[block], k1, Nk, Nm, m_QtB[block], nconv);

        std::vector<complex_t> TDa1(m_Nm);
        for (int j = k1; j < Nk; ++j) {
          TDa1[j] = m_Ht2B[block][index(j, j)];
        }

        // here small eigenvalue problem should be checked for confirmation
        std::vector<complex_t> Xt(Nm * Nm);

        eigenvector_Ht(Xt, m_Ht2B[block], Nk, Nm);
        mult_Qt(m_YtB[block], m_QtB[block], Xt, Nk, Nm);
        check_eigen_Ht(m_HtB[block], TDa1, m_YtB[block], Nk, Nm);
        // now m_Yt is eigenvectors of Ht

        for (int j = k1; j < Nk; ++j) {
          TDa[j + m_Nm * block] = TDa1[j];
        }
      }
      if(ith == 0){
        Kdis       = 0;
        Kthreshold = 0;
        Kalmost    = 0;
      }
#pragma omp barrier

      // convergenece test
      vout.general(m_vl, "                        TDa            "
                         "        abs           diff2\n");

      for (int j = 0; j < k1; ++j) {

        m_fopr->mult_dd(m_v, m_B[j]);

#pragma omp barrier

        if(ith == 0){
          for(int block = 0; block < m_num_blocks; ++block){
            m_bvec1[block] = TDa[j + m_Nm * block];
          }
	}
#pragma omp barrier

	block_axpy(m_v, &m_bvec1[0], m_B[j], real_t(-1.0),
                   *m_index_block);

        block_norm2(&m_rvec1[0], m_B[j], *m_index_block);
        block_norm2(&m_rvec2[0], m_v,    *m_index_block);

#pragma omp barrier

        if(ith == 0){
          real_t vv = 0.0;
          for(int block = 0; block < m_num_blocks; ++block){
            m_rvec2[block] = m_rvec2[block]/m_rvec1[block];
            vv += m_rvec2[block];
	  }

          complex_t TDaF = TDa[j + m_Nm * blockF];
          vout.general(m_vl, " %4d  %4d  (%12.8f, %12.8f)  %12.8f  %12.4e\n",
                       j, blockF, real(TDaF), imag(TDaF), abs(TDaF), vv);

        }
#pragma omp barrier

      } // j-loop

      for (int j = k1; j < Nk; ++j) {

        m_B[j].set(0.0);

#pragma omp barrier

        for (int k = 0; k < Nk; ++k) {

          if(ith == 0){
            for(int block = 0; block < m_num_blocks; ++block){
              m_bvec1[block] = m_YtB[block][index(k, j)];
            }
          }
#pragma omp barrier

          block_axpy(m_B[j], &m_bvec1[0], vk[k], real_t(1.0),
                     *m_index_block);
#pragma omp barrier
        }

        m_fopr->mult_dd(m_v, m_B[j]);

        if(ith == 0){
          for(int block = 0; block < m_num_blocks; ++block){
            m_bvec1[block] = TDa[j + m_Nm * block];
          }
        }
#pragma omp barrier

	block_axpy(m_v, &m_bvec1[0], m_B[j], real_t(-1.0),
                   *m_index_block);

        block_norm2(&m_rvec1[0], m_B[j], *m_index_block);
        block_norm2(&m_rvec2[0], m_v,    *m_index_block);

#pragma omp barrier

        if(ith == 0){
          real_t vv = 0.0;
          for(int block = 0; block < m_num_blocks; ++block){
            m_rvec2[block] = m_rvec2[block]/m_rvec1[block];
            vv += m_rvec2[block];
	  }

          complex_t TDaF = TDa[j + m_Nm * blockF];
          vout.general(m_vl, " %4d  %4d  (%12.8f, %12.8f)  %12.8f  %12.4e\n",
                     j, blockF, real(TDaF), imag(TDaF), abs(TDaF), vv);
        }
#pragma omp barrier

        if(ith == 0){
          for(int block = 0; block < m_num_blocks; ++block){

            m_TDa2[block][j] = TDa[j + m_Nm * block];

            if (m_rvec2[block] < Enorm_eigen) {
              m_Iconv[block][Kdis] = j;
              ++m_Kdis[block];
              if (!m_sorter->comp(m_TDa2[block][j], Vthreshold)) {
                ++m_Kthreshold[block];
              }
            } else if ((m_rvec2[block] < Enorm_almost) &&
                       m_sorter->comp(m_TDa2[block][j], Vthreshold)) {
              ++m_Kalmost[block];
            }
          }
        }
#pragma omp barrier

      }   // j-loop end

      if(ith == 0){
        Kdis       = 0;
        Kthreshold = 0;
        Kalmost    = 0;
        for(int block = 0; block < m_num_blocks; ++block){
          Kdis       += m_Kdis[block];
          if(m_Kthreshold[block] > 0) ++Kthreshold;
          if(m_Kalmost[block] > 0)    ++Kalmost;
        }
      }

      vout.detailed(m_vl, " #modes already converged: %4d\n", k1);
      vout.detailed(m_vl, " #modes newly converged:   %4d\n", Kdis);
      vout.detailed(m_vl, " #modes almost converged:  %4d\n", Kalmost);

      if ((Kthreshold == m_num_blocks) && (Kalmost == m_num_blocks)) {

        if(ith == 0){
          Nsbt  = 0;
          for(int block = 0; block < m_num_blocks; ++block){
            // (there is a converged eigenvalue larger than Vthreshold.)

            Nsbt  += m_Kdis[block] - m_Kthreshold[block];
            Nconv = iter;
          }
        }
#pragma omp barrier

      //break;
	goto break_iteration;
        //iter = Niter_eigen;
      }else if(iter == Niter_eigen){
        //break;
	goto break_iteration;
      }

      if(ith == 0){
        for(int block = 0; block < m_num_blocks; ++block){
          m_HtB[block][index(k2, k2 - 1)] = beta[block];
        }
      }
#pragma omp barrier

      // if not converged, Krylov subspace is extended

      for (int k = k2 - 1; k < Nm; ++k) {
        step(Nm, k, vk, m_f);
      }

      for(int block = iblock; block < nblock; ++block){

        for (int j = 0; j < Nm; ++j) {
          for (int k = 0; k < Nm; ++k) {
            m_Ht2B[block][index(k, j)] = m_HtB[block][index(k, j)];
          }
        }

        //- getting eigenvalues
        setUnit_Qt(Nm, m_QtB[block]);

        // vout.detailed(m_vl, "QR transformation\n");

        tqri(m_Ht2B[block], k1, Nm, Nm, m_QtB[block], nconv);

        for (int j = 0; j < Nm; ++j) {
          m_TDa2[block][j] = m_Ht2B[block][index(j, j)];
        }

        // here small eigenvalue problem should be checked for confirmation
        std::vector<complex_t> Xt(Nm * Nm);

        eigenvector_Ht(Xt, m_Ht2B[block], Nm, Nm);
        mult_Qt(m_YtB[block], m_QtB[block], Xt, Nm, Nm);
        check_eigen_Ht(m_HtB[block], m_TDa2[block], m_YtB[block], Nm, Nm);
        // now m_Yt is eigenvectors of Ht

        //- sorting
        m_sorter->sort(m_TDa2[block], Nm);
        // CAUTION: if already converged eigenvalues become the
        //   i >= k2 range, the following algorithm may fail.
        //   Some care is necessary.

        //m_Ht2 = m_Ht;

        setUnit_Qt(Nm, m_QtB[block]);

        //- implicitly shifted QR transformations
        for (int ip = k2; ip < Nm; ++ip) {
          complex_t Dsh  = m_TDa2[block][ip];
          int       kmin = k1 + 1;
          int       kmax = Nm;
          // vout.general(m_vl, "ip = %d  Dsh = (%f, %f)\n",
          //              ip, real(Dsh), imag(Dsh));
          qrtrf(m_HtB[block], Nm, Nm, m_QtB[block], Dsh, kmin, kmax);
          // check_Qt(Nm, Nm, m_Qt, m_Ht, m_Ht2);
        }
      }
#pragma omp barrier

      for (int j = k1; j <= k2; ++j) {
        m_B[j].set(0.0);
#pragma omp barrier

        //for (int k = 0; k < Nm; ++k) {
        for (int k = k1; k < Nm; ++k) {

	  if(ith == 0){
            for(int block = 0; block < m_num_blocks; ++block){
              m_bvec1[block] = m_QtB[block][index(k, j)];
            }
	  }
#pragma omp barrier

          block_axpy(m_B[j], &m_bvec1[0], vk[k], real_t(1.0),
                     *m_index_block);
#pragma omp barrier
        }
      }

      for (int j = k1; j <= k2; ++j) {
        copy(vk[j], m_B[j]);
      }
#pragma omp barrier

      schmidt(k2, vk);

      //- Compressed vector f and beta(k2)
      copy(m_f, m_B[k2]);
#pragma omp barrier

      if(ith == 0){
        for(int block = 0; block < m_num_blocks; ++block){
          m_bvec1[block] = m_HtB[block][index(k2, k2 - 1)];
        }
      }
#pragma omp barrier

      block_scal(m_f, &m_bvec1[0], *m_index_block);
#pragma omp barrier

      block_norm2(&m_rvec1[0], m_f, *m_index_block);
#pragma omp barrier

      if(ith == 0){
        real_t beta_k = 0.0;
        for(int block = 0; block < m_num_blocks; ++block){
          beta_k += m_rvec1[block] * m_rvec1[block];
        }
        beta_k = sqrt(beta_k);
        vout.detailed(m_vl, " beta(k) = %20.14f\n", beta_k);
      }
#pragma omp barrier

      // reconstructing Hessenberg matrix
      reconst_Ht(m_Ht2B, beta, vk, Nk);

      if(ith == 0){
        real_t diff = 0.0;
        for(int block = 0; block < m_num_blocks; ++block){
	  for (int i = 0; i < Nk; ++i) {
            for (int j = 0; j < Nk; ++j) {
              real_t diff1 = abs(m_HtB[block][index(i, j)]
				 - m_Ht2B[block][index(i, j)]);
              diff += diff1 * diff1;
            }
          }
        }
        diff = diff / real_t(Nk * Nk * m_num_blocks);
        vout.general(m_vl, " difference of Hessenberg = %e\n", diff);

        m_HtB = m_Ht2B;
      }
#pragma omp barrier

    } // end of iter loop

    break_iteration:
      vout.general(m_vl, "iteration loop was broken\n");

    // not sorted yet
    for (int i = 0; i < m_Nk; ++i) {
      copy(vk[i], m_B[i]);
    }
#pragma omp barrier

    // sort
    sort_eigenvectors(TDa, vk);

    // check of eigenvectors
    vout.general(m_vl, "    ik      diff2(ave)      diff2(max)  "
                       "              TDa[%d]\n", blockF);

    for(int ik = 0; ik < Nk; ++ik){

      m_fopr->mult_dd(m_v, vk[ik]);

      if(ith == 0){
        for(int block = 0; block < m_num_blocks; ++block){
          m_bvec1[block] = TDa[ik + m_Nm * block];
        }
      }
#pragma omp barrier

      block_axpy(m_v, &m_bvec1[0], vk[ik], real_t(-1.0), *m_index_block);

      block_norm2(&m_rvec1[0], m_v, *m_index_block);

#pragma omp barrier

      if(ith == 0){
        real_t anorm2_ave = 0.0;
        real_t anorm2_max = 0.0;
        for(int block = 0; block < m_num_blocks; ++block){
          anorm2_ave += m_rvec1[block];
          anorm2_max = std::max(m_rvec1[block], anorm2_max);
        }
        anorm2_ave = Communicator::reduce_sum(anorm2_ave);
        anorm2_max = Communicator::reduce_max(anorm2_max);
        int NPE = CommonParameters::NPE();
        anorm2_ave = anorm2_ave/real_t(m_num_blocks * NPE);
        vout.general(m_vl, "  %4d   %14.6e  %14.6e   (%14.6e, %14.6e)\n",
                     ik, anorm2_ave, anorm2_max,
                   real(m_bvec1[blockF]), imag(m_bvec1[blockF]));
      }
    } // ik loop

  } // omp parallel

  for(int ik = 0; ik < Nm; ++ik){
    for(int block = 0; block < m_num_blocks; ++block){
      TDk[block + m_num_blocks * ik] = TDa[ik + m_Nm * block];
    }
  }

  if (Nconv == -1) {
    vout.crucial(m_vl, "%s: NOT converged.\n", class_name.c_str());
    //exit(EXIT_FAILURE);
  } else {
    vout.general(m_vl, "\n Converged:\n");
    vout.general(m_vl, "  Nconv   = %d\n", Nconv);
    vout.general(m_vl, "  Kdis    = %d\n", Kdis);
    vout.general(m_vl, "  Nsbt    = %d\n", Nsbt);
  }
  vout.general(m_vl, "%s: solve finished.\n\n", class_name.c_str());
}

//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRArnoldi_block<FIELD, FOPR>::sort_eigenvectors(
                                     std::vector<complex_t>& TDa,
                                     std::vector<FIELD>& vk)
{
  int ith = ThreadManager::get_thread_id();

  if(ith == 0){
    for(int block = 0; block < m_num_blocks; ++block){

      std::vector<int>  reorder1(m_Nk);
      for(int ik = 0; ik < m_Nk; ++ik){
        reorder1[ik] = ik;
      }

      for(int ik = 0; ik < m_Nk; ++ik){
        int ikp = reorder1[ik];
        complex_t TD_tmp  = TDa[ikp + m_Nm * block];
        int       ik_tmp  = ikp;
        int       ik2_tmp = ik;

        for(int ik2 = ik+1; ik2 < m_Nk; ++ik2){
          int ik2p = reorder1[ik2];
          if (m_sorter->comp(TDa[ik2p + m_Nm * block], TD_tmp)) {
            TD_tmp  = TDa[ik2p + m_Nm * block];
            ik_tmp  = ik2p;
            ik2_tmp = ik2;
          }
        }

        if(ik_tmp != reorder1[ik]){
          for(int ik2 = ik2_tmp; ik2 > ik; --ik2){
            reorder1[ik2] = reorder1[ik2-1];
          }
          reorder1[ik] = ik_tmp;
        }
      }
      m_reorder[block] = reorder1;
    }

    for(int ik = 0; ik < m_Nk; ++ik){
      for(int block = 0; block < m_num_blocks; ++block){
        m_TDa2[block][ik] = TDa[m_reorder[block][ik] + m_Nm * block];
      }
    }

    for(int ik = 0; ik < m_Nk; ++ik){
      for(int block = 0; block < m_num_blocks; ++block){
        TDa[ik + m_Nm * block] = m_TDa2[block][ik];
      }
    }

  } // if(ith == 0)
#pragma omp barrier

 for(int ik = 0; ik < m_Nk; ++ik){

    vk[ik].set(real_t(0.0));
#pragma omp barrier

    for(int ik2 = 0; ik2 < m_Nk; ++ik2){

      if(ith == 0){
        for(int block = 0; block < m_num_blocks; ++block){
          m_rvec1[block] = real_t(0.0);
          if(m_reorder[block][ik] == ik2) m_rvec1[block] = 1.0;
        }
      }
#pragma omp barrier

      int add_any = 0;
      for(int block = 0; block < m_num_blocks; ++block){
        if(m_reorder[block][ik] == ik2) ++add_any;
      }

      if(add_any > 0){
        block_axpy(vk[ik], &m_rvec1[0], m_B[ik2], real_t(1.0),
                   *m_index_block);
      }
#pragma omp barrier
    }

    block_norm2(&m_rvec1[0], vk[ik], *m_index_block);

#pragma omp barrier

    if(ith == 0){
      for(int block = 0; block < m_num_blocks; ++block){
        m_rvec1[block] = real_t(1.0)/sqrt(m_rvec1[block]);
      }
    }
#pragma omp barrier

    block_scal(vk[ik], &m_rvec1[0], *m_index_block);

#pragma omp barrier
 } // ik loop

}
//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRArnoldi_block<FIELD, FOPR>::step(
                                             int Nm, int k,
                                             std::vector<FIELD>& vk,
                                             FIELD& w)
{
  if (k >= Nm) {
    vout.crucial(m_vl, "Error at %s: k is larger than Nm.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  int ith = ThreadManager::get_thread_id();

  m_fopr->mult_dd(w, vk[k]);

  for (int j = 0; j <= k; ++j) {
    block_dotc(&m_bvec1[0], vk[j], w, *m_index_block);
    block_axpy(w, &m_bvec1[0], vk[j], real_t(-1.0),
               *m_index_block);

#pragma omp barrier

    if(ith == 0){
      for(int block = 0; block < m_num_blocks; ++block){
        m_HtB[block][index(j, k)] = m_bvec1[block];
      }
    }
#pragma omp barrier
  }

  block_norm2(&m_rvec1[0], w, *m_index_block);  // m_rvec1 <= beta

#pragma omp barrier

  if(ith == 0){
    for(int block = 0; block < m_num_blocks; ++block){
      real_t beta = m_rvec1[block];
      beta = sqrt(beta);
      m_HtB[block][index(k + 1, k)] = cmplx(beta, real_t(0.0));
      m_rvec1[block] = 1.0 / beta;
    }
  }
#pragma omp barrier

  if (k < Nm - 1) {

    copy(vk[k + 1], w);
#pragma omp barrier

    block_scal(vk[k + 1], &m_rvec1[0], *m_index_block);
#pragma omp barrier

  }

}


//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRArnoldi_block<FIELD, FOPR>::schmidt(
                                            int Nk,
                                            std::vector<FIELD>& vk)
{
  // block orthogonalization

  int ith = ThreadManager::get_thread_id();

#pragma omp barrier

  for (int j = 0; j < Nk; ++j) {
    for (int i = 0; i < j; ++i) {
      block_dotc(&m_bvec1[0], vk[i], vk[j], *m_index_block);
      block_axpy(vk[j], &m_bvec1[0], vk[i], real_t(-1.0),
                 *m_index_block);
    }
    block_norm2(&m_rvec1[0], vk[j], *m_index_block);

#pragma omp barrier

    if(ith == 0){
      for(int block = 0; block < m_num_blocks; ++block){
        m_rvec1[block] = 1.0/sqrt(m_rvec1[block]);
      }
    }
#pragma omp barrier

    block_scal(vk[j], &m_rvec1[0], *m_index_block);
  }

#pragma omp barrier

}

//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRArnoldi_block<FIELD, FOPR>::setUnit_Qt(
  int Nm, std::vector<complex_t>& Qt)
{
  for (int i = 0; i < Qt.size(); ++i) {
    Qt[i] = cmplx(0.0, 0.0);
  }

  for (int k = 0; k < Nm; ++k) {
    Qt[index(k, k)] = cmplx(1.0, 0.0);
  }
}


//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRArnoldi_block<FIELD, FOPR>::tqri(
              std::vector<complex_t>& Ht, int k1,
              int Nk, int Nm, std::vector<complex_t>& Qt, int& nconv)
{
  int Niter = 100 * Nm;
  int kmin  = k1 + 1;
  int kmax  = Nk;
  // (these parameters should be tuned)
  // vout.general("tqri: k1 = %d\n", k1);

  nconv = -1;

  std::vector<complex_t> Hr(Nm * Nm);
  for (int i = 0; i < Nm; ++i) {
    for (int j = 0; j < Nm; ++j) {
      Hr[index(i, j)] = Ht[(index(i, j))];
    }
  }

  for (int iter = 0; iter < Niter; ++iter) {
    complex_t Dsh;  // shift
    if (kmax > kmin) {
      shift_wilkinson(Dsh, Ht[index(kmax - 2, kmax - 2)],
                      Ht[index(kmax - 2, kmax - 1)],
                      Ht[index(kmax - 1, kmax - 2)],
                      Ht[index(kmax - 1, kmax - 1)]);
    } else {
      Dsh = Ht[index(kmax - 1, kmax - 1)];
    }


    //- transformation
    qrtrf(Ht, Nk, Nm, Qt, Dsh, kmin, kmax);

    check_Qt(Nk, Nm, Qt, Ht, Hr);

    //- Convergence criterion (redef of kmin and kmax)
    for (int j = kmax - 1; j >= kmin; --j) {
      real_t dds = abs(Ht[index(j - 1, j - 1)]) + abs(Ht[index(j, j)]);
      if (abs(Ht[index(j, j - 1)]) + dds > dds) {
        kmax = j + 1;
        for (int i = 0; i < kmax - 1; ++i) {
          real_t dds = abs(Ht[index(j, j)]) + abs(Ht[index(j + 1, j + 1)]);
          if (abs(Ht[index(i + 1, i)]) + dds > dds) {
            kmin = i + 1;
            break;
          }
        }
        break;
      }

      if (j == kmin) {
        nconv = iter;
        vout.paranoiac(m_vl, "  tqri: converged at iter = %d\n", nconv);
        reunit_Qt(Qt, Nk); // reunitarization of Qt
        return;
      }
    } // end of j loop
  }   // end of iter loop

  if (nconv == -1) {
    vout.crucial(m_vl, "Error at %s: QR method NOT converged, Niter = %d.\n",
                 class_name.c_str(), Niter);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRArnoldi_block<FIELD, FOPR>::qrtrf(
                                          std::vector<complex_t>& Ht,
                                          int Nk, int Nm,
                                          std::vector<complex_t>& Qt,
                                          complex_t Dsh,
                                          int kmin, int kmax)
{
  int kmin1 = kmin;
  int kmax1 = kmax;

  for(int k2 = 0; k2 < Nk; ++k2){

    kmax1 = kmax;
    for(int k = kmin1-1; k < kmax; ++k){
      if(abs(Ht[index(k+1,k)]) == 0.0){
        kmax1 = k+1;
        break;
      }
    }
    if(kmax1 == kmin1){
      ++kmin1;
      continue;
    }

    if(kmax1 != kmax || kmin1 != kmin){
      vout.paranoiac(m_vl, "qrtrf: kmin1 = %d  kmax1 = %d\n", kmin1, kmax1);
      int k = kmax1 - 1;
      vout.paranoiac(m_vl, "(1) Ht[%d, %d] = (%e, %e)\n", k+1, k,
                     real(Ht[index(k+1, k)]), imag(Ht[index(k+1, k)]));
    }

    qrtrf1(Ht, Nk, Nm, Qt, Dsh, kmin1, kmax1);

    kmin1 = kmax1 + 1;
    if(kmin1 >= kmax) break;
  }

}

//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRArnoldi_block<FIELD, FOPR>::qrtrf1(
                                          std::vector<complex_t>& Ht,
                                          int Nk, int Nm,
                                          std::vector<complex_t>& Qt,
                                          complex_t Dsh,
                                          int kmin, int kmax)
{
  int k = kmin - 1;

  real_t f1 = abs(Ht[index(k, k)] - Dsh);
  real_t f2 = abs(Ht[index(k + 1, k)]);

  //stabilized version
  complex_t c, s, beta;
  if(f2 == 0.0){
    vout.paranoiac(m_vl, "f2 = 0.0\n");
    c = cmplx(1.0, 0.0);
    s = cmplx(0.0, 0.0);
    beta = cmplx(1.0, 0.0);
  }else if(f1 == 0.0){
    vout.paranoiac(m_vl, "f1 = 0.0\n");
    c    = cmplx(0.0, 0.0);
    s    = cmplx(1.0, 0.0);
    beta = cmplx(1.0, 0.0);
  }else{
    real_t Fden = 1.0/pythag(f1, f2);
    beta =  (Ht[index(k,k)]-Dsh)/f1;
    complex_t tmp1 = cmplx(real_t(1.0)/f1, real_t(0.0));
    beta =  (Ht[index(k,k)] - Dsh) * tmp1;
    c = cmplx( f1 * Fden, real_t(0.0));
    s = cmplx(-real(Ht[index(k + 1, k)]) * Fden, real_t(0.0));
  }

  for (int j = 0; j < Nk; ++j) {
    complex_t tmp1 = Ht[index(k, j)];
    complex_t tmp2 = Ht[index(k + 1, j)];
    Ht[index(k, j)]     = c * tmp1 - beta * s * tmp2;
    Ht[index(k + 1, j)] = s * tmp1 + beta * c * tmp2;
  }

  vout.paranoiac(m_vl, "Ht[%d, %d] = (%e, %e)\n", k+1, k,
                 real(Ht[index(k+1, k)]), imag(Ht[index(k+1, k)]));

  for (int j = 0; j < Nk; ++j) {
    complex_t tmp1 = Ht[index(j, k)];
    complex_t tmp2 = Ht[index(j, k + 1)];
    Ht[index(j, k)]     = tmp1 * c - tmp2 * conj(beta) * s;
    Ht[index(j, k + 1)] = tmp1 * s + tmp2 * conj(beta) * c;
  }

  vout.paranoiac(m_vl, "Ht[index(%d, %d) = (%e, %e)\n", k+1, k,
                 real(Ht[index(k+1, k)]), imag(Ht[index(k+1, k)]));

  for (int j = 0; j < Nk; ++j) {
    complex_t tmp1 = Qt[index(j, k)];
    complex_t tmp2 = Qt[index(j, k + 1)];
    Qt[index(j, k)]     = tmp1 * c - tmp2 * conj(beta) * s;
    Qt[index(j, k + 1)] = tmp1 * s + tmp2 * conj(beta) * c;
  }

  //- Givens transformations
  for (int k = kmin; k < kmax - 1; ++k) {

    // stabilized version
    real_t f1 = abs(Ht[index(k,  k-1)]);
    real_t f2 = abs(Ht[index(k+1,k-1)]);
    complex_t c, s, alph, beta;
    if(f1 + f2 == 0.0){
      vout.paranoiac(m_vl, "f1 + f2 = 0\n");
      continue;
    }else if(f2 == 0.0){
      vout.paranoiac(m_vl, "f2 = 0\n");
      c = cmplx(1.0, 0.0);
      s = cmplx(0.0, 0.0);
      complex_t tmp1 = cmplx(real_t(1.0)/f1, real_t(0.0));
      alph = conj(Ht[index(k, k-1)]) * tmp1;
      beta = cmplx(1.0, 0.0);
    }else if(f1 == 0.0){
      vout.paranoiac(m_vl, "f1 = 0\n");
      c = cmplx(0.0, 0.0);
      s = cmplx(1.0, 0.0);
      alph = cmplx(1.0, 0.0);
      complex_t tmp2 = cmplx(real_t(1.0)/f2, real_t(0.0));
      beta = conj(Ht[index(k+1,k-1)]) * tmp2;
    }else{
      real_t Fden = 1.0/pythag(f1, f2);
      complex_t tmp1 = cmplx(real_t(1.0)/f1, real_t(0.0));
      alph = conj(Ht[index(k,  k-1)]) * tmp1;
      complex_t tmp2 = cmplx(real_t(1.0)/f2, real_t(0.0));
      beta = conj(Ht[index(k+1,k-1)]) * tmp2;
      c = cmplx( f1 * Fden, real_t(0.0));
      s = cmplx(-f2 * Fden, real_t(0.0));
    }

    for (int j = 0; j < Nk; ++j) {
      complex_t tmp1 = Ht[index(k, j)];
      complex_t tmp2 = Ht[index(k + 1, j)];
      Ht[index(k, j)]     = alph * c * tmp1 - beta * s * tmp2;
      Ht[index(k + 1, j)] = alph * s * tmp1 + beta * c * tmp2;
    }
    for (int j = 0; j < Nk; ++j) {
      complex_t tmp1  = Ht[index(j, k)];
      complex_t tmp2  = Ht[index(j, k + 1)];
      complex_t alphc = conj(alph);
      complex_t betac = conj(beta);
      Ht[index(j, k)]     = tmp1 * alphc * c - tmp2 * betac * s;
      Ht[index(j, k + 1)] = tmp1 * alphc * s + tmp2 * betac * c;
    }

    if (k < kmax - 2) {
      Ht[index(k + 1, k - 1)] = cmplx(real_t(0.0), real_t(0.0));
    }

    for (int j = 0; j < Nk; ++j) {
      complex_t tmp1 = Qt[index(j, k)];
      complex_t tmp2 = Qt[index(j, k + 1)];
      Qt[index(j, k)]     = tmp1 * conj(alph) * c - tmp2 * conj(beta) * s;
      Qt[index(j, k + 1)] = tmp1 * conj(alph) * s + tmp2 * conj(beta) * c;
    }

  }

  k = kmax - 1;
  complex_t tmp1 = Ht[index(k, k - 1)];
  complex_t tmpd = abs(Ht[index(k, k)]);

  if (abs(tmp1) + abs(tmpd) != abs(tmpd)) {
    complex_t tmp2 = cmplx(real_t(1.0)/abs(tmp1), real_t(0.0));
    complex_t alph = conj(tmp1) * tmp2;

    for (int j = 0; j < Nk; ++j) {
      Ht[index(k, j)] *= alph;
    }

    for (int j = 0; j < Nk; ++j) {
      Ht[index(j, k)] *= conj(alph);
    }

    for (int j = 0; j < Nk; ++j) {
      Qt[index(j, k)] *= conj(alph);
    }
  }

  // enforce Hessenberg form to Ht
  for (int j = 0; j < Nk; ++j) {
    for (int k = j + 2; k < Nk; ++k) {
      Ht[index(k, j)] = cmplx(0.0, 0.0);
    }
  }
  for (int j = 0; j < Nk - 1; ++j) {
    int k = j + 1;
    Ht[index(k, j)] = cmplx(real(Ht[index(k, j)]), real_t(0.0));
  }

  // reunitarization of Qt
  reunit_Qt(Qt, Nk);
}

//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRArnoldi_block<FIELD, FOPR>::check_Qt(
  const int Nk, const int Nm,
  std::vector<complex_t>& Qt,
  std::vector<complex_t>& Ht,
  std::vector<complex_t>& At)
{  // confirm At = Qt * Ht * Qt^\dag
  std::vector<complex_t> Bt(Nm * Nm);
  std::vector<complex_t> Ct(Nm * Nm);

  for (int i = 0; i < Nk; ++i) {
    for (int j = 0; j < Nk; ++j) {
      Bt[index(i, j)] = cmplx(0.0, 0.0);
      for (int k = 0; k < Nk; ++k) {
        Bt[index(i, j)] += Ht[index(i, k)] * conj(Qt[index(j, k)]);
      }
    }
  }

  for (int i = 0; i < Nk; ++i) {
    for (int j = 0; j < Nk; ++j) {
      Ct[index(i, j)] = cmplx(0.0, 0.0);
      for (int k = 0; k < Nk; ++k) {
        Ct[index(i, j)] += Qt[index(i, k)] * Bt[index(k, j)];
      }
    }
  }

  real_t res = 0.0;
  for (int i = 0; i < Nk; ++i) {
    for (int j = 0; j < Nk; ++j) {
      Ct[index(i, j)] -= At[index(i, j)];
      res             += abs(Ct[index(i, j)]) * abs(Ct[index(i, j)]);
    }
  }

  // vout.detailed(m_vl, "  check A = Q^dag H Q: res = %e\n", res);
  if (res > 1.e-24) {
    vout.crucial(m_vl, "%s: check A = Q^d H Q: too large res = %e\n",
                 class_name.c_str(), res);
  }
}


//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRArnoldi_block<FIELD, FOPR>::shift_wilkinson(
  complex_t& kappa,
  const complex_t a, const complex_t b,
  const complex_t c, const complex_t d)
{
  real_t s = abs(a) + abs(b) + abs(c) + abs(d);

  complex_t p = cmplx(real_t(0.5), real_t(0.0)) * (a - d);
  complex_t q = b * c;
  complex_t r = sqrt(p * p + q);

  real_t abst1 = abs(p) * abs(p) + abs(r) * abs(r);
  real_t abst2 = 2.0 * real(p * conj(r));

  real_t absr1 = abst1 + abst2;
  real_t absr2 = abst1 - abst2;

  complex_t rmax, rmin;
  if (absr1 > absr2) {
    rmax = p + r;
  } else {
    rmax = p - r;
  }
  rmin = -q / rmax;

  kappa = rmin + d;

  // check
  complex_t res1 = rmax * rmax - (a - d) * rmax - b * c;
  complex_t res2 = rmin * rmin - (a - d) * rmin - b * c;
  real_t    res  = (abs(res1) + abs(res2)) / s;

  //vout.detailed(m_vl, " rmax  = (%e, %e)\n", real(rmax), imag(rmax));
  //vout.detailed(m_vl, " rmim  = (%e, %e)\n", real(rmin), imag(rmin));
  //vout.detailed(m_vl, " res   = %e  s = %e\n", res, s);
  //vout.detailed(m_vl, " shift = (%e, %e)\n", real(kappa), imag(kappa));

}


//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRArnoldi_block<FIELD, FOPR>::eigenvector_Ht(
                                          std::vector<complex_t>& Yt,
                                          std::vector<complex_t>& Ht,
                                          int km, int Nm)
{
  int i = 0;
  Yt[index(0, i)] = cmplx(1.0, 0.0);
  for (int j = 1; j < Nm; ++j) {
    Yt[index(j, i)] = cmplx(0.0, 0.0);
  }

  for (int i = 1; i < km; ++i) {
    complex_t lambda = Ht[index(i, i)];

    Yt[index(i, i)] = cmplx(1.0, 0.0);

    for (int j = i + 1; j < Nm; ++j) {
      Yt[index(j, i)] = cmplx(0.0, 0.0);
    }

    Yt[index(i - 1, i)] = - Ht[index(i-1, i)]
                           / (Ht[index(i-1, i-1)] - lambda);

    for (int j = i - 2; j >= 0; --j) {
      Yt[index(j, i)] = -Ht[index(j, i)];
      for (int k = j + 1; k < i; ++k) {
        Yt[index(j, i)] += -Ht[index(j, k)] * Yt[index(k, i)];
      }
      Yt[index(j, i)] *= real_t(1.0) / (Ht[index(j, j)] - lambda);
    }
  }
}


//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRArnoldi_block<FIELD, FOPR>::mult_Qt(
  std::vector<complex_t>& Yt,
  std::vector<complex_t>& Qt,
  std::vector<complex_t>& Xt,
  int km, int Nm)
{  //  Yt = Qt * Xt.
  for (int j = 0; j < km; ++j) {
    for (int i = 0; i < Nm; ++i) {
      Yt[index(i, j)] = cmplx(0.0, 0.0);
      for (int k = 0; k < Nm; ++k) {
        Yt[index(i, j)] += Qt[index(i, k)] * Xt[index(k, j)];
      }
    }
  }
}


//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRArnoldi_block<FIELD, FOPR>::check_eigen_Ht(
  std::vector<complex_t>& Ht,
  std::vector<complex_t>& TDa,
  std::vector<complex_t>& Xt,
  int km, int Nm)
{
  real_t crit = 1.0e-20;

  real_t diff = 0.0;
  for (int i = 0; i < km; ++i) {
    real_t diff1 = 0.0;
    for (int j = 0; j < Nm; ++j) {
      complex_t diff2 = -TDa[i] * Xt[index(j, i)];
      for (int k = 0; k < km; ++k) {
        diff2 += Ht[index(j, k)] * Xt[index(k, i)];
      }
      diff1 += abs(diff2) * abs(diff2);
    }

    diff += diff1;
  }

  diff = diff / real_t(km * Nm);

  //vout.detailed(m_vl, "  eigenrelation of Ht: diff2 = %e\n", diff);
  vout.paranoiac(m_vl, "  eigenrelation of Ht: diff2 = %e\n", diff);
  if (diff > crit) {
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRArnoldi_block<FIELD, FOPR>::reunit_Qt(
  std::vector<complex_t>& Qt,
  int Nk)
{
  for (int j = 0; j < Nk; ++j) {
    for (int i = 0; i < j; ++i) {
      complex_t qq = cmplx(real_t(0.0), real_t(0.0));
      for (int k = 0; k < Nk; ++k) {
        //qq += conj(Qt[index(k,j)]) * Qt[index(k,i)];
        qq += conj(Qt[index(j, k)]) * Qt[index(i, k)];
      }

      for (int k = 0; k < Nk; ++k) {
        //Qt[index(k,i)] -= qq * Qt[index(k,j)];
        Qt[index(i, k)] -= qq * Qt[index(j, k)];
      }
    }

    real_t qq = 0.0;
    for (int k = 0; k < Nk; ++k) {
      //real_t qq1 = abs(Qt[index(k,j)]);
      real_t qq1 = abs(Qt[index(j, k)]);
      qq += qq1 * qq1;
    }
    qq = 1.0 / sqrt(qq);
    for (int k = 0; k < Nk; ++k) {
      //Qt[index(k,j)] *= qq;
      Qt[index(j, k)] *= qq;
    }
  }
}


//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRArnoldi_block<FIELD, FOPR>::reconst_Ht(
                          std::vector<std::vector<complex_t> >& Ht,
                          std::vector<complex_t>& beta,
                          std::vector<FIELD>& vk,
                          int Nk)
{
  int ith = ThreadManager::get_thread_id();

  std::vector<complex_t> u(Nk);
  std::vector<complex_t> x(Nk);

  for (int j = 0; j < Nk; ++j) {

    m_fopr->mult_dd(m_v, vk[j]);

    for (int i = 0; i < Nk; ++i) {

      block_dotc(&m_bvec1[0], vk[i], m_v, *m_index_block);
#pragma omp barrier

      if(ith == 0){
        for(int block = 0; block < m_num_blocks; ++block){
          Ht[block][index(i, j)] = m_bvec1[block];
        }
      }
#pragma omp barrier
    }
  }
  // now m_v = A * vk[Nk-1]

  int j = Nk - 1;
  for (int i = 0; i < Nk; ++i) {

    if(ith == 0){
      for(int block = 0; block < m_num_blocks; ++block){
        m_bvec1[block] = -Ht[block][index(i, j)];
      }
    }
#pragma omp barrier

    block_axpy(m_v, &m_bvec1[0], vk[i], real_t(1.0), *m_index_block);
#pragma omp barrier
  }

  block_norm2(&m_rvec1[0], m_v, *m_index_block);

#pragma omp barrier

  if(ith == 0){
    for(int block = 0; block < m_num_blocks; ++block){
      beta[block] = cmplx(sqrt(m_rvec1[block]), real_t(0.0));

      // enforce Hessenberg form to Ht
      for (int j = 0; j < Nk; ++j) {
        for (int k = j + 2; k < Nk; ++k) {
          Ht[block][index(k, j)] = cmplx(real_t(0.0), real_t(0.0));
        }
      }

      for (int j = 0; j < Nk - 1; ++j) {
        int k = j + 1;
        Ht[block][index(k, j)]
                = cmplx(real(Ht[block][index(k, j)]), real_t(0.0));
      }
    }
  }

#pragma omp barrier

}

//============================================================END=====
