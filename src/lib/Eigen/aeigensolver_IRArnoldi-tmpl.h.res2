/*!
        @file    aeigensolver_IRArnoldi-tmpl.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "lib/Eigen/aeigensolver_IRArnoldi.h"

#ifdef __PGI
#pragma global opt=1
// This is becaouse for PGI compiler on some environment,
// compilation with optimization level higher than O2 fails.
//                                  [12 Feb 2021 H.Matsufuru]
#endif

using std::string;

template<typename FIELD, typename FOPR>
const std::string AEigensolver_IRArnoldi<FIELD, FOPR>::
class_name = "AEigensolver_IRArnoldi";

//====================================================================
template<typename FIELD, typename FOPR>
AEigensolver_IRArnoldi<FIELD, FOPR>::~AEigensolver_IRArnoldi()
{
  if (m_sorter) delete m_sorter;
}


//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRArnoldi<FIELD, FOPR>::set_parameters(
  const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  string str_sortfield_type;
  int    Nk, Np;
  int    Niter_eigen;
  double Enorm_eigen, Vthreshold;

  int err = 0;
  err += params.fetch_string("eigensolver_mode", str_sortfield_type);
  err += params.fetch_int("number_of_wanted_eigenvectors", Nk);
  err += params.fetch_int("number_of_working_eigenvectors", Np);
  err += params.fetch_int("maximum_number_of_iteration", Niter_eigen);
  err += params.fetch_double("convergence_criterion_squared", Enorm_eigen);
  err += params.fetch_double("threshold_value", Vthreshold);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(str_sortfield_type, Nk, Np, Niter_eigen, Enorm_eigen,
                 Vthreshold);
}


//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRArnoldi<FIELD, FOPR>::get_parameters(Parameters& params) const
{
  params.set_string("eigensolver_mode", m_sort_type);
  params.set_int("number_of_wanted_eigenvectors", m_Nk);
  params.set_int("number_of_working_eigenvectors", m_Np);
  params.set_int("maximum_number_of_iteration", m_Niter_eigen);
  params.set_double("convergence_criterion_squared", m_Enorm_eigen);
  params.set_double("threshold_value", m_Vthreshold);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRArnoldi<FIELD, FOPR>::set_parameters(
  const std::string& sort_type,
  int Nk, int Np,
  int Niter_eigen,
  double Enorm_eigen,
  double Vthreshold)
{
  if (m_sorter) delete m_sorter;

  m_sort_type = sort_type;
  m_sorter    = new Sorter<complex_t>(sort_type);

  set_parameters(Nk, Np, Niter_eigen, Enorm_eigen, Vthreshold);
}


//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRArnoldi<FIELD, FOPR>::set_parameters(
  int Nk, int Np,
  int Niter_eigen,
  double Enorm_eigen,
  double Vthreshold)
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

  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Nk          = %d\n", m_Nk);
  vout.general(m_vl, "  Np          = %d\n", m_Np);
  vout.general(m_vl, "  Niter_eigen = %d\n", m_Niter_eigen);
  vout.general(m_vl, "  Enorm_eigen = %16.8e\n", m_Enorm_eigen);
  vout.general(m_vl, "  Vthreshold  = %16.8e\n", m_Vthreshold);

  int Nm = m_Nk + m_Np;

  m_TDa2.resize(Nm);

  m_Qt.resize(Nm * Nm);
  m_Iconv.resize(Nm);

  m_Ht.resize((Nm + 1) * Nm);  // (Nm+1) * Nm matrix
  m_Ht2.resize((Nm + 1) * Nm);
  m_Yt.resize(Nm * Nm);

  int Nin  = m_fopr->field_nin();
  int Nvol = m_fopr->field_nvol();
  int Nex  = m_fopr->field_nex();

  m_B.resize(Nm);
  for (int k = 0; k < Nm; ++k) {
    m_B[k].reset(Nin, Nvol, Nex);
  }

  m_f.reset(Nin, Nvol, Nex);
  m_v.reset(Nin, Nvol, Nex);
}


//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRArnoldi<FIELD, FOPR>::solve(
  std::vector<complex_t>& TDa,
  std::vector<FIELD>& vk,
  int& Nsbt, int& Nconv, const FIELD& b)
{
  int    Nk          = m_Nk;
  int    Np          = m_Np;
  int    Niter_eigen = m_Niter_eigen;
  real_t Enorm_eigen = m_Enorm_eigen;
  real_t Vthreshold  = m_Vthreshold;

  real_t Enorm_almost = 1.e-12; // this is an experimental value

  int Nm = Nk + Np;

  if (Nk + Np > TDa.size()) {
    vout.crucial(m_vl, "Error at %s: std::vector TDa is too small.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  } else if (Nk + Np > vk.size()) {
    vout.crucial(m_vl, "Error at %s: std::vector vk is too small.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  vout.general(m_vl, "Implicitly Restarted Arnoldi algorithm start\n");
  vout.general(m_vl, "  Nk = %d  Np = %d\n", Nk, Np);
  vout.general(m_vl, "  Nm = %d\n", Nm);
  vout.general(m_vl, "  size of TDa  = %d\n", TDa.size());
  vout.general(m_vl, "  size of vk   = %d\n", vk.size());
  vout.general(m_vl, "  Enorm_eigen  = %e\n", Enorm_eigen);
  vout.general(m_vl, "  Enorm_almost = %e\n", Enorm_almost);

  Nconv = -1;
  Nsbt  = 0;

  int k1 = 0;
  int k2 = Nk;

  int Kdis       = 0;
  int Kthreshold = 0;
  int Kalmost    = 0;

  int nconv;

#pragma omp parallel
  {
    //- set initial vector
    real_t bnorm2 = b.norm2();
    if (bnorm2 > 1.e-12) { // b is regarded as a nonzero vector
      copy(vk[0], b);
    } else {
      vk[0].set(1.0);  // start with a uniform vector
    }

    real_t vnorm = vk[0].norm2();
    vnorm = 1.0 / sqrt(vnorm);
    scal(vk[0], vnorm);

#pragma omp master
    {
      for (int j = 0; j < Nm; ++j) {
        for (int k = 0; k < Nm + 1; ++k) {
          m_Ht[index(k, j)]  = cmplx(real_t(0.0), real_t(0.0));
          m_Ht2[index(k, j)] = cmplx(real_t(0.0), real_t(0.0));
        }
      }
    }
#pragma omp barrier

    //- initial Nk steps
    for (int k = 0; k < k2; ++k) {
      step(Nm, k, vk, m_f);
    }

#pragma omp barrier

    //- restarting loop begins
    for (int iter = 0; iter < Niter_eigen; ++iter) {
      vout.detailed(m_vl, "\n iter = %d\n", iter);

      // Convergence test

      complex_t beta = cmplx(real_t(0.0), real_t(0.0));

#pragma omp master
      {
        // here size Nk eigenproblem is solved.
        beta = m_Ht[index(k2, k2 - 1)];
        m_Ht[index(k2, k2 - 1)] = cmplx(real_t(0.0), real_t(0.0));

        m_Ht2 = m_Ht;
        setUnit_Qt(Nm, m_Qt);

        tqri(m_Ht2, k1, Nk, Nm, m_Qt, nconv);

        for (int j = k1; j < Nk; ++j) {
          TDa[j] = m_Ht2[index(j, j)];
        }

        // here small eigenvalue problem should be checked for confirmation
        std::vector<complex_t> Xt(Nm * Nm);

        eigenvector_Ht(Xt, m_Ht2, Nk, Nm);
        mult_Qt(m_Yt, m_Qt, Xt, Nk, Nm);
        check_eigen_Ht(m_Ht, TDa, m_Yt, Nk, Nm);
        // now m_Yt is eigenvectors of Ht

        Kdis       = 0;
        Kthreshold = 0;
        Kalmost    = 0;
      }
#pragma omp barrier

      // convergenece test
      vout.general(m_vl, "                   TDa            "
                         "        abs           diff2\n");

      for (int j = 0; j < k1; ++j) {
        m_fopr->mult(m_v, m_B[j]);
        axpy(m_v, -TDa[j], m_B[j]);
        real_t bnorm = m_B[j].norm2();
        real_t vv    = m_v.norm2() / bnorm;
        vout.general(m_vl, " %4d* (%12.8f, %12.8f)  %12.8f  %12.4e\n",
                     j, real(TDa[j]), imag(TDa[j]), abs(TDa[j]), vv);
      }

      for (int j = k1; j < Nk; ++j) {
        m_B[j].set(0.0);
        for (int k = 0; k < Nk; ++k) {
          axpy(m_B[j], m_Yt[index(k, j)], vk[k]);
        }
        m_fopr->mult(m_v, m_B[j]);
        axpy(m_v, -TDa[j], m_B[j]);
        real_t bnorm = m_B[j].norm2();
        real_t vv    = m_v.norm2() / bnorm;
        vout.general(m_vl, " %4d  (%12.8f, %12.8f)  %12.8f  %12.4e\n",
                     j, real(TDa[j]), imag(TDa[j]), abs(TDa[j]), vv);

#pragma omp master
        {
          m_TDa2[j] = TDa[j];

          if (vv < Enorm_eigen) {
            m_Iconv[Kdis] = j;
            ++Kdis;
            if (!m_sorter->comp(m_TDa2[j], Vthreshold)) {
              ++Kthreshold;
            }
          } else if ((vv < Enorm_almost) &&
                     m_sorter->comp(m_TDa2[j], Vthreshold)) {
            ++Kalmost;
          }
        } // omp master
#pragma omp barrier
      }   // j-loop end

      vout.detailed(m_vl, " #modes already converged: %4d\n", k1);
      vout.detailed(m_vl, " #modes newly converged:   %4d\n", Kdis);
      vout.detailed(m_vl, " #modes almost converged:  %4d\n", Kalmost);

      if ((Kthreshold > 0) && (Kalmost == 0)) {
#pragma omp master
        {
          // (there is a converged eigenvalue larger than Vthreshold.)

          //- Sorting
          for (int i = k1; i < Kdis; ++i) {
            TDa[i] = m_TDa2[m_Iconv[i]];
          }
          Nsbt  = Kdis - Kthreshold;
          Nconv = iter;
        }
#pragma omp barrier

        // break;
        iter = Niter_eigen;
      }

      // deflation - this does not work yet

      /*
      if(Kdis > 0){
        deflation(k1, k2, Kdis, TDa, vk, beta);
        k1 += Kdis;
      }
      */

#pragma omp master
      {
        m_Ht[index(k2, k2 - 1)] = beta;
      }
#pragma omp barrier

      // if not converged, Krylov subspace is extended


      for (int k = k2 - 1; k < Nm; ++k) {
        step(Nm, k, vk, m_f);
      }

#pragma omp master
      {
        for (int j = 0; j < Nm; ++j) {
          for (int k = 0; k < Nm; ++k) {
            m_Ht2[index(k, j)] = m_Ht[index(k, j)];
          }
        }

        //- getting eigenvalues
        setUnit_Qt(Nm, m_Qt);

        // vout.detailed(m_vl, "QR transformation\n");

        tqri(m_Ht2, k1, Nm, Nm, m_Qt, nconv);

        for (int j = 0; j < Nm; ++j) {
          m_TDa2[j] = m_Ht2[index(j, j)];
        }

        // here small eigenvalue problem should be checked for confirmation
        std::vector<complex_t> Xt(Nm * Nm);

        eigenvector_Ht(Xt, m_Ht2, Nm, Nm);
        mult_Qt(m_Yt, m_Qt, Xt, Nm, Nm);
        check_eigen_Ht(m_Ht, m_TDa2, m_Yt, Nm, Nm);
        // now m_Yt is eigenvectors of Ht

        //- sorting
        m_sorter->sort(m_TDa2, Nm);
        // CAUTION: if already converged eigenvalues become the
        //   i >= k2 range, the following algorithm may fail.
        //   Some care is necessary.

        // m_Ht2 = m_Ht;

        setUnit_Qt(Nm, m_Qt);

        //- implicitly shifted QR transformations
        for (int ip = k2; ip < Nm; ++ip) {
          complex_t Dsh  = m_TDa2[ip];
          int       kmin = k1 + 1;
          int       kmax = Nm;
          // vout.general(m_vl, "ip = %d  Dsh = (%f, %f)\n",
          //              ip, real(Dsh), imag(Dsh));
          qrtrf(m_Ht, Nm, Nm, m_Qt, Dsh, kmin, kmax);
          // check_Qt(Nm, Nm, m_Qt, m_Ht, m_Ht2);
        }
      }
#pragma omp barrier

      for (int j = k1; j <= k2; ++j) {
        m_B[j].set(0.0);
        //for (int k = 0; k < Nm; ++k) {
        for (int k = k1; k < Nm; ++k) {
          axpy(m_B[j], m_Qt[index(k, j)], vk[k]);
        }
      }

      for (int j = k1; j <= k2; ++j) {
        copy(vk[j], m_B[j]);
      }

      schmidt(k2, vk);

      //- Compressed vector f and beta(k2)
      copy(m_f, m_B[k2]);

      real_t ff = m_f.norm2();
      // vout.detailed(m_vl, " ff = %20.14f\n", ff);

      scal(m_f, m_Ht[index(k2, k2 - 1)]);

      real_t beta_k = m_f.norm2();
      beta_k = sqrt(beta_k);

      vout.detailed(m_vl, " beta(k) = %20.14f\n", beta_k);

      // reconstructing Hessenberg matrix
      reconst_Ht(m_Ht2, m_Qt, beta, vk, Nk);

#pragma omp master
      {
        real_t diff = 0.0;
        for (int i = 0; i < Nk; ++i) {
          for (int j = 0; j < Nk; ++j) {
            real_t diff1 = abs(m_Ht[index(i, j)] - m_Ht2[index(i, j)]);
            diff += diff1 * diff1;
          }
        }
        diff = diff / real_t(Nk * Nk);
        vout.general(m_vl, " difference of Hessenberg = %e\n", diff);

        /*
        if(diff > 1.0e-16){
          for(int i1 = 0; i1 < Nk; ++i1){
           for(int i2 = 0; i2 < Nk; ++i2){
             vout.general(m_vl, " Ht[%d,%d] = (%e,%e)  (%e,%e)\n", i2, i1,
                       real(m_Ht[index(i2,i1)]), imag(m_Ht[index(i2,i1)]),
                       real(m_Ht2[index(i2,i1)]), imag(m_Ht2[index(i2,i1)]));
           }
          }
        }
        */

        m_Ht = m_Ht2;
      }
#pragma omp barrier
    } // end of iter loop
  }   // omp parallel

  std::vector<int> idx = m_sorter->sort_index(TDa, Kdis);
  for (int i = 0; i < Kdis; ++i) {
    copy(vk[i], m_B[m_Iconv[idx[i]]]);

    /*
    real_t vv = vk[i].norm2();
    vv = 1.0/sqrt(vv);
    scal(vk[i], vv);
    vv = vk[i].norm2();
    vout.paranoiac(m_vl, " vk[%d]:  norm2 = %22.14e\n", i, vv);
    */
  }

  if (Nconv == -1) {
    vout.crucial(m_vl, "Error at %s: NOT converged.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  } else {
    vout.general(m_vl, "\n Converged:\n");
    vout.general(m_vl, "  Nconv   = %d\n", Nconv);
    //  vout.general(m_vl, "  beta(k) = %20.14e\n", beta_k);
    vout.general(m_vl, "  Kdis    = %d\n", Kdis);
    vout.general(m_vl, "  Nsbt    = %d\n", Nsbt);
    vout.general(m_vl, "%s: solve finished.\n\n", class_name.c_str());
  }
}


//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRArnoldi<FIELD, FOPR>::deflation(
  int k1, int k2, int Kdis,
  std::vector<complex_t>& TDa,
  std::vector<FIELD>& vk,
  complex_t& beta)
{
  for (int i = k2 - 1; i >= k1 + Kdis; --i) {
    copy(vk[i], vk[i - Kdis]);
  }

  for (int i = 0; i < Kdis; ++i) {
    copy(vk[k1 + i], m_B[m_Iconv[i]]);
    if ((k1 + i) != m_Iconv[i]) copy(m_B[k1 + i], m_B[m_Iconv[i]]);
#pragma omp master
    {
      TDa[k1 + i] = m_TDa2[m_Iconv[i]];
    }
#pragma omp barrier
  }

  schmidt(k2, vk);
  reconst_Ht(m_Ht, m_Qt, beta, vk, k2);

#pragma omp master
  {
    for (int i = 0; i < k1 + Kdis; ++i) {
      m_Ht[index(i + 1, i)] = cmplx(0.0, 0.0);
    }
  }
#pragma omp barrier
}


//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRArnoldi<FIELD, FOPR>::step(int Nm, int k,
                                               std::vector<FIELD>& vk,
                                               FIELD& w)
{
  if (k >= Nm) {
    vout.crucial(m_vl, "Error at %s: k is larger than Nm.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  m_fopr->mult(w, vk[k]);

  for (int j = 0; j <= k; ++j) {
    complex_t alph = dotc(vk[j], w);
    axpy(w, -alph, vk[j]);
#pragma omp master
    m_Ht[index(j, k)] = alph;
#pragma omp barrier
  }

  real_t beta = w.norm2();
  beta = sqrt(beta);

  if (k < Nm - 1) {
    real_t beta_r = 1.0 / beta;
    copy(vk[k + 1], w);
    scal(vk[k + 1], beta_r);
  }

#pragma omp master
  m_Ht[index(k + 1, k)] = cmplx(beta, real_t(0.0));
#pragma omp barrier
}


//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRArnoldi<FIELD, FOPR>::schmidt(
  int Nk,
  std::vector<FIELD>& vk)
{
  for (int j = 0; j < Nk; ++j) {
    for (int i = 0; i < j; ++i) {
      complex_t vv = dotc(vk[i], vk[j]);
      axpy(vk[j], -vv, vk[i]);
    }
    real_t vv = vk[j].norm2();
    vv = 1.0 / sqrt(vv);
    scal(vk[j], vv);
  }
}


//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRArnoldi<FIELD, FOPR>::setUnit_Qt(
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
void AEigensolver_IRArnoldi<FIELD, FOPR>::tqri(
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
        vout.detailed(m_vl, "  tqri: converged at iter = %d\n", nconv);
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
void AEigensolver_IRArnoldi<FIELD, FOPR>::qrtrf(
  std::vector<complex_t>& Ht,
  int Nk, int Nm, std::vector<complex_t>& Qt,
  complex_t Dsh, int kmin, int kmax)
{
  int k = kmin - 1;

  real_t    f1   = abs(Ht[index(k, k)] - Dsh);
  real_t    f2   = abs(Ht[index(k + 1, k)]);
  real_t    Fden = 1.0 / sqrt(f1 * f1 + f2 * f2);
  complex_t c    = cmplx(f1 * Fden, real_t(0.0));
  complex_t s    = cmplx(-f2 * Fden, real_t(0.0));
  complex_t beta = (Ht[index(k, k)] - Dsh) / f1;

  for (int j = 0; j < Nk; ++j) {
    complex_t tmp1 = Ht[index(k, j)];
    complex_t tmp2 = Ht[index(k + 1, j)];
    Ht[index(k, j)]     = c * tmp1 - beta * s * tmp2;
    Ht[index(k + 1, j)] = s * tmp1 + beta * c * tmp2;
  }
  for (int j = 0; j < Nk; ++j) {
    complex_t tmp1 = Ht[index(j, k)];
    complex_t tmp2 = Ht[index(j, k + 1)];
    Ht[index(j, k)]     = tmp1 * c - tmp2 * conj(beta) * s;
    Ht[index(j, k + 1)] = tmp1 * s + tmp2 * conj(beta) * c;
  }

  for (int j = 0; j < Nk; ++j) {
    complex_t tmp1 = Qt[index(j, k)];
    complex_t tmp2 = Qt[index(j, k + 1)];
    Qt[index(j, k)]     = c * tmp1 - conj(beta) * s * tmp2;
    Qt[index(j, k + 1)] = s * tmp1 + conj(beta) * c * tmp2;
  }

  //- Givens transformations
  for (int k = kmin; k < kmax - 1; ++k) {
    real_t    f1   = abs(Ht[index(k, k - 1)]);
    real_t    f2   = abs(Ht[index(k + 1, k - 1)]);
    real_t    Fden = 1.0 / sqrt(f1 * f1 + f2 * f2);
    complex_t alph = conj(Ht[index(k, k - 1)]) / f1;
    complex_t beta = conj(Ht[index(k + 1, k - 1)]) / f2;
    complex_t c    = cmplx(f1 * Fden, real_t(0.0));
    complex_t s    = cmplx(-f2 * Fden, real_t(0.0));

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
      Qt[index(j, k)]     = conj(alph) * c * tmp1 - conj(beta) * s * tmp2;
      Qt[index(j, k + 1)] = conj(alph) * s * tmp1 + conj(beta) * c * tmp2;
    }
  }

  k = kmax - 1;
  complex_t tmp1 = Ht[index(k, k - 1)];
  complex_t tmpd = abs(Ht[index(k, k)]);

  if (abs(tmp1) + abs(tmpd) != abs(tmpd)) {
    complex_t tmp2 = abs(tmp1);
    complex_t alph = conj(tmp1) / tmp2;

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

  // enforce Hessenberg form to Ht
  for (int j = 0; j < Nk - 1; ++j) {
    int k = j + 1;
    Ht[index(k, j)] = cmplx(real(Ht[index(k, j)]), real_t(0.0));
  }

  // reunitarization of Qt
  reunit_Qt(Qt, Nk);
}


//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRArnoldi<FIELD, FOPR>::check_Qt(
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
void AEigensolver_IRArnoldi<FIELD, FOPR>::shift_wilkinson(
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

  /*
  vout.detailed(m_vl, " rmax  = (%e, %e)\n", real(rmax), imag(rmax));
  vout.detailed(m_vl, " rmim  = (%e, %e)\n", real(rmin), imag(rmin));
  vout.detailed(m_vl, " res   = %e  s = %e\n", res, s);
  vout.detailed(m_vl, " shift = (%e, %e)\n", real(kappa), imag(kappa));
  */
}


//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRArnoldi<FIELD, FOPR>::eigenvector_Ht(
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

    Yt[index(i - 1, i)] = -Ht[index(i - 1, i)] / (Ht[index(i - 1, i - 1)] - lambda);

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
void AEigensolver_IRArnoldi<FIELD, FOPR>::mult_Qt(
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
void AEigensolver_IRArnoldi<FIELD, FOPR>::check_eigen_Ht(
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

  vout.detailed(m_vl, "  eigenrelation of Ht: diff2 = %e\n", diff);
  if (diff > crit) {
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRArnoldi<FIELD, FOPR>::reunit_Qt(
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
void AEigensolver_IRArnoldi<FIELD, FOPR>::reconst_Ht(
  std::vector<complex_t>& Ht,
  std::vector<complex_t>& Qt,
  complex_t& beta,
  std::vector<FIELD>& vk,
  int Nk)
{
  std::vector<complex_t> u(Nk);
  std::vector<complex_t> x(Nk);

  for (int j = 0; j < Nk; ++j) {
    m_fopr->mult(m_v, vk[j]);
    for (int i = 0; i < Nk; ++i) {
      complex_t Htmp = dotc(vk[i], m_v);
#pragma omp master
      Ht[index(i, j)] = Htmp;
#pragma omp barrier
    }
  }
  // now m_v = A * vk[Nk-1]
  int j = Nk - 1;
  for (int i = 0; i < Nk; ++i) {
    axpy(m_v, -Ht[index(i, j)], vk[i]);
  }
  real_t vv = m_v.norm2();
#pragma omp master
  {
    beta = cmplx(real_t(sqrt(vv)), real_t(0.0));

    // enforce Hessenberg form to Ht
    for (int j = 0; j < Nk; ++j) {
      for (int k = j + 2; k < Nk; ++k) {
        Ht[index(k, j)] = cmplx(real_t(0.0), real_t(0.0));
      }
    }

    for (int j = 0; j < Nk - 1; ++j) {
      int k = j + 1;
      Ht[index(k, j)] = cmplx(real(Ht[index(k, j)]), real_t(0.0));
    }
  }
#pragma omp barrier

  // The following code intends to apply Householder transformation
  // to enforce the Hessenberg form to Ht.
  // This code is not used in current version of algorithm,
  // since it seems not appropriate for present situation.
  // Nevertheless the code is retained for potential use in future.
  // Note that it may contain unfound bugs.
  //                                       [24 Feb 2020 H.Matsufuru]

  // Householder

  /*
  for(int j = 0; j < Nk-1; ++j){

    u[j] = cmplx(0.0, 0.0);
    x[j] = cmplx(0.0, 0.0);

    complex_t asd  = Ht[index(j+1,j)];
    real_t    asd2 = abs(asd) * abs(asd);

    real_t x2 = 0.0;
    for(int k = j+2; k < Nk; ++k){
      real_t x1 = abs(Ht[index(k,j)]);
      x2 += x1 * x1;
      u[k] = Ht[index(k,j)];
    }

    real_t sqx2 = 0.0;
    if(x2 > 0.0) sqx2 = sqrt(x2);

    if(Ht[index(j,j)] + sqx2 == Ht[index(j,j)]){
      // in this case, deflation is better.
      continue;
    }

    real_t    zabs = abs(asd) - sqrt(asd2 + x2);
    complex_t ph   = cmplx(1.0, 0.0);
    if(abs(asd) > 1.e-12) ph = asd/abs(asd);
    complex_t zeta = ph * cmplx(zabs, 0.0);
    real_t   uabs2 = zabs * zabs + x2;

    vout.general(m_vl, " j = %d  uabs2 = %e  x2 = %e\n",j, uabs2, x2);

    u[j+1] = zeta;
    Ht[index(j+1,j)] = asd - zeta;

    for(int k = j+2; k < Nk; ++k){
      Ht[index(k,j)] = cmplx(0.0, 0.0);
    }

    // mult H(u) to H from left
    for(int i = j+1; i < Nk; ++i){

      complex_t prod = cmplx(0.0, 0.0);
      for(int k = j+1; k < Nk; ++k){
        prod += conj(u[k]) * Ht[index(k,i)];
      }
      prod = prod * 2.0/uabs2;

      for(int k = j+1; k < Nk; ++k){
        Ht[index(k,i)] -= prod * u[k];
      }

    }

    // mult H(u) to H from right
    for(int i = 0; i < Nk; ++i){

      complex_t prod = cmplx(0.0, 0.0);
      for(int k = j+1; k < Nk; ++k){
        prod += conj(Ht[index(i,k)]) * u[k];
      }
      prod = prod * 2.0/uabs2;

      for(int k = j+1; k < Nk; ++k){
        Ht[index(i,k)] -= prod * conj(u[k]);
      }

    }

    // mult H(u) to Q from left
    for(int i = j+1; i < Nk; ++i){

      complex_t prod = cmplx(0.0, 0.0);
      for(int k = j+1; k < Nk; ++k){
        prod += conj(u[k]) * Qt[index(k,i)];
      }
      prod = prod * 2.0/uabs2;
      for(int k = j+1; k < Nk; ++k){
        Qt[index(k,i)] -= prod * u[k];
      }

    }

    // mult phase factor to H
    Ht[index(j+1,j)] *= conj(ph);
    for(int k = j+1; k < Nk; ++k){
      Ht[index(j+1,k)] *= ph;
    }

    // mult phase factor to Q
    for(int k = 0; k < Nk; ++k){
      Qt[index(j,k)] *= ph;
    }

  } // j-loop

  // mult Qt^d to V from right
  for(int j = 0; j < Nk; ++j){
    m_B[j].set(0.0);
    for (int k = 0; k < Nk; ++k) {
      axpy(m_B[j], conj(Qt[index(j,k)]), vk[k]);
    }
  }
  */
}


//====================================================================
//============================================================END=====
