/*!
        @file    aeigensolver_IRLanczos-tmpl.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate:: 2022-05-07 23:30:01 #$

        @version $LastChangedRevision: 2373 $
*/

#include "aeigensolver_IRLanczos.h"

#ifdef __PGI
#pragma global opt=1
// This is becaouse for PGI compiler on some environment,
// compilation with optimization level higher than O2 fails.
//                                  [12 Feb 2021 H.Matsufuru]
#endif

using std::string;

template<typename FIELD, typename FOPR>
const std::string AEigensolver_IRLanczos<FIELD, FOPR>::
class_name = "AEigensolver_IRLanczos";

//====================================================================
template<typename FIELD, typename FOPR>
AEigensolver_IRLanczos<FIELD, FOPR>::~AEigensolver_IRLanczos()
{
  if (m_sorter) delete m_sorter;
}


//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRLanczos<FIELD, FOPR>::set_parameters(
                                            const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  string str_sortfield_type;
  int    Nk, Np;
  int    Nk_conv;
  int    Niter_eigen;
  double Enorm_eigen, Vthreshold;

  int err = 0;
  err += params.fetch_string("eigensolver_mode", str_sortfield_type);
  err += params.fetch_int("number_of_wanted_eigenvectors", Nk);
  err += params.fetch_int("number_of_working_eigenvectors", Np);
  // err += params.fetch_int("number_of_converged_eigenvectors", Nk_conv);
  err += params.fetch_int("maximum_number_of_iteration", Niter_eigen);
  err += params.fetch_double("convergence_criterion_squared", Enorm_eigen);
  err += params.fetch_double("threshold_value", Vthreshold);

  if(params.fetch_int("number_of_converged_eigenvectors", Nk_conv)){
    Nk_conv = 0;  // Nk_conv is not given as parameter
  }

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(str_sortfield_type, Nk, Np, Nk_conv, Niter_eigen, Enorm_eigen,
                 Vthreshold);

}

//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRLanczos<FIELD, FOPR>::set_parameters(
                                       const std::string& sort_type,
                                       int Nk, int Np,
                                       int Nk_conv,
                                       int Niter_eigen,
                                       double Enorm_eigen,
                                       double Vthreshold)
{
  if (m_sorter) delete m_sorter;

  m_sorter = new Sorter<real_t>(sort_type);

  set_parameters(Nk, Np, Nk_conv, Niter_eigen, Enorm_eigen, Vthreshold);

}

//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRLanczos<FIELD, FOPR>::set_parameters(
                                             int Nk, int Np,
                                             int Nk_conv,
                                             int Niter_eigen,
                                             double Enorm_eigen,
                                             double Vthreshold)
{
  //- range check
  int err = 0;
  err += ParameterCheck::non_negative(Nk);
  err += ParameterCheck::non_negative(Np);
  err += ParameterCheck::non_negative(Nk_conv);
  err += ParameterCheck::non_negative(Niter_eigen);
  err += ParameterCheck::square_non_zero(Enorm_eigen);
  // NB. Vthreshold == 0 is allowed.

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  //- store values
  m_Nk_conv     = Nk_conv;
  m_Nk          = Nk;
  m_Np          = Np;
  m_Niter_eigen = Niter_eigen;
  m_Enorm_eigen = real_t(Enorm_eigen);
  m_Vthreshold  = real_t(Vthreshold);

  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Nk          = %d\n", m_Nk);
  vout.general(m_vl, "  Np          = %d\n", m_Np);
  vout.general(m_vl, "  Nk_conv     = %d\n", m_Nk_conv);
  vout.general(m_vl, "  Niter_eigen = %d\n", m_Niter_eigen);
  vout.general(m_vl, "  Enorm_eigen = %16.8e\n", m_Enorm_eigen);
  vout.general(m_vl, "  Vthreshold  = %16.8e\n", m_Vthreshold);

  int Nm = m_Nk + m_Np;
  m_TDb.resize(Nm);
  m_TDa2.resize(Nm);
  m_TDb2.resize(Nm);
  m_Qt.resize(Nm * Nm);
  m_Iconv.resize(Nm);

  int Nin  = m_fopr->field_nin();
  int Nvol = m_fopr->field_nvol();
  int Nex  = m_fopr->field_nex();

  m_B.resize(Nm);
  for (int k = 0; k < Nm; ++k) {
    m_B[k].reset(Nin, Nvol, Nex);
  }

  m_f.reset(Nin, Nvol, Nex);
  m_v.reset(Nin, Nvol, Nex);

  // do not use Ritz pairs wiht large Ritz value in checking convergence
  m_Enorm_ignore_factor = 1.0e5 * 1.e5;

}

//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRLanczos<FIELD, FOPR>::solve(
                                   std::vector<complex_t>& TDac,
                                   std::vector<FIELD>& vk,
                                   int& Nsbt, int& Nconv, const FIELD& b)
{
  int size = TDac.size();
  std::vector<real_t> TDa(size);

  solve(TDa, vk, Nsbt, Nconv, b);

  for(int i = 0; i < size; ++i){
    TDac[i] = cmplx(TDa[i], real_t(0.0));
  }


}

//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRLanczos<FIELD, FOPR>::solve(
                                   std::vector<real_t>& TDa,
                                   std::vector<FIELD>& vk,
                                   int& Nsbt, int& Nconv, const FIELD& b)
{
  int    Nk          = m_Nk;
  int    Np          = m_Np;
  int    Niter_eigen = m_Niter_eigen;
  real_t Enorm_eigen = m_Enorm_eigen;
  real_t Vthreshold  = m_Vthreshold;
  real_t Enorm_ignore_factor = m_Enorm_ignore_factor;
  std::vector<real_t> Enorm(Nk);

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

  vout.general(m_vl, "  Nk = %d  Np = %d\n", Nk, Np);
  vout.general(m_vl, "  Nm = %d\n", Nm);
  vout.general(m_vl, "  size of TDa = %d\n", TDa.size());
  vout.general(m_vl, "  size of vk  = %d\n", vk.size());
  vout.general(m_vl, "  Enorm_eigen  = %e\n", Enorm_eigen);
  vout.general(m_vl, "  Enorm_ignore_factor  = %e\n", Enorm_ignore_factor);

  Nconv = -1;
  Nsbt  =  0;

  int k1 = 1;
  int k2 = Nk;
  //  int kredundant = (20<(Np/2))?20:(Np/2);
  //  int k2 = Nk + kredundant;
  vout.general(m_vl, "  k2 = %d\n", k2);
  int Kdis = 0;
  int Kthreshold = 0;
  bool all_converged=false;

  int nconv;

#pragma omp parallel
 {

  //- set initial vector
  real_t bnorm2 = b.norm2();
  if(bnorm2 > 1.e-12){ // b is regarded as a nonzero vector
    copy(vk[0], b);
  }else{
    vk[0].set(1.0);    // start with a uniform vector
  }

  //- normalizing the initial vector
  real_t vnorm = vk[0].norm2();
  vnorm = 1.0/sqrt(vnorm);
  scal(vk[0], vnorm);

  //- initial Nk steps
  for (int k = 0; k < k2; ++k) {
    step(Nm, k, TDa, m_TDb, vk, m_f);
  }

#pragma omp barrier

  //- restarting loop begins
  for (int iter = 0; iter < Niter_eigen; ++iter) {

    vout.detailed(m_vl, "\n iter=%d\n", iter);

    for (int k = k2; k < Nm; ++k) {
      step(Nm, k, TDa, m_TDb, vk, m_f);
    }

    scal(m_f, m_TDb[Nm - 1]);

#pragma omp master
  {
    for (int k = 0; k < Nm; ++k) {
      m_TDa2[k] =   TDa[k + k1 - 1];
      m_TDb2[k] = m_TDb[k + k1 - 1];
    }

    //- getting eigenvalues
    setUnit_Qt(Nm, m_Qt);
    tqri(m_TDa2, m_TDb2, Nm, Nm, m_Qt, nconv);
  }
#pragma omp barrier

    if(nconv == -1){
      vout.crucial(m_vl, "Error at %s: QR method NOT converged.\n",
                   class_name.c_str());
      exit(EXIT_FAILURE);
    }else{
      vout.paranoiac(m_vl, "  tqri: converged at iter = %d\n", nconv);
    }

#pragma omp master
 {
    //- sorting
    m_sorter->sort(m_TDa2, Nm);

    //- implicitly shifted QR transformations
    setUnit_Qt(Nm, m_Qt);
    for (int ip = k2; ip < Nm; ++ip) {
      real_t Dsh  = m_TDa2[ip];
      int    kmin = k1;
      int    kmax = Nm;
      qrtrf(TDa, m_TDb, Nm, Nm, m_Qt, Dsh, kmin, kmax);
    }
 }
#pragma omp barrier

    for (int i = 0; i < (Nk + 1); ++i) {
      m_B[i].set(0.0);
    }

    for (int j = k1 - 1; j < k2 + 1; ++j) {
      for (int k = 0; k < Nm; ++k) {
        axpy(m_B[j], m_Qt[k + Nm * j], vk[k]);
      }
    }

    for (int j = k1 - 1; j < k2 + 1; ++j) {
      copy(vk[j], m_B[j]);
    }

    //- Compressed vector f and beta(k2)
    scal(m_f, m_Qt[Nm - 1 + Nm * (k2 - 1)]);

    axpy(m_f, m_TDb[k2 - 1], vk[k2]);

    real_t beta_k;
    beta_k = m_f.norm2();
    beta_k = sqrt(beta_k);

    vout.detailed(m_vl, " beta(k) = %20.14f\n", beta_k);

    real_t beta_r = 1.0 / beta_k;  
    copy(vk[k2], m_f);
    scal(vk[k2], beta_r);

    //- Convergence test

#pragma omp master
  {
    m_TDb[k2 - 1] = beta_k;

    m_TDa2 = TDa;
    m_TDb2 = m_TDb;

    setUnit_Qt(Nm, m_Qt);
    //    int nconv = -1;
    tqri(m_TDa2, m_TDb2, Nk, Nm, m_Qt, nconv);
  }
#pragma omp barrier

    if(nconv == -1){
      vout.crucial(m_vl, "Error at %s: QR method NOT converged.\n",
                   class_name.c_str());
      exit(EXIT_FAILURE);
    }else{
      vout.paranoiac(m_vl, "  tqri: converged at iter = %d\n", nconv);
    }

    for (int k = 0; k < Nk; ++k) {
      m_B[k].set(0.0);
    }

    for (int j = 0; j < Nk; ++j) {
      for (int k = 0; k < Nk; ++k) {
        axpy(m_B[j], m_Qt[k + j * Nm], vk[k]);
      }
    }

#pragma omp master
 {
    Kdis       = 0;
    Kthreshold = 0;
 }
    for (int i = 0; i < Nk; ++i) {

      m_fopr->mult(m_v, m_B[i]);
      real_t vnum = dot(m_B[i], m_v);
      real_t vden = dot(m_B[i], m_B[i]);
      // real_t vden = m_B[i].norm2();
      real_t vnorm = dot(m_v, m_v);
      //real_t vnorm = m_v.norm2();
      real_t TDa2_tmp = vnum / vden;
      axpy(m_v, -TDa2_tmp, m_B[i]);

     //real_t vv = m_v.norm2();
     real_t vv = dot(m_v, m_v);
     vout.detailed(m_vl, "  %4d  %18.14f  %18.14e\n", i, TDa2_tmp, vv);

#pragma omp master
     {
      m_TDa2[i] = TDa2_tmp;
      Enorm[i] = vv;
      if (vv < Enorm_eigen) {
        m_Iconv[Kdis] = i;
        ++Kdis;
        if (!m_sorter->comp(m_TDa2[i], Vthreshold)) {
          ++Kthreshold;
        }
      }
     } // omp master
#pragma omp barrier

    }  // i-loop end

    vout.detailed(m_vl, " #modes converged: %d\n", Kdis);
    if (Kthreshold > 0) {

#pragma omp master
     {
      // (there is a converged eigenvalue larger than Vthreshold.)
       std::vector<bool> converged(Nk, false);
       for (int i = 0; i < Kdis; ++i) {
         converged[m_Iconv[i]] = true;
       }
       std::vector<real_t> tmp_TDa(m_TDa2);
       std::vector<int> idx = m_sorter->sort_index(tmp_TDa, Nk);
       all_converged=true;
       for(int i = 0; i < m_Nk_conv; ++i){
         if(m_sorter->comp(m_TDa2[idx[i]], Vthreshold) && !converged[idx[i]]) {
           if ( Enorm[idx[i]] < Enorm_eigen * Enorm_ignore_factor ){
             all_converged=false;
           }
         }
       }
       if(all_converged){
         Nsbt = Kdis - Kthreshold;
       }
     }
#pragma omp barrier
    } else if(Kdis >= m_Nk_conv){
#pragma omp master
      { // check if all needed modes are converged
        std::vector<bool> converged(Nk, false);
        for (int i = 0; i < Kdis; ++i) {
          converged[m_Iconv[i]] = true;
        }
        std::vector<real_t> tmp_TDa(m_TDa2);
        std::vector<int> idx = m_sorter->sort_index(tmp_TDa, Nk);
        all_converged=false;
        int num_converged=0;
        for(int i=0; i<Nk; ++i){
          if(converged[idx[i]]){
            num_converged++;
          } else if ( Enorm[idx[i]] < Enorm_eigen * Enorm_ignore_factor ){
            break;
          }
          if(num_converged >= m_Nk_conv){
            all_converged=true;
            break;
          }
        }
        if(all_converged){ // no missing modes
          Nsbt = m_Nk_conv;
        } else {
          vout.general("%s: %d modes are convrged but seems some are still missing.\n",
                       class_name.c_str(), Kdis, m_Nk_conv);
        }
      } // omp master
#pragma omp barrier
    }
    if(all_converged){ // break
#pragma omp master
      {
        for (int i = 0; i < Kdis; ++i) {
          TDa[i] = m_TDa2[m_Iconv[i]];
        }
Nconv = iter;
      }
#pragma omp barrier
      iter = Niter_eigen;
    }
    fflush(stdout);
#pragma omp barrier
  } // end of iter loop
} // omp parallel

  std::vector<int> idx = m_sorter->sort_index(TDa, Kdis);
  for (int i = 0; i < Kdis; ++i) {
    copy(vk[i], m_B[m_Iconv[idx[i]]]);
    real_t vv = vk[i].norm2();
    vv = 1.0/sqrt(vv);
    scal(vk[i], vv);
    vv = vk[i].norm2();
    vout.paranoiac(m_vl, " vk[%d]:  norm2 = %22.14e\n", i, vv);
  }

  if (Nconv == -1) {
    vout.crucial(m_vl, "Error at %s: NOT converged.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }else{
    vout.general(m_vl, "\n Converged:\n");
    vout.general(m_vl, "  Nconv   = %d\n", Nconv);
    //  vout.general(m_vl, "  beta(k) = %20.14e\n", beta_k);
    vout.general(m_vl, "  Kdis    = %d\n", Kdis);
    vout.general(m_vl, "  Nsbt    = %d\n", Nsbt);
  }

  return;

}

//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRLanczos<FIELD, FOPR>::step(
                          int Nm, int k, std::vector<real_t>& TDa,
                          std::vector<real_t>& TDb, std::vector<FIELD>& vk,
                          FIELD& w)
{
  if (k >= Nm) {
    vout.crucial(m_vl, "Error at %s: k is larger than Nm.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);

  } else if (k == 0) {  // Initial step

    m_fopr->mult(w, vk[k]);

    real_t alph = dot(vk[k], w);

    axpy(w, -alph, vk[k]);

    real_t beta = dot(w, w);
    beta = sqrt(beta);
    real_t beta_r = 1.0 / beta;
    copy(vk[k + 1], w);
    scal(vk[k + 1], beta_r);

    //  vout.paranoiac(m_vl, "  step %d  alph = %e  beta = %e\n",
    //                k, alph, beta);
#pragma omp master
   {
    TDa[k] = alph;
    TDb[k] = beta;
   }

  } else {   // Iteration step

    m_fopr->mult(w, vk[k]);

    axpy(w, -TDb[k - 1], vk[k - 1]);

    real_t alph = dot(vk[k], w);

    axpy(w, -alph, vk[k]);

    real_t beta = dot(w, w);
    beta = sqrt(beta);
    real_t beta_r = 1.0 / beta;
    scal(w, beta_r);

    //  vout.paranoiac(m_vl, "  step %d  alph = %e  beta = %e\n",
    //                k, alph, beta);

#pragma omp master
   {
    TDa[k] = alph;
    TDb[k] = beta;
   }

    schmidt_orthogonalization(w, vk, k);

    if (k < Nm - 1) copy(vk[k + 1], w);
  }

}


//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRLanczos<FIELD, FOPR>::schmidt_orthogonalization(
                                             FIELD& w,
                                             std::vector<FIELD>& vk,
                                             int k)
{
  for (int j = 0; j < k; ++j) {
    // dcomplex prod = dotc(vk[j], w);
    complex_t prod = dotc(vk[j], w);
    prod *= cmplx(-1.0, 0.0);
    axpy(w, prod, vk[j]);

    // vout.paranoiac(m_vl, "Schmitd (%d, %d): %e  %e\n",
    //                j, k, real(prod), imag(prod) );
  }
}


//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRLanczos<FIELD, FOPR>::setUnit_Qt(
                                    int Nm, std::vector<real_t>& Qt)
{
  for (int i = 0; i < Qt.size(); ++i) {
    Qt[i] = 0.0;
  }

  for (int k = 0; k < Nm; ++k) {
    Qt[k + k * Nm] = 1.0;
  }
}


//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRLanczos<FIELD, FOPR>::tqri(
                 std::vector<real_t>& TDa, std::vector<real_t>& TDb,
                 int Nk, int Nm, std::vector<real_t>& Qt, int& nconv)
{
  int Niter = 100 * Nm;
  int kmin  = 1;
  int kmax  = Nk;
  // (these parameters should be tuned)

  nconv = -1;

  for (int iter = 0; iter < Niter; ++iter) {
    //- determination of 2x2 leading submatrix
    real_t dsub = TDa[kmax - 1] - TDa[kmax - 2];
    real_t dd   = sqrt(dsub * dsub + 4.0 * TDb[kmax - 2] * TDb[kmax - 2]);
    real_t Dsh  = 0.5 * (TDa[kmax - 2] + TDa[kmax - 1]
                         + fabs(dd) * (dsub / fabs(dsub)));
    // (Dsh: shift)

    //- transformation
    qrtrf(TDa, TDb, Nk, Nm, Qt, Dsh, kmin, kmax);

    //- Convergence criterion (redef of kmin and kmax)
    for (int j = kmax - 1; j >= kmin; --j) {
      real_t dds = fabs(TDa[j - 1]) + fabs(TDa[j]);
      if (fabs(TDb[j - 1]) + dds > dds) {
        kmax = j + 1;

        for (int j = 0; j < kmax - 1; ++j) {
          real_t dds = fabs(TDa[j]) + fabs(TDa[j + 1]);

          if (fabs(TDb[j]) + dds > dds) {
            kmin = j + 1;

            break;
          }
        }

        break;
      }

      if (j == kmin) {
        nconv = iter;
	//  vout.paranoiac(m_vl, "  tqri: converged at iter = %d\n", Nconv);
        return;
      }
    } // end of j loop
  }   // end of iter loop

  /*
  if (Nconv == -1) {
    vout.crucial(m_vl, "Error at %s: QR method NOT converged, Niter = %d.\n",
                 class_name.c_str(), Niter);
    exit(EXIT_FAILURE);
  }
  */

}

//====================================================================
template<typename FIELD, typename FOPR>
void AEigensolver_IRLanczos<FIELD, FOPR>::qrtrf(
                                  std::vector<real_t>& TDa,
                                  std::vector<real_t>& TDb,
                                  int Nk, int Nm, std::vector<real_t>& Qt,
                                  real_t Dsh, int kmin, int kmax)
{
  int    k = kmin - 1;
  real_t x;

  real_t Fden = 1.0 / sqrt((TDa[k] - Dsh) * (TDa[k] - Dsh)
                           + TDb[k] * TDb[k]);
  real_t c = (TDa[k] - Dsh) * Fden;
  real_t s = -TDb[k] * Fden;

  real_t tmpa1 = TDa[k];
  real_t tmpa2 = TDa[k + 1];
  real_t tmpb  = TDb[k];

  TDa[k]     = c * c * tmpa1 + s * s * tmpa2 - 2.0 * c * s * tmpb;
  TDa[k + 1] = s * s * tmpa1 + c * c * tmpa2 + 2.0 * c * s * tmpb;
  TDb[k]     = c * s * (tmpa1 - tmpa2) + (c * c - s * s) * tmpb;
  x          = -s * TDb[k + 1];
  TDb[k + 1] = c * TDb[k + 1];

  for (int i = 0; i < Nk; ++i) {
    real_t Qtmp1 = Qt[i + Nm * k];
    real_t Qtmp2 = Qt[i + Nm * (k + 1)];
    Qt[i + Nm * k]       = c * Qtmp1 - s * Qtmp2;
    Qt[i + Nm * (k + 1)] = s * Qtmp1 + c * Qtmp2;
  }

  //- Givens transformations
  for (int k = kmin; k < kmax - 1; ++k) {
    real_t Fden = 1.0 / sqrt(x * x + TDb[k - 1] * TDb[k - 1]);
    real_t c    = TDb[k - 1] * Fden;
    real_t s    = -x * Fden;

    real_t tmpa1 = TDa[k];
    real_t tmpa2 = TDa[k + 1];
    real_t tmpb  = TDb[k];
    TDa[k]     = c * c * tmpa1 + s * s * tmpa2 - 2.0 * c * s * tmpb;
    TDa[k + 1] = s * s * tmpa1 + c * c * tmpa2 + 2.0 * c * s * tmpb;
    TDb[k]     = c * s * (tmpa1 - tmpa2) + (c * c - s * s) * tmpb;
    TDb[k - 1] = c * TDb[k - 1] - s * x;
    if (k != kmax - 2) {
      x          = -s * TDb[k + 1];
      TDb[k + 1] = c * TDb[k + 1];
    }

    for (int i = 0; i < Nk; ++i) {
      real_t Qtmp1 = Qt[i + Nm * k];
      real_t Qtmp2 = Qt[i + Nm * (k + 1)];
      Qt[i + Nm * k]       = c * Qtmp1 - s * Qtmp2;
      Qt[i + Nm * (k + 1)] = s * Qtmp1 + c * Qtmp2;
    }
  }

}

//====================================================================
//============================================================END=====
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
