/*!
        @file    aeigensolver_IRArnoldi.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate:: 2022-05-07 23:30:01 #$

        @version $LastChangedRevision: 2373 $
*/

#ifndef AEIGENSOLVER_IRARNOLDI_INCLUDED
#define AEIGENSOLVER_IRARNOLDI_INCLUDED

#include <vector>

#include "Eigen/aeigensolver.h"
#include "Tools/sorter_alt.h"
#include "bridge_complex.h"
#include "complexTraits.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Eigenvalue solver with Implicitly Restarted Arnoldi algorithm.

/*!
    This template class determines eigenvalues and eigenvectors for
    a given non-hermitian fermion operator.
                                          [24 Feb 2020 H.Matsufuru]
 */

template<typename FIELD, typename FOPR>
  class AEigensolver_IRArnoldi : public AEigensolver<FIELD,FOPR>
{
 public:
  typedef typename FIELD::real_t real_t;
  typedef typename ComplexTraits<real_t>::complex_t complex_t;

  static const std::string class_name;
  Bridge::VerboseLevel m_vl;

 private:
  int    m_Nk;
  int    m_Np;
  int    m_Nk_conv;
  int    m_Nm;  //!< Nm = Nk + Np
  int    m_Niter_eigen;
  real_t m_Enorm_eigen;
  real_t m_Vthreshold;  //!< given as an absolute value
  real_t m_Enorm_ignore_factor;

  FOPR   *m_fopr;
  Sorter<complex_t> *m_sorter;

  std::vector<complex_t> m_TDa2;
  std::vector<complex_t> m_Qt;

  std::vector<complex_t> m_Ht;  //!< Hessenberg matrix
  std::vector<complex_t> m_Ht2; //!< temporary Hessenberg matrix
  std::vector<complex_t> m_Yt;

  std::vector<int> m_Iconv;

  std::vector<FIELD> m_B;
  FIELD m_f;
  FIELD m_v;

 public:
  AEigensolver_IRArnoldi(FOPR *fopr)
    : m_vl(CommonParameters::Vlevel()),
    m_fopr(fopr), m_sorter(0) {}

  AEigensolver_IRArnoldi(FOPR *fopr, const Parameters& params)
    : m_vl(CommonParameters::Vlevel()),
    m_fopr(fopr), m_sorter(0)
    {
      set_parameters(params);
    }

  ~AEigensolver_IRArnoldi();

  void set_parameters(const Parameters& params);
  void set_parameters(int Nk, int Np, int Nk_conv, int Niter_eigen,
                      double Enorm_eigen, double Vthreshold);
  void set_parameters(const std::string& sort_type,
                      int Nk, int Np, int Nk_conf, int Niter_eigen,
                      double Enorm_eigen, double Vthreshold);

  void solve(std::vector<complex_t>& TDa, std::vector<FIELD>& vk,
             int& Nsbt, int& Nconv, const FIELD& b);

 private:

  void step(int Nm, int k,
            std::vector<FIELD>& vk, FIELD& f);

  //! deflation of converged eigenvectors.
  void deflation(int k1, int k2, int Kdis,
                 std::vector<complex_t>& TDa,
                 std::vector<FIELD>& vk,
                 complex_t& beta);


  void qrtrf(std::vector<complex_t>& Ht,
             int Nk, int Nm, std::vector<complex_t>& Qt,
             complex_t Dsh, int kmin, int kmax);

  // void tqri(std::vector<complex_t>& Ht,
  //           int Nk, int Nm, std::vector<complex_t>& Qt, int& nconv);
  void tqri(std::vector<complex_t>& Ht, int k1,
            int Nk, int Nm, std::vector<complex_t>& Qt, int& nconv);

  void setUnit_Qt(int Nm, std::vector<complex_t>& Qt);

  void schmidt(int Nk, std::vector<FIELD>& vk);

  void shift_wilkinson(complex_t& kappa,
                       const complex_t a, const complex_t b,
                       const complex_t c, const complex_t d);

  void check_Qt(const int Nk, const int Nm,
                std::vector<complex_t>& Qt,
                std::vector<complex_t>& Ht,
                std::vector<complex_t>& At);

  void eigenvector_Ht(std::vector<complex_t>& Yt,
                      std::vector<complex_t>& St,
                      int km, int Nm);

  //! Yt = Qt * Xt.
  void mult_Qt(std::vector<complex_t>& Yt,
               std::vector<complex_t>& Qt,
               std::vector<complex_t>& Xt,
               int km, int Nm);

  void check_eigen_Ht(std::vector<complex_t>& Ht,
                      std::vector<complex_t>& TDa,
                      std::vector<complex_t>& Xt,
		      int km, int Nm);

  void reconst_Ht(std::vector<complex_t>& Ht,
                  std::vector<complex_t>& Qt,
                  complex_t& beta,
                  std::vector<FIELD>& vk, int Nk);

  void reunit_Qt(std::vector<complex_t>& Qt, int Nk);

  inline int index(int i, int j){ return j + i * m_Nm; }
  // note that i can be Nm. max: i = Nm, j = Nm-1, 
  //  j+i*Nm = Nm-1 + Nm*Nm = (Nm+1)*Nm-1, size = (Nm+1)*Nm.

#ifdef USE_FACTORY
 private:
  static AEigensolver<FIELD,FOPR> *create_object(FOPR *fopr)
    { return new AEigensolver_IRArnoldi<FIELD,FOPR>(fopr); }


  static AEigensolver<FIELD,FOPR> *create_object_with_params(
                            FOPR *fopr, const Parameters& params)
    { return new AEigensolver_IRArnoldi<FIELD,FOPR>(fopr, params); }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= AEigensolver<FIELD, FOPR>::Factory_fopr::Register(
                               "IRArnoldi", create_object);
    init &= AEigensolver<FIELD, FOPR>::Factory_fopr_params::Register(
                               "IRArnoldi", create_object_with_params);
    return init;
  }
#endif

};
#endif
