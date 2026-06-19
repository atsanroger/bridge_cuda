/*!
        @file    aeigensolver_IRLanczos.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate:: 2022-01-30 02:14:34 #$

        @version $LastChangedRevision: 2344 $
*/

#ifndef AEIGENSOLVER_IRLANCZOS_INCLUDED
#define AEIGENSOLVER_IRLANCZOS_INCLUDED

#include <vector>

#include "Eigen/aeigensolver.h"
#include "Tools/sorter_alt.h"
#include "bridge_complex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Eigenvalue solver with Implicitly Restarted Lanczos algorithm.

/*!
    This class determines eigenvalues and eigenvectors for a given
    fermion operator.
    Low- or high-lying eigenmodes are determined by chaning
    SortField class object.
                                 [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.         [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                 [21 Mar 2015 Y.Namekawa]
    Refactored as a template class.
                                 [29 May 2018 H.Matsufuru]
 */

template<typename FIELD, typename FOPR>
class AEigensolver_IRLanczos : public AEigensolver<FIELD,FOPR>
{
 public:
  typedef typename FIELD::real_t real_t;
  typedef std::complex<real_t> complex_t;

  static const std::string class_name;
  Bridge::VerboseLevel m_vl;

 private:
  int    m_Nk;
  int    m_Np;
  int    m_Nk_conv;
  int    m_Niter_eigen;
  real_t m_Enorm_eigen;
  real_t m_Vthreshold;
  real_t m_Enorm_ignore_factor;

  FOPR   *m_fopr;
  Sorter<real_t> *m_sorter;

  std::vector<real_t> m_TDb;
  std::vector<real_t> m_TDa2;
  std::vector<real_t> m_TDb2;
  std::vector<real_t> m_Qt;

  std::vector<int> m_Iconv;

  std::vector<FIELD> m_B;
  FIELD m_f;
  FIELD m_v;

 public:
  AEigensolver_IRLanczos(FOPR *fopr)
    : m_vl(CommonParameters::Vlevel()), m_fopr(fopr), m_sorter(0) {}

  AEigensolver_IRLanczos(FOPR *fopr, const Parameters& params)
    : m_vl(CommonParameters::Vlevel()), m_fopr(fopr), m_sorter(0)
    {
      set_parameters(params);
    }

  ~AEigensolver_IRLanczos();

  void set_parameters(const Parameters& params);
  void set_parameters(int Nk, int Np, int Nk_conv, int Niter_eigen, double Enorm_eigen,
                      double Vthreshold);
  void set_parameters(const std::string& sort_type,
                      int Nk, int Np, int Nk_conv, int Niter_eigen, double Enorm_eigen,
                      double Vthreshold);

  void solve(std::vector<real_t>& TDa, std::vector<FIELD>& vk,
             int& Nsbt, int& Nconv, const FIELD& b);

  void solve(std::vector<complex_t>& TDa, std::vector<FIELD>& vk,
             int& Nsbt, int& Nconv, const FIELD& b);

 private:

  void step(int Nm, int k,
            std::vector<real_t>& TDa,
            std::vector<real_t>& TDb,
            std::vector<FIELD>& vk, FIELD& f);

  void qrtrf(std::vector<real_t>& TDa, std::vector<real_t>& TDb,
             int Nk, int Nm, std::vector<real_t>& Qt,
             real_t Dsh, int kmin, int kmax);

  void tqri(std::vector<real_t>& TDa, std::vector<real_t>& TDb,
            int Nk, int Nm, std::vector<real_t>& Qt, int& nconv);

  void setUnit_Qt(int Nm, std::vector<real_t>& Qt);

  void schmidt_orthogonalization(FIELD& w, std::vector<FIELD>& vk, int k);

#ifdef USE_FACTORY
 private:
  static AEigensolver<FIELD,FOPR> *create_object(FOPR *fopr)
    { return new AEigensolver_IRLanczos<FIELD,FOPR>(fopr); }

  static AEigensolver<FIELD,FOPR> *create_object_with_params(
                           FOPR *fopr, const Parameters& params)
    { return new AEigensolver_IRLanczos<FIELD,FOPR>(fopr, params); }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= AEigensolver<FIELD,FOPR>::Factory_fopr::Register(
                               "IRLanczos", create_object);
    init &= AEigensolver<FIELD,FOPR>::Factory_fopr_params::Register(
                               "IRLanczos", create_object_with_params);
    return init;
  }
#endif

};
#endif
