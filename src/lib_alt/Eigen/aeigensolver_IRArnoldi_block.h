/*!
        @file    aeigensolver_IRArnoldi_block.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2024-04-05 23:37:58 #$

        @version $LastChangedRevision: 2599 $
*/

#ifndef AEIGENSOLVER_IRARNOLDI_BLOCK_INCLUDED
#define AEIGENSOLVER_IRARNOLDI_BLOCK_INCLUDED

#include <vector>

#include "bridge_complex.h"
#include "complexTraits.h"

#include "lib/Eigen/aeigensolver.h"
#include "lib/Tools/sorter_alt.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib_alt/Field/aindex_block_lex_base.h"
#include "lib_alt/Fopr/afopr_dd.h"

//! Blocked version of the Implicitly Restarted Arnoldi algorithm.

/*!
    This template class implements the blocked version of the
    Implicitly Restarted Arnoldi algorithm of eigenvalue solver.
    The eigenproblem is solved in units of block for the
    domain-decomposed fermion operator.
                                          [27 Mar 2024 H.Matsufuru]
 */

template<typename FIELD, typename FOPR>
class AEigensolver_IRArnoldi_block : public AEigensolver<FIELD, FOPR>
{
 public:
  typedef typename FIELD::real_t real_t;
  typedef typename ComplexTraits<real_t>::complex_t complex_t;
  typedef AIndex_block_lex<real_t, FIELD::IMPL> AIndex_block;
  typedef AFopr_dd<FIELD> FOPR_DD;

  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;
  int m_Ndim;

  int m_Nk;
  int m_Np;
  int m_Nm;     //!< Nm = Nk + Np
  int m_Niter_eigen;
  real_t m_Enorm_eigen;
  real_t m_Vthreshold;   //!< given as an absolute value

  std::vector<int> m_fine_lattice;    //!< fine lattice size
  std::vector<int> m_coarse_lattice;  //!< coarse lattice size
  int m_num_blocks;                   //!< number of blocks
  AIndex_block *m_index_block;

  real_t m_ecrit;     //!< criteria of too small number

  std::string m_sort_type;

  FOPR_DD *m_fopr;
  Sorter<complex_t> *m_sorter;

  std::vector<std::vector<complex_t> > m_TDa2;
  std::vector<std::vector<complex_t> > m_HtB;  //!< Hessenberg matrix
  std::vector<std::vector<complex_t> > m_Ht2B;
  std::vector<std::vector<complex_t> > m_YtB;
  std::vector<std::vector<complex_t> > m_QtB;

  std::vector<std::vector<int> > m_Iconv;
  std::vector<int> m_Kdis;
  std::vector<int> m_Kthreshold;
  std::vector<int> m_Kalmost;

  std::vector<complex_t> m_bvec1;

  std::vector<real_t> m_rvec1, m_rvec2;

  std::vector<std::vector<int> > m_reorder;

  std::vector<FIELD> m_B;
  FIELD m_f;
  FIELD m_v;

 public:
  AEigensolver_IRArnoldi_block(FOPR *fopr, const Parameters& params):
    m_fopr((FOPR_DD*)fopr), m_sorter(0)
  { init(params);  }

  ~AEigensolver_IRArnoldi_block(){ tidyup(); }

  void set_parameters(const Parameters& params);

  void set_parameters(const std::string& sort_type,
                      int Nk, int Np, int Niter_eigen,
                      double Enorm_eigen, double Vthreshold,
                      std::vector<int> fine_lattice,
                      std::vector<int> coarse_lattice);

  void set_parameters(int Nk, int Np, int Niter_eigen,
                      double Enorm_eigen, double Vthreshold,
                      std::vector<int> fine_lattice,
                      std::vector<int> coarse_lattice);

  void get_parameters(Parameters& params) const;

  void solve(std::vector<complex_t>& TDa, std::vector<FIELD>& vk,
             int& Nsbt, int& Nconv, const FIELD& b);

 private:

  void init(const Parameters& params);

  void tidyup();

  void step(int Nm, int k,
            std::vector<FIELD>& vk, FIELD& f);

  void sort_eigenvectors(std::vector<complex_t>& TDa,
                         std::vector<FIELD>& vk);

  //! deflation of converged eigenvectors.
  void deflation(int k1, int k2, int Kdis,
                 std::vector<complex_t>& TDa,
                 std::vector<FIELD>& vk,
                 complex_t& beta);


  void qrtrf(std::vector<complex_t>& Ht,
             int Nk, int Nm, std::vector<complex_t>& Qt,
             complex_t Dsh, int kmin, int kmax);

  void qrtrf1(std::vector<complex_t>& Ht,
              int Nk, int Nm, std::vector<complex_t>& Qt,
              complex_t Dsh, int kmin, int kmax);

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

  void reconst_Ht(std::vector<std::vector<complex_t> >& Ht,
                  std::vector<complex_t>& beta,
                  std::vector<FIELD>& vk, int Nk);

  void reunit_Qt(std::vector<complex_t>& Qt, int Nk);

  inline int index(int i, int j) { return j + i * m_Nm; }
  // note that i can be Nm. max: i = Nm, j = Nm-1,
  //  j+i*Nm = Nm-1 + Nm*Nm = (Nm+1)*Nm-1, size = (Nm+1)*Nm.

  inline int index_block(int i, int j, int block)
  { return j + i * m_Nm + m_Nm * m_Nm * block; }

#ifdef USE_FACTORY
 private:
  static AEigensolver<FIELD, FOPR> *create_object_with_params(
                                     FOPR *fopr, const Parameters& params)
  { return new AEigensolver_IRArnoldi_block<FIELD, FOPR>(fopr, params); }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= AEigensolver<FIELD, FOPR>::Factory_fopr_params::Register(
                           "IRArnoldi_block", create_object_with_params);
    return init;
  }
#endif
};
#endif
