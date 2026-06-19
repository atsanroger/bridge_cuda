/*!
        @file    asolver_SAP_dw2.h
        @brief   SAP solver (another trial version for Domainwall fermion)
        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate: 2024-06-01 20:47:37 +0900 (2024年06月01日 (土)) $
        @version $LastChangedRevision: 2604 $
*/

#ifndef ASOLVER_SAP_DW2_H
#define ASOLVER_SAP_DW2_H

/*!
  Multiplicative SAP solver
 */

#include <cstdio>
#include <cstdlib>
#include <vector>
using std::vector;
#include <string>
using std::string;

#include "lib_alt/Solver/asolver.h"
#include "lib_alt/Fopr/afopr_dd.h"
#include "lib_alt/alt_impl.h"

template<typename REAL_TYPE, Impl IMPL>
class AIndex_block_lex;

template<typename AFIELD>
class ASolver_SAP_dw2 : public ASolver<AFIELD>
{
 public:
  typedef typename AFIELD::real_t                  real_t;
  typedef AIndex_block_lex<real_t, AFIELD::IMPL>   block_index_t;
  using ASolver<AFIELD>::m_vl;
  static const std::string class_name;

 protected:

  const block_index_t *m_block_index;
  AFopr_dd<AFIELD> *m_fopr;   // given from outside
  ASolver<AFIELD> *m_solver_in_block;
  // PV operator
  ASolver<AFIELD> *m_solver_PV;
  AFopr_dd<AFIELD> *m_fopr_PV;

  int m_Niter;           //!< maximum iteration number.
  real_t m_Stop_cond;    //!< stopping criterion (squared).
  int m_Nconv;           //!< iteratoin number to calculate flop

  //  const int m_min_res_iter = 6;
  //const int m_min_res_iter = 20;

  //! to remember convergence iteration to provide flop count.
  int m_nconv;

  //! calling constructor without fermion operator is forbidden.
  ASolver_SAP_dw2() { }

  //! working vectors.
  AFIELD m_b, m_r, m_p, m_v, m_x;


 public:
  //! constructor.
  //  ASolver_SAP_dw(AFopr<AFIELD>* fopr, const block_index_t *block_index)
  ASolver_SAP_dw2(AFopr_dd<AFIELD> *fopr, const block_index_t *block_index,
                  ASolver<AFIELD> *solver_in_block,
                  ASolver<AFIELD> *solver_PV)
    : m_Niter(0), m_Stop_cond(0.0L)
  {
    m_fopr        = fopr;
    m_block_index = block_index;
    m_solver_in_block = solver_in_block;
    m_solver_PV = solver_PV;
    this->init();
  }

  //! destructor.
  ~ASolver_SAP_dw2() { this->tidyup(); }

  //! setting parameters by a Parameter object.
  void set_parameters(const Parameters& params);

  //! setting parameters.
  void set_parameters(const int Niter, const real_t Stop_cond);

  //! solver main.
  void solve(AFIELD& xq, const AFIELD& b, int& nconv, real_t& diff);

  //! returns the pointer to the fermion operator.
  AFopr<AFIELD> *get_fopr() { return m_fopr; }

  //! returns the floating point operation count.
  double flop_count();

 protected:

  void init(void);

  void tidyup(void);
};

#endif // ASOLVER_SAP_DW2_H
//============================================================END=====
