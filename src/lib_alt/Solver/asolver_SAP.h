/*!
        @file    asolver_SAP.h
        @brief   SAP solver (Alt-version)
        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate: 2023-03-20 10:52:44 +0900 (Mon, 20 Mar 2023) $
        @version $LastChangedRevision: 2499 $
*/

#ifndef ASOLVER_SAP_H
#define ASOLVER_SAP_H

/*!
  Multiplicative SAP solver with MINRES.
 */

#include <cstdio>
#include <cstdlib>
#include <vector>
using std::vector;
#include <string>
using std::string;

#include "lib_alt/Solver/asolver.h"
#include "lib_alt/Fopr/afopr_dd.h"
#include "asolver_SAP_MINRES.h"

template<typename AFIELD>
class ASolver_SAP : public ASolver<AFIELD>
{
 public:
  typedef typename AFIELD::real_t                  real_t;
  typedef AIndex_block_lex<real_t, AFIELD::IMPL>   block_index_t;
  using ASolver<AFIELD>::m_vl;
  static const std::string class_name;

 protected:

  const block_index_t *m_block_index;
  AFopr_dd<AFIELD> *m_fopr;   // given from outside
  unique_ptr<ASolver_SAP_MINRES<AFIELD> > m_sap_minres;

  int m_Niter;           //!< maximum iteration number.
  real_t m_Stop_cond;    //!< stopping criterion (squared).
  int m_Nconv;           //!< iteratoin number to calculate flop

  const int m_min_res_iter = 6;
  //const int m_min_res_iter = 20;

  //! to remember convergence iteration to provide flop count.
  int m_nconv;

  //! calling constructor without fermion operator is forbidden.
  ASolver_SAP() { }

  //! working vectors.
  AFIELD m_x, m_r, m_p;


 public:
  //! constructor.
  //  ASolver_SAP(AFopr<AFIELD>* fopr, const block_index_t *block_index)
  ASolver_SAP(AFopr_dd<AFIELD> *fopr, const block_index_t *block_index)
    : m_Niter(0), m_Stop_cond(0.0L)
  {
    m_fopr        = fopr;
    m_block_index = block_index;
    this->init();
  }

  //! destructor.
  ~ASolver_SAP() { this->tidyup(); }

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

#endif // ASOLVER_SAP_H
//============================================================END=====
