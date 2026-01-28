/*!
        @file    asolver_SAP_QWS.h
        @brief   SAP solver (qws version)
        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate: 2023-03-20 10:52:44 +0900 (Mon, 20 Mar 2023) $
        @version $LastChangedRevision: 2499 $
*/

#ifndef ASOLVER_SAP_QWS_H
#define ASOLVER_SAP_QWS_H

/*!
  Multiplicative SAP solver implemented in qws
 */

#include <cstdio>
#include <cstdlib>
#include <vector>
using std::vector;
#include <string>
using std::string;

#include "lib_alt/Solver/asolver.h"
#include "lib_alt_QXS/Fopr/afopr_Clover_QWS_dd.h"

#ifdef USE_QWSLIB
#include "lib_alt_QXS/extra/qws_lib.h"
#endif


template<typename AFIELD>
class ASolver_SAP_QWS : public ASolver<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  using ASolver<AFIELD>::m_vl;
  static const std::string class_name;

 protected:

  //  const block_index_t *m_block_index;
  AFopr_Clover_QWS_dd<AFIELD> *m_fopr; // given from outside

  int m_Niter;                         //!< maximum iteration number.
  real_t m_Stop_cond;                  //!< stopping criterion (squared).
  int m_Nconv;                         //!< iteratoin number to calculate flop
  double m_flop;                       //!< flop count
  const int m_nm = 10;                 //!< fixted iteration for the inner jacobi iteration

  //! to remember convergence iteration to provide flop count.
  int m_nconv;

  //! calling constructor without fermion operator is forbidden.
  ASolver_SAP_QWS() { }

  //! working vectors.
  AFIELD m_x, m_r, m_b;

#ifdef USE_QWSLIB
  scs_t *m_s, *m_q;
#endif

 public:
  //! constructor.
  ASolver_SAP_QWS(AFopr_dd<AFIELD> *fopr)
    : m_Niter(0), m_Stop_cond(0.0L)
  {
#ifndef USE_QWSLIB
    vout.crucial("%s: USE_QWSLIB is not defined\n", class_name.c_str());
    exit(EXIT_FAILURE);
#endif

    m_fopr = dynamic_cast<AFopr_Clover_QWS_dd<AFIELD> *>(fopr);
    if (m_fopr == nullptr) {
      vout.crucial(m_vl, "asolver_SAP_QWS: bad fopr is given, must be AFopr_Clover_QWS_dd\n");
      exit(EXIT_FAILURE);
    }
    this->init();
  }

  //! destructor.
  ~ASolver_SAP_QWS() { this->tidyup(); }

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

#endif // ASOLVER_SAP_QWS_H
//============================================================END=====
