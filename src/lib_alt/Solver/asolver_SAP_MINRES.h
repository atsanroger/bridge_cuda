/*!
        @file    afopr_SAP_MINRES.h
        @brief   MinRes solver inside a SAP solver (Alt-version)
        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate::  $
        @version $LastChangedRevision: 2499 $
 */

#ifndef ASOLVER_SAP_MINRES_H
#define ASOLVER_SAP_MINRES_H

#include <cstdio>
#include <cstdlib>
#include <vector>
using std::vector;
#include <string>
using std::string;

#include "lib_alt/Solver/asolver.h"
#include "lib_alt/Fopr/afopr_dd.h"
#include "lib_alt/Field/aindex_block_lex_base.h"


template<typename AFIELD>
class ASolver_SAP_MINRES : public ASolver<AFIELD>
{
 public:
  typedef typename AFIELD::real_t                  real_t;
  typedef AIndex_block_lex<real_t, AFIELD::IMPL>   block_index_t;
  using ASolver<AFIELD>::m_vl;
  static const std::string class_name;

 protected:

  const block_index_t *m_block_index;

  AFopr_dd<AFIELD> *m_fopr; //!< need mult_sap()

  int m_Niter;              //!< maximum iteration number.
  real_t m_Stop_cond;       //!< stopping criterion (squared).

  //! to remember convergence iteration to provide flop count.
  int m_nconv;

  //! calling constructor without fermion operator is forbidden.
  ASolver_SAP_MINRES() { }

  //! working vectors.
  AFIELD m_r, m_p;

  using complex_t = typename AFIELD::complex_t;

  //  aligned_vector<real_t> m_r2_block;
  //  aligned_vector<real_t> m_p2_block;
  //  aligned_vector<complex_t> m_alpha_block;
  std::vector<real_t> m_r2_block;
  std::vector<real_t> m_p2_block;
  std::vector<complex_t> m_alpha_block;


 public:
  //! constructor.
  ASolver_SAP_MINRES(AFopr_dd<AFIELD> *fopr, const block_index_t *block_index)
    : m_Niter(0), m_Stop_cond(-1.0)
  {
    m_fopr        = fopr;
    m_block_index = block_index;
    this->init();
  }

  ASolver_SAP_MINRES(AFopr<AFIELD> *fopr)
    : m_Niter(0), m_Stop_cond(-1.0)
  {
    vout.crucial(m_vl, "partial construct is not allowed for Asolver_SAP_MINRES\n");
    abort();
    m_fopr        = nullptr;
    m_block_index = nullptr;
  }

  //! destructor.
  ~ASolver_SAP_MINRES() { this->tidyup(); }

  //! setting block.
  void set_block(block_index_t *block_index)
  {
    m_block_index = block_index;
    this->init();
  }

  //! setting parameters by a Parameter object.
  void set_parameters(const Parameters& params);

  //! setting parameters.
  void set_parameters(const int Niter, const real_t Stop_cond);

  //! solver main.
  void solve(AFIELD& xq, const AFIELD& b, int& nconv, real_t& diff, const int eo);

  //! returns the pointer to the fermion operator.
  AFopr<AFIELD> *get_afopr() { return m_fopr; }

  //! returns the floating point operation count.
  double flop_count();


 protected:

  void init(void);

  void tidyup(void);

 private:

  double m_flop_each;
  void calc_flop_each();

#ifdef USE_FACTORY
 private:
  static ASolver<AFIELD> *create_object(AFopr<AFIELD> *afopr)
  {
    return new ASolver_SAP_MINRES<AFIELD>(afopr);
  }

 public:
  static bool register_factory()
  {
    return ASolver<AFIELD>::Factory_fopr::Register("SAP_MINRES", create_object);
  }
#endif
};

#endif // ASOLVER_SAP_MINRES_H

//============================================================END=====
