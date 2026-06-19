/*!
      @file    aprecond_Mixedprec-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#include <cassert>

#include "lib_alt/Solver/aprecond_Mixedprec.h"


template<typename AFIELD, typename AFIELD2>
const std::string APrecond_Mixedprec<AFIELD, AFIELD2>::class_name
  = "APrecond_Mixedprec";
//====================================================================
template<typename AFIELD, typename AFIELD2>
void APrecond_Mixedprec<AFIELD, AFIELD2>::init()
{
  ThreadManager::assert_single_thread(class_name);

  int nin  = m_solver->get_fopr()->field_nin();
  int nvol = m_solver->get_fopr()->field_nvol();
  int nex  = m_solver->get_fopr()->field_nex();

  m_v2.reset(nin, nvol, nex);
  m_w2.reset(nin, nvol, nex);

  int Nvol = CommonParameters::Nvol();
  if (nvol == Nvol) {
    m_field_type = LEXICAL;
  } else if (nvol == Nvol / 2) {
    m_field_type = EVEN_ODD;
  } else {
    vout.crucial("%s: unsupported field_type.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  m_accum_flop = 0.0;
}


//====================================================================
template<typename AFIELD, typename AFIELD2>
void APrecond_Mixedprec<AFIELD, AFIELD2>::tidyup()
{
  // ThreadManager::assert_single_thread(class_name);
  //  delete m_solver;
}


//====================================================================
template<typename AFIELD, typename AFIELD2>
void APrecond_Mixedprec<AFIELD, AFIELD2>::mult(AFIELD& v,
                                               const AFIELD& w)
{
  int nin  = v.nin();
  int nvol = v.nvol();
  int nex  = v.nex();
  assert(w.check_size(nin, nvol, nex));
  assert(m_v2.check_size(nin, nvol, nex));
  assert(m_w2.check_size(nin, nvol, nex));

  int     nconv = -1;
  real_t2 diff;

#pragma omp barrier

  if (m_field_type == LEXICAL) {
    AIndex_lex<double, AFIELD::IMPL> index_d;
    AIndex_lex<float, AFIELD2::IMPL> index_f;
    convert(index_f, m_w2, index_d, w);
  } else if (m_field_type == EVEN_ODD) {
    AIndex_eo<double, AFIELD::IMPL> index_d;
    AIndex_eo<float, AFIELD2::IMPL> index_f;
    convert_h(index_f, m_w2, index_d, w);
  } else {
    vout.crucial("%s: unsupported field_type.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

#pragma omp barrier

  m_solver->solve(m_v2, m_w2, nconv, diff);
#pragma omp barrier

#pragma omp master
  {
    m_accum_flop += m_solver->flop_count();
  }

  vout.detailed("%s: nconv = %d  diff = %e\n",
                class_name.c_str(), nconv, diff);

  if (m_field_type == LEXICAL) {
    AIndex_lex<double, AFIELD::IMPL> index_d;
    AIndex_lex<float, AFIELD2::IMPL> index_f;
    convert(index_d, v, index_f, m_v2);
  } else if (m_field_type == EVEN_ODD) {
    AIndex_eo<double, AFIELD::IMPL> index_d;
    AIndex_eo<float, AFIELD2::IMPL> index_f;
    convert_h(index_d, v, index_f, m_v2);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD, typename AFIELD2>
void APrecond_Mixedprec<AFIELD, AFIELD2>::reset_flop_count()
{
#pragma omp master
  {
    m_accum_flop = 0.0;
  }
}


//====================================================================
template<typename AFIELD, typename AFIELD2>
double APrecond_Mixedprec<AFIELD, AFIELD2>::flop_count()
{
  return m_accum_flop;
}


//============================================================END=====
