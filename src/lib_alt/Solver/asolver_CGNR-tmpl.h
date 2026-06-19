/*!
        @file    solver_CGNR.cpp
        @brief
        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate: 2013-04-27 12:28:50 #$
        @version $LastChangedRevision: 2499 $
*/

#include "lib_alt/Solver/asolver_CGNR.h"
using Bridge::vout;

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = ASolver_CGNR::register_factory();
}
#endif

template<typename AFIELD>
const std::string ASolver_CGNR<AFIELD>::class_name = "ASolver_CGNR";

//====================================================================
template<typename AFIELD>
void ASolver_CGNR<AFIELD>::init()
{
  int nin  = m_fopr->field_nin();
  int nvol = m_fopr->field_nvol();
  int nex  = m_fopr->field_nex();

  m_b2.reset(nin, nvol, nex);
}


//====================================================================
template<typename AFIELD>
void ASolver_CGNR<AFIELD>::set_parameters(const Parameters& params)
{
  std::string str_vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);

  vout.general(m_vl, "%s:\n", class_name.c_str());

  this->ASolver_CG<AFIELD>::set_parameters(params);
}


//====================================================================
template<typename AFIELD>
void ASolver_CGNR<AFIELD>::solve(AFIELD& xq, const AFIELD& b,
                                 int& Nconv, real_t& diff)
{
  int ith = ThreadManager::get_thread_id();

  const std::string mode_prev = m_fopr->get_mode();

  vout.detailed(m_vl, "%s: solver starts\n", class_name.c_str());

  if ((mode_prev == "DdagD") || (mode_prev == "DdagD_prec")) {
    vout.detailed(m_vl, "%s: calling CG, mode=%s\n",
                  class_name.c_str(), mode_prev.c_str());
    this->ASolver_CG<AFIELD>::solve(xq, b, Nconv, diff);
    if (ith == 0) {
      m_mode_fopr = mode_prev;
      m_flop      = this->ASolver_CG<AFIELD>::flop_count();
    }

    return;
  }

#pragma omp barrier

  if ((mode_prev == "D") || (mode_prev == "Ddag") ||
      (mode_prev == "D_prec") || (mode_prev == "Ddag_prec")) {
    vout.detailed(m_vl, "  mode_prev = %s\n", mode_prev.c_str());
  } else {
    vout.crucial(m_vl, "Error at %s: unsupported mode for fopr %s.",
                 class_name.c_str(), mode_prev.c_str());
    exit(EXIT_FAILURE);
  }

#pragma omp barrier

  if (ith == 0) {
    if ((m_b2.nin() != b.nin()) || (m_b2.nvol() != b.nvol()) ||
        (m_b2.nex() != b.nex())) {
      vout.detailed(m_vl, "  b2 size is reset.\n");
      m_b2.reset(b.nin(), b.nvol(), b.nex());
    }
  }

#pragma omp barrier

  m_fopr->mult_dag(m_b2, b);

  if (ith == 0) {
    m_flop = m_fopr->flop_count();
  }

#pragma omp barrier

  if (ith == 0) {
    if (mode_prev == "D") {
      m_mode_fopr = "DdagD";
    } else if (mode_prev == "Ddag") {
      m_mode_fopr = "DDdag";
    } else if (mode_prev == "D_prec") {
      m_mode_fopr = "DdagD_prec";
    } else if (mode_prev == "Ddag_prec") {
      m_mode_fopr = "DDdag_prec";
    } else {
      vout.crucial(m_vl, "%s: unknown mode for fopr: %s\n",
                   class_name.c_str(), mode_prev.c_str());
    }
  }

  m_fopr->set_mode(m_mode_fopr);

  this->ASolver_CG<AFIELD>::solve(xq, m_b2, Nconv, diff);
  if (ith == 0) {
    m_flop += this->ASolver_CG<AFIELD>::flop_count();
  }

  m_fopr->set_mode(mode_prev);

  // #pragma omp barrier  // set_mode has a barrier
}


//====================================================================
template<typename AFIELD>
double ASolver_CGNR<AFIELD>::flop_count()
{
  double flop_fopr = m_fopr->flop_count();
  if (flop_fopr < CommonParameters::epsilon_criterion()) {
    vout.crucial(m_vl, "Warning at %s: no fopr->flop_count() is"
                       "available, setting flop = 0\n", class_name.c_str());
    return 0.0;
  }

  return m_flop;

  /*
  int NPE = CommonParameters::NPE();

  std::string mode_keep = this->ASolver_CG<AFIELD>::get_fopr()->get_mode();
  this->ASolver_CG<AFIELD>::get_fopr()->set_mode(m_mode_fopr);

  double gflop_fopr = this->ASolver_CG<AFIELD>::get_fopr()->flop_count();
  if (gflop_fopr < CommonParameters::epsilon_criterion()) {
    vout.crucial(m_vl, "Warning at %s: no fopr->flop_count() is"
                 "available, setting flop = 0\n", class_name.c_str());
    return 0.0;
  }

  double gflop_solver = this->ASolver_CG<AFIELD>::flop_count();
  this->ASolver_CG<AFIELD>::get_fopr()->set_mode(mode_keep);

  double gflop = gflop_solver + gflop_fopr;

  return gflop;
  */
}


//====================================================================
//============================================================END=====
