/*!
        @file    action_F_alt_Standard_lex.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-10-17 13:51:32 #$
        @version $LastChangedRevision: 2550 $
*/

#include "Action/Fermion/action_F_alt_Standard_lex.h"

#include "lib/Tools/timer.h"


template<typename AFIELD>
const std::string Action_F_alt_Standard_lex<AFIELD>::class_name
                                       = "Action_F_alt_Standard_lex";

//====================================================================
template<typename AFIELD>
void Action_F_alt_Standard_lex<AFIELD>::init()
{
  m_vl = CommonParameters::Vlevel();

}

//====================================================================
template<typename AFIELD>
void Action_F_alt_Standard_lex<AFIELD>::tidyup()
{
  // do nothing.
}

//====================================================================
template<typename AFIELD>
void Action_F_alt_Standard_lex<AFIELD>::set_parameters(
                                           const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);

}

//====================================================================
template<typename AFIELD>
void Action_F_alt_Standard_lex<AFIELD>::get_parameters(
                                         Parameters& params) const
{
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}

//====================================================================
template<typename AFIELD>
void Action_F_alt_Standard_lex<AFIELD>::set_config(Field *U)
{
  m_U = U;

  m_fopr->set_config(U);

  m_fopr_force->set_config(U);

}

//====================================================================
template<typename AFIELD>
double Action_F_alt_Standard_lex<AFIELD>::langevin(RandomNumbers *rand)
{
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  int NinF     = m_fopr->field_nin();
  int NvolF    = m_fopr->field_nvol();
  int NexF     = m_fopr->field_nex();
  int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  assert(NvolF == Nvol);
  m_psf.reset(NinF, NvolF, NexF);

  Timer timer;
  timer.start();

  vout.general(m_vl, "%s: %s\n", class_name.c_str(), m_label.c_str());
  vout.increase_indent();

  Field xi(NinF, NvolF, NexF);
  rand->gauss_lex_global(xi);

  AFIELD xiA(NinF, NvolF, NexF);

  double H_psf;

#pragma omp parallel
 {
  if(m_fopr->needs_convert()){
    m_fopr->convert(xiA, xi);
  }else{
    AIndex_lex<real_t, AFIELD::IMPL> index_alt;
    convert_spinor(index_alt, xiA, xi);
  }

  m_fopr->set_config(m_U);

  m_fopr->set_mode("Ddag");

  m_fopr->mult(m_psf, xiA);
  
  double H_psf1 = xi.norm2();

  int ith = ThreadManager::get_thread_id();
  if(ith == 0) H_psf = H_psf1;
 }

  timer.stop();
  double elapsed_time = timer.elapsed_sec();

  vout.general(m_vl, "H_F     = %18.8f\n", H_psf);
  vout.general(m_vl, "H_F/dof = %18.8f\n", H_psf / size_psf);
  vout.general(m_vl, "elapsed time: %12.6f [sec]\n", elapsed_time);
  vout.decrease_indent();

  return H_psf;

}

//====================================================================
template<typename AFIELD>
double Action_F_alt_Standard_lex<AFIELD>::calcH()
{
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  int NinF     = m_fopr->field_nin();
  int NvolF    = m_fopr->field_nvol();
  int NexF     = m_fopr->field_nex();
  int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  Timer timer;
  timer.start();

  vout.general(m_vl, "%s: %s\n", class_name.c_str(), m_label.c_str());
  vout.increase_indent();

  AFIELD v1(NinF, NvolF, NexF);
  int    Nconv;
  double diff;

  m_fprop_H->set_config(m_U);
  m_fprop_H->set_mode("DdagD");
  m_fprop_H->invert(v1, m_psf, Nconv, diff);

  vout.general(m_vl, "Fprop_H: %6d %18.15e\n", Nconv, diff);

  double H_psf;

#pragma omp parallel
 {
  double H_psf1 = dot(v1, m_psf);

  int ith = ThreadManager::get_thread_id();
  if(ith == 0) H_psf = H_psf1;
 }

  timer.stop();
  double elapsed_time = timer.elapsed_sec();

  vout.general(m_vl, "H_F     = %18.8f\n", H_psf);
  vout.general(m_vl, "H_F/dof = %18.8f\n", H_psf / size_psf);
  vout.general(m_vl, "elapsed time: %12.6f [sec]\n", elapsed_time);
  vout.decrease_indent();

  return H_psf;

}

//====================================================================
template<typename AFIELD>
void Action_F_alt_Standard_lex<AFIELD>::force(Field& force)
{
  const int Nin  = m_U->nin();
  const int Nvol = m_U->nvol();
  const int Nex  = m_U->nex();
  const int Nc   = CommonParameters::Nc();
  const int Ndim = CommonParameters::Ndim();

  assert(force.nin() == Nin);
  assert(force.nvol() == Nvol);
  assert(force.nex() == Nex);

  AFIELD force1(Nin, Nvol, Nex);

  const int   NinF  = m_fopr->field_nin();
  const int   NvolF = m_fopr->field_nvol();
  const int   NexF  = m_fopr->field_nex();

  Timer timer;
  timer.start();

  vout.general(m_vl, "%s: %s\n", class_name.c_str(), m_label.c_str());
  vout.increase_indent();

  AFIELD eta(NinF, NvolF, NexF);
  int    Nconv;
  double diff;

  m_fprop_MD->set_config(m_U);
  m_fprop_MD->set_mode("DdagD");
  m_fprop_MD->invert(eta, m_psf, Nconv, diff);

  vout.general(m_vl, "Fprop_MD: %6d %18.15e\n", Nconv, diff);

#pragma omp parallel
 {
  m_fopr->set_config(m_U);
  m_fopr_force->set_config(m_U);

  m_fopr_force->force_core(force1, eta);

  AIndex_lex<real_t, AFIELD::IMPL> index_lex;
  reverse_gauge(index_lex, force, force1);

  double Fave, Fmax, Fdev;
  force.stat(Fave, Fmax, Fdev);
  vout.general(m_vl, "Fave = %12.6f  Fmax = %12.6f  Fdev = %12.6f\n",
               Fave, Fmax, Fdev);
 }

 timer.stop();
  double elapsed_time = timer.elapsed_sec();
  vout.general(m_vl, "elapsed time: %12.6f [sec]\n", elapsed_time);
  vout.decrease_indent();

}

//============================================================END=====
