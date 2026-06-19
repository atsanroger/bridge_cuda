/*!
        @file    action_F_alt_Standard_eo.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-07-03 23:32:09 #$
        @version $LastChangedRevision: 2529 $
*/

#include "Action/Fermion/action_F_alt_Standard_eo.h"

#include "lib/Tools/timer.h"
#include "lib/Field/index_eo.h"


template<typename AFIELD>
const std::string Action_F_alt_Standard_eo<AFIELD>::class_name
                                       = "Action_F_alt_Standard_eo";

//====================================================================
template<typename AFIELD>
void Action_F_alt_Standard_eo<AFIELD>::init()
{
  m_vl = CommonParameters::Vlevel();

}

//====================================================================
template<typename AFIELD>
void Action_F_alt_Standard_eo<AFIELD>::tidyup()
{
  // do nothing.
}

//====================================================================
template<typename AFIELD>
void Action_F_alt_Standard_eo<AFIELD>::set_parameters(
                                           const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);

}

//====================================================================
template<typename AFIELD>
void Action_F_alt_Standard_eo<AFIELD>::get_parameters(
                                         Parameters& params) const
{
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}

//====================================================================
template<typename AFIELD>
void Action_F_alt_Standard_eo<AFIELD>::set_config(Field *U)
{
  m_U = U;

  m_fopr->set_config(U);

  m_fopr_force->set_config(U);

}

//====================================================================
template<typename AFIELD>
double Action_F_alt_Standard_eo<AFIELD>::langevin(RandomNumbers *rand)
{
  int Nx   = CommonParameters::Nx();
  int Ny   = CommonParameters::Ny();
  int Nz   = CommonParameters::Nz();
  int Nt   = CommonParameters::Nt();
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();
  int Nx2 = Nx/2;

  int NinF     = m_fopr->field_nin();
  int NvolF    = m_fopr->field_nvol();
  int NexF     = m_fopr->field_nex();
  int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  Timer timer;
  timer.start();

  assert((2 * NvolF) == Nvol);
  m_psf.reset(NinF, NvolF, NexF);

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  Field xi_lex(NinF, Nvol, NexF);
  rand->gauss_lex_global(xi_lex);

  Field xi(NinF, NvolF, NexF);
  Index_eo idx_eo;

  AFIELD xiA(NinF, NvolF, NexF);

  double H_psf;
  int ith = ThreadManager::get_thread_id();

#pragma omp parallel
 {
  idx_eo.convertField(xi, xi_lex, 0);

  if(m_fopr->needs_convert()){
      m_fopr->convert(xiA, xi);
    }else{
      AIndex_lex<real_t, AFIELD::IMPL> index_lex(Nx2, Ny, Nz, Nt);
      convert_spinor(index_lex, xiA, xi);
    }

    m_fopr->set_config(m_U);

    m_fopr->set_mode("Ddag");

    m_fopr->mult(m_psf, xiA);

    double H_psf1 = xiA.norm2();
    if(ith == 0) H_psf = H_psf1;
  }

  timer.stop();
  double elapsed_time = timer.elapsed_sec();

  vout.general(m_vl, "    H_F      = %18.8f\n", H_psf);
  vout.general(m_vl, "    H_F/dof  = %18.8f\n", H_psf / size_psf);
  vout.general(m_vl, "    elapsed time: %12.6f [sec]\n", elapsed_time);

  return H_psf;

}

//====================================================================
template<typename AFIELD>
double Action_F_alt_Standard_eo<AFIELD>::calcH()
{
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  int NinF     = m_fopr->field_nin();
  int NvolF    = m_fopr->field_nvol();
  int NexF     = m_fopr->field_nex();
  int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  Timer timer;
  timer.start();

  AFIELD v1(NinF, NvolF, NexF);
  int    Nconv;
  double diff;

  m_fprop_H->set_config(m_U);
  m_fprop_H->set_mode("DdagD_even");
  m_fprop_H->invert(v1, m_psf, Nconv, diff);

  vout.general(m_vl, "    Fprop_H: %6d %18.15e\n", Nconv, diff);

  int ith = ThreadManager::get_thread_id();
  double H_psf;

#pragma omp parallel
 {
  double H_psf1 = dot(v1, m_psf); 
  if(ith == 0) H_psf = H_psf1;
 }

  timer.stop();
  double elapsed_time = timer.elapsed_sec();

  vout.general(m_vl, "    H_Fstandard  = %18.8f\n", H_psf);
  vout.general(m_vl, "    H_F/dof      = %18.8f\n", H_psf / size_psf);
  vout.general(m_vl, "    elapsed time: %12.6f [sec]\n", elapsed_time);

  return H_psf;

}

//====================================================================
template<typename AFIELD>
void Action_F_alt_Standard_eo<AFIELD>::force(Field& force)
{
  int Nin  = m_U->nin();
  int Nvol = m_U->nvol();
  int Nex  = m_U->nex();
  int Nc   = CommonParameters::Nc();
  int Ndim = CommonParameters::Ndim();

  assert(force.nin() == Nin);
  assert(force.nvol() == Nvol);
  assert(force.nex() == Nex);

  AFIELD force1(Nin, Nvol, Nex);
  AFIELD force2(Nin, Nvol, Nex);

  int NinF   = m_fopr->field_nin();
  int NvolF  = m_fopr->field_nvol();
  int NexF   = m_fopr->field_nex();

  Timer timer;
  timer.start();

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  AFIELD eta(NinF, NvolF, NexF);

  int    Nconv;
  double diff;

  m_fprop_MD->set_config(m_U);
  m_fprop_MD->set_mode("DdagD_even");
  m_fprop_MD->invert(eta, m_psf, Nconv, diff);

  vout.general(m_vl, "    Fprop_MD: %6d %18.15e\n", Nconv, diff);

#pragma omp parallel
{
  m_fopr->set_config(m_U);

  m_fopr_force->set_config(m_U);
  m_fopr_force->force_core(force1, eta);

  AIndex_lex<real_t, AFIELD::IMPL> index_lex;
  reverse_gauge(index_lex, force, force1);

  double Fave, Fmax, Fdev;
  force.stat(Fave, Fmax, Fdev);
  vout.general(m_vl, "    Fave = %12.6f  Fmax = %12.6f  Fdev = %12.6f\n",
               Fave, Fmax, Fdev);
 }

  timer.stop();
  double elapsed_time = timer.elapsed_sec();
  vout.general(m_vl, "    elapsed time: %12.6f [sec]\n", elapsed_time);

}

//============================================================END=====
