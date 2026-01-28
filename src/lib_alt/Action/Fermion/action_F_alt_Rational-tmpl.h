/*!
      @file    action_F_alt_Rational-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-10-17 13:51:32 #$
      @version $LastChangedRevision: 2550 $
*/

#include "lib/Tools/timer.h"
#include "lib/Field/index_eo.h"

template<typename AFIELD>
const std::string Action_F_alt_Rational<AFIELD>::class_name
                                           = "Action_F_alt_Rational";

//====================================================================
template<typename AFIELD>
void Action_F_alt_Rational<AFIELD>::init()
{
  m_vl = CommonParameters::Vlevel();

  const int NinF     = m_fopr_langev->field_nin();
  const int NvolF    = m_fopr_langev->field_nvol();
  const int NexF     = m_fopr_langev->field_nex();

  m_psf.reset(NinF, NvolF, NexF);

}

//====================================================================
template<typename AFIELD>
void Action_F_alt_Rational<AFIELD>::tidyup()
{
  // do nothing.
}

//====================================================================
template<typename AFIELD>
void Action_F_alt_Rational<AFIELD>::set_parameters(
                                            const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);
}

//====================================================================
template<typename AFIELD>
void Action_F_alt_Rational<AFIELD>::get_parameters(
                                            Parameters& params) const
{
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
template<typename AFIELD>
void Action_F_alt_Rational<AFIELD>::set_config(Field *U)
{
  m_U = U;

  m_fopr_langev->set_config(U);

  m_fopr_H->set_config(U);

  m_forceF_MD->set_config(U);
}

//====================================================================
template<typename AFIELD>
double Action_F_alt_Rational<AFIELD>::langevin(RandomNumbers *rand)
{
  const int NinF     = m_fopr_langev->field_nin();
  const int NvolF    = m_fopr_langev->field_nvol();
  const int NexF     = m_fopr_langev->field_nex();
  const int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  Timer timer;
  timer.start();

  vout.general(m_vl, "%s: %s initial langevin step\n",
               class_name.c_str(), m_label.c_str());
  vout.increase_indent();

  const int Nvol = CommonParameters::Nvol();

  Field xi(NinF, NvolF, NexF);

  AFIELD xiA(NinF, NvolF, NexF);

  if(NvolF == Nvol){

    vout.detailed(m_vl, "  Gaussian filed with lexical action.\n");
    rand->gauss_lex_global(xi);

  }else if(NvolF == Nvol/2){

    vout.detailed(m_vl, "Gaussian filed with even-odd action.\n");

    Field xi_lex(NinF, Nvol, NexF);
    rand->gauss_lex_global(xi_lex);

    Index_eo idx_eo;
    idx_eo.convertField(xi, xi_lex, 0);

  }else{

     vout.crucial(m_vl, "Unsupported field volume.\n");
     exit(EXIT_FAILURE);

  }    

  double H_psf;

  #pragma omp parallel
  {
    if(m_fopr_langev->needs_convert()){
      m_fopr_langev->convert(xiA, xi);
    }else{
      AIndex_lex<real_t,AFIELD::IMPL> index_alt;
      convert(index_alt, xiA, xi);
    }

    m_fopr_langev->set_config(m_U);

    m_fopr_langev->mult(m_psf, xiA);

    double xi2 = xi.norm2();

    int ith = ThreadManager::get_thread_id();
    if(ith == 0) H_psf = xi2;
  }

  timer.stop();
  double elapsed_time = timer.elapsed_sec();

  vout.general(m_vl, "H_F     = %18.8f\n", H_psf);
  vout.general(m_vl, "H_F/dof = %18.8f\n", H_psf / size_psf);
  vout.general(m_vl, "elapsed time = %18.8f\n", elapsed_time);
  vout.decrease_indent();

  return H_psf;

}

//====================================================================
template<typename AFIELD>
double Action_F_alt_Rational<AFIELD>::calcH()
{
  const int NinF     = m_fopr_H->field_nin();
  const int NvolF    = m_fopr_H->field_nvol();
  const int NexF     = m_fopr_H->field_nex();
  const int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  vout.general(m_vl, "%s: %s\n", class_name.c_str(), m_label.c_str());
  vout.increase_indent();

  Timer timer;
  timer.start();

  AFIELD v1(NinF, NvolF, NexF);

  double H_psf;

#pragma omp parallel
 {
  m_fopr_H->set_config(m_U);
  m_fopr_H->mult(v1, m_psf);

  double H_psf1 = dot(v1, m_psf);

  int ith = ThreadManager::get_thread_id();
  if(ith == 0) H_psf = H_psf1;
 }

  timer.stop();
  double elapsed_time = timer.elapsed_sec();

  vout.general(m_vl, "H_F     = %18.8f\n", H_psf);
  vout.general(m_vl, "H_F/dof = %18.8f\n", H_psf / size_psf);
  vout.general(m_vl, "elapsed time = %18.8f\n", elapsed_time);
  vout.decrease_indent();

  return H_psf;

}

//====================================================================
template<typename AFIELD>
void Action_F_alt_Rational<AFIELD>::force(Field& force)
{
  vout.general(m_vl, "%s: %s\n", class_name.c_str(), m_label.c_str());
  vout.increase_indent();

  Timer timer;
  timer.start();

  const int Nin  = m_U->nin();
  const int Nvol = m_U->nvol();
  const int Nex  = m_U->nex();

  assert(force.nin() == Nin);
  assert(force.nvol() == Nvol);
  assert(force.nex() == Nex);

  AFIELD force1(Nin, Nvol, Nex);

#pragma omp parallel
 {
  m_forceF_MD->set_config(m_U);

  m_forceF_MD->force_core(force1, m_psf);

  AIndex_lex<real_t,AFIELD::IMPL> index_lex;
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
