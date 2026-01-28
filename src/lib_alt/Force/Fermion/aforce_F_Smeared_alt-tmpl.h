/*!
        @file    aforce_F_Smeared_alt-tmpl.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-07-03 23:32:09 #$
        @version $LastChangedRevision: 2529 $
*/

#include "lib_alt/Force/Fermion/aforce_F_Smeared_alt.h"

template<typename AFIELD>
const std::string AForce_F_Smeared_alt<AFIELD>::class_name
                                    = "AForce_F_Smeared_alt<AFIELD>";

//====================================================================
template<typename AFIELD>
void AForce_F_Smeared_alt<AFIELD>::init()
{
  int Nc = CommonParameters::Nc();
  int NinG = 2 * Nc * Nc;
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  m_force1.reset(NinG, Nvol, Ndim);   // AFIELD

  m_forceG1.reset(Nvol, Ndim);        // Field_G
  m_forceG2.reset(Nvol, Ndim);        // Field_G

}

//====================================================================
template<typename AFIELD>
void AForce_F_Smeared_alt<AFIELD>::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

}

//====================================================================
template<typename AFIELD>
void AForce_F_Smeared_alt<AFIELD>::set_config(Field *U)
{
  int ith = ThreadManager::get_thread_id();
  int nth = ThreadManager::get_num_threads();
  int nth_max = ThreadManager::get_num_threads_available();

  if(ith == 0) m_U = (Field_G *)U;
#pragma omp barrier

  if (nth > 1 || nth_max == 1) {
    set_config_impl(U);
  } else {
    set_config_omp(U);
  }

  m_director_smear->set_config(U);

  m_force->set_config(m_director_smear->get_config());

}

//====================================================================
template<typename AFIELD>
void AForce_F_Smeared_alt<AFIELD>::set_config_omp(Field *U)
{
#pragma omp parallel
  { set_config_impl(U); }
}

//====================================================================
template<typename AFIELD>
void AForce_F_Smeared_alt<AFIELD>::set_config_impl(Field *U)
{
  AIndex_lex<real_t,AFIELD::IMPL> index_lex;
  convert_gauge(index_lex, m_Ucp, *U);
}

//====================================================================
template<typename AFIELD>
void AForce_F_Smeared_alt<AFIELD>::force_udiv(AFIELD& force_,
                                              const AFIELD& eta)
{
  int Nsmear = m_director_smear->get_Nsmear();

  if (Nsmear == 0) {

    m_force->force_udiv(force_, eta);

  } else {

    AIndex_lex<real_t,AFIELD::IMPL> index_lex;

    Field_G* Uptr = m_director_smear->get_config();

    m_force->set_config(Uptr);

    m_force->force_udiv(m_force1, eta);

    reverse_gauge(index_lex, m_forceG2, m_force1);

    mult_jacobian(m_forceG2);

    convert_gauge(index_lex, force_, m_forceG2);

  }

}

//====================================================================
template<typename AFIELD>
void AForce_F_Smeared_alt<AFIELD>::force_udiv1(AFIELD& force_,
                                               const AFIELD& zeta,
                                               const AFIELD& eta)
{
  int Nsmear = m_director_smear->get_Nsmear();

  if (Nsmear == 0) {

    m_force->force_udiv1(force_, zeta, eta);

  } else {

    AIndex_lex<real_t,AFIELD::IMPL> index_lex;

    Field_G *Uptr = m_director_smear->get_config();

    m_force->set_config(Uptr);

    m_force->force_udiv1(m_force1, zeta, eta);

    reverse_gauge(index_lex, m_forceG2, m_force1);

    mult_jacobian(m_forceG2);

    convert_gauge(index_lex, force_, m_forceG2);

  }

}

//====================================================================
template<typename AFIELD>
void AForce_F_Smeared_alt<AFIELD>::mult_jacobian(Field_G& force)
{  // this function is not alt-coded.

  const int Nsmear = m_director_smear->get_Nsmear();

  copy(m_forceG1, force);

  for(int ismear = Nsmear - 1; ismear >= 0; --ismear) {

    Field *Uptr = m_director_smear->get_config(ismear);
    m_director_smear->force_udiv(force, m_forceG1, *Uptr);

    if (ismear > 0) copy(m_forceG1, force);

  }

}

//============================================================END=====
