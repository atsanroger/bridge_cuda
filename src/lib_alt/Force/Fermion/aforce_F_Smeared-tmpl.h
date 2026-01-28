/*!
        @file    aforce_F_Smeared-tmpl.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-04-17 11:32:37 #$
        @version $LastChangedRevision: 2513 $
*/

//#include "lib_alt/Force/Fermion/aforce_F_Smeared.h"

template<typename AFIELD>
const std::string AForce_F_Smeared<AFIELD>::class_name
                                      = "AForce_F_Smeared<AFIELD>";

//====================================================================
template<typename AFIELD>
void AForce_F_Smeared<AFIELD>::init()
{
  //  set_parameters(params);
}

//====================================================================
template<typename AFIELD>
void AForce_F_Smeared<AFIELD>::set_config(Field *U)
{
  m_U = (Field_G *)U;

  AIndex_lex<real_t,AFIELD::IMPL> index_lex;

#pragma omp parallel
  { convert_gauge(index_lex, m_Ucp, *U); }

  m_director_smear->set_config(U);

  m_force->set_config(m_director_smear->get_config());

}

//====================================================================
template<typename AFIELD>
void AForce_F_Smeared<AFIELD>::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);
}


//====================================================================
template<typename AFIELD>
void AForce_F_Smeared<AFIELD>::force_udiv(AFIELD& force_,
                                          const AFIELD& eta)
{
  int Nc = CommonParameters::Nc();
  int NinG = 2 * Nc * Nc;
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  int Nsmear = m_director_smear->get_Nsmear();

  AFIELD force1(NinG, Nvol, Ndim);
  Field_G force2(Nvol, Ndim);

  if (Nsmear == 0) {

    m_force->force_udiv(force_, eta);

  } else {

    AIndex_lex<real_t,AFIELD::IMPL> index_lex;

    Field_G *Uptr = m_director_smear->get_config();

    m_force->set_config(Uptr);

    m_force->force_udiv(force1, eta);

    #pragma omp parallel
    { reverse_gauge(index_lex, force2, force1); }

    mult_jacobian(force2);

    #pragma omp parallel
    { convert_gauge(index_lex, force_, force2); }

  }

  //  copy(force_, force1);

}

//====================================================================
template<typename AFIELD>
void AForce_F_Smeared<AFIELD>::force_udiv1(AFIELD& force_,
                               const AFIELD& zeta, const AFIELD& eta)
{
  int Nc = CommonParameters::Nc();
  int NinG = 2 * Nc * Nc;
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  int Nsmear = m_director_smear->get_Nsmear();

  AFIELD  force1(NinG, Nvol, Ndim);
  Field_G force2(Nvol, Ndim);

  if (Nsmear == 0) {

    m_force->force_udiv1(force_, zeta, eta);

  } else {

    AIndex_lex<real_t,AFIELD::IMPL> index_lex;

    Field_G *Uptr = m_director_smear->get_config();

    m_force->set_config(Uptr);

    m_force->force_udiv1(force1, zeta, eta);

    #pragma omp parallel
    { reverse_gauge(index_lex, force2, force1); }

    mult_jacobian(force2);

    #pragma omp parallel
    { convert_gauge(index_lex, force_, force2); }

  }

  //  copy(force_, force); // force_ = force;

}

//====================================================================
template<typename AFIELD>
void AForce_F_Smeared<AFIELD>::mult_jacobian(Field_G& force)
{  // this function is not alt-coded.

  const int Nsmear = m_director_smear->get_Nsmear();

  Field_G f_tmp(force);  // copy to temporal field.

  for(int ismear = Nsmear - 1; ismear >= 0; --ismear) {

    Field *Uptr = m_director_smear->get_config(ismear);

    m_force_smear->force_udiv(force, f_tmp, *Uptr);

    if (ismear > 0) copy(f_tmp, force);  // ftmp = force;
  }
}

//============================================================END=====
