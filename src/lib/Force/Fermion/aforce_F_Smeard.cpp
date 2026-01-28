/*!
      @file    aforce_F_Smeared.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-04-17 11:32:37 #$
      @version $LastChangedRevision: 2513 $
*/

#ifdef USE_ALT_CODE

#include "lib_alt/Force/Fermion/aforce_F_Smeared.h"

#include "lib/ResourceManager/threadManager.h"
#include "lib/Fopr/afopr_Smeared.h"

#include "lib/Field/field.h"


template<>
const std::string AForce_F_Smeared<Field>::class_name
                                      = "AForce_F_Smeared<Field>";

//====================================================================
template<>
void AForce_F_Smeared<Field>::init()
{

}

//====================================================================
template<>
void AForce_F_Smeared<Field>::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);
}

//====================================================================
template<>
void AForce_F_Smeared<Field>::set_config(Field *U)
{
  m_U = (Field_G *)U;
  m_director_smear->set_config(U);
  m_force->set_config(m_director_smear->get_config());
}

//====================================================================
template<>
void AForce_F_Smeared<Field>::mult_jacobian(Field_G& force)
{
  const int Nsmear = m_director_smear->get_Nsmear();

  for (int i_smear = Nsmear - 1; i_smear >= 0; --i_smear) {
    Field_G *Uptr = m_director_smear->get_config(i_smear);

    Field_G f_tmp(force);  // copy to temporal field.
    // m_director_smear->force_udiv(force, f_tmp, *Uptr);
    m_force_smear->force_udiv(force, f_tmp, *Uptr);

    if (i_smear > 0) copy(f_tmp, force);  // ftmp = force;
  }
}

//====================================================================
template<>
void AForce_F_Smeared<Field>::force_udiv(Field& force_, const Field& eta)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  const int Nsmear = m_director_smear->get_Nsmear();

  Field_G force(Nvol, Ndim);

  if (Nsmear == 0) {
    m_force->force_udiv(force, eta);
  } else {
    Field_G *Uptr = m_director_smear->get_config();

    m_force->set_config(Uptr);
    m_force->force_udiv(force, eta);

    mult_jacobian(force);
  }

  copy(force_, force); // force_ = force;
}


//====================================================================
template<>
void AForce_F_Smeared<Field>::force_udiv1(Field& force_, const Field& zeta, const Field& eta)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  const int Nsmear = m_director_smear->get_Nsmear();

  Field_G force(Nvol, Ndim);

  if (Nsmear == 0) {
    m_force->force_udiv1(force, zeta, eta);
  } else {
    Field_G *Uptr = m_director_smear->get_config();

    m_force->set_config(Uptr);
    m_force->force_udiv1(force, zeta, eta);

    mult_jacobian(force);
  }

  copy(force_, force); // force_ = force;
}


//====================================================================

//#include "lib_alt/Force/Fermion/aforce_F_Smeared-tmpl.h"


template class AForce_F_Smeared<Field>;

#endif
//============================================================END=====
