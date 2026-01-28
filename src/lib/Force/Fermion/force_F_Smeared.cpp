/*!
        @file    force_F_Smeared.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-04-17 11:32:37 #$

        @version $LastChangedRevision: 2513 $
*/

#include "force_F_Smeared.h"

const std::string Force_F_Smeared::class_name = "Force_F_Smeared";

//====================================================================
void Force_F_Smeared::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }
}


//====================================================================
void Force_F_Smeared::get_parameters(Parameters& params) const
{
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Force_F_Smeared::force_udiv(Field& force_, const Field& eta)
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
void Force_F_Smeared::force_udiv1(Field& force_, const Field& zeta, const Field& eta)
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
void Force_F_Smeared::mult_jacobian(Field_G& force)
{
  const int Nsmear = m_director_smear->get_Nsmear();

  for (int i_smear = Nsmear - 1; i_smear >= 0; --i_smear) {
    Field_G *Uptr = m_director_smear->get_config(i_smear);

    Field_G f_tmp(force);  // copy to temporal field.
    m_force_smear->force_udiv(force, f_tmp, *Uptr);

    if (i_smear > 0) copy(f_tmp, force);  // ftmp = force;
  }
}


//============================================================END=====
