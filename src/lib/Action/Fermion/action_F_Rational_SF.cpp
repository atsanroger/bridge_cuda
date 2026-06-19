/*!
        @file    action_F_Rational_SF.cpp

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "action_F_Rational_SF.h"

const std::string Action_F_Rational_SF::class_name = "Action_F_Rational_SF";

//====================================================================
void Action_F_Rational_SF::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }
}


//====================================================================
void Action_F_Rational_SF::get_parameters(Parameters& params) const
{
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
double Action_F_Rational_SF::langevin(RandomNumbers *rand)
{
  const int NinF     = m_fopr_langev->field_nin();
  const int NvolF    = m_fopr_langev->field_nvol();
  const int NexF     = m_fopr_langev->field_nex();
  const int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  m_psf.reset(NinF, NvolF, NexF);

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  Field xi(NinF, NvolF, NexF);
  rand->gauss_lex_global(xi);

  m_fopr_langev->set_config(m_U);
  m_fopr_langev->mult(m_psf, xi);

  Field_SF::set_boundary_zero(xi);

  const double xi2   = xi.norm();
  const double H_psf = xi2 * xi2;

  vout.general(m_vl, "    H_Frational  = %18.8f\n", H_psf);
  vout.general(m_vl, "    H_F/dof      = %18.8f\n", H_psf / size_psf);

  return H_psf;
}


//====================================================================
double Action_F_Rational_SF::calcH()
{
  const int NinF     = m_fopr_H->field_nin();
  const int NvolF    = m_fopr_H->field_nvol();
  const int NexF     = m_fopr_H->field_nex();
  const int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  Field v1(NinF, NvolF, NexF);
  m_fopr_H->set_config(m_U);
  m_fopr_H->mult(v1, m_psf);

  const double H_psf = dot(v1, m_psf);

  vout.general(m_vl, "    H_Frational  = %18.8f\n", H_psf);
  vout.general(m_vl, "    H_F/dof      = %18.8f\n", H_psf / size_psf);

  return H_psf;
}


//====================================================================
void Action_F_Rational_SF::force(Field& force)
{
  const int Nvol = m_U->nvol();
  const int Nex  = m_U->nex();

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  Field_G force1(Nvol, Nex);
  m_fopr_force_MD->set_config(m_U);
  m_fopr_force_MD->force_core(force1, m_psf);

  Field_G Fboundary(force1);
  Field_SF::set_boundary_spatial_link_zero(Fboundary);

  force = (Field)Fboundary;

  double Fave, Fmax, Fdev;
  force.stat(Fave, Fmax, Fdev);
  vout.general(m_vl, "    Frational_ave = %12.6f  Frational_max = %12.6f  Fratio\
nal_dev = %12.6f\n",
               Fave, Fmax, Fdev);
}


//====================================================================
//============================================================END=====
