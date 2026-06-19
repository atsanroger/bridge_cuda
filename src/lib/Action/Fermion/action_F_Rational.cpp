/*!
        @file    action_F_Rational.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-10-17 13:51:32 #$

        @version $LastChangedRevision: 2550 $
*/

#include "action_F_Rational.h"

const std::string Action_F_Rational::class_name = "Action_F_Rational";

//====================================================================
void Action_F_Rational::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }
}


//====================================================================
void Action_F_Rational::get_parameters(Parameters& params) const
{
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
double Action_F_Rational::langevin(RandomNumbers *rand)
{
  const int NinF     = m_fopr_langev->field_nin();
  const int NvolF    = m_fopr_langev->field_nvol();
  const int NexF     = m_fopr_langev->field_nex();
  const int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  m_psf.reset(NinF, NvolF, NexF);

  vout.general(m_vl, "%s: %s\n", class_name.c_str(), m_label.c_str());

  Field xi(NinF, NvolF, NexF);
  rand->gauss_lex_global(xi);

  m_fopr_langev->set_config(m_U);
  m_fopr_langev->mult(m_psf, xi);

  const double xi2   = xi.norm();
  const double H_psf = xi2 * xi2;

  vout.general(m_vl, "  H_Frational  = %18.8f\n", H_psf);
  vout.general(m_vl, "  H_F/dof      = %18.8f\n", H_psf / size_psf);

  return H_psf;
}


//====================================================================
double Action_F_Rational::calcH()
{
  const int NinF     = m_fopr_H->field_nin();
  const int NvolF    = m_fopr_H->field_nvol();
  const int NexF     = m_fopr_H->field_nex();
  const int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  vout.general(m_vl, "%s: %s\n", class_name.c_str(), m_label.c_str());

  Field v1(NinF, NvolF, NexF);
  m_fopr_H->set_config(m_U);
  m_fopr_H->mult(v1, m_psf);

  const double H_psf = dot(v1, m_psf);

  vout.general(m_vl, "  H_Frational  = %18.8f\n", H_psf);
  vout.general(m_vl, "  H_F/dof      = %18.8f\n", H_psf / size_psf);

  return H_psf;
}


//====================================================================
void Action_F_Rational::force(Field& force)
{
  const int Nin  = m_U->nin();
  const int Nvol = m_U->nvol();
  const int Nex  = m_U->nex();

  assert(force.nin() == Nin);
  assert(force.nvol() == Nvol);
  assert(force.nex() == Nex);

  vout.general(m_vl, "%s: %s\n", class_name.c_str(), m_label.c_str());

  m_fopr_force_MD->set_config(m_U);
  m_fopr_force_MD->force_core(force, m_psf);

  double Fave, Fmax, Fdev;
  force.stat(Fave, Fmax, Fdev);
  vout.general(m_vl, "  Fave = %12.6f  Fmax = %12.6f  Fdev = %12.6f\n",
               Fave, Fmax, Fdev);
}


//============================================================END=====
