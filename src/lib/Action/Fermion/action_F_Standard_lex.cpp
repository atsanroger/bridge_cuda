/*!
        @file    action_F_Standard_lex.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-10-17 13:51:32 #$

        @version $LastChangedRevision: 2550 $
*/

#include "action_F_Standard_lex.h"

const std::string Action_F_Standard_lex::class_name = "Action_F_Standard_lex";

//====================================================================
void Action_F_Standard_lex::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }
}


//====================================================================
void Action_F_Standard_lex::get_parameters(Parameters& params) const
{
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Action_F_Standard_lex::set_config(Field *U)
{
  m_U = U;

  m_fopr->set_config(U);
  m_fopr_force->set_config(U);
}


//====================================================================
double Action_F_Standard_lex::langevin(RandomNumbers *rand)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  const int NinF     = m_fopr->field_nin();
  const int NvolF    = m_fopr->field_nvol();
  const int NexF     = m_fopr->field_nex();
  const int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  assert(NvolF == Nvol);
  m_psf.reset(NinF, NvolF, NexF);

  vout.general(m_vl, "%s: %s\n", class_name.c_str(), m_label.c_str());

  Field xi(NinF, NvolF, NexF);
  rand->gauss_lex_global(xi);

  m_fopr->set_config(m_U);
  m_fopr->set_mode("Ddag");

  m_fopr->mult(m_psf, xi);

  const double xi2   = xi.norm();
  const double H_psf = xi2 * xi2;

  vout.general(m_vl, "  H_Fstandard  = %18.8f\n", H_psf);
  vout.general(m_vl, "  H_F/dof      = %18.8f\n", H_psf / size_psf);

  return H_psf;
}


//====================================================================
double Action_F_Standard_lex::calcH()
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  const int NinF     = m_fopr->field_nin();
  const int NvolF    = m_fopr->field_nvol();
  const int NexF     = m_fopr->field_nex();
  const int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  vout.general(m_vl, "%s: %s\n", class_name.c_str(), m_label.c_str());

  Field  v1(NinF, NvolF, NexF);
  int    Nconv;
  double diff;
  m_fprop_H->set_config(m_U);
  m_fprop_H->invert_DdagD(v1, m_psf, Nconv, diff);

  vout.general(m_vl, "  Fprop_H: %6d %18.15e\n", Nconv, diff);

  const double H_psf = dot(v1, m_psf);

  vout.general(m_vl, "  H_Fstandard  = %18.8f\n", H_psf);
  vout.general(m_vl, "  H_F/dof      = %18.8f\n", H_psf / size_psf);

  return H_psf;
}


//====================================================================
void Action_F_Standard_lex::force(Field& force)
{
  const int Nin  = m_U->nin();
  const int Nvol = m_U->nvol();
  const int Nex  = m_U->nex();
  const int Nc   = CommonParameters::Nc();
  const int Ndim = CommonParameters::Ndim();

  assert(force.nin() == Nin);
  assert(force.nvol() == Nvol);
  assert(force.nex() == Nex);

  const int NinF  = m_fopr->field_nin();
  const int NvolF = m_fopr->field_nvol();
  const int NexF  = m_fopr->field_nex();

  vout.general(m_vl, "%s: %s\n", class_name.c_str(), m_label.c_str());

  Field  eta(NinF, NvolF, NexF);
  int    Nconv;
  double diff;
  m_fprop_MD->set_config(m_U);
  m_fprop_MD->invert_DdagD(eta, m_psf, Nconv, diff);

  vout.detailed(m_vl, "  Fprop_MD: %6d %18.15e\n", Nconv, diff);

  m_fopr_force->set_config(m_U);

  m_fopr_force->force_core(force, eta);

  double Fave, Fmax, Fdev;
  force.stat(Fave, Fmax, Fdev);
  vout.general(m_vl, "  Fave = %12.6f  Fmax = %12.6f  Fdev = %12.6f\n",
               Fave, Fmax, Fdev);
}


//============================================================END=====
