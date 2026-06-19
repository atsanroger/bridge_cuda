/*!
        @file    action_F_Overlap_Nf2.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-09-08 06:55:16 #$

        @version $LastChangedRevision: 2314 $
*/

#include "action_F_Overlap_Nf2.h"

const std::string Action_F_Overlap_Nf2::class_name = "Action_F_Overlap_Nf2";

//====================================================================
void Action_F_Overlap_Nf2::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }
}


//====================================================================
void Action_F_Overlap_Nf2::get_parameters(Parameters& params) const
{
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));

  return;
}


//====================================================================
double Action_F_Overlap_Nf2::langevin(RandomNumbers *rand)
{
  const int Nc = CommonParameters::Nc();
  const int Nd = CommonParameters::Nd();

  const int Nvol = CommonParameters::Nvol();
  const int NPE  = CommonParameters::NPE();

  const int NinF  = m_fopr->field_nin();
  const int NvolF = m_fopr->field_nvol();
  const int NexF  = m_fopr->field_nex();

  assert(NvolF == Nvol);
  m_psf.reset(NinF, NvolF, NexF);

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  Field xi(NinF, NvolF, NexF);
  rand->gauss_lex_global(xi);

  m_fopr->set_config(m_U);
  // m_fopr->set_mode("Ddag");
  m_fopr->set_mode("H");

  m_fopr->mult(m_psf, xi);

  const double xi2   = xi.norm();
  const double H_psf = xi2 * xi2;

  vout.general(m_vl, "    H_Foverlap = %18.8f\n", H_psf);
  vout.general(m_vl, "    H_F/dof    = %18.8f\n", H_psf / Nvol / NPE / Nc / Nd);

  return H_psf;
}


//====================================================================
double Action_F_Overlap_Nf2::calcH()
{
  const int Nc = CommonParameters::Nc();
  const int Nd = CommonParameters::Nd();

  const int Nvol = CommonParameters::Nvol();
  const int NPE  = CommonParameters::NPE();

  const int NinF  = m_fopr->field_nin();
  const int NvolF = m_fopr->field_nvol();
  const int NexF  = m_fopr->field_nex();

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  Field  v1(NinF, NvolF, NexF);
  int    Nconv;
  double diff;
  m_fprop_H->set_config(m_U);
  m_fprop_H->invert_DdagD(v1, m_psf, Nconv, diff);

  vout.general(m_vl, "    Fprop_H: %6d %18.15e\n", Nconv, diff);

  const double H_psf = dot(v1, m_psf);

  vout.general(m_vl, "    H_F_overlap = %18.8f\n", H_psf);
  vout.general(m_vl, "    H_F/dof     = %18.8f\n", H_psf / Nvol / NPE / Nc / Nd);

  return H_psf;
}


//====================================================================
void Action_F_Overlap_Nf2::force(Field& force)
{
  const int NinF  = m_fopr->field_nin();
  const int NvolF = m_fopr->field_nvol();
  const int NexF  = m_fopr->field_nex();

  const int Nin  = m_U->nin();
  const int Nvol = m_U->nvol();
  const int Nex  = m_U->nex();

  assert(force.nin() == Nin);
  assert(force.nvol() == Nvol);
  assert(force.nex() == Nex);

  Field  psi(NinF, NvolF, NexF);
  int    Nconv;
  double diff;
  m_fprop_MD->set_config(m_U);
  m_fprop_MD->invert_DdagD(psi, m_psf, Nconv, diff);

  vout.general(m_vl, "  Solver(ov):  Nconv = %d  diff  = %.8e\n", Nconv, diff);

  m_fopr_force->set_config(m_U);
  m_fopr_force->force_core(force, psi);

  double Fave, Fmax, Fdev;
  force.stat(Fave, Fmax, Fdev);
  vout.general(m_vl,
               "    Foverlap_ave = %12.6f  Foverlap_max = %12.6f  Foverlap_dev = %12.6f\n",
               Fave, Fmax, Fdev);
}


//====================================================================
//===========================================================END======
