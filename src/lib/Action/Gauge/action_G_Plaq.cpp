/*!
        @file    action_G_Plaq.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-10-17 13:51:32 #$

        @version $LastChangedRevision: 2550 $
*/

#include "action_G_Plaq.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Action_G_Plaq::register_factory();
}
#endif

const std::string Action_G_Plaq::class_name = "Action_G_Plaq";

//====================================================================
void Action_G_Plaq::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  double beta;

  int err = 0;
  err += params.fetch_double("beta", beta);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(beta);

  //- post-process
  m_force_G->set_parameters(params);
}


//====================================================================
void Action_G_Plaq::get_parameters(Parameters& params) const
{
  params.set_double("beta", m_beta);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Action_G_Plaq::set_parameters(const double beta)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  beta = %8.4f\n", beta);

  //- range check
  // NB. beta == 0 is allowed.

  //- store values
  m_beta = beta;
}


//====================================================================
double Action_G_Plaq::langevin(RandomNumbers *rand)
{
  const double H_U = calcH(); // calculate action H_U=beta*(1-Plaq)*Lvol*6 (SA)

  return H_U;
}


//====================================================================
double Action_G_Plaq::calcH()
{
  const int Ndim  = CommonParameters::Ndim();
  const int Ndim2 = Ndim * (Ndim - 1) / 2;

  const int Nvol = CommonParameters::Nvol();
  const int NPE  = CommonParameters::NPE();

  const double plaq = m_staple.plaquette(*m_U);                          // calculate plaquette (SA)
  const double H_U  = m_beta * (1.0 - plaq) * Nvol * NPE * Ndim2;        // action (SA)

  vout.general(m_vl, "H_Gplaq    = %18.8f\n", H_U);                      // total action (SA)
  vout.general(m_vl, "H_G/dof    = %18.8f\n", H_U / Nvol / NPE / Ndim2); // action per dof (SA)

  return H_U;
}


//====================================================================
void Action_G_Plaq::force(Field& force)
{
  //- check of argument type
  assert(force.nin() == m_U->nin());
  assert(force.nvol() == m_U->nvol());
  assert(force.nex() == m_U->nex());

  vout.general(m_vl, "%s:  %s\n", class_name.c_str(), m_label.c_str());

  force.set(0.0);

  m_force_G->force_core(force, m_U);
}


//====================================================================
//============================================================END=====
