/*!
        @file    force_G_Plaq.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "force_G_Plaq.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Force_G_Plaq::register_factory();
}
#endif

const std::string Force_G_Plaq::class_name = "Force_G_Plaq";

//====================================================================
void Force_G_Plaq::set_parameters(const Parameters& params)
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
}


//====================================================================
void Force_G_Plaq::get_parameters(Parameters& params) const
{
  params.set_double("beta", m_beta);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Force_G_Plaq::set_parameters(const double beta)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  beta = %8.4f\n", beta);

  //- range check
  // NB. kappa == 0 is allowed.

  //- store values
  m_beta = beta;
}


//====================================================================
void Force_G_Plaq::force_core(Field& force)
{
  const int Nc   = CommonParameters::Nc();
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  assert(m_U->nin() == Nc * Nc * 2);
  assert(m_U->nvol() == Nvol);
  assert(m_U->nex() == Ndim);

  assert(force.nin() == Nc * Nc * 2);
  assert(force.nvol() == Nvol);
  assert(force.nex() == Ndim);

  for (int mu = 0; mu < Ndim; ++mu) {
    Field_G st;
    m_staple.staple(st, *m_U, mu);
    // Calculate staple for m_U(all site,\mu) (SA)

    /* -->--
       |     |
       |     |
       U_mu
    */
    Field_G force1;

    for (int site = 0; site < Nvol; ++site) {
      Mat_SU_N ut(Nc);                            // Nc x Nc complex matrix (SA)
      ut = m_U->mat(site, mu) * st.mat_dag(site); // U_\mu * (staple)^\dagger (SA)
      ut.at();                                    // anti-hermitian and traceless (SA)
      //ut *= -beta/Nc;
      // force = -beta*{U_\mu *staple^\dagger}_{traceless & anti-hermitian) (SA)
      force1.set_mat(site, 0, ut);
    }

    axpy(force, mu, -(m_beta / Nc), force1, 0);
  }

  double Fave, Fmax, Fdev;
  force.stat(Fave, Fmax, Fdev);

  //- calculate average, max, deviation of force over (site, mu)
  vout.general(m_vl, "    Fave = %12.6f  Fmax = %12.6f  Fdev = %12.6f\n",
               Fave, Fmax, Fdev);
}


//====================================================================
//============================================================END=====
