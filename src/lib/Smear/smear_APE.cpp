/*!
        @file    smear_APE.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "smear_APE.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Smear_APE::register_factory();
}
#endif

const std::string Smear_APE::class_name = "Smear_APE";

//====================================================================
void Smear_APE::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  double rho1;

  int err = 0;
  err += params.fetch_double("rho_uniform", rho1);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(rho1);
}


//====================================================================
void Smear_APE::get_parameters(Parameters& params) const
{
  std::vector<double> rho(m_rho.size());
  // std::copy(rho.begin(), m_rho.begin(), m_rho.end());
  for (size_t i = 0; i < m_rho.size(); ++i) {
    rho[i] = m_rho[i];
  }
  params.set_double_vector("rho", rho);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Smear_APE::set_parameters(const double rho1)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  rho = %8.4f\n", rho1);

  //- range check
  // NB. rho == 0 is allowed.

  //- store values
  // m_rho.resize(m_Ndim * m_Ndim);  // already resized in init.
  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      m_rho[mu + nu * m_Ndim] = rho1;
    }
  }
}


//====================================================================
void Smear_APE::set_parameters(const std::vector<double>& rho)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  rho[%d] = %8.4f\n", mu, rho[mu]);
  }

  // range check
  // NB. rho == 0 is allowed.
  assert(rho.size() == m_Ndim * m_Ndim);

  //- store values
  // m_rho.resize(m_Ndim * m_Ndim);  // already resized in init.
  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      m_rho[mu + nu * m_Ndim] = rho[mu + nu * m_Ndim];
    }
  }
}


//====================================================================
void Smear_APE::smear(Field_G& Usmear, const Field_G& U)
{
  const int Nvol = CommonParameters::Nvol();

  assert(U.nvol() == Nvol);
  assert(U.nex() == m_Ndim);
  assert(Usmear.nvol() == Nvol);
  assert(Usmear.nex() == m_Ndim);

  Usmear.set(0.0);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    Field_G c_tmp;
    c_tmp.set(0.0);

    Field_G u_tmp;
    u_tmp.setpart_ex(0, U, mu);

    Field_G u_tmp2;

    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu != mu) {
        Staple_lex staple;

        double rho = m_rho[mu + m_Ndim * nu];
        staple.upper(u_tmp2, U, mu, nu);
        c_tmp.addpart_ex(0, u_tmp2, 0, rho);

        staple.lower(u_tmp2, U, mu, nu);
        c_tmp.addpart_ex(0, u_tmp2, 0, rho);
      }
    }

    double rho0 = m_rho[mu + m_Ndim * mu];
    m_proj->project(u_tmp2, rho0, c_tmp, u_tmp);
    Usmear.setpart_ex(mu, u_tmp2, 0);
  }
}


//====================================================================
//============================================================END=====
