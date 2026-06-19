/*!
        @file    aforceSmear_APE-tmpl.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-07-03 23:32:09 #$

        @version $LastChangedRevision: 2529 $
*/

#include "lib_alt_QXS/Smear/aforceSmear_APE.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = ForceSmear_APE::register_factory();
}
#endif

template<typename AFIELD>
const std::string AForceSmear_APE<AFIELD>::class_name
                                      = "AForceSmear_APE<AFIELD>";

//====================================================================
template<typename AFIELD>
void AForceSmear_APE<AFIELD>::init()
{
  m_Ndim = CommonParameters::Ndim();
  m_Nvol = CommonParameters::Nvol();

  m_rho.resize(m_Ndim * m_Ndim);
  m_U.resize(m_Ndim);
  m_iTheta.resize(m_Ndim);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      m_rho[mu + nu * m_Ndim] = 0.0;
    }
  }

  int Nc = CommonParameters::Nc();
  int Nin = 2 * Nc * Nc;
  m_shift = new ShiftField_lex(Nin);

}


//====================================================================
template<typename AFIELD>
void AForceSmear_APE<AFIELD>::tidyup()
{
  delete m_shift;
}


//====================================================================
template<typename AFIELD>
void AForceSmear_APE<AFIELD>::set_parameters(const Parameters& params)
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
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(rho1);
}


//====================================================================
template<typename AFIELD>
void AForceSmear_APE<AFIELD>::get_parameters(Parameters& params) const
{
  params.set_double_vector("rho", m_rho);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
template<typename AFIELD>
void AForceSmear_APE<AFIELD>::set_parameters(const double rho1)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  rho = %8.4f\n", rho1);

  //- store values
  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      m_rho[mu + nu * m_Ndim] = rho1;
    }
  }
}


//====================================================================
template<typename AFIELD>
void AForceSmear_APE<AFIELD>::set_parameters(const std::vector<double>& rho)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  rho[%d] = %8.4f\n", mu, rho[mu]);
  }

  //- range check
  // NB. rho == 0 is allowed.
  assert(rho.size() == m_Ndim * m_Ndim);

  //- store values
  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      m_rho[mu + nu * m_Ndim] = rho[mu + nu * m_Ndim];
    }
  }
}


//====================================================================
template<typename AFIELD>
void AForceSmear_APE<AFIELD>::force_udiv(Field_G& Sigma,
                                         const Field_G& Sigmap,
                                         const Field_G& U)
{
  const int Nc = CommonParameters::Nc();

  assert(Sigmap.nin() == (2 * Nc * Nc));
  assert(Sigmap.nvol() == m_Nvol);
  assert(Sigmap.nex() == m_Ndim);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    copy(m_U[mu], 0, U, mu);
  }
#pragma omp barrier

  for (int mu = 0; mu < m_Ndim; ++mu) {
    m_C.set(0.0);
#pragma omp barrier

    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;
      double rho = m_rho[mu + m_Ndim * nu];

      staple(m_Ctmp, m_U[mu], m_U[nu], mu, nu);
      axpy(m_C, 0, rho, m_Ctmp, 0);
#pragma omp barrier
    }

    copy(m_sigmap_tmp, 0, Sigmap, mu);
#pragma omp barrier

    double alpha = m_rho[mu + m_Ndim * mu];

    m_proj->force_recursive(m_Xi, m_iTheta[mu],
                            alpha, m_sigmap_tmp, m_C, m_U[mu]);
    copy(Sigma, mu, m_Xi, 0);
#pragma omp barrier
  }

  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;
      double rho = m_rho[mu + m_Ndim * nu];

      force_each(m_sigma_tmp, m_U[mu], m_U[nu],
                 m_iTheta[mu], m_iTheta[nu], mu, nu);
      axpy(Sigma, mu, rho, m_sigma_tmp, 0);
#pragma omp barrier
    }
  }
}


//====================================================================
template<typename AFIELD>
void AForceSmear_APE<AFIELD>::force_each(Field_G& Sigma_mu,
                                         const Field_G& V_mu,
                                         const Field_G& V_nu,
                                         const Field_G& iTheta_mu,
                                         const Field_G& iTheta_nu,
                                         const int mu, const int nu)
{
#pragma omp barrier

  Sigma_mu.set(0.0);
#pragma omp barrier

  m_shift->backward(m_vt1, V_nu, mu);
  m_shift->backward(m_vt2, V_mu, nu);

  mult_Field_Gnd(m_vt3, 0, m_vt1, 0, m_vt2, 0);
#pragma omp barrier

  multadd_Field_Gnd(Sigma_mu, 0, m_vt3, 0, iTheta_nu, 0, 1.0);
#pragma omp barrier

  mult_Field_Gdn(m_vt3, 0, iTheta_mu, 0, V_nu, 0);
#pragma omp barrier

  mult_Field_Gdn(m_vt2, 0, m_vt1, 0, m_vt3, 0);
#pragma omp barrier

  m_shift->forward(m_vt3, m_vt2, nu);
  axpy(Sigma_mu, 1.0, m_vt3);
#pragma omp barrier

  mult_Field_Gdn(m_vt3, 0, V_mu, 0, iTheta_nu, 0);
#pragma omp barrier

  mult_Field_Gdn(m_vt2, 0, m_vt1, 0, m_vt3, 0);
#pragma omp barrier

  m_shift->forward(m_vt3, m_vt2, nu);
  axpy(Sigma_mu, 1.0, m_vt3);
#pragma omp barrier

  m_shift->backward(m_vt1, iTheta_nu, mu);
  m_shift->backward(m_vt2, V_mu, nu);

  mult_Field_Gnd(m_vt3, 0, m_vt1, 0, m_vt2, 0);
#pragma omp barrier

  multadd_Field_Gnd(Sigma_mu, 0, m_vt3, 0, V_nu, 0, 1.0);
#pragma omp barrier

  mult_Field_Gdd(m_vt2, 0, m_vt1, 0, V_mu, 0);
#pragma omp barrier

  mult_Field_Gnn(m_vt3, 0, m_vt2, 0, V_nu, 0);
#pragma omp barrier

  m_shift->forward(m_vt2, m_vt3, nu);

  axpy(Sigma_mu, 1.0, m_vt2);
#pragma omp barrier

  m_shift->backward(m_vt1, V_nu, mu);
  m_shift->backward(m_vt2, iTheta_mu, nu);

  mult_Field_Gnd(m_vt3, 0, m_vt1, 0, m_vt2, 0);
#pragma omp barrier

  multadd_Field_Gnd(Sigma_mu, 0, m_vt3, 0, V_nu, 0, 1.0);
#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AForceSmear_APE<AFIELD>::staple(Field_G& c,
                            const Field_G& u_mu, const Field_G& u_nu,
                            const int mu, const int nu)
{
#pragma omp barrier

  //- upper direction
  m_shift->backward(m_vt1, u_mu, nu);

  mult_Field_Gnn(m_vt2, 0, u_nu, 0, m_vt1, 0);
#pragma omp barrier

  m_shift->backward(m_vt1, u_nu, mu);

  mult_Field_Gnd(c, 0, m_vt2, 0, m_vt1, 0);
#pragma omp barrier

  //- lower direction
  m_shift->backward(m_vt2, u_nu, mu);

  mult_Field_Gnn(m_vt1, 0, u_mu, 0, m_vt2, 0);
#pragma omp barrier

  mult_Field_Gdn(m_vt2, 0, u_nu, 0, m_vt1, 0);
#pragma omp barrier

  m_shift->forward(m_vt1, m_vt2, nu);

  axpy(c, 0, real_t(1.0), m_vt1, 0);
#pragma omp barrier
}


//============================================================END=====
