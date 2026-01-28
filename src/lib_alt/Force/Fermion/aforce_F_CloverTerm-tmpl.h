/*!
        @file    aforce_F_CloverTerm-tmpl.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-10-17 13:51:32 #$

        @version $LastChangedRevision: 2550 $
*/

#include "lib_alt/Force/Fermion/aforce_F_CloverTerm.h"

template<typename AFIELD>
const std::string AForce_F_CloverTerm<AFIELD>::class_name
                                 = "AForce_F_CloverTerm<AFIELD>";

//====================================================================
template<typename AFIELD>
void AForce_F_CloverTerm<AFIELD>::init(const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  } else {
    m_vl = CommonParameters::Vlevel();
  }

  vout.general(m_vl, "%s: construction\n", class_name.c_str());
  vout.increase_indent();

  m_fopr_csw = new AFopr_CloverTerm<AFIELD>(params);

  //m_forceF_ct = new Force_F_CloverTerm(params);

  int Nc = CommonParameters::Nc();
  int Nd = CommonParameters::Nd();
  int NinF = Nc * Nd * 2;
  int NinG = Nc * Nc * 2;
  int Nvol = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();

  m_shift  = new ShiftAField_lex<AFIELD>(NinF);

  m_staple = new AStaple_lex<AFIELD>;

  m_force1.reset(NinG, Nvol, m_Ndim);

  m_Cup.reset(NinG, Nvol, m_Ndim * (m_Ndim-1));
  m_Cdn.reset(NinG, Nvol, m_Ndim * (m_Ndim-1));

  m_Utmp1.reset(NinG, Nvol, 1);
  m_Utmp2.reset(NinG, Nvol, 1);

  m_zeta.reset(NinF, Nvol, 1);
  m_eta2.reset(NinF, Nvol, 1);
  m_eta3.reset(NinF, Nvol, 1);
  m_zeta_mu.reset(NinF, Nvol, 1);

  m_vt1.reset(NinF, Nvol, 1);
  m_vt2.reset(NinF, Nvol, 1);
  m_vt3.reset(NinF, Nvol, 1);

  set_parameters(params);

  vout.decrease_indent();

  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());
}

//====================================================================
template<typename AFIELD>
void AForce_F_CloverTerm<AFIELD>::tidyup()
{
  delete m_shift;
  delete m_staple;

  //delete m_forceF_ct;
  delete m_fopr_csw;
}

//====================================================================
template<typename AFIELD>
void AForce_F_CloverTerm<AFIELD>::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  std::string repr = params.get_string("gamma_matrix_type");
  m_repr   = repr;

  double kappa, cSW;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_double("hopping_parameter", kappa);
  err += params.fetch_double("clover_coefficient", cSW);
  err += params.fetch_int_vector("boundary_condition", bc);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(kappa, cSW, bc);

  m_fopr_csw->set_parameters(params);
}


//====================================================================
template<typename AFIELD>
void AForce_F_CloverTerm<AFIELD>::set_parameters(
                                        const double kappa,
                                        const double cSW,
                                        const std::vector<int> bc)
{
  assert(bc.size() == m_Ndim);

  int ith = ThreadManager::get_thread_id();

  if (ith == 0) {
    m_kappa = kappa;
    m_cSW   = cSW;

    m_boundary.resize(m_Ndim);
    m_boundary = bc;
  }

  //- print input parameters
  vout.general(m_vl, "%s: set parameters\n", class_name.c_str());
  vout.general(m_vl, "  gamma-matrix type = %s\n", m_repr.c_str());
  vout.general(m_vl, "  kappa = %12.8f\n", m_kappa);
  vout.general(m_vl, "  cSW   = %12.8f\n", cSW);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }

}


//====================================================================
template<typename AFIELD>
void AForce_F_CloverTerm<AFIELD>::get_parameters(Parameters& params) const
{
  params.set_double("hopping_parameter", m_kappa);
  params.set_double("clover_coefficient", m_cSW);
  params.set_int_vector("boundary_condition", m_boundary);
  params.set_string("gamma_matrix_type", m_repr);
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
template<typename AFIELD>
void AForce_F_CloverTerm<AFIELD>::set_config(Field* U)
{
  m_fopr_csw->set_config(U);
  //m_forceF_ct->set_config(U);

  int nth = ThreadManager::get_num_threads();

  if (nth > 1) {
    set_config_impl(U);
  } else {
    set_config_omp(U);
  }

}

//====================================================================
template<typename AFIELD>
void AForce_F_CloverTerm<AFIELD>::set_config_omp(Field* U)
{
#pragma omp parallel
 {
   set_config_impl(U);
 }

}

//====================================================================
template<typename AFIELD>
void AForce_F_CloverTerm<AFIELD>::set_config_impl(Field* U)
{
  AIndex_lex<real_t, AFIELD::IMPL> index;
  convert_gauge(index, m_Ucp, *U);

  set_component();
}

//====================================================================
template<typename AFIELD>
void AForce_F_CloverTerm<AFIELD>::force_udiv(AFIELD& force,
                                             const AFIELD& eta)
{
  m_fopr_csw->set_mode("H");
  m_fopr_csw->mult(m_zeta, eta);

  force_udiv1_impl(m_force1, m_zeta, eta);
  copy(force, m_force1);

  force_udiv1_impl(m_force1, eta, m_zeta);
  axpy(force, 1.0, m_force1);
}


//====================================================================
template<typename AFIELD>
void AForce_F_CloverTerm<AFIELD>::force_udiv1(AFIELD& force,
                                              const AFIELD& zeta,
                                              const AFIELD& eta)
{
  force_udiv1_impl(force, zeta, eta);
}


//====================================================================
template<typename AFIELD>
void AForce_F_CloverTerm<AFIELD>::set_component()
{
#pragma omp barrier

  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;

      m_staple->upper(m_Utmp1, m_Ucp, mu, nu);
      copy(m_Cup, index_dir(mu, nu), m_Utmp1, 0);

      m_staple->lower(m_Utmp2, m_Ucp, mu, nu);
      copy(m_Cdn, index_dir(mu, nu), m_Utmp2, 0);
    }
  }

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AForce_F_CloverTerm<AFIELD>::force_udiv1_impl(AFIELD& force,
                                                   const AFIELD& zeta,
                                                   const AFIELD& eta)
{
#pragma omp barrier

  force_udiv1_impl2(force, zeta, eta);

#pragma omp barrier
}

//====================================================================
/*
template<typename AFIELD>
void AForce_F_CloverTerm<AFIELD>::force_udiv1_impl1(AFIELD& force,
                                                   const AFIELD& zeta,
                                                   const AFIELD& eta)
{
  ThreadManager::assert_single_thread(class_name);

  const int Nc   = CommonParameters::Nc();
  const int Nd   = CommonParameters::Nd();
  const int NinF = Nc * Nd * 2;
  const int NinG = Nc * Nc * 2;
  const int Nvol = CommonParameters::Nvol();

  Field zeta_c(NinF, Nvol, 1);
  Field eta_c(NinF, Nvol, 1);
  Field force_c(NinG, Nvol, m_Ndim);

  AIndex_lex<real_t, AFIELD::IMPL> index_alt;

  if (m_fopr_csw->needs_convert()) {
    vout.detailed(m_vl, "convert required.\n");
    m_fopr_csw->reverse(zeta_c, zeta);
    m_fopr_csw->reverse(eta_c,  eta);
  } else {
    vout.detailed(m_vl, "convert not required.\n");
    reverse(index_alt, zeta_c, zeta);
    reverse(index_alt, eta_c,  eta);
  }

  m_forceF_ct->force_udiv1(force_c, zeta_c, eta_c);

  convert_gauge(index_alt, force, force_c);

}
*/

//====================================================================
template<typename AFIELD>
int AForce_F_CloverTerm<AFIELD>::index_dir(const int mu, const int nu)
{
  int nu2 = nu;
  if(nu > mu) nu2 = nu - 1;

  return mu + m_Ndim * nu2;
}

//====================================================================
template<typename AFIELD>
void AForce_F_CloverTerm<AFIELD>::force_udiv1_impl2(AFIELD& force,
                                                   const AFIELD& zeta,
                                                   const AFIELD& eta)
{
#pragma omp barrier

  force.set(0.0);

#pragma omp barrier

  real_t fac = -m_kappa * m_cSW / 8.0;

  m_fopr_csw->mult_gm5(m_eta2, eta);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;

      m_fopr_csw->mult_isigma(m_eta3, m_eta2, mu, nu);

      copy(m_Utmp1, 0, m_Cup, index_dir(mu, nu));
      axpy(m_Utmp1, 0, -1.0, m_Cdn, index_dir(mu, nu));

      // R(1) and R(5)
      Alt_Spinor::mult_Gdv(m_vt1, 0, m_Utmp1, 0, m_eta3, 0);
      Alt_Spinor::tensorProd(m_Utmp2, 0, zeta, m_vt1);
      axpy(force, mu, fac, m_Utmp2, 0);

      // R(4) and R(8)
      m_shift->backward(m_vt1, m_eta3, mu);
      m_shift->backward(m_zeta_mu, zeta, mu);
      Alt_Spinor::mult_Gnv(m_vt2, 0, m_Utmp1, 0, m_zeta_mu, 0);
      Alt_Spinor::tensorProd(m_Utmp2, 0, m_vt2, m_vt1);
      axpy(force, mu, fac, m_Utmp2, 0);

      // R(2)
      m_shift->backward(m_vt1, m_eta3, nu);
      Alt_Spinor::mult_Gnv(m_vt2, 0, m_Cup, index_dir(nu, mu), m_vt1, 0);
      Alt_Spinor::mult_Gdv(m_vt1, 0, m_Ucp, mu, m_vt2, 0);
      m_shift->backward(m_vt2, zeta, nu);
      Alt_Spinor::mult_Gnv(m_vt3, 0, m_Ucp, nu, m_vt2, 0);
      Alt_Spinor::tensorProd(m_Utmp2, 0, m_vt3, m_vt1);
      axpy(force, mu, fac, m_Utmp2, 0);

      // R(3)
      m_shift->backward(m_vt1, m_eta3, nu);
      Alt_Spinor::mult_Gnv(m_vt2, 0, m_Ucp, nu, m_vt1, 0);
      m_shift->backward(m_vt1, m_vt2, mu);
      Alt_Spinor::mult_Gnv(m_vt3, 0, m_Ucp, mu, m_zeta_mu, 0);
      m_shift->backward(m_vt2, m_vt3, nu);
      Alt_Spinor::mult_Gnv(m_vt3, 0, m_Ucp, nu, m_vt2, 0);
      Alt_Spinor::tensorProd(m_Utmp2, 0, m_vt3, m_vt1);
      axpy(force, mu, fac, m_Utmp2, 0);

      // R(6)
      Alt_Spinor::mult_Gdv(m_vt1, 0, m_Cup, index_dir(nu, mu), m_eta3, 0);
      m_shift->forward(m_vt2, m_vt1, nu);
      Alt_Spinor::mult_Gdv(m_vt1, 0, m_Ucp, mu, m_vt2, 0);
      Alt_Spinor::mult_Gdv(m_vt2, 0, m_Ucp, nu, zeta, 0);
      m_shift->forward(m_vt3, m_vt2, nu);
      Alt_Spinor::tensorProd(m_Utmp2, 0, m_vt3, m_vt1);
      axpy(force, mu, -fac, m_Utmp2, 0);

      // R(7)
      Alt_Spinor::mult_Gdv(m_vt1, 0, m_Ucp, nu, m_eta3, 0);
      m_shift->backward(m_vt2, m_vt1, mu);
      m_shift->forward(m_vt1, m_vt2, nu);
      Alt_Spinor::mult_Gnv(m_vt2, 0, m_Ucp, mu, m_zeta_mu, 0);
      Alt_Spinor::mult_Gdv(m_vt3, 0, m_Ucp, nu, m_vt2, 0);
      m_shift->forward(m_vt2, m_vt3, nu);
      Alt_Spinor::tensorProd(m_Utmp2, 0, m_vt2, m_vt1);
      axpy(force, mu, -fac, m_Utmp2, 0);
    }
  }

#pragma omp barrier
}

//============================================================END=====
