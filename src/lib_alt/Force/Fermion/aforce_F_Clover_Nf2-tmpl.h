/*!
        @file    force_F_Clover_Nf2.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-10-17 13:51:32 #$

        @version $LastChangedRevision: 2550 $
*/

#include "lib_alt/Force/Fermion/aforce_F_Clover_Nf2.h"

template<typename AFIELD>
const std::string AForce_F_Clover_Nf2<AFIELD>::class_name
                                 = "AForce_F_Clover_Nf2<AFIELD>";

//====================================================================
template<typename AFIELD>
void AForce_F_Clover_Nf2<AFIELD>::init(const Parameters& params)
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

  m_fopr_w = new AFopr_Clover<AFIELD>(params);

  //m_forceF = new Force_F_Clover_Nf2(params);
  //m_forceF_w  = new Force_F_Wilson_Nf2(params);
  m_forceF_w  = new AForce_F_Wilson_Nf2<AFIELD>(params);
  m_forceF_ct = new AForce_F_CloverTerm<AFIELD>(params);

  int Nc = CommonParameters::Nc();
  int Nd = CommonParameters::Nd();
  int NinF = Nc * Nd * 2;
  int NinG = Nc * Nc * 2;
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  m_zeta.reset(NinF, Nvol, 1);
  m_eta2.reset(NinF, Nvol, 1);
  m_eta3.reset(NinF, Nvol, 1);

  //  m_zeta_c.reset(NinF, Nvol, 1);
  //  m_eta_c.reset(NinF, Nvol, 1);
  //  m_force_c.reset(NinG, Nvol, Ndim);

  m_U.reset(NinG, Nvol, Ndim);
  m_force1.reset(NinG, Nvol, Ndim);
  m_force2.reset(NinG, Nvol, Ndim);

  set_parameters(params);

  vout.decrease_indent();

  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());
}

//====================================================================
template<typename AFIELD>
void AForce_F_Clover_Nf2<AFIELD>::tidyup()
{
  delete m_forceF_ct;
  delete m_forceF_w;
  delete m_fopr_w;
}

//====================================================================
template<typename AFIELD>
void AForce_F_Clover_Nf2<AFIELD>::set_parameters(const Parameters& params)
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

  m_fopr_w->set_parameters(params);
}


//====================================================================
template<typename AFIELD>
void AForce_F_Clover_Nf2<AFIELD>::set_parameters(
                                        const double kappa,
                                        const double cSW,
                                        const std::vector<int> bc)
{
  const int Ndim = CommonParameters::Ndim();
  assert(bc.size() == Ndim);

  int ith = ThreadManager::get_thread_id();

  if (ith == 0) {
    m_kappa = kappa;
    m_cSW   = cSW;

    m_boundary.resize(Ndim);
    m_boundary = bc;
  }

  //- print input parameters
  vout.general(m_vl, "%s: set parameters\n", class_name.c_str());
  vout.general(m_vl, "  gamma-matrix type = %s\n", m_repr.c_str());
  vout.general(m_vl, "  kappa = %12.8f\n", m_kappa);
  vout.general(m_vl, "  cSW   = %12.8f\n", cSW);
  for (int mu = 0; mu < Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }

}


//====================================================================
template<typename AFIELD>
void AForce_F_Clover_Nf2<AFIELD>::get_parameters(Parameters& params) const
{
  params.set_double("hopping_parameter", m_kappa);
  params.set_double("clover_coefficient", m_cSW);
  params.set_int_vector("boundary_condition", m_boundary);
  params.set_string("gamma_matrix_type", m_repr);
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
template<typename AFIELD>
void AForce_F_Clover_Nf2<AFIELD>::set_config(Field* U)
{
  m_fopr_w->set_config(U);
  m_forceF_w->set_config(U);
  m_forceF_ct->set_config(U);

  AIndex_lex<real_t, AFIELD::IMPL> index;
  convert_gauge(index, m_Ucp, *U);
}

//====================================================================
template<typename AFIELD>
void AForce_F_Clover_Nf2<AFIELD>::force_udiv(AFIELD& force,
                                             const AFIELD& eta)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  m_fopr_w->set_mode("H");
  m_fopr_w->mult(m_zeta, eta);

  force_udiv1_impl(m_force1, m_zeta, eta);
  copy(force, m_force1);

  force_udiv1_impl(m_force1, eta, m_zeta);
  axpy(force, 1.0, m_force1);
}


//====================================================================
template<typename AFIELD>
void AForce_F_Clover_Nf2<AFIELD>::force_udiv1(AFIELD& force,
                                              const AFIELD& zeta,
                                              const AFIELD& eta)
{
  force_udiv1_impl(force, zeta, eta);
}


//====================================================================
template<typename AFIELD>
void AForce_F_Clover_Nf2<AFIELD>::force_udiv1_impl(AFIELD& force,
                                                   const AFIELD& zeta,
                                                   const AFIELD& eta)
{
#pragma omp barrier

  m_forceF_ct->force_udiv1(force, zeta, eta);

  m_forceF_w->force_udiv1(m_force2, zeta, eta);

  axpy(force, 1.0, m_force2);

#pragma omp barrier
}

//============================================================END=====
