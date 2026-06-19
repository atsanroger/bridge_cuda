/*!
        @file    force_F_Wilson_eo_Nf2.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-07-03 23:32:09 #$

        @version $LastChangedRevision: 2529 $
*/

#include "lib_alt/Force/Fermion/aforce_F_Wilson_eo_Nf2.h"

template<typename AFIELD>
const std::string AForce_F_Wilson_eo_Nf2<AFIELD>::class_name
                                 = "AForce_F_Wilson_eo_Nf2<AFIELD>";

//====================================================================
template<typename AFIELD>
void AForce_F_Wilson_eo_Nf2<AFIELD>::init(const Parameters& params)
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

  m_fopr_w = new AFopr_Wilson_eo<AFIELD>(params);

  m_force_w = new AForce_F_Wilson_Nf2<AFIELD>(params);

  int Nc = CommonParameters::Nc();
  //int Nd = CommonParameters::Nd();
  //int NinF = Nc * Nd * 2;
  int NinG = Nc * Nc * 2;
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  //m_eta2.reset(NinF, Nvol, 1);
  //m_eta3.reset(NinF, Nvol, 1);

  int NinF   = m_fopr_w->field_nin();
  int NvolF  = m_fopr_w->field_nvol();
  int NexF   = m_fopr_w->field_nex();

  if(NvolF != Nvol/2){
    vout.crucial(m_vl, "%s: wrong volume size\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  //  int NvolF = Nvol/2;
  m_zeta.reset(NinF, NvolF, NexF);
  m_v1.reset(NinF, NvolF, NexF);

  m_eta_lex.reset(NinF, Nvol, NexF);
  m_zeta_lex.reset(NinF, Nvol, NexF);

  m_force1.reset(NinG, Nvol, Ndim);
  //  m_force2.reset(NinG, Nvol, Ndim);

  m_at_construct = true;
  set_parameters(params);
  m_at_construct = false;

  vout.decrease_indent();

  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());
}

//====================================================================
template<typename AFIELD>
void AForce_F_Wilson_eo_Nf2<AFIELD>::tidyup()
{
  delete m_force_w;
  delete m_fopr_w;
}

//====================================================================
template<typename AFIELD>
void AForce_F_Wilson_eo_Nf2<AFIELD>::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  std::string repr = params.get_string("gamma_matrix_type");
  m_repr   = repr;

  double kappa;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_double("hopping_parameter", kappa);
  err += params.fetch_int_vector("boundary_condition", bc);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(kappa, bc);

  if(!m_at_construct){
    m_fopr_w->set_parameters(params);
    m_force_w->set_parameters(params);
  }

}

//====================================================================
template<typename AFIELD>
void AForce_F_Wilson_eo_Nf2<AFIELD>::set_parameters(const double kappa,
                                        const std::vector<int> bc)
{
  const int Ndim = CommonParameters::Ndim();
  assert(bc.size() == Ndim);

  int ith = ThreadManager::get_thread_id();

  if (ith == 0) {
    m_kappa = kappa;

    m_boundary.resize(Ndim);
    m_boundary = bc;
  }
#pragma omp barrier

  //- print input parameters
  vout.general(m_vl, "%s: set parameters\n", class_name.c_str());
  vout.general(m_vl, "  gamma-matrix type = %s\n", m_repr.c_str());
  vout.general(m_vl, "  kappa = %12.8f\n", m_kappa);
  for (int mu = 0; mu < Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }

}


//====================================================================
template<typename AFIELD>
void AForce_F_Wilson_eo_Nf2<AFIELD>::get_parameters(Parameters& params) const
{
  params.set_double("hopping_parameter", m_kappa);
  params.set_int_vector("boundary_condition", m_boundary);
  params.set_string("gamma_matrix_type", m_repr);
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
template<typename AFIELD>
void AForce_F_Wilson_eo_Nf2<AFIELD>::set_config(Field* U)
{
  m_fopr_w->set_config(U);
  m_force_w->set_config(U);

  AIndex_lex<real_t, AFIELD::IMPL> index;
  convert_gauge(index, m_Ucp, *U);
}

//====================================================================
template<typename AFIELD>
void AForce_F_Wilson_eo_Nf2<AFIELD>::force_udiv(AFIELD& force,
                                             const AFIELD& eta)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  int Nin  = force.nin();
  //  int Nvol = force.nvol();
  int Nex  = force.nex();
  //  AFIELD force1(Nin, Nvol, Nex);

  int NinF   = m_fopr_w->field_nin();
  int NvolF  = m_fopr_w->field_nvol();
  int NexF   = m_fopr_w->field_nex();

  // AFIELD zeta(NinF, NvolF, NexF);

  m_fopr_w->set_mode("H");
  m_fopr_w->mult(m_zeta, eta);

  force_udiv1_impl(force, m_zeta, eta);
  force_udiv1_impl(m_force1, eta, m_zeta);
  axpy(force, real_t(1.0), m_force1);
  //  scal(force, real_t(-1.0));

}


//====================================================================
template<typename AFIELD>
void AForce_F_Wilson_eo_Nf2<AFIELD>::force_udiv1(AFIELD& force,
                                              const AFIELD& zeta,
                                              const AFIELD& eta)
{
  force_udiv1_impl(force, zeta, eta);
  //  scal(force, real_t(-1.0));
}


//====================================================================
template<typename AFIELD>
void AForce_F_Wilson_eo_Nf2<AFIELD>::force_udiv1_impl(AFIELD& force,
                                                   const AFIELD& zeta,
                                                   const AFIELD& eta)
{
#pragma omp barrier

  // const int Nvol = CommonParameters::Nvol();
  //const int Ndim = CommonParameters::Ndim();

  /*
  int Nin  = force.nin();
  //  int Nvol = force.nvol();
  int Nex  = force.nex();
  //AFIELD force1(Nin, Nvol, Nex);

  int NinF   = m_fopr_w->field_nin();
  int NvolF  = m_fopr_w->field_nvol();
  int NexF   = m_fopr_w->field_nex();

  //AFIELD v1(NinF, NvolF, NexF);

  //AFIELD eta_lex(NinF, Nvol, NexF);
  //AFIELD zeta_lex(NinF, Nvol, NexF);
  */
  AIndex_eo<real_t, AFIELD::IMPL> index_eo;

  force.set(0.0);
#pragma omp barrier

  m_fopr_w->mult(m_v1, eta, "Doe");
  index_eo.merge(m_eta_lex, eta, m_v1);
 
  m_fopr_w->mult(m_v1, zeta, "Doe");
  index_eo.merge(m_zeta_lex, zeta, m_v1);

  scal(m_zeta_lex, real_t(-1.0));
#pragma omp barrier

  m_force_w->force_udiv1(force, m_zeta_lex, m_eta_lex);

  /*
  for (int mu = 0; mu < Ndim; ++mu) {
    m_eta3.set(0.0);
    m_fopr_w->mult_up(mu, m_eta3, eta);
    m_fopr_w->mult_gm5(m_eta2, m_eta3);

    Alt_Spinor::mult_Gdv(m_eta3, 0, m_Ucp, mu, m_eta2, 0);

    Alt_Spinor::tensorProd(force, mu, zeta, m_eta3);
  }

  scal(force, -m_kappa); // force *= -m_kappa;
  */
#pragma omp barrier

}

//============================================================END=====
