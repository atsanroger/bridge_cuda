/*!
        @file    aforce_F_Wilson_Nf2-tmpl.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
        @version $LastChangedRevision: 2160 $
*/

//#include "Force/Fermion/aforce_F_Wilson_Nf2.h"

template<typename AFIELD>
const std::string AForce_F_Wilson_Nf2<AFIELD>::class_name
                                  = "AForce_F_Wilson_Nf2<AFIELD>";

//====================================================================
template<typename AFIELD>
void AForce_F_Wilson_Nf2<AFIELD>::init(const Parameters& params)
{
  set_parameters(params);

  int Nc   = CommonParameters::Nc();
  int NinG = 2 * Nc * Nc;
  int NinF = m_fopr_w->field_nin();
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  m_force1.reset(NinG, Nvol, Ndim);

  m_zeta.reset(NinF, Nvol, 1);
  m_eta2.reset(NinF, Nvol, 1);
  m_eta3.reset(NinF, Nvol, 1);

}

//====================================================================
template<typename AFIELD>
void AForce_F_Wilson_Nf2<AFIELD>::tidyup()
{
  // do nothing.
}

//====================================================================
template<typename AFIELD>
void AForce_F_Wilson_Nf2<AFIELD>::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
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

}

//====================================================================
template<typename AFIELD>
void AForce_F_Wilson_Nf2<AFIELD>::set_parameters(
                                          const double kappa,
                                          const std::vector<int> bc)
{
  const int Ndim = CommonParameters::Ndim();

  //- range check
  // NB. kappa == 0 is allowed.
  assert(bc.size() == Ndim);

  m_kappa = kappa;

  m_boundary.resize(Ndim);
  m_boundary = bc;

  // m_fopr_w->set_parameters(m_kappa, m_boundary);
  // - Note that m_fopr_w is given outside this class.

  //- print input parameters
  vout.general(m_vl, "Parameters of %s:\n", class_name.c_str());
  vout.general(m_vl, "  kappa = %12.8f\n", m_kappa);
  for (int mu = 0; mu < Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }

}

//====================================================================
template<typename AFIELD>
void AForce_F_Wilson_Nf2<AFIELD>::set_config(Field *U)
{
  ThreadManager_OpenMP::assert_single_thread(class_name);

  m_U = (Field_G*)U;

  Index_lex_alt<real_t,AFIELD::IMPL> index_lex;

#pragma omp parallel
{ convert_gauge(index_lex, m_Ucp, *U); }

  // m_fopr_w->set_config(U);

}

//====================================================================
template<typename AFIELD>
void AForce_F_Wilson_Nf2<AFIELD>::force_udiv(AFIELD& force,
                                             const AFIELD& eta)
{
  m_fopr_w->set_mode("H");
  m_fopr_w->mult(m_zeta, eta);

  force_udiv1_impl(m_force1, m_zeta, eta);
#pragma omp parallel
{ copy(force, m_force1); }

  force_udiv1_impl(m_force1, eta, m_zeta);
#pragma omp parallel
{ axpy(force, 1.0, m_force1); }

}

//====================================================================
template<typename AFIELD>
void AForce_F_Wilson_Nf2<AFIELD>::force_udiv1(AFIELD& force,
                              const AFIELD& zeta, const AFIELD& eta)
{
  force_udiv1_impl(force, zeta, eta);
}

//====================================================================
template<typename AFIELD>
void AForce_F_Wilson_Nf2<AFIELD>::force_udiv1_impl(AFIELD& force,
                              const AFIELD& zeta, const AFIELD& eta)
{
  int Ndim = CommonParameters::Ndim();

  force.set(0.0);

#pragma omp parallel
 {
  for (int mu = 0; mu < Ndim; ++mu) {
    m_eta3.set(0.0);
    m_fopr_w->mult_up(mu, m_eta3, eta);
    m_fopr_w->mult_gm5(m_eta2, m_eta3);

    mult_Gd(m_eta3, 0, m_Ucp, mu, m_eta2, 0);

    tensorProd_4spinor(force, mu, zeta, 0, m_eta3, 0);

  }
  #pragma omp barrier
  scal(force, -m_kappa);
 }

}

//====================================================================
//============================================================END=====
