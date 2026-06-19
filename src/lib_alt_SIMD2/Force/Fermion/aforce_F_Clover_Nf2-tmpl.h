/*!
        @file    aforce_F_Clover_Nf2-tmpl.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
        @version $LastChangedRevision: 2160 $
*/


template<typename AFIELD>
const std::string AForce_F_Clover_Nf2<AFIELD>::class_name
                                  = "AForce_F_Clover_Nf2<AFIELD>";

//====================================================================
template<typename AFIELD>
void AForce_F_Clover_Nf2<AFIELD>::init(const Parameters& params)
{
  double kappa, cSW;
  std::vector<int> bc;

  extract_parameters(kappa, cSW, bc, params);

  set_parameters(kappa, cSW, bc);

  m_force_w   = new AForce_F_Wilson_Nf2<AFIELD>(m_fopr, params);
  m_force_csw = new AForce_F_CloverTerm<AFIELD>(m_fopr, params);

  int Nc   = CommonParameters::Nc();
  int NinG = 2 * Nc * Nc;
  int NinF = m_fopr->field_nin();
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  m_force1.reset(NinG, Nvol, Ndim);
  m_force2.reset(NinG, Nvol, Ndim);

  m_zeta.reset(NinF, Nvol, 1);

}

//====================================================================
template<typename AFIELD>
void AForce_F_Clover_Nf2<AFIELD>::tidyup()
{
  delete m_force_w;
  delete m_force_csw;
}

//====================================================================
template<typename AFIELD>
void AForce_F_Clover_Nf2<AFIELD>::extract_parameters(
                                          double& kappa,
                                          double& cSW,
                                          std::vector<int>& bc,
                                          const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);

  int err = 0;
  err += params.fetch_double("hopping_parameter", kappa);
  err += params.fetch_double("clover_coefficient", cSW);
  err += params.fetch_int_vector("boundary_condition", bc);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

}

//====================================================================
template<typename AFIELD>
void AForce_F_Clover_Nf2<AFIELD>::set_parameters(const Parameters& params)
{
  double kappa, cSW;
  std::vector<int> bc;

  extract_parameters(kappa, cSW, bc, params);
  set_parameters(kappa, cSW, bc);

  m_force_w->set_parameters(params);
  m_force_csw->set_parameters(params);

}

//====================================================================
template<typename AFIELD>
void AForce_F_Clover_Nf2<AFIELD>::set_parameters(
                                          const double kappa,
                                          const double cSW,
                                          const std::vector<int> bc)
{
  int Ndim = CommonParameters::Ndim();

  m_kappa = real_t(kappa);
  m_cSW   = real_t(cSW);

  assert(bc.size() == Ndim);
  m_boundary.resize(Ndim);
  m_boundary = bc;

  //- print input parameters
  vout.general(m_vl, "Parameters of %s:\n", class_name.c_str());
  vout.general(m_vl, "  kappa = %12.8f\n", m_kappa);
  vout.general(m_vl, "  cSW   = %12.8f\n", m_cSW);
  for (int mu = 0; mu < Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }

}

//====================================================================
template<typename AFIELD>
void AForce_F_Clover_Nf2<AFIELD>::set_config(Field *U)
{
  ThreadManager_OpenMP::assert_single_thread(class_name);

  m_U = (Field_G*)U;

  Index_lex_alt<real_t,AFIELD::IMPL> index_lex;

#pragma omp parallel
{ convert_gauge(index_lex, m_Ucp, *U); }

  // m_fopr->set_config(U);
  m_force_w->set_config(U);
  m_force_csw->set_config(U);

}

//====================================================================
template<typename AFIELD>
void AForce_F_Clover_Nf2<AFIELD>::force_udiv(AFIELD& force,
                                             const AFIELD& eta)
{
  m_fopr->set_mode("H");
  m_fopr->mult(m_zeta, eta);

  force_udiv1_impl(force, m_zeta, eta);

  force_udiv1_impl(m_force1, eta, m_zeta);
  axpy(force, 1.0, m_force1);

}

//====================================================================
template<typename AFIELD>
void AForce_F_Clover_Nf2<AFIELD>::force_udiv1(AFIELD& force,
                              const AFIELD& zeta, const AFIELD& eta)
{
  force_udiv1_impl(force, zeta, eta);
}

//====================================================================
template<typename AFIELD>
void AForce_F_Clover_Nf2<AFIELD>::force_udiv1_impl(AFIELD& force,
                              const AFIELD& zeta, const AFIELD& eta)
{
  m_force_w->force_udiv1(force, zeta, eta);

  m_force_csw->force_udiv1(m_force2, zeta, eta);
  axpy(force, 1.0, m_force2);

}

//====================================================================
//============================================================END=====
