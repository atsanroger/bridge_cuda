/*!
        @file    aforce_F_CloverTerm-tmpl.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
        @version $LastChangedRevision: 2160 $
*/

#include "lib_alt_SIMD2/Force/Fermion/aforce_F_CloverTerm.h"

template<typename AFIELD>
const std::string AForce_F_CloverTerm<AFIELD>::class_name
                                  = "AForce_F_CloverTerm<AFIELD>";

//====================================================================
template<typename AFIELD>
void AForce_F_CloverTerm<AFIELD>::init(const Parameters& params)
{
  set_parameters(params);

  int Nc   = CommonParameters::Nc();
  int NinG = 2 * Nc * Nc;
  int NinF = m_fopr->field_nin();
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();
  m_Ndim = Ndim;

  m_shiftF = new ShiftAField_lex<AFIELD>(NinF);
  m_shiftG = new ShiftAField_lex<AFIELD>(NinG);
  m_staple = new AStaple_lex<AFIELD>;

  m_Cud.reset(NinG, Nvol, Ndim * Ndim);

  m_force1.reset(NinG, Nvol, Ndim);

  m_zeta.reset(NinF, Nvol, 1);

  m_SG.resize(Ndim * Ndim);
  unique_ptr<GammaMatrixSet> gmset(GammaMatrixSet::New(m_repr));

  m_SG[index_dir(0, 1)] = gmset->get_GM(gmset->SIGMA12);
  m_SG[index_dir(1, 2)] = gmset->get_GM(gmset->SIGMA23);
  m_SG[index_dir(2, 0)] = gmset->get_GM(gmset->SIGMA31);
  m_SG[index_dir(3, 0)] = gmset->get_GM(gmset->SIGMA41);
  m_SG[index_dir(3, 1)] = gmset->get_GM(gmset->SIGMA42);
  m_SG[index_dir(3, 2)] = gmset->get_GM(gmset->SIGMA43);

  m_SG[index_dir(1, 0)] = m_SG[index_dir(0, 1)].mult(-1);
  m_SG[index_dir(2, 1)] = m_SG[index_dir(1, 2)].mult(-1);
  m_SG[index_dir(0, 2)] = m_SG[index_dir(2, 0)].mult(-1);
  m_SG[index_dir(0, 3)] = m_SG[index_dir(3, 0)].mult(-1);
  m_SG[index_dir(1, 3)] = m_SG[index_dir(3, 1)].mult(-1);
  m_SG[index_dir(2, 3)] = m_SG[index_dir(3, 2)].mult(-1);

  m_SG[index_dir(0, 0)] = gmset->get_GM(gmset->UNITY);
  m_SG[index_dir(1, 1)] = gmset->get_GM(gmset->UNITY);
  m_SG[index_dir(2, 2)] = gmset->get_GM(gmset->UNITY);
  m_SG[index_dir(3, 3)] = gmset->get_GM(gmset->UNITY);

}

//====================================================================
template<typename AFIELD>
void AForce_F_CloverTerm<AFIELD>::tidyup()
{
  delete m_shiftF;
  delete m_shiftG;
  delete m_staple;
}

//====================================================================
template<typename AFIELD>
void AForce_F_CloverTerm<AFIELD>::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  double kappa;
  double cSW;
  std::vector<int> bc;
  std::string repr;

  int err = 0;
  err += params.fetch_double("hopping_parameter", kappa);
  err += params.fetch_double("clover_coefficient", cSW);
  err += params.fetch_int_vector("boundary_condition", bc);
  err += params.fetch_string("gamma_matrix_type", repr);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  if(repr != "Dirac"){
    vout.crucial(m_vl,
                 "Error at %s: gamma matrix type = %s is unsupported.\n",
                 class_name.c_str(), repr.c_str());
    exit(EXIT_FAILURE);
  }

  m_repr = repr;

  set_parameters(kappa, cSW, bc);

}

//====================================================================
template<typename AFIELD>
void AForce_F_CloverTerm<AFIELD>::set_parameters(
                                          const double kappa,
                                          const double cSW,
                                          const std::vector<int> bc)
{
  const int Ndim = CommonParameters::Ndim();

  assert(bc.size() == Ndim);

  m_kappa = real_t(kappa);
  m_cSW   = real_t(cSW);

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
void AForce_F_CloverTerm<AFIELD>::set_config(Field *U)
{
  ThreadManager_OpenMP::assert_single_thread(class_name);

  m_U = (Field_G*)U;

  Index_lex_alt<real_t,AFIELD::IMPL> index_lex;

#pragma omp parallel
{ convert_gauge(index_lex, m_Ucp, *U); }

 set_component();

  // m_fopr->set_config(U);
  // this should be done outside this class

}

//====================================================================
template<typename AFIELD>
void AForce_F_CloverTerm<AFIELD>::force_udiv(AFIELD& force,
                                             const AFIELD& eta)
{
  vout.crucial(m_vl, "Error at %s: force_udiv() is not available.\n",
               class_name.c_str());
  exit(EXIT_FAILURE);
}

//====================================================================
template<typename AFIELD>
void AForce_F_CloverTerm<AFIELD>::force_udiv1(AFIELD& force,
                              const AFIELD& zeta, const AFIELD& eta)
{
  force_udiv1_impl(force, zeta, eta);
}

//====================================================================
template<typename AFIELD>
void AForce_F_CloverTerm<AFIELD>::force_udiv1_impl(AFIELD& force,
                              const AFIELD& zeta, const AFIELD& eta)
{
  int Nc   = CommonParameters::Nc();
  int Nd   = CommonParameters::Nd();
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();
  int NinG = 2 * Nc * Nc;
  int NinF = 2 * Nc * Nd;

  AFIELD eta2(NinF, Nvol, 1), eta3(NinF, Nvol, 1);
  AFIELD vt1(NinF, Nvol, 1), vt2(NinF, Nvol, 1);
  AFIELD vt3(NinF, Nvol, 1), vt4(NinF, Nvol, 1);
  AFIELD zeta_mu(NinF, Nvol, 1);

  AFIELD Umu(NinG, Nvol, 1), Unu(NinG, Nvol, 1);
  AFIELD force1(NinG, Nvol, 1), force2(NinG, Nvol, 1);
  AFIELD Utmp(NinG, Nvol, 1), Utmp2(NinG, Nvol, 1);

#pragma omp parallel
 {
  force.set(0.0);
  #pragma omp barrier

  m_fopr->mult_gm5(eta2, eta);
  #pragma omp barrier

  for (int mu = 0; mu < Ndim; ++mu) {
    for (int nu = 0; nu < Ndim; ++nu) {
      if (nu == mu) continue;

      mult_isigma(eta3, eta2, mu, nu);

      copy(Umu, 0, m_Ucp, mu);
      copy(Unu, 0, m_Ucp, nu);
      #pragma omp barrier


      // R(1) and R(5)
      mult_Gd(vt1, 0, m_Cud, index_dir(mu, nu), eta3, 0);
      tensorProd_4spinor(force1, 0, zeta, 0, vt1, 0);
      copy(force2, force1);
      #pragma omp barrier

      // R(2)
      mult_Gd(vt3, 0, Umu, 0, eta3, 0);
      m_shiftF->backward(vt1, vt3, nu);
      m_shiftF->backward(vt2, zeta, nu);
      m_shiftG->backward(Utmp, Unu, mu);
      mult_Gn(vt3, 0, Utmp, 0, vt1, 0);
      mult_Gn(vt4, 0, Unu, 0, vt2, 0);
      tensorProd_4spinor(force1, 0, vt4, 0, vt3, 0);
      axpy(force2, 1.0, force1);
      #pragma omp barrier

      // R(4) and R(8)
      m_shiftF->backward(vt1, eta3, mu);
      m_shiftF->backward(zeta_mu, zeta, mu);
      mult_Gn(vt4, 0, m_Cud, index_dir(mu, nu), zeta_mu, 0);
      tensorProd_4spinor(force1, 0, vt4, 0, vt1, 0);
      axpy(force2, 1.0, force1);
      #pragma omp barrier

      // R(3)
      m_shiftF->backward(vt1, eta3, nu);
      mult_Gn(vt3, 0, Unu, 0, vt1, 0);
      mult_Gn(vt4, 0, Umu, 0, zeta_mu, 0);
      m_shiftF->backward(vt1, vt3, mu);
      m_shiftF->backward(vt2, vt4, nu);
      mult_Gn(vt4, 0, Unu, 0, vt2, 0);
      tensorProd_4spinor(force1, 0, vt4, 0, vt1, 0);
      axpy(force2, 1.0, force1);
      #pragma omp barrier

      // R(6)
      m_shiftG->backward(Utmp, Unu, mu);
      mult_Gdd(Utmp2, 0, Utmp, 0, Umu, 0);
      mult_Gn(vt1, 0, Utmp2, 0, eta3, 0);
      mult_Gd(vt2, 0, Unu, 0, zeta, 0);
      m_shiftF->forward(vt3, vt1, nu);
      m_shiftF->forward(vt4, vt2, nu);
      tensorProd_4spinor(force1, 0, vt4, 0, vt3, 0);
      axpy(force2, -1.0, force1);
      #pragma omp barrier

      // R(7)
      mult_Gd(vt1, 0, Unu, 0, eta3, 0);
      mult_Gn(vt2, 0, Umu, 0, zeta_mu, 0);
      m_shiftF->backward(vt3, vt1, mu);
      m_shiftF->forward(vt1, vt3, nu);
      mult_Gd(vt4, 0, Unu, 0, vt2, 0);
      m_shiftF->forward(vt2, vt4, nu);
      tensorProd_4spinor(force1, 0, vt2, 0, vt1, 0);
      axpy(force2, -1.0, force1);
      #pragma omp barrier

      double fac = -m_kappa * m_cSW / 8.0;
      scal(force2, fac);
      #pragma omp barrier

      axpy(force, mu, 1.0, force2, 0);
      #pragma omp barrier
    }
  }
 }

}

//====================================================================
template<typename AFIELD>
void AForce_F_CloverTerm<AFIELD>::set_component()
{
  int Nc   = CommonParameters::Nc();
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();
  int NinG = 2 * Nc * Nc;

  AFIELD Cmu_ud1(NinG, Nvol, 1);
  AFIELD Cmu_ud2(NinG, Nvol, 1);

#pragma omp parallel
 {
  for (int mu = 0; mu < Ndim; ++mu) {
    for (int nu = 0; nu < Ndim; ++nu) {
      if (nu == mu) continue;
      m_staple->upper_th(Cmu_ud1, m_Ucp, mu, nu);
      m_staple->lower_th(Cmu_ud2, m_Ucp, mu, nu);
      axpy(Cmu_ud1, -1.0, Cmu_ud2);
      copy(m_Cud, index_dir(mu, nu), Cmu_ud1, 0);
    }
  }
 }

  /*
  for (int mu = 0; mu < Ndim; ++mu) {
    for (int nu = 0; nu < Ndim; ++nu) {
      if (nu == mu) continue;
      m_staple->upper(Cmu_ud1, m_Ucp, mu, nu);
      m_staple->lower(Cmu_ud2, m_Ucp, mu, nu);
      axpy(Cmu_ud1, -1.0, Cmu_ud2);
      copy(m_Cud, index_dir(mu, nu), Cmu_ud1, 0);
    }
  }
  */

}

//====================================================================
template<typename AFIELD>
void AForce_F_CloverTerm<AFIELD>::mult_isigma(AFIELD& v,
                                   const AFIELD& w, int mu, int nu)
{
  // this function may be called inside parallel region.

  assert(mu != nu);
  mult_iGM(v, 0, m_SG[index_dir(mu, nu)], w, 0);

}

//============================================================END=====
