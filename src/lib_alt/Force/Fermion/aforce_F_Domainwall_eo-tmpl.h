/*!
        @File    aforce_F_Domainwall_eo-tmpl.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-07-03 23:32:09 #$

        @version $LastChangedRevision: 2529 $
*/

#include "lib_alt/Force/Fermion/aforce_F_Domainwall_eo.h"

template<typename AFIELD>
const std::string AForce_F_Domainwall_eo<AFIELD>::class_name
                               = "AForce_F_Domainwall_eo<AFIELD>";
//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_eo<AFIELD>::init(const Parameters& params)
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

  int Nc   = CommonParameters::Nc();
  int Nd   = CommonParameters::Nd();
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();
  int NinF = Nc * Nd * 2;
  int NinG = Nc * Nc * 2;

  m_boundary.resize(Ndim);

  set_parameters(params);

  m_fopr_dw = new AFopr_Domainwall_eo<AFIELD>(params);

  int err = 0;
  err += params.fetch_string("kernel_type", m_kernel_type);
  if (err > 0) {
    vout.crucial(m_vl, "%s: Error: kernel_type is not specified.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  Parameters params_kernel = params;
  double kappa = 1.0 / (8.0 - 2.0 * m_M0);
  params_kernel.set_double("hopping_parameter", kappa);

  if(m_kernel_type == "Wilson"){
    m_fopr_w  = new AFopr_Wilson_eo<AFIELD>(params_kernel);
  }else{
    vout.crucial(m_vl, "%s: Error: unsupported kernel_type: %s\n",
                 class_name.c_str(), m_kernel_type.c_str());
    exit(EXIT_FAILURE);
  }

  m_force_dw = new AForce_F_Domainwall<AFIELD>(params);

  int Nvol2 = Nvol/2;

  m_v1.reset(NinF, Nvol2, m_Ns);
  m_v2.reset(NinF, Nvol2, m_Ns);
  m_v3.reset(NinF, Nvol2, m_Ns);

  m_eta_lex.reset( NinF, Nvol, m_Ns);
  m_zeta_lex.reset(NinF, Nvol, m_Ns);

  m_z1.reset(NinF, Nvol2, m_Ns);

  m_force1.reset(NinG, Nvol, Ndim);
  m_force2.reset(NinG, Nvol, Ndim);

  m_index_eo = new AIndex_eo<real_t, AFIELD::IMPL>();

  vout.decrease_indent();

  vout.general(m_vl, "%s: construction finished.\n",
	       class_name.c_str());

}

//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_eo<AFIELD>::tidyup()
{
  delete m_force_dw;
  delete m_fopr_dw;
  delete m_fopr_w;
  delete m_index_eo;
}

//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_eo<AFIELD>::set_parameters(const Parameters& params)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    if(ith == 0) m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  double mq, M0, b, c;
  int Ns;
  std::vector<int> bc;
  std::string repr;

  int err = 0;
  err += params.fetch_double("quark_mass", mq);
  err += params.fetch_double("domain_wall_height", M0);
  err += params.fetch_int("extent_of_5th_dimension", Ns);
  err += params.fetch_int_vector("boundary_condition", bc);
  err += params.fetch_string("gamma_matrix_type", repr);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  if(ith == 0) m_repr = repr;

  int err2 = 0;
  err2 += params.fetch_double("coefficient_b", b);
  err2 += params.fetch_double("coefficient_c", c);

  if (err2) {
    vout.general(m_vl, "  coefficients b, c are not provided:"
                       " set to Shamir's form.\n");
    b = 1.0;
    c = 0.0;
  }

  set_parameters(real_t(mq), real_t(M0), Ns, bc, real_t(b), real_t(c));

  // setting parameters for kernel objects are necessary
  //  const real_t kappa = 1.0 / (8.0 - 2.0 * m_M0);
  //  m_fopr_w->set_parameters(kappa, m_boundary);
  //  m_force_w->set_parameters(kappa, m_boundary);

#pragma omp barrier

}


//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_eo<AFIELD>::set_parameters(
                                          const real_t mq,
                                          const real_t M0,
                                          const int Ns,
                                          const std::vector<int> bc,
                                          const real_t b,
                                          const real_t c)
{
  const int Ndim = CommonParameters::Ndim();

  //- range check
  int err = 0;
  err += ParameterCheck::non_zero(mq);
  err += ParameterCheck::non_zero(M0);
  err += ParameterCheck::non_zero(Ns);

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  assert(bc.size() == Ndim);

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) {
    m_mq = mq;
    m_M0 = M0;
    m_Ns = Ns;

    m_boundary.resize(Ndim);
    m_boundary = bc;

    if (m_b.size() != m_Ns) {
      m_b.resize(m_Ns);
      m_c.resize(m_Ns);
    }
    for (int is = 0; is < m_Ns; ++is) {
      m_b[is] = b;
      m_c[is] = c;
    }
  }

  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  gamma-matrix type = %s\n", m_repr.c_str());
  vout.general(m_vl, "  mq   = %12.8f\n", mq);
  vout.general(m_vl, "  M0   = %12.8f\n", M0);
  vout.general(m_vl, "  Ns   = %4d\n", Ns);
  for (int mu = 0; mu < Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, bc[mu]);
  }
  vout.general(m_vl, "  coefficients:\n");
  for (int is = 0; is < m_Ns; ++is) {
    vout.general(m_vl, "  b[%2d] = %16.10f  c[%2d] = %16.10f\n",
                 is, m_b[is], is, m_c[is]);
  }

  //-- post-process
  //- Domain-wall operator
  //m_fopr_dw->set_parameters(m_mq, m_M0, m_Ns, m_boundary, b, c);

}

//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_eo<AFIELD>::get_parameters(Parameters& params) const
{
  std::vector<double> b(m_Ns);
  std::vector<double> c(m_Ns);
  for (int is = 0; is < m_Ns; ++is) {
    b[is] = double(m_b[is]);
    c[is] = double(m_c[is]);
  }

  params.set_double("quark_mass", double(m_mq));
  params.set_double("domain_wall_height", double(m_M0));
  params.set_int("extent_of_5th_dimension", m_Ns);
  params.set_int_vector("boundary_condition", m_boundary);
  params.set_double_vector("coefficient_b", b);
  params.set_double_vector("coefficient_c", c);
  params.set_string("gamma_matrix_type", m_repr);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_eo<AFIELD>::set_config(Field* U)
{
  int nth = ThreadManager::get_num_threads();
  if(nth > 1){
    set_config_impl(U);
  } else {
    set_config_omp(U);
  }

  m_fopr_w->set_config(U);
  m_fopr_dw->set_config(U);
  m_force_dw->set_config(U);
}

//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_eo<AFIELD>::set_config_omp(Field* U)
{
#pragma omp parallel
 { set_config_impl(U); }
}

//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_eo<AFIELD>::set_config_impl(Field* U)
{
#pragma omp barrier

  AIndex_lex<real_t, AFIELD::IMPL> index;
  convert_gauge(index, m_Ucp, *U);

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_eo<AFIELD>::set_mode(const std::string& mode)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_mode = mode;

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_eo<AFIELD>::force_udiv(AFIELD& force,
                                             const AFIELD& eta)
{
#pragma omp barrier

  m_fopr_dw->set_mode("H");
  m_fopr_dw->mult(m_z1, eta);

  set_mode("H");
  force_udiv1(force, m_z1, eta);

  set_mode("Hdag");
  force_udiv1(m_force2, eta, m_z1);
  axpy(force, real_t(1.0), m_force2);

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_eo<AFIELD>::force_udiv1(AFIELD& force,
                                                 const AFIELD& zeta,
                                                 const AFIELD& eta)
{
#pragma omp barrier

  force.set(0.0);
#pragma omp barrier

  if(m_mode == "H"){

    m_fopr_dw->mult(m_v1, eta, "Doe");
    m_fopr_dw->mult(m_v2, m_v1,  "Doo_inv");
    m_index_eo->merge(m_eta_lex, eta, m_v2);

    m_fopr_dw->mult(m_v3, zeta, "Dee_inv");
    m_fopr_dw->mult(m_v2, m_v3,   "Doe");
    m_fopr_dw->mult(m_v1, m_v2,   "Doo_inv");
    m_index_eo->merge(m_zeta_lex, m_v3, m_v1);

    m_force_dw->set_mode("H");
    m_force_dw->force_udiv1(force, m_zeta_lex, m_eta_lex);
    scal(force, real_t(-1.0));

  }else if(m_mode == "Hdag"){

    m_fopr_dw->mult(m_v3, eta, "Dee_inv");
    m_fopr_dw->mult(m_v1,  m_v3, "Doe");
    m_fopr_dw->mult(m_v2,  m_v1, "Doo_inv");
    m_index_eo->merge(m_eta_lex, m_v3, m_v2);

    m_fopr_dw->mult(m_v2, zeta, "Doe");
    m_fopr_dw->mult(m_v1, m_v2,   "Doo_inv");
    m_index_eo->merge(m_zeta_lex, zeta, m_v1);

    m_force_dw->set_mode("Hdag");
    m_force_dw->force_udiv1(force, m_zeta_lex, m_eta_lex);
    scal(force, real_t(-1.0));

  } else {
    vout.crucial(m_vl, "%s: mode undefined.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

#pragma omp barrier

}

//============================================================END=====
