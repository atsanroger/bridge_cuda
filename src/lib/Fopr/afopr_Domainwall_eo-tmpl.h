/*!
        @file    afopr_Domainwall_eo-tmpl.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2025-01-19 19:38:52 #$
        @version $LastChangedRevision: 2633 $
*/

#include "lib/Fopr/afopr_Domainwall_eo.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
using namespace std;

#include "lib/Parameters/commonParameters.h"
#include "lib/Communicator/communicator.h"
#include "lib/ResourceManager/threadManager.h"

template<typename AFIELD>
const std::string AFopr_Domainwall_eo<AFIELD>::class_name
  = "AFopr_Domainwall_eo";

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_eo<AFIELD>::init(const Parameters& params)
{
  int Nc = CommonParameters::Nc();
  int Nd = CommonParameters::Nd();
  m_NinF = 2 * Nc * Nd;

  m_Nvol  = CommonParameters::Nvol();
  m_Nvol2 = m_Nvol / 2;
  m_Ndim  = CommonParameters::Ndim();

  // setup verbose level
  string vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(vlevel);

  vout.general(m_vl, "%s: Initialization start\n", class_name.c_str());
  int err = 0;

  // setup kernel operator
  std::string kernel_type;
  err += params.fetch_string("kernel_type", kernel_type);
  if (err > 0) {
    vout.crucial(m_vl, "Error at %s: kernel_type is not specified.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
  kernel_type += "_eo";

  Parameters params_kernel = params;
  double     M0;
  err += params.fetch_double("domain_wall_height", M0);
  if (err > 0) {
    vout.crucial(m_vl, "Error at %s: domain_wall_height is not specified.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
  m_M0 = real_t(M0);

  double kappa = 1.0 / (8.0 - 2.0 * M0);
  params_kernel.set_double("hopping_parameter", kappa);

  // Factory is assumed to work
  m_foprw = AFopr<AFIELD>::New(kernel_type, params_kernel);
  m_foprw->set_mode("D");

  m_Ns = 0;

  set_parameters(params);

  m_w4.reset(m_NinF, m_Nvol2, 1);
  m_v4.reset(m_NinF, m_Nvol2, 1);
  m_y4.reset(m_NinF, m_Nvol2, 1);
  m_t4.reset(m_NinF, m_Nvol2, 1);

  if (needs_convert()) {
    m_w4lex.reset(m_NinF, m_Nvol, 1);
    m_v4lex.reset(m_NinF, m_Nvol, 1);
  }

  m_index_eo = new Index_eo_Domainwall<AFIELD>;

}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_eo<AFIELD>::tidyup()
{
  delete m_foprw;
  delete m_index_eo;
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_eo<AFIELD>::set_parameters(
  const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  string           gmset_type;
  double           mq, M0;
  int              Ns;
  std::vector<int> bc;
  double           b, c;
  double           alpha;

  int err_optional = 0;
  err_optional += params.fetch_string("gamma_matrix_type", gmset_type);

  int err = 0;
  err += params.fetch_double("quark_mass", mq);
  err += params.fetch_double("domain_wall_height", M0);
  err += params.fetch_int("extent_of_5th_dimension", Ns);
  err += params.fetch_int_vector("boundary_condition", bc);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  std::string repr;
  if (!params.fetch_string("gamma_matrix_type", repr)) {
    m_repr = repr;
  } else {
    m_repr = "Dirac";  // default
    vout.general(m_vl, "gamma_matrix_type is not given: defalt = %s\n",
                 m_repr.c_str());
  }

  int err2 = 0;
  err2 += params.fetch_double("coefficient_b", b);
  err2 += params.fetch_double("coefficient_c", c);
  if (err2) {
    vout.general(m_vl, "  coefficients b, c are not provided:"
                       " set to Shamir's form.\n");
    b = 1.0;
    c = 0.0;
  }

  int err3 = params.fetch_double("parameter_alpha", alpha);
  if (err3) {
    vout.general(m_vl, "  parameter alpha is not provided: set to 1.0.\n");
    alpha = 1.0;
  }

  set_parameters(real_t(mq), real_t(M0), Ns, bc,
                 real_t(b), real_t(c), real_t(alpha));

  if (real_t(M0) != m_M0) set_kernel_parameters(params);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_eo<AFIELD>::get_parameters(
                                            Parameters& params) const
{
  params.set_string("kernel_type", m_kernel_type);
  params.set_double("quark_mass", double(m_mq));
  params.set_double("domain_wall_height", double(m_M0));
  params.set_int("extent_of_5th_dimension", m_Ns);
  params.set_int_vector("boundary_condition", m_boundary);
  params.set_double("coefficient_b", double(m_b[0]));
  params.set_double("coefficient_c", double(m_c[0]));
  params.set_double("parameter_alpha", double(m_alpha));
  params.set_string("gamma_matrix_type", m_repr);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_eo<AFIELD>::set_parameters(
  const real_t mq,
  const real_t M0,
  const int Ns,
  const std::vector<int> bc,
  const real_t b,
  const real_t c,
  const real_t alpha)
{
  int ith = ThreadManager::get_thread_id();
  if (ith == 0) {
    m_M0 = real_t(M0);
    m_mq = real_t(mq);
    m_Ns = Ns;
    m_alpha = alpha;

    m_boundary.resize(m_Ndim);
    assert(bc.size() == m_Ndim);
    for (int mu = 0; mu < m_Ndim; ++mu) {
      m_boundary[mu] = bc[mu];
    }

    if (m_b.size() != m_Ns) {
      m_b.resize(m_Ns);
      m_c.resize(m_Ns);
    }
    for (int is = 0; is < m_Ns; ++is) {
      m_b[is] = real_t(b);
      m_c[is] = real_t(c);
    }
  }
#pragma omp barrier

  vout.general(m_vl, "Parameters of %s:\n", class_name.c_str());
  vout.general(m_vl, "   mq   = %8.4f\n", m_mq);
  vout.general(m_vl, "   M0   = %8.4f\n", m_M0);
  vout.general(m_vl, "   Ns   = %4d\n", m_Ns);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "   boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }
  //vout.general(m_vl, "  coefficients b = %16.10f  c = %16.10f\n",
  //                                                m_b[0], m_c[0]);
  vout.general(m_vl, "  coefficients:\n");
  for (int is = 0; is < m_Ns; ++is) {
    vout.general(m_vl, "  b[%2d] = %16.10f  c[%2d] = %16.10f\n",
                 is, m_b[is], is, m_c[is]);
  }
  vout.general(m_vl, "  alpha = %8.4f\n", m_alpha);

  set_precond_parameters();

  // working 5d vectors.
  if (m_w1.nex() != Ns) {
    m_w1.reset(m_NinF, m_Nvol2, m_Ns);
    m_v1.reset(m_NinF, m_Nvol2, m_Ns);
    m_v2.reset(m_NinF, m_Nvol2, m_Ns);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_eo<AFIELD>::set_kernel_parameters(
  const Parameters& params)
{
  Parameters params_kernel = params;

  double M0;
  params.fetch_double("domain_wall_height", M0);

  double kappa = 1.0 / (8.0 - 2.0 * M0);
  params_kernel.set_double("hopping_parameter", kappa);

  m_foprw->set_parameters(params_kernel);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_eo<AFIELD>::set_precond_parameters()
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) {
    if (m_dp.size() != m_Ns) {
      m_dp.resize(m_Ns);
      m_dm.resize(m_Ns);
      m_e.resize(m_Ns - 1);
      m_f.resize(m_Ns - 1);
    }

    for (int is = 0; is < m_Ns; ++is) {
      m_dp[is] = m_alpha * (1.0 + m_b[is] * (4.0 - m_M0));
      m_dm[is] = m_alpha * (1.0 - m_c[is] * (4.0 - m_M0));
    }

    m_e[0] = m_mq * m_dm[m_Ns - 1] / m_dp[0];
    m_f[0] = m_mq * m_dm[0]/m_alpha;
    for (int is = 1; is < m_Ns - 1; ++is) {
      m_e[is] = m_e[is - 1] * m_dm[is - 1] / m_dp[is];
      m_f[is] = m_f[is - 1] * m_dm[is] / m_dp[is - 1];
    }

    m_g = m_e[m_Ns - 2] * m_dm[m_Ns - 2];
  }
#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_eo<AFIELD>::set_coefficients(
  const std::vector<real_t> vec_b,
  const std::vector<real_t> vec_c)
{
  if ((vec_b.size() != m_Ns) || (vec_c.size() != m_Ns)) {
    vout.crucial(m_vl, "%s: size of coefficient vectors incorrect.\n",
                 class_name.c_str());
  }

  vout.general(m_vl, "%s: coefficient vectors are set:\n",
               class_name.c_str());

  for (int is = 0; is < m_Ns; ++is) {
    m_b[is] = vec_b[is];
    m_c[is] = vec_c[is];
    vout.general(m_vl, "b[%2d] = %16.10f  c[%2d] = %16.10f\n",
                 is, m_b[is], is, m_c[is]);
  }

  set_precond_parameters();
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_eo<AFIELD>::convert(AFIELD& v, const Field& w)
{
  if (!needs_convert()) {
    vout.crucial(m_vl, "%s: convert is not necessary.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

#pragma omp barrier

  int Nex = w.nex();
  for (int ex = 0; ex < Nex; ++ex) {
    copy(m_w4lex, 0, w, ex);
    m_foprw->convert(m_v4lex, m_w4lex);
    copy(v, ex, m_v4lex, 0);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_eo<AFIELD>::reverse(Field& v, const AFIELD& w)
{
  if (!needs_convert()) {
    vout.crucial(m_vl, "%s: convert is not necessary.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

#pragma omp barrier

  int Nex = w.nex();
  for (int ex = 0; ex < Nex; ++ex) {
    copy(m_v4lex, 0, w, ex);
    m_foprw->reverse(m_w4lex, m_v4lex);
    copy(v, ex, m_w4lex, 0);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_eo<AFIELD>::set_mode(std::string mode)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_mode = mode;
  vout.paranoiac(m_vl, "  mode is set to %s\n", mode.c_str());

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_eo<AFIELD>::mult(AFIELD& v, const AFIELD& w)
{
  if (m_mode == "D") {
    D(v, w);
  } else if (m_mode == "Ddag") {
    Ddag(v, w);
  } else if (m_mode == "DdagD") {
    DdagD(v, w);
  } else if (m_mode == "DDdag") {
    Ddag(m_w1, w);
    D(v, m_w1);
  } else if (m_mode == "H") {
    H(v, w);
  } else {
    vout.crucial(m_vl, "mode undeifined in %s.\n", class_name.c_str());
    vout.crucial(m_vl, "in mult, mode=%s.\n", m_mode.c_str());
    abort();
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_eo<AFIELD>::mult_dag(AFIELD& v, const AFIELD& w)
{
  if (m_mode == "D") {
    Ddag(v, w);
  } else if (m_mode == "Ddag") {
    D(v, w);
  } else if (m_mode == "DdagD") {
    DdagD(v, w);
  } else if (m_mode == "DDdag") {
    Ddag(m_w1, w);
    D(v, m_w1);
  } else if (m_mode == "H") {
    Hdag(v, w);
  } else {
    vout.crucial(m_vl, "mode undeifined in %s.\n", class_name.c_str());
    vout.crucial(m_vl, "in mult_dag, mode=%s.\n", m_mode.c_str());
    abort();
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_eo<AFIELD>::mult(AFIELD& v, const AFIELD& w,
                                       std::string mode)
{
  assert(w.check_size(m_NinF, m_Nvol2, m_Ns));
  assert(v.check_size(m_NinF, m_Nvol2, m_Ns));

  if (mode == "Deo") {
    D_eo(v, w, 0);
  } else if (mode == "Doe") {
    D_eo(v, w, 1);
  } else if (mode == "Dee") {
    D_ee(v, w, 0);
  } else if (mode == "Doo") {
    D_ee(v, w, 1);
  } else if (mode == "Dee_inv") {
    L_inv(m_v1, w);
    U_inv(v, m_v1);
  } else if (mode == "Doo_inv") {
    L_inv(m_v1, w);
    U_inv(v, m_v1);
  } else {
    std::cout << "mode undeifined in AFopr_Domainwall_eo.\n";
    abort();
  }
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_eo<AFIELD>::preProp(AFIELD& Be, AFIELD& bo, const AFIELD& b){
  // bo = Doo_inv b[odd]
  // Be = Dee_inv (b[even] - Deo Doo_inv b[odd])
  //   or
  // bo = Doo_dag_inv b[odd]
  // Be = b[even] - Doe_dag Doo_dag_inv b[odd]
  m_index_eo->split(Be, m_w1, b);
#pragma omp barrier
  if(m_mode == "D"){ // D=(1- Dee_inv Deo Doo_inv Doo)
    L_inv(m_v1, m_w1); //   mult(bo, m_w1, "Doo_inv");
    U_inv(bo, m_v1);
    D_eo(m_w1, bo, 0); //   mult(m_w1, bo, "Deo");
    axpy(Be, real_t(-1), m_w1);
    L_inv(m_v1, Be); //   mult(Be, Be, "Dee_inv");
    U_inv(Be, m_v1);
  } else {
    Udag_inv(m_v1, m_w1); // Doo_dag_inv
    Ldag_inv(bo, m_v1);
    Ddag_eo(m_w1, bo, 0);
    axpy(Be, real_t(-1), m_w1);
  }
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_eo<AFIELD>::postProp(AFIELD& x, const AFIELD& xe, const AFIELD& bo)
  {
    // x[even] = xe
    // x[odd]  = Doo_inv b[odd] - Doo_inv Doe x[even]
    //         = bo - Doo_inv Doe xe
    //   or
    // x[even] =Dee_dag_inv xe
    // x[odd]  = Doo_dag_inv b[odd] - Doo_dag_inv Deo_dag x[even]
    //         = bo - Doo_dag_inv Doe_dag xe

  if(m_mode == "D"){  // D=(1- Dee_inv Deo Doo_inv Doo)
    D_eo(m_w1, xe, 1); //   mult(m_w1, xe, "Doe");
    L_inv(m_v1, m_w1); //   mult(m_w1, m_w1, "Doo_inv");
    U_inv(m_w1, m_v1);
    aypx(real_t(-1.0), m_w1, bo);
#pragma omp barrier
    m_index_eo->merge(x, xe, m_w1);
  } else {
    Udag_inv(m_v1, xe);   // Dee_dag_inv
    Ldag_inv(m_v2, m_v1);
    Ddag_eo(m_w1, m_v2, 1);
    Udag_inv(m_v1, m_w1); // Doo_dag_inv
    Ldag_inv(m_w1, m_v1);
    aypx(real_t(-1.0), m_w1, bo);
#pragma omp barrier
    m_index_eo->merge(x, m_v2, m_w1);
  }



  }

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_eo<AFIELD>::DdagD(AFIELD& v, const AFIELD& w)
{
  D_eo(m_v1, w, 1);
  L_inv(m_v2, m_v1);
  U_inv(m_v1, m_v2);
  D_eo(m_v2, m_v1, 0);
  L_inv(m_v1, m_v2);
  U_inv(m_v2, m_v1);

  copy(m_v1, w);
  axpy(m_v1, real_t(-1.0), m_v2);
#pragma omp barrier

  Udag_inv(v, m_v1);
  Ldag_inv(m_v2, v);
  Ddag_eo(v, m_v2, 1);
  Udag_inv(m_v2, v);
  Ldag_inv(v, m_v2);
  Ddag_eo(m_v2, v, 0);

  copy(v, m_v1);
  axpy(v, real_t(-1.0), m_v2);
#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_eo<AFIELD>::D(AFIELD& v, const AFIELD& w)
{
  D_eo(m_v1, w, 1);
  L_inv(m_v2, m_v1);
  U_inv(m_v1, m_v2);
  D_eo(m_v2, m_v1, 0);
  L_inv(m_v1, m_v2);
  U_inv(m_v2, m_v1);

  copy(v, w);
  axpy(v, real_t(-1.0), m_v2);
#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_eo<AFIELD>::H(AFIELD& v, const AFIELD& w)
{
  D_eo(m_v1, w, 1);
  L_inv(m_v2, m_v1);
  U_inv(m_v1, m_v2);
  D_eo(m_v2, m_v1, 0);
  L_inv(m_v1, m_v2);
  U_inv(m_v2, m_v1);

  copy(m_v1, w);
  axpy(m_v1, real_t(-1.0), m_v2);
#pragma omp barrier

  mult_gm5R(v, m_v1);

}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_eo<AFIELD>::Ddag(AFIELD& v, const AFIELD& w)
{
  assert(w.check_size(m_NinF, m_Nvol2, m_Ns));
  assert(v.check_size(m_NinF, m_Nvol2, m_Ns));

  Udag_inv(m_v1, w);
  Ldag_inv(m_v2, m_v1);
  Ddag_eo(m_v1, m_v2, 1);
  Udag_inv(m_v2, m_v1);
  Ldag_inv(m_v1, m_v2);
  Ddag_eo(m_v2, m_v1, 0);

  copy(v, w);
  axpy(v, real_t(-1.0), m_v2);
#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_eo<AFIELD>::Hdag(AFIELD& v, const AFIELD& w)
{
  assert(w.check_size(m_NinF, m_Nvol2, m_Ns));
  assert(v.check_size(m_NinF, m_Nvol2, m_Ns));

  mult_gm5R(v, w);

  Udag_inv(m_v1, v);
  Ldag_inv(m_v2, m_v1);
  Ddag_eo(m_v1, m_v2, 1);
  Udag_inv(m_v2, m_v1);
  Ldag_inv(m_v1, m_v2);
  Ddag_eo(m_v2, m_v1, 0);

  //copy(v, w);
  axpy(v, real_t(-1.0), m_v2);
#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_eo<AFIELD>::mult_gm5(AFIELD& v, const AFIELD& w)
{
  int Nex = v.nex();
  assert(Nex == w.nex());

  // omp barrier at the beginning and end of m_foprw->mult_gm5 is
  // assumed.
  if (Nex == 1) {
    mult_gm5_4d(v, w);
  } else {
#pragma omp barrier
    for (int ex = 0; ex < Nex; ++ex) {
      copy(m_w4, 0, w, ex);
      m_foprw->mult_gm5(m_v4, m_w4);
      copy(v, ex, m_v4, 0);
#pragma omp barrier
    }
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_eo<AFIELD>::mult_gm5_4d(AFIELD& v,
                                              const AFIELD& w)
{
  m_foprw->mult_gm5(v, w);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_eo<AFIELD>::mult_R(AFIELD& v, const AFIELD& w)
{
  assert(w.check_size(m_NinF, m_Nvol2, m_Ns));
  assert(v.check_size(m_NinF, m_Nvol2, m_Ns));

#pragma omp barrier

  for (int is = 0; is < m_Ns; ++is) {
    copy(v, m_Ns - 1 - is, w, is);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_eo<AFIELD>::mult_gm5R(AFIELD& v, const AFIELD& w)
{
  assert(w.check_size(m_NinF, m_Nvol2, m_Ns));
  assert(v.check_size(m_NinF, m_Nvol2, m_Ns));

#pragma omp barrier

  for (int is = 0; is < m_Ns; ++is) {
    copy(m_w4, 0, w, is);
    mult_gm5_4d(m_v4, m_w4);
    copy(v, m_Ns - 1 - is, m_v4, 0);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_eo<AFIELD>::D_eo(AFIELD& v, const AFIELD& w,
                                       const int ieo)
{  // ieo = 0: even < odd, ieo = 1: odd <- even
#pragma omp barrier

  for (int is = 0; is < m_Ns; ++is) {
    m_v4.set(0.0);
    int    is_up = (is + 1) % m_Ns;
    real_t Fup   = 0.5 * m_alpha;
    if (is == m_Ns-1) Fup = -0.5 * m_mq;
    copy(m_y4, 0, w, is_up);
    mult_gm5_4d(m_t4, m_y4);
    axpy(m_y4, real_t(-1.0), m_t4);
    axpy(m_v4, 0, Fup, m_y4, 0);  // m_v4 = - P_ w(is+1)

    int    is_dn = (is - 1 + m_Ns) % m_Ns;
    real_t Fdn   = 0.5 * m_alpha;
    if (is == 0) Fdn = -0.5 * m_mq;
    copy(m_y4, 0, w, is_dn);
    mult_gm5_4d(m_t4, m_y4);
    axpy(m_y4, real_t(1.0), m_t4);
    axpy(m_v4, 0, Fdn, m_y4, 0);  // m_v4 += - P+ w(is-1)

    copy(m_w4, 0, w, is);

    if(is == 0){
      real_t fac1 = 0.5 * ( 1.0 + m_alpha);
      real_t fac2 = 0.5 * (-1.0 + m_alpha);
      if(fac2 != 0.0){
        mult_gm5_4d(m_t4, m_w4);
        scal(m_w4, fac1);
        axpy(m_w4, fac2, m_t4);
      }
    }else if(is == m_Ns-1){
      real_t fac1 = 0.5 * (1.0 + m_alpha);
      real_t fac2 = 0.5 * (1.0 - m_alpha);
      if(fac2 != 0.0){
        mult_gm5_4d(m_t4, m_w4);
        scal(m_w4, fac1);
        axpy(m_w4, fac2, m_t4);
      }
    }else{
      scal(m_w4, m_alpha);
    }
 
    scal(m_w4, m_b[is]);
    axpy(m_w4, m_c[is], m_v4);

    if (ieo == 0) {
      m_foprw->mult(m_v4, m_w4, "Deo");
    } else {
      m_foprw->mult(m_v4, m_w4, "Doe");
    }

    scal(m_v4, real_t(4.0 - m_M0));

    copy(v, is, m_v4, 0);

#pragma omp barrier
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_eo<AFIELD>::Ddag_eo(AFIELD& v, const AFIELD& w,
                                          const int ieo)
{  // ieo = 0: even < odd, ieo = 1: odd <- even
#pragma omp barrier

  v.set(0.0);

#pragma omp barrier

  for (int is = 0; is < m_Ns; ++is) {
    copy(m_w4, 0, w, is);

    // m_foprw->Mdageo(m_v4, m_w4, ieo);
    mult_gm5_4d(m_v4, m_w4);
    if (ieo == 0) {
      m_foprw->mult(m_y4, m_v4, "Deo");
    } else {
      m_foprw->mult(m_y4, m_v4, "Doe");
    }
    mult_gm5_4d(m_v4, m_y4);

    scal(m_v4, real_t(4.0 - m_M0));

   copy(m_w4, m_v4);

    if(is == 0){
      real_t fac1 = 0.5 * ( 1.0 + m_alpha);
      real_t fac2 = 0.5 * (-1.0 + m_alpha);
      if(fac2 != 0.0){
        mult_gm5_4d(m_t4, m_v4);
        scal(m_v4, fac1);
        axpy(m_v4, fac2, m_t4);
      }
    }else if(is == m_Ns-1){
      real_t fac1 = 0.5 * (1.0 + m_alpha);
      real_t fac2 = 0.5 * (1.0 - m_alpha);
      if(fac2 != 0.0){
        mult_gm5_4d(m_t4, m_v4);
        scal(m_v4, fac1);
        axpy(m_v4, fac2, m_t4);
      }
    }else{
      scal(m_v4, m_alpha);
    }

    axpy(v, is, m_b[is], m_v4, 0);

    scal(m_w4, m_c[is]);

    int    is_up = (is + 1) % m_Ns;
    real_t Fup   = 0.5 * m_alpha;
    if (is_up == 0) Fup = -0.5 * m_mq;
    mult_gm5_4d(m_y4, m_w4);
    axpy(m_y4, real_t(-1.0), m_w4); // m_y4 = - (1 - gm5) m_v4
    axpy(v, is_up, -Fup, m_y4, 0);  // v += Fup * (1 - gm5) m_y4

    int    is_dn = (is - 1 + m_Ns) % m_Ns;
    real_t Fdn   = 0.5 * m_alpha;
    if (is_dn == m_Ns - 1) Fdn = -0.5 * m_mq;
    mult_gm5_4d(m_y4, m_w4);
    axpy(m_y4, real_t(1.0), m_w4);  // m_y4 = (1 + gm5) m_v4
    axpy(v, is_dn, Fdn, m_y4, 0);   // v += Fdn * (1 + gm5) m_y4

#pragma omp barrier
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_eo<AFIELD>::D_ee(AFIELD& v, const AFIELD& w,
                                       const int ieo)
{  // ieo = 0: even, ieo = 1: odd (but they are identical)
#pragma omp barrier

  for (int is = 0; is < m_Ns; ++is) {

    m_v4.set(0.0);

    int    is_up = (is + 1) % m_Ns;
    real_t Fup   = 0.5 * m_alpha;
    if (is == m_Ns - 1) Fup = -0.5 * m_mq;
    copy(m_y4, 0, w, is_up);
    mult_gm5_4d(m_t4, m_y4);
    axpy(m_y4, real_t(-1.0), m_t4);
    axpy(m_v4, 0, Fup, m_y4, 0);

    int    is_dn = (is - 1 + m_Ns) % m_Ns;
    real_t Fdn   = 0.5 * m_alpha;
    if (is == 0) Fdn = -0.5 * m_mq;
    copy(m_y4, 0, w, is_dn);
    mult_gm5_4d(m_t4, m_y4);
    axpy(m_y4, real_t(1.0), m_t4);
    axpy(m_v4, 0, Fdn, m_y4, 0);

    copy(m_w4, 0, w, is);

    if(is == 0){
      real_t fac1 = 0.5 * ( 1.0 + m_alpha);
      real_t fac2 = 0.5 * (-1.0 + m_alpha);
      if(fac2 != 0.0){
        mult_gm5_4d(m_t4, m_w4);
        scal(m_w4, fac1);
        axpy(m_w4, fac2, m_t4);
      }
    }else if(is == m_Ns-1){
      real_t fac1 = 0.5 * (1.0 + m_alpha);
      real_t fac2 = 0.5 * (1.0 - m_alpha);
      if(fac2 != 0.0){
        mult_gm5_4d(m_t4, m_w4);
        scal(m_w4, fac1);
        axpy(m_w4, fac2, m_t4);
      }
    }else{
      scal(m_w4, m_alpha);
    }

    real_t F1 = m_b[is] * (4.0 - m_M0) + 1.0;
    scal(m_w4, F1);

    real_t F2 = m_c[is] * (4.0 - m_M0) - 1.0;
    axpy(m_w4, F2, m_v4);

    copy(v, is, m_w4, 0);

#pragma omp barrier
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_eo<AFIELD>::Ddag_ee(AFIELD& v, const AFIELD& w,
                                       const int ieo)
{  // ieo = 0: even, ieo = 1: odd (but they are identical)
#pragma omp barrier

  v.set(0.0);
#pragma omp barrier

  for (int is = 0; is < m_Ns; ++is) {

    copy(m_w4, 0, w, is);
    real_t F1 = m_b[is] * (4.0 - m_M0) + 1.0;
    scal(m_w4, F1);

    if(is == 0){
      real_t fac1 = 0.5 * ( 1.0 + m_alpha);
      real_t fac2 = 0.5 * (-1.0 + m_alpha);
      if(fac2 != 0.0){
        mult_gm5_4d(m_t4, m_w4);
        scal(m_w4, fac1);
	axpy(m_w4, fac2, m_t4);
      }
    }else if(is == m_Ns-1){
      real_t fac1 = 0.5 * (1.0 + m_alpha);
      real_t fac2 = 0.5 * (1.0 - m_alpha);
      if(fac2 != 0.0){
        mult_gm5_4d(m_t4, m_w4);
        scal(m_w4, fac1);
        axpy(m_w4, fac2, m_t4);
      }
    }else{
      scal(m_w4, m_alpha);
    }

    axpy(v, is, real_t(1.0), m_w4, 0);

    copy(m_y4, 0, w, is);
    real_t F2 = m_c[is] * (4.0 - m_M0) - 1.0;
    scal(m_y4, F2);

    int    is_up = (is + 1) % m_Ns;
    real_t Fup   = 0.5 * m_alpha;
    if (is_up == 0) Fup = -0.5 * m_mq;
    mult_gm5_4d(m_t4, m_y4);
    axpy(m_t4, real_t(-1.0), m_y4); // m_t4 = - (1 - gm5) * m_y4
    axpy(v, is_up, -Fup, m_t4, 0);  // v+= Fup * (1 - gm5) * m_t4

    int    is_dn = (is - 1 + m_Ns) % m_Ns;
    real_t Fdn   = 0.5 * m_alpha;
    if (is_dn == m_Ns-1) Fdn = -0.5 * m_mq;
    mult_gm5_4d(m_t4, m_y4);
    axpy(m_t4, real_t(1.0), m_y4);  // m_t4 = (1 + gm5) * m_y4
    axpy(v, is_dn, Fdn, m_t4, 0);   // v+= Fdn * (1 + gm5) * m_t4
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_eo<AFIELD>::L_inv(AFIELD& v, const AFIELD& w)
{
#pragma omp barrier

  copy(v, 0, w, 0);
  copy(m_y4, 0, w, 0);
  scal(m_y4, m_e[0]);

#pragma omp barrier

  for (int is = 1; is < m_Ns - 1; ++is) {
    copy(v, is, w, is);

    copy(m_v4, 0, v, is - 1);
    mult_gm5_4d(m_w4, m_v4);
    axpy(m_v4, real_t(1.0), m_w4);
    scal(m_v4, real_t(0.5) * m_dm[is] / m_dp[is - 1]);

    axpy(v, is, real_t(1.0), m_v4, 0);
    copy(m_w4, 0, v, is);
    axpy(m_y4, m_e[is], m_w4);

#pragma omp barrier
  }

  int is = m_Ns - 1;
  copy(v, is, w, is);
  copy(m_v4, 0, v, is - 1);
  mult_gm5_4d(m_w4, m_v4);
  axpy(m_v4, real_t(1.0), m_w4);
  scal(m_v4, real_t(0.5) * m_dm[is] / m_dp[is - 1]);

  axpy(v, is, real_t(1.0), m_v4, 0);

  mult_gm5_4d(m_w4, m_y4);
  axpy(m_y4, real_t(-1.0), m_w4);
  scal(m_y4, real_t(-0.5));
  axpy(v, is, real_t(1.0), m_y4, 0);

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_eo<AFIELD>::U_inv(AFIELD& v, const AFIELD& w)
{
#pragma omp barrier

  int is = m_Ns - 1;
  copy(m_y4, 0, w, is);

  // multiply (alpha * P+ + P-)
  real_t fac1 = 0.5 * ( 1.0 + m_alpha);
  real_t fac2 = 0.5 * (-1.0 + m_alpha);
  if(fac2 != 0.0){
    mult_gm5_4d(m_t4, m_y4);
    scal(m_y4, fac1);
    axpy(m_y4, fac2, m_t4);
  }
  scal(m_y4, real_t(1.0) / (m_dp[is] + m_g));
  copy(v, is, m_y4, 0);

  mult_gm5_4d(m_w4, m_y4);
  axpy(m_y4, real_t(1.0), m_w4);
  scal(m_y4, real_t(0.5));

#pragma omp barrier

  for (int is = m_Ns - 2; is >= 0; --is) {
    copy(m_v4, 0, w, is);

    copy(m_w4, 0, v, is + 1);
    mult_gm5_4d(m_t4, m_w4);
    axpy(m_w4, real_t(-1.0), m_t4);

    axpy(m_v4, real_t(0.5) * m_dm[is], m_w4);

    axpy(m_v4, -m_f[is], m_y4);

    scal(m_v4, real_t(1.0) / m_dp[is]);

    if(is == 0){  // multiply (alpha * P- + P+)
      real_t fac1 = 0.5 * (1.0 + m_alpha);
      real_t fac2 = 0.5 * (1.0 - m_alpha);
      if(fac2 != 0.0){
        mult_gm5_4d(m_t4, m_v4);
        scal(m_v4, fac1);
	axpy(m_v4, fac2, m_t4);
      }
    }

    copy(v, is, m_v4, 0);

#pragma omp barrier
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_eo<AFIELD>::Udag_inv(AFIELD& v, const AFIELD& w)
{
#pragma omp barrier

  copy(m_v4, 0, w, 0);

  // multiply (alpha * P- + P+)
  real_t fac1 = 0.5 * (1.0 + m_alpha);
  real_t fac2 = 0.5 * (1.0 - m_alpha);
  if(fac2 != 0.0){
    mult_gm5_4d(m_t4, m_v4);
    scal(m_v4, fac1);
    axpy(m_v4, fac2, m_t4);
  }
  scal(m_v4, real_t(1.0) / m_dp[0]);
  copy(v, 0, m_v4, 0);

  copy(m_y4, m_v4);
  scal(m_y4, m_f[0]);

#pragma omp barrier

  for (int is = 1; is < m_Ns - 1; ++is) {
    copy(m_t4, 0, w, is);

    copy(m_v4, 0, v, is - 1);
    mult_gm5_4d(m_w4, m_v4);
    axpy(m_v4, real_t(-1.0), m_w4);
    axpy(m_t4, real_t(0.5) * m_dm[is - 1], m_v4);

    scal(m_t4, real_t(1.0) / m_dp[is]);
    copy(v, is, m_t4, 0);

    axpy(m_y4, m_f[is], m_t4);

#pragma omp barrier
  }

  int is = m_Ns - 1;

  copy(m_t4, 0, w, is);

  copy(m_v4, 0, v, is - 1);
  mult_gm5_4d(m_w4, m_v4);
  axpy(m_v4, real_t(-1.0), m_w4);
  axpy(m_t4, real_t(0.5) * m_dm[is - 1], m_v4);

  mult_gm5_4d(m_w4, m_y4);
  axpy(m_y4, real_t(1.0), m_w4);
  scal(m_y4, real_t(0.5));

  axpy(m_t4, real_t(-1.0), m_y4);

  scal(m_t4, real_t(1.0) / (m_dp[is] + m_g));
  //copy(v, is, m_t4, 0);

  // multiply (alpha * P+ + P-)
  fac1 = 0.5 * ( 1.0 + m_alpha);
  fac2 = 0.5 * (-1.0 + m_alpha);
  if(fac2 != 0.0){
    mult_gm5_4d(m_y4, m_t4);
    scal(m_t4, fac1);
    axpy(m_t4, fac2, m_y4);
  }
  copy(v, is, m_t4, 0);

#pragma omp barrier

}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_eo<AFIELD>::Ldag_inv(AFIELD& v, const AFIELD& w)
{
#pragma omp barrier

  int is = m_Ns - 1;
  copy(v, is, w, is);

  copy(m_y4, 0, w, is);
  mult_gm5_4d(m_t4, m_y4);
  axpy(m_y4, real_t(-1.0), m_t4);
  scal(m_y4, real_t(0.5));

#pragma omp barrier

  for (int is = m_Ns - 2; is >= 0; --is) {
    copy(v, is, w, is);

    copy(m_v4, 0, v, is + 1);
    mult_gm5_4d(m_t4, m_v4);
    axpy(m_v4, real_t(1.0), m_t4);
    scal(m_v4, real_t(0.5) * m_dm[is + 1] / m_dp[is]);

    axpy(m_v4, -m_e[is], m_y4);

    axpy(v, is, real_t(1.0), m_v4, 0);

#pragma omp barrier
  }
}


//====================================================================
template<typename AFIELD>
double AFopr_Domainwall_eo<AFIELD>::flop_count(std::string mode)
{
  int    Lvol2 = CommonParameters::Lvol() / 2;
  double vsite = static_cast<double>(Lvol2);
  double vNs   = static_cast<double>(m_Ns);

  double flop_Wilson = m_foprw->flop_count("Deo");
  //  double flop_Wilson = m_foprw->flop_count("Meo");
  //  double flop_Wilson = m_foprw->flop_count();

  double axpy1 = static_cast<double>(2 * m_NinF) * vsite;
  double scal1 = static_cast<double>(1 * m_NinF) * vsite;

  double flop_Deo = (flop_Wilson + 5.0 * axpy1 + 2.0 * scal1) * vNs;

  double flop_Dee = vNs * (7.0 * axpy1 + scal1);

  double flop_LU_inv =
    2.0 * ((3.0 * axpy1 + scal1) * (vNs - 1.0) + axpy1 + 2.0 * scal1);

  double flop = 0.0;
  if ((mode == "Meo") || (mode == "Moe")) {
    flop = flop_Deo + flop_LU_inv;
  } else if ((mode == "Dee_inv") || (mode == "Doo_inv")) {
    flop = flop_LU_inv;
  } else if ((mode == "Deo") || (mode == "Doe")) {
    flop = flop_Deo;
  } else if ((mode == "Dee") || (mode == "Doo")) {
    flop = flop_Dee;
  } else if ((mode == "D") || (mode == "Ddag")) {
    flop = 2.0 * (flop_LU_inv + flop_Deo) + vNs * axpy1;
  } else if (mode == "DdagD") {
    flop = 2.0 * (2.0 * (flop_LU_inv + flop_Deo) + vNs * axpy1);
  } else {
    vout.crucial(m_vl, "Error at %s: input mode %s is undefined.\n",
                 class_name.c_str(), mode.c_str());
    exit(EXIT_FAILURE);
  }

  return flop;
}


//============================================================END=====
