/*!
      @file    afopr_Domainwall-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2513 $
*/

#include "lib/Fopr/afopr_Domainwall.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
using namespace std;

#include "lib/ResourceManager/threadManager.h"
#include "lib/Parameters/commonParameters.h"
#include "lib/Communicator/communicator.h"


template<typename AFIELD>
const std::string AFopr_Domainwall<AFIELD>::class_name
  = "AFopr_Domainwall";
//====================================================================
template<typename AFIELD>
void AFopr_Domainwall<AFIELD>::init(const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: construction\n", class_name.c_str());
  vout.increase_indent();

  int Nc = CommonParameters::Nc();
  int Nd = CommonParameters::Nd();
  m_NinF = 2 * Nc * Nd;

  m_Nvol = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();

  // setup verbose level
  string vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(vlevel);

  int err = 0;

  // setup kernel operator
  err += params.fetch_string("kernel_type", m_kernel_type);
  if (err > 0) {
    vout.crucial(m_vl, "%s: Error: kernel_type is not specified.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

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
  m_foprw          = AFopr<AFIELD>::New(m_kernel_type, params_kernel);
  m_kernel_created = true;

  m_foprw->set_mode("D");

  m_Ns = 0;

  set_parameters(params);

  m_w4.reset(m_NinF, m_Nvol, 1);
  m_v4.reset(m_NinF, m_Nvol, 1);
  m_t4.reset(m_NinF, m_Nvol, 1);
  m_y4.reset(m_NinF, m_Nvol, 1);

  if (needs_convert()) {
    m_w4lex.reset(m_NinF, m_Nvol, 1);
    m_v4lex.reset(m_NinF, m_Nvol, 1);
  }

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall<AFIELD>::init(AFOPR *fopr,
                                    const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: construction\n", class_name.c_str());
  vout.increase_indent();

  int Nc = CommonParameters::Nc();
  int Nd = CommonParameters::Nd();
  m_NinF = 2 * Nc * Nd;

  m_Nvol = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();

  // setup verbose level
  string vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(vlevel);

  vout.general(m_vl, "%s: Initialization start\n", class_name.c_str());
  int err = 0;

  m_foprw          = fopr;
  m_kernel_type    = "unknown";
  m_kernel_created = false;

  m_foprw->set_mode("D");

  m_Ns = 0;

  set_parameters(params);

  m_w4.reset(m_NinF, m_Nvol, 1);
  m_v4.reset(m_NinF, m_Nvol, 1);
  m_t4.reset(m_NinF, m_Nvol, 1);
  m_y4.reset(m_NinF, m_Nvol, 1);

  if (needs_convert()) {
    m_w4lex.reset(m_NinF, m_Nvol, 1);
    m_v4lex.reset(m_NinF, m_Nvol, 1);
  }

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall<AFIELD>::init(AFOPR *fopr)
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: construction\n", class_name.c_str());
  vout.increase_indent();

  int Nc = CommonParameters::Nc();
  int Nd = CommonParameters::Nd();
  m_NinF = 2 * Nc * Nd;

  m_Nvol = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();

  vout.general(m_vl, "%s: Initialization start\n", class_name.c_str());
  int err = 0;

  m_foprw          = fopr;
  m_kernel_type    = "unknown";
  m_kernel_created = false;

  m_foprw->set_mode("D");

  m_Ns = 0;

  m_w4.reset(m_NinF, m_Nvol, 1);
  m_v4.reset(m_NinF, m_Nvol, 1);
  m_t4.reset(m_NinF, m_Nvol, 1);
  m_y4.reset(m_NinF, m_Nvol, 1);

  if (needs_convert()) {
    m_w4lex.reset(m_NinF, m_Nvol, 1);
    m_v4lex.reset(m_NinF, m_Nvol, 1);
  }

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall<AFIELD>::tidyup()
{
  if (m_kernel_created == true) delete m_foprw;
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall<AFIELD>::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  string           gmset_type;
  double           mq, M0;
  int              Ns;
  std::vector<int> bc;
  double           b, c;

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

  int err2 = 0;
  err2 += params.fetch_double("coefficient_b", b);
  err2 += params.fetch_double("coefficient_c", c);

  if (err2) {
    vout.general(m_vl, "  coefficients b, c are not provided:"
                       " set to Shamir's form.\n");
    b = 1.0;
    c = 0.0;
  }

  set_parameters(real_t(mq), real_t(M0), Ns, bc,
                 real_t(b), real_t(c));

  if (real_t(M0) != m_M0) set_kernel_parameters(params);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall<AFIELD>::get_parameters(Parameters& params) const
{
  params.set_string("kernel_type", m_kernel_type);
  params.set_double("quark_mass", double(m_mq));
  params.set_double("domain_wall_height", double(m_M0));
  params.set_int("extent_of_5th_dimension", m_Ns);
  params.set_int_vector("boundary_condition", m_boundary);
  params.set_double("coefficient_b", double(m_b[0]));
  params.set_double("coefficient_c", double(m_c[0]));
  params.set_string("gamma_matrix_type", m_repr);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall<AFIELD>::set_parameters(
  const real_t mq,
  const real_t M0,
  const int Ns,
  const std::vector<int> bc,
  const real_t b,
  const real_t c)
{
  int ith = ThreadManager::get_thread_id();

  if (ith == 0) {
    m_M0 = M0;
    m_mq = mq;
    m_Ns = Ns;

    assert(bc.size() == m_Ndim);
    if (m_boundary.size() != m_Ndim) m_boundary.resize(m_Ndim);
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

  vout.general(m_vl, "%s: input parameters\n", class_name.c_str());
  vout.general(m_vl, "  mq   = %8.4f\n", m_mq);
  vout.general(m_vl, "  M0   = %8.4f\n", m_M0);
  vout.general(m_vl, "  Ns   = %4d\n", m_Ns);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }
  vout.general(m_vl, "  coefficients:\n");
  for (int is = 0; is < m_Ns; ++is) {
    vout.general(m_vl, "  b[%2d] = %16.10f  c[%2d] = %16.10f\n",
                 is, m_b[is], is, m_c[is]);
  }

  set_precond_parameters();

  // working 5d vectors.
  if (m_w1.nex() != Ns) {
    m_w1.reset(m_NinF, m_Nvol, m_Ns);
    m_v1.reset(m_NinF, m_Nvol, m_Ns);
    m_v2.reset(m_NinF, m_Nvol, m_Ns);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall<AFIELD>::set_parameters(
  const real_t mq,
  const real_t M0,
  const int Ns,
  const std::vector<int> bc,
  const std::vector<real_t> vec_b,
  const std::vector<real_t> vec_c)
{
  int ith = ThreadManager::get_thread_id();

  if (ith == 0) {
    m_M0 = M0;
    m_mq = mq;
    m_Ns = Ns;

    assert(bc.size() == m_Ndim);
    if (m_boundary.size() != m_Ndim) m_boundary.resize(m_Ndim);
    for (int mu = 0; mu < m_Ndim; ++mu) {
      m_boundary[mu] = bc[mu];
    }

    if (m_b.size() != m_Ns) {
      m_b.resize(m_Ns);
      m_c.resize(m_Ns);
    }
    for (int is = 0; is < m_Ns; ++is) {
      m_b[is] = real_t(vec_b[is]);
      m_c[is] = real_t(vec_c[is]);
    }
  }

  
  vout.general(m_vl, "%s: parameters\n", class_name.c_str());
  vout.general(m_vl, "  mq   = %8.4f\n", m_mq);
  vout.general(m_vl, "  M0   = %8.4f\n", m_M0);
  vout.general(m_vl, "  Ns   = %4d\n", m_Ns);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }
  vout.general(m_vl, "  coefficients:\n");
  for (int is = 0; is < m_Ns; ++is) {
    vout.general(m_vl, "  b[%2d] = %16.10f  c[%2d] = %16.10f\n",
                 is, m_b[is], is, m_c[is]);
  }

  set_precond_parameters();

  // working 5d vectors.
  if (m_w1.nex() != Ns) {
    m_w1.reset(m_NinF, m_Nvol, m_Ns);
    m_v1.reset(m_NinF, m_Nvol, m_Ns);
    m_v2.reset(m_NinF, m_Nvol, m_Ns);
  }

}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall<AFIELD>::set_coefficients(
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
    m_b[is] = real_t(vec_b[is]);
    m_c[is] = real_t(vec_c[is]);
    vout.general(m_vl, "b[%2d] = %16.10f  c[%2d] = %16.10f\n",
                 is, m_b[is], is, m_c[is]);
  }

  set_precond_parameters();
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall<AFIELD>::set_kernel_parameters(
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
void AFopr_Domainwall<AFIELD>::set_precond_parameters()
{
  if (m_dp.size() != m_Ns) {
    m_dp.resize(m_Ns);
    m_dm.resize(m_Ns);
    m_e.resize(m_Ns - 1);
    m_f.resize(m_Ns - 1);
  }

  for (int is = 0; is < m_Ns; ++is) {
    m_dp[is] = 1.0 + m_b[is] * (4.0 - m_M0);
    m_dm[is] = 1.0 - m_c[is] * (4.0 - m_M0);
  }

  m_e[0] = m_mq * m_dm[m_Ns - 1] / m_dp[0];
  m_f[0] = m_mq * m_dm[0];
  for (int is = 1; is < m_Ns - 1; ++is) {
    m_e[is] = m_e[is - 1] * m_dm[is - 1] / m_dp[is];
    m_f[is] = m_f[is - 1] * m_dm[is] / m_dp[is - 1];
  }

  m_g = m_e[m_Ns - 2] * m_dm[m_Ns - 2];
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall<AFIELD>::convert(AFIELD& v, const Field& w)
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
void AFopr_Domainwall<AFIELD>::reverse(Field& v, const AFIELD& w)
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
void AFopr_Domainwall<AFIELD>::set_mode(std::string mode)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_mode = mode;

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall<AFIELD>::mult(AFIELD& v, const AFIELD& w)
{
  if (m_mode == "D") {
    D(v, w);
  } else if (m_mode == "Ddag") {
    Ddag(v, w);
  } else if (m_mode == "DdagD") {
    DdagD(v, w);
  } else if (m_mode == "DDdag") {
    DDdag(v, w);
  } else if (m_mode == "H") {
    H(v, w);
  } else if (m_mode == "Hdag") {
    Hdag(v, w);
  } else if (m_mode == "D_prec") {
    D_prec(v, w);
  } else if (m_mode == "Ddag_prec") {
    Ddag_prec(v, w);
  } else if (m_mode == "DdagD_prec") {
    DdagD_prec(v, w);
  } else if (m_mode == "Prec") {
    Prec(v, w);
  } else {
    vout.crucial(m_vl, "mode undeifined in %s.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall<AFIELD>::mult_dag(AFIELD& v, const AFIELD& w)
{
  if (m_mode == "D") {
    Ddag(v, w);
  } else if (m_mode == "Ddag") {
    D(v, w);
  } else if (m_mode == "DdagD") {
    DdagD(v, w);
  } else if (m_mode == "DDdag") {
    DDdag(v, w);
  } else if (m_mode == "H") {
    Hdag(v, w);
  } else if (m_mode == "Hdag") {
    H(v, w);
  } else if (m_mode == "D_prec") {
    Ddag_prec(v, w);
  } else if (m_mode == "Ddag_prec") {
    D_prec(v, w);
  } else if (m_mode == "DdagD_prec") {
    DdagD_prec(v, w);
  } else if (m_mode == "Prec") {
    Precdag(v, w);
  } else {
    vout.crucial(m_vl, "mode undeifined in %s.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall<AFIELD>::mult(AFIELD& v, const AFIELD& w,
                                    std::string mode)
{
  assert(w.check_size(m_NinF, m_Nvol, m_Ns));
  assert(v.check_size(m_NinF, m_Nvol, m_Ns));

  if (mode == "Prec") {
    Prec(v, w);
  } else if (mode == "Precdag") {
    Precdag(v, w);
  } else if (mode == "D") {
    D(v, w);
  } else if (mode == "Ddag") {
    Ddag(v, w);
  } else if (mode == "DdagD") {
    DdagD(v, w);
  } else if (mode == "DDdag") {
    DDdag(v, w);
  } else if (mode == "D_prec") {
    D_prec(v, w);
  } else if (mode == "Ddag_prec") {
    Ddag_prec(v, w);
  } else if (mode == "DdagD_prec") {
    DdagD_prec(v, w);
  } else {
    vout.crucial(m_vl, "%s: undefined mode = %s\n",
                 class_name.c_str(), mode.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall<AFIELD>::mult_dag(AFIELD& v, const AFIELD& w,
                                        std::string mode)
{
  assert(w.check_size(m_NinF, m_Nvol, m_Ns));
  assert(v.check_size(m_NinF, m_Nvol, m_Ns));

  if (mode == "Prec") {
    Precdag(v, w);
  } else if (mode == "Precdag") {
    Prec(v, w);
  } else if (mode == "D") {
    Ddag(v, w);
  } else if (mode == "Ddag") {
    D(v, w);
  } else if (mode == "DdagD") {
    DdagD(v, w);
  } else if (mode == "DDdag") {
    DDdag(v, w);
  } else if (mode == "D_prec") {
    Ddag_prec(v, w);
  } else if (mode == "Ddag_prec") {
    D_prec(v, w);
  } else if (mode == "DdagD_prec") {
    DdagD_prec(v, w);
  } else {
    std::cout << "mode undeifined in AFopr_Domainwall.\n";
    abort();
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall<AFIELD>::mult_gm5(AFIELD& v, const AFIELD& w)
{
  int Nex = w.nex();
  assert(w.check_size(m_NinF, m_Nvol, Nex));
  assert(v.check_size(m_NinF, m_Nvol, Nex));

  // omp barrier at the beginning and end of m_foprw->mult_gm5 is
  // assumed.
  if (Nex == 1) {
    m_foprw->mult_gm5(v, w);
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
void AFopr_Domainwall<AFIELD>::mult_chproj_4d(AFIELD& v,
                                              const AFIELD& w,
                                              const int ipm)
{
#pragma omp barrier

  copy(v, w);

  m_foprw->mult_gm5(m_w4, w);

  axpy(v, real_t(ipm), m_w4);
  scal(v, real_t(0.5));

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall<AFIELD>::mult_gm5_4d(AFIELD& v, const AFIELD& w)
{
  m_foprw->mult_gm5(v, w);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall<AFIELD>::DdagD(AFIELD& v, const AFIELD& w)
{
  assert(w.check_size(m_NinF, m_Nvol, m_Ns));
  assert(v.check_size(m_NinF, m_Nvol, m_Ns));

  D(m_v1, w);
  Ddag(v, m_v1);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall<AFIELD>::DDdag(AFIELD& v, const AFIELD& w)
{
  assert(w.check_size(m_NinF, m_Nvol, m_Ns));
  assert(v.check_size(m_NinF, m_Nvol, m_Ns));

  Ddag(m_v1, w);
  D(v, m_v1);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall<AFIELD>::DdagD_prec(AFIELD& v, const AFIELD& w)
{
  assert(w.check_size(m_NinF, m_Nvol, m_Ns));
  assert(v.check_size(m_NinF, m_Nvol, m_Ns));

  D_prec(m_v1, w);
  Ddag_prec(v, m_v1);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall<AFIELD>::D_prec(AFIELD& v, const AFIELD& w)
{
  L_inv(v, w);
  U_inv(m_v2, v);
  D(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall<AFIELD>::Ddag_prec(AFIELD& v, const AFIELD& w)
{
  Ddag(v, w);
  Udag_inv(m_v2, v);
  Ldag_inv(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall<AFIELD>::Prec(AFIELD& v, const AFIELD& w)
{
  L_inv(m_v2, w);
  U_inv(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall<AFIELD>::Precdag(AFIELD& v, const AFIELD& w)
{
  Udag_inv(m_v2, w);
  Ldag_inv(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall<AFIELD>::H(AFIELD& v, const AFIELD& w)
{
  D(m_v2, w);
  mult_gm5R(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall<AFIELD>::Hdag(AFIELD& v, const AFIELD& w)
{
  mult_gm5R(m_v2, w);
  Ddag(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall<AFIELD>::mult_gm5R(AFIELD& v, const AFIELD& w)
{
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
void AFopr_Domainwall<AFIELD>::mult_R(AFIELD& v, const AFIELD& w)
{
  assert(w.check_size(m_NinF, m_Nvol, m_Ns));
  assert(v.check_size(m_NinF, m_Nvol, m_Ns));

#pragma omp barrier

  for (int is = 0; is < m_Ns; ++is) {
    copy(v, m_Ns - 1 - is, w, is);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall<AFIELD>::D(AFIELD& v, const AFIELD& w)
{
#pragma omp barrier

  for (int is = 0; is < m_Ns; ++is) {
    m_y4.set(0.0);
    int    is_up = (is + 1) % m_Ns;
    real_t Fup   = -0.5;
    if (is == m_Ns - 1) Fup = 0.5 * m_mq;
    copy(m_v4, 0, w, is_up);
    mult_gm5_4d(m_t4, m_v4);
    axpy(m_v4, real_t(-1.0), m_t4);
    axpy(m_y4, 0, Fup, m_v4, 0);

    int    is_dn = (is - 1 + m_Ns) % m_Ns;
    real_t Fdn   = -0.5;
    if (is == 0) Fdn = 0.5 * m_mq;
    copy(m_v4, 0, w, is_dn);
    mult_gm5_4d(m_t4, m_v4);
    axpy(m_v4, real_t(1.0), m_t4);
    axpy(m_y4, 0, Fdn, m_v4, 0);

    copy(m_w4, 0, w, is);

    copy(v, is, m_w4, 0);
    axpy(v, is, real_t(1.0), m_y4, 0);

    scal(m_w4, m_b[is]);
    axpy(m_w4, -m_c[is], m_y4);

    m_foprw->mult(m_v4, m_w4);

    axpy(v, is, real_t(4.0) - m_M0, m_v4, 0);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall<AFIELD>::Ddag(AFIELD& v, const AFIELD& w)
{
#pragma omp barrier

  v.set(0.0);
#pragma omp barrier

  for (int is = 0; is < m_Ns; ++is) {
    copy(m_w4, 0, w, is);
    copy(m_y4, m_w4);
    m_foprw->mult_dag(m_v4, m_w4);

    axpy(m_w4, m_b[is] * (real_t(4.0) - m_M0), m_v4);
    axpy(v, is, real_t(1.0), m_w4, 0);

    axpy(m_y4, -m_c[is] * (real_t(4.0) - m_M0), m_v4);

    int    is_up = (is + 1) % m_Ns;
    real_t Fup   = -0.5;
    if (is_up == 0) Fup = 0.5 * m_mq;
    mult_gm5_4d(m_t4, m_y4);
    axpy(m_t4, real_t(-1.0), m_y4);   // m_t4 = - (1 - gm5) * m_y4
    axpy(v, is_up, -Fup, m_t4, 0);

    int    is_dn = (is - 1 + m_Ns) % m_Ns;
    real_t Fdn   = -0.5;
    if (is_dn == m_Ns - 1) Fdn = 0.5 * m_mq;
    mult_gm5_4d(m_t4, m_y4);
    axpy(m_t4, real_t(1.0), m_y4);  // m_t4 = (1 + gm5) * m_y4
    axpy(v, is_dn, Fdn, m_t4, 0);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall<AFIELD>::L_inv(AFIELD& v, const AFIELD& w)
{
#pragma omp barrier

  copy(v, 0, w, 0);
  copy(m_y4, 0, w, 0);
  scal(m_y4, m_e[0]);

  for (int is = 1; is < m_Ns - 1; ++is) {
    copy(v, is, w, is);

    copy(m_v4, 0, v, is - 1);
    mult_gm5_4d(m_t4, m_v4);
    axpy(m_v4, real_t(1.0), m_t4);
    scal(m_v4, real_t(0.5) * m_dm[is] / m_dp[is - 1]);

    axpy(v, is, real_t(1.0), m_v4, 0);
    copy(m_t4, 0, v, is);
    axpy(m_y4, m_e[is], m_t4);
  }

  int is = m_Ns - 1;
  copy(v, is, w, is);
  copy(m_v4, 0, v, is - 1);
  mult_gm5_4d(m_t4, m_v4);
  axpy(m_v4, real_t(1.0), m_t4);
  scal(m_v4, real_t(0.5) * m_dm[is] / m_dp[is - 1]);

  axpy(v, is, real_t(1.0), m_v4, 0);

  mult_gm5_4d(m_t4, m_y4);
  axpy(m_y4, real_t(-1.0), m_t4);
  scal(m_y4, real_t(-0.5));
  axpy(v, is, real_t(1.0), m_y4, 0);

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall<AFIELD>::U_inv(AFIELD& v, const AFIELD& w)
{
  v.set(0.0);

#pragma omp barrier

  int is = m_Ns - 1;
  copy(m_y4, 0, w, is);
  scal(m_y4, real_t(1.0) / (m_dp[is] + m_g));
  copy(v, is, m_y4, 0);

  mult_gm5_4d(m_w4, m_y4);
  axpy(m_y4, real_t(1.0), m_w4);
  scal(m_y4, real_t(0.5));

  for (int is = m_Ns - 2; is >= 0; --is) {
    copy(m_v4, 0, w, is);

    copy(m_w4, 0, v, is + 1);
    mult_gm5_4d(m_t4, m_w4);
    axpy(m_w4, real_t(-1.0), m_t4);

    axpy(m_v4, real_t(0.5) * m_dm[is], m_w4);

    axpy(m_v4, -m_f[is], m_y4);

    scal(m_v4, real_t(1.0) / m_dp[is]);

    copy(v, is, m_v4, 0);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall<AFIELD>::Udag_inv(AFIELD& v, const AFIELD& w)
{
#pragma omp barrier

  copy(m_v4, 0, w, 0);
  scal(m_v4, real_t(1.0) / m_dp[0]);
  copy(v, 0, m_v4, 0);

  copy(m_y4, m_v4);
  scal(m_y4, m_f[0]);

  for (int is = 1; is < m_Ns - 1; ++is) {
    copy(m_t4, 0, w, is);

    copy(m_v4, 0, v, is - 1);
    mult_gm5_4d(m_w4, m_v4);
    axpy(m_v4, real_t(-1.0), m_w4);
    axpy(m_t4, real_t(0.5) * m_dm[is - 1], m_v4);

    scal(m_t4, real_t(1.0) / m_dp[is]);
    copy(v, is, m_t4, 0);

    axpy(m_y4, m_f[is], m_t4);
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
  copy(v, is, m_t4, 0);

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall<AFIELD>::Ldag_inv(AFIELD& v, const AFIELD& w)
{
  int is = m_Ns - 1;

#pragma omp barrier

  copy(v, is, w, is);

  copy(m_y4, 0, w, is);
  mult_gm5_4d(m_w4, m_y4);
  axpy(m_y4, real_t(-1.0), m_w4);
  scal(m_y4, real_t(0.5));

  for (int is = m_Ns - 2; is >= 0; --is) {
    copy(v, is, w, is);

    copy(m_v4, 0, v, is + 1);
    mult_gm5_4d(m_w4, m_v4);
    axpy(m_v4, real_t(1.0), m_w4);
    scal(m_v4, real_t(0.5) * m_dm[is + 1] / m_dp[is]);

    axpy(m_v4, -m_e[is], m_y4);

    axpy(v, is, real_t(1.0), m_v4, 0);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
double AFopr_Domainwall<AFIELD>::flop_count(std::string mode)
{
  int    Lvol  = CommonParameters::Lvol();
  double vsite = static_cast<double>(Lvol);
  double vNs   = static_cast<double>(m_Ns);

  // double flop_Wilson = m_foprw->flop_count("D");
  double flop_Wilson = m_foprw->flop_count();

  double axpy1 = static_cast<double>(2 * m_NinF);
  double scal1 = static_cast<double>(1 * m_NinF);

  double flop_DW = vNs * (flop_Wilson + vsite * (6 * axpy1 + 2 * scal1));
  // In Ddag case, flop_Wilson + 7 axpy which equals flop_DW.

  double flop_LU_inv = 2.0 * vsite *
                       ((3.0 * axpy1 + scal1) * (vNs - 1.0) + axpy1 + 2.0 * scal1);

  double flop = 0.0;
  if (mode == "Prec") {
    flop = flop_LU_inv;
  } else if ((mode == "D") || (mode == "Ddag")) {
    flop = flop_DW;
  } else if (mode == "DdagD") {
    flop = 2.0 * flop_DW;
  } else if ((mode == "D_prec") || (mode == "Ddag_prec")) {
    flop = flop_LU_inv + flop_DW;
  } else if (mode == "DdagD_prec") {
    flop = 2.0 * (flop_LU_inv + flop_DW);
  } else {
    vout.crucial(m_vl, "Error at %s: input repr is undefined.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  return flop;
}


//============================================================END=====
