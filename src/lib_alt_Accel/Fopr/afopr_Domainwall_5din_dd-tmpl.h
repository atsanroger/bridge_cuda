/*!
      @file    afopr_Domainwall_5din_dd-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2590 $
*/

#include "lib_alt_Accel/Fopr/afopr_Domainwall_5din_dd.h"

#include "lib_alt_Accel/Fopr/afopr_common_th-inc.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
using namespace std;

#include "lib/ResourceManager/threadManager.h"
#include "lib/Parameters/commonParameters.h"
#include "lib/Communicator/communicator.h"


template<typename AFIELD>
const std::string AFopr_Domainwall_5din_dd<AFIELD>::class_name
                                         = "AFopr_Domainwall_5din_dd";
//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::init(const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  //int req_comm = 1;  // set 1 if communication forced any time
  int req_comm = 0;  // set 1 if communication forced any time

  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  } else {
    m_vl = CommonParameters::Vlevel();
  }

  vout.general(m_vl, "%s: construction\n", class_name.c_str());


  m_repr = "Dirac";  // now only the Dirac repr is available.

  m_alpha = real_t(1.0);  // Neff alpha; overridden in set_parameters if given.

  std::string repr;
  if (!params.fetch_string("gamma_matrix_type", repr)) {
    if (repr != "Dirac") {
      vout.crucial("  Error at %s: unsupported gamma-matrix type: %s\n",
                   class_name.c_str(), repr.c_str());
      exit(EXIT_FAILURE);
    }
  }

  int Nc = CommonParameters::Nc();
  if (Nc != 3) {
    vout.crucial("%s: only applicable to Nc = 3\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  int Nd = CommonParameters::Nd();
  m_Nvcd = 2 * Nc * Nd;
  m_Ndf  = 2 * Nc * Nc;

  m_Nvol = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();
  m_Nx   = CommonParameters::Nx();
  m_Ny   = CommonParameters::Ny();
  m_Nz   = CommonParameters::Nz();
  m_Nt   = CommonParameters::Nt();

  m_Nsize[0] = m_Nx;
  m_Nsize[1] = m_Ny;
  m_Nsize[2] = m_Nz;
  m_Nsize[3] = m_Nt;

  do_comm_any = 0;
  for (int mu = 0; mu < m_Ndim; ++mu) {
    do_comm[mu] = 1;
    if ((req_comm == 0) && (Communicator::npe(mu) == 1)) do_comm[mu] = 0;
    do_comm_any += do_comm[mu];
    vout.general("  do_comm[%d] = %d\n", mu, do_comm[mu]);
  }

  vout.increase_indent();

  m_Ns = 0; // temporary set
  set_parameters(params);

  vout.decrease_indent();

  // gauge configuration.
  m_U.reset(m_Ndf, m_Nvol, m_Ndim);

  // gauge configuration with block condition.
  m_Ublock.reset(NDF, m_Nvol, m_Ndim);  // inside block
  m_Uhop.reset(  NDF, m_Nvol, m_Ndim);  // hopping part: U - Ublock

  m_Nbdsize.resize(m_Ndim);
  int Nvst = (m_Nvcd / 2) * m_Ns * m_Nvol;
  m_Nbdsize[0] = Nvst / m_Nx;
  m_Nbdsize[1] = Nvst / m_Ny;
  m_Nbdsize[2] = Nvst / m_Nz;
  m_Nbdsize[3] = Nvst / m_Nt;

  setup_channels();

  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::tidyup()
{
  ThreadManager::assert_single_thread(class_name);

  // openacc device memory clean up
  for(int mu = 0; mu < m_Ndim; ++mu){

#ifdef USE_MPI
    real_t* buf_dn1 = (real_t*)chsend_dn[mu].ptr();
    real_t* buf_dn2 = (real_t*)chrecv_dn[mu].ptr();
    real_t* buf_up1 = (real_t*)chrecv_up[mu].ptr();
    real_t* buf_up2 = (real_t*)chsend_up[mu].ptr();
    BridgeACC::afield_tidyup(buf_dn1, m_Nbdsize[mu]);
    BridgeACC::afield_tidyup(buf_dn2, m_Nbdsize[mu]);
    BridgeACC::afield_tidyup(buf_up1, m_Nbdsize[mu]);
    BridgeACC::afield_tidyup(buf_up2, m_Nbdsize[mu]);
#else
    real_t* buf_up = (real_t*)chsend_dn[mu].ptr();
    real_t* buf_dn = (real_t*)chsend_up[mu].ptr();
    BridgeACC::afield_tidyup(buf_up, m_Nbdsize[mu]);
    BridgeACC::afield_tidyup(buf_dn, m_Nbdsize[mu]);
#endif

  }

#ifdef USE_ACCEL_CUDA
  if (m_coeff_on_device) {
    BridgeACC::afield_tidyup(m_e.data(),     m_Ns - 1);
    BridgeACC::afield_tidyup(m_f.data(),     m_Ns - 1);
    BridgeACC::afield_tidyup(m_dpinv.data(), m_Ns);
    BridgeACC::afield_tidyup(m_dm.data(),    m_Ns);
    m_coeff_on_device = false;
  }
#endif

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::setup_channels()
{
  ThreadManager::assert_single_thread(class_name);

  chsend_up.resize(m_Ndim);
  chrecv_up.resize(m_Ndim);
  chsend_dn.resize(m_Ndim);
  chrecv_dn.resize(m_Ndim);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    size_t Nvsize = m_Nbdsize[mu] * sizeof(real_t);

    chsend_dn[mu].send_init(Nvsize, mu, -1);
    chsend_up[mu].send_init(Nvsize, mu, 1);
#ifdef USE_MPI
    chrecv_up[mu].recv_init(Nvsize, mu, 1);
    chrecv_dn[mu].recv_init(Nvsize, mu, -1);
#else
    void *buf_up = (void *)chsend_dn[mu].ptr();
    chrecv_up[mu].recv_init(Nvsize, mu, 1, buf_up);
    void *buf_dn = (void *)chsend_up[mu].ptr();
    chrecv_dn[mu].recv_init(Nvsize, mu, -1, buf_dn);
#endif

    if (do_comm[mu] == 1) {
      chset_send.append(chsend_up[mu]);
      chset_send.append(chsend_dn[mu]);
      chset_recv.append(chrecv_up[mu]);
      chset_recv.append(chrecv_dn[mu]);
    }

    // openacc device memory allocation
#ifdef USE_MPI
    real_t* buf_dn1 = (real_t*)chsend_dn[mu].ptr();
    real_t* buf_dn2 = (real_t*)chrecv_dn[mu].ptr();
    real_t* buf_up1 = (real_t*)chrecv_up[mu].ptr();
    real_t* buf_up2 = (real_t*)chsend_up[mu].ptr();
    BridgeACC::afield_init(buf_dn1, m_Nbdsize[mu]);
    BridgeACC::afield_init(buf_dn2, m_Nbdsize[mu]);
    BridgeACC::afield_init(buf_up1, m_Nbdsize[mu]);
    BridgeACC::afield_init(buf_up2, m_Nbdsize[mu]);
#else
    BridgeACC::afield_init((real_t*)buf_up, m_Nbdsize[mu]);
    BridgeACC::afield_init((real_t*)buf_dn, m_Nbdsize[mu]);
#endif

  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::set_parameters(const Parameters& params)
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
  std::vector<int> block_size;

  int err_optional = 0;
  err_optional += params.fetch_string("gamma_matrix_type", gmset_type);

  int err = 0;
  err += params.fetch_double("quark_mass", mq);
  err += params.fetch_double("domain_wall_height", M0);
  err += params.fetch_int("extent_of_5th_dimension", Ns);
  err += params.fetch_int_vector("boundary_condition", bc);
  err += params.fetch_int_vector("block_size", block_size);

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

  // Neff alpha (bulk-vs-boundary 5D scaling); default 1.0 = standard Moebius.
  // Must be set before set_precond_parameters().
  double alpha = 1.0;
  if (params.fetch_double("parameter_alpha", alpha)) {
    vout.general(m_vl, "  parameter alpha is not provided: set to 1.0.\n");
    alpha = 1.0;
  }
  if (ThreadManager::get_thread_id() == 0) m_alpha = real_t(alpha);
#pragma omp barrier

  set_parameters(real_t(mq), real_t(M0), Ns, bc,
                 real_t(b), real_t(c), block_size);

  //  if (real_t(M0) != m_M0) set_kernel_parameters(params);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::get_parameters(
                                            Parameters& params) const
{
  params.set_string("kernel_type", m_kernel_type);
  params.set_double("quark_mass", double(m_mq));
  params.set_double("domain_wall_height", double(m_M0));
  params.set_int("extent_of_5th_dimension", m_Ns);
  params.set_int_vector("boundary_condition", m_boundary);
  params.set_double("coefficient_b", double(m_b[0]));
  params.set_double("coefficient_c", double(m_c[0]));
  params.set_string("gamma_matrix_type", m_repr);
  params.set_int_vector("block_size", m_block_size);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::set_parameters(
                                   const real_t mq,
                                   const real_t M0,
                                   const int Ns,
                                   const std::vector<int> bc,
                                   const real_t b,
                                   const real_t c,
                                   const std::vector<int> block_size)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  if (ith == 0) {
    m_M0   = real_t(M0);
    m_mq   = real_t(mq);
    m_Ns   = Ns;
    m_NinF = m_Nvcd * m_Ns;

    assert(bc.size() == m_Ndim);
    if (m_boundary.size()   != m_Ndim) m_boundary.resize(m_Ndim);
    if (m_block_size.size() != m_Ndim) m_block_size.resize(m_Ndim);

    for (int mu = 0; mu < m_Ndim; ++mu) {
      m_boundary[mu] = bc[mu];
      m_bc[mu]  = 1;
      if(do_comm[mu] > 0){ // do communication
        if(Communicator::ipe(mu) == 0) m_bc[mu] = m_boundary[mu];
        m_bc2[mu] = 0;
      }else{  // no communication
        m_bc[mu]  = 0;    // for boundar part (dummy)
	m_bc2[mu] = m_boundary[mu];  // for bulk part
      }
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

  for (int mu = 0; mu < m_Ndim; ++mu) {
    m_block_size[mu] = block_size[mu];
  }

#pragma omp barrier

  vout.general(m_vl, "%s: input parameters\n", class_name.c_str());
  vout.general(m_vl, "  mq   = %8.4f\n", m_mq);
  vout.general(m_vl, "  M0   = %8.4f\n", m_M0);
  vout.general(m_vl, "  Ns   = %4d\n", m_Ns);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }
  vout.general(m_vl, "  coefficients:\n");
  for (int is = 0; is < m_Ns; ++is) {
    vout.general(m_vl, "    b[%2d] = %16.10f  c[%2d] = %16.10f\n",
                 is, m_b[is], is, m_c[is]);
  }
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  block_size[%d] = %2d\n", mu, m_block_size[mu]);
  }

  set_precond_parameters();

  int NBx = m_Nx/m_block_size[0];
  int NBy = m_Ny/m_block_size[1];
  int NBz = m_Nz/m_block_size[2];
  int NBt = m_Nt/m_block_size[3];
  int ipex = Communicator::ipe(0);
  int ipey = Communicator::ipe(1);
  int ipez = Communicator::ipe(2);
  int ipet = Communicator::ipe(3);
  m_Ieo = (NBx * ipex + NBy * ipey + NBz * ipez + NBt * ipet) % 2;

  // working 5d vectors.
  if (m_w1.nin() != m_NinF) {
    if (ith == 0) {
      m_w1.reset(m_NinF, m_Nvol, 1);
      m_w2.reset(m_NinF, m_Nvol, 1);
      m_w3.reset(m_NinF, m_Nvol, 1);
      m_w4.reset(m_NinF, m_Nvol, 1);
      m_v1.reset(m_NinF, m_Nvol, 1);
      m_v2.reset(m_NinF, m_Nvol, 1);
    }
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::set_coefficients(
                                    const std::vector<real_t> vec_b,
                                    const std::vector<real_t> vec_c)
{
  if ((vec_b.size() != m_Ns) || (vec_c.size() != m_Ns)) {
    vout.crucial(m_vl, "%s: size of coefficient vectors incorrect.\n",
                 class_name.c_str());
  }

  vout.general(m_vl, "%s: coefficient vectors are set:\n",
               class_name.c_str());

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) {
    for (int is = 0; is < m_Ns; ++is) {
      m_b[is] = real_t(vec_b[is]);
      m_c[is] = real_t(vec_c[is]);
      vout.general(m_vl, "b[%2d] = %16.10f  c[%2d] = %16.10f\n",
                   is, m_b[is], is, m_c[is]);
    }
  }

  set_precond_parameters();
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::set_precond_parameters()
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) {

    if (m_dp.size() != m_Ns) {
      m_dp.resize(m_Ns);
      m_dm.resize(m_Ns);
      m_dpinv.resize(m_Ns);
      m_e.resize(m_Ns - 1);
      m_f.resize(m_Ns - 1);
    }

    // Neff alpha congruence on the 5D bulk diagonal (same as the non-eo 5din):
    // m_dp/m_dm scale by alpha; alpha cancels in the m_e/m_f LU factors because
    // m_f[0] is divided by alpha. alpha = 1 reproduces the standard coefficients.
    for (int is = 0; is < m_Ns; ++is) {
      m_dp[is] = m_alpha * (1.0 + m_b[is] * (4.0 - m_M0));
      m_dm[is] = m_alpha * (1.0 - m_c[is] * (4.0 - m_M0));
    }

    m_e[0] = m_mq * m_dm[m_Ns - 1] / m_dp[0];
    m_f[0] = m_mq * m_dm[0] / m_alpha;
    for (int is = 1; is < m_Ns - 1; ++is) {
      m_e[is] = m_e[is - 1] * m_dm[is - 1] / m_dp[is];
      m_f[is] = m_f[is - 1] * m_dm[is] / m_dp[is - 1];
    }

    m_g = m_e[m_Ns - 2] * m_dm[m_Ns - 2];

    for (int is = 0; is < m_Ns - 1; ++is) {
      m_dpinv[is] = 1.0 / m_dp[is];
    }
    m_dpinv[m_Ns - 1] = 1.0 / (m_dp[m_Ns - 1] + m_g);

#ifdef USE_ACCEL_CUDA
    // The dd SAP kernels (mult_domainwall_5din_dd_5dir/5dirdag_dirac) read the
    // Moebius b/c from __constant__ memory. Upload them here (b/c only; mass-
    // independent, so it will not clobber any eo operator's mass-dependent
    // const_e/f/dpinv living in the same process).
    BridgeACC::setDomainwallConstantBC(m_b.data(), m_c.data(), m_Ns);

    // Mass-dependent LU coefficients (e/f/dpinv/dm) are per-operator, so they
    // are kept as this operator's own device copy (mapped once, re-uploaded on
    // rebuild) and read by the L/U-inverse kernels via dev_ptr -- NOT shared
    // __constant__ memory (which a coexisting D/D_PV would clobber).
    if (!m_coeff_on_device) {
      BridgeACC::afield_init(m_e.data(),     m_Ns - 1);
      BridgeACC::afield_init(m_f.data(),     m_Ns - 1);
      BridgeACC::afield_init(m_dpinv.data(), m_Ns);
      BridgeACC::afield_init(m_dm.data(),    m_Ns);
      m_coeff_on_device = true;
    }
    BridgeACC::copy_to_device(m_e.data(),     m_Ns - 1);
    BridgeACC::copy_to_device(m_f.data(),     m_Ns - 1);
    BridgeACC::copy_to_device(m_dpinv.data(), m_Ns);
    BridgeACC::copy_to_device(m_dm.data(),    m_Ns);
#endif
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::set_config(Field *u)
{
  int nth = ThreadManager::get_num_threads();

  vout.detailed(m_vl, "%s: set_config is called: num_threads = %d\n",
                class_name.c_str(), nth);

  if (nth > 1) {
    set_config_impl(u);
  } else {
    set_config_omp(u);
  }

  vout.detailed(m_vl, "%s: set_config finished\n", class_name.c_str());
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::set_config_omp(Field *u)
{
  vout.detailed(m_vl, "  set_config_omp is called.\n");

#pragma omp parallel
  {
    set_config_impl(u);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::set_config_impl(Field *u)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  //  if (ith == 0) m_conf = u;

  AIndex_lex<real_t,AFIELD::IMPL> index_lex;

  convert_gauge(index_lex, m_U, *u);

  copy(m_Ublock, m_U);

  if (ith == 0){
    BridgeACC::set_block_config(m_Ublock.ptr(0), m_Nsize, &m_block_size[0]);
  }

  copy(m_Uhop, m_U);
  m_Uhop.axpy(real_t(-1.0), m_Ublock);

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::convert(AFIELD& v, const Field& w)
{
#pragma omp barrier

  int ith, nth, isite, nsite;
  set_threadtask_afopr(ith, nth, isite, nsite, m_Nvol);

  AIndex_lex<real_t, AFIELD::IMPL> index;

  for (int site = isite; site < nsite; ++site) {
    for (int is = 0; is < m_Ns; ++is) {
      for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
         int in_alt = ivcd + NVCD * is;
         real_t vt = real_t(w.cmp(ivcd, site, is));
         v.set_host(index.idx(in_alt, m_NinF, site, 0), vt);
      }
    }
  } // site-llop

  v.update_device();

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::reverse(Field& v, const AFIELD& w)
{
#pragma omp barrier

  int ith, nth, isite, nsite;
  set_threadtask_afopr(ith, nth, isite, nsite, m_Nvol);

  w.update_host();

  AIndex_lex<real_t, AFIELD::IMPL> index;

  for (int site = isite; site < nsite; ++site) {
    for (int is = 0; is < m_Ns; ++is) {
      for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
        int in_alt = ivcd + NVCD * is;
        double vt = double(w.cmp_host(index.idx(in_alt, m_NinF, site, 0)));
        v.set(ivcd, site, is, vt);
      }
    }
  } // site-llop

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::set_mode(std::string mode)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_mode = mode;

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::mult(AFIELD& v, const AFIELD& w)
{
  if (m_mode == "D") {
    //D(v, w);
    D_alt(v, w);
  } else if (m_mode == "Ddag") {
    //Ddag(v, w);
    Ddag_alt(v, w);
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
void AFopr_Domainwall_5din_dd<AFIELD>::mult_dag(AFIELD& v, const AFIELD& w)
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
void AFopr_Domainwall_5din_dd<AFIELD>::mult(AFIELD& v, const AFIELD& w,
                                    std::string mode)
{
  assert(w.check_size(m_NinF, m_Nvol, 1));
  assert(v.check_size(m_NinF, m_Nvol, 1));

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
void AFopr_Domainwall_5din_dd<AFIELD>::mult_dag(AFIELD& v, const AFIELD& w,
                                        std::string mode)
{
  assert(w.check_size(m_NinF, m_Nvol, 1));
  assert(v.check_size(m_NinF, m_Nvol, 1));

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
    std::cout << "mode undeifined in AFopr_Domainwall_5din_dd.\n";
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::mult_sap(AFIELD& v,
                                                const AFIELD& w,
                                                const int ieo)
{
  if (m_mode == "D") {
    if(ieo < 0){
      mult_D_dd(v, w);
    }else{
      mult_D_sap(v, w, ieo);
    }
  } else if (m_mode == "Ddag") {
    if(ieo < 0){
      mult_Ddag_dd(v, w);
    }else{
      mult_Ddag_sap(v, w, ieo);
    }
  } else if (m_mode == "DdagD") {
    if(ieo < 0){
      //mult_D_dd(m_v1, w);
      //mult_Ddag_dd(v, m_v1);
      mult_D_sap(m_v1, w, ieo);
      mult_Ddag_sap(v, m_v1, ieo);
    }else{
      mult_D_sap(m_v1, w, ieo);
      mult_Ddag_sap(v, m_v1, ieo);
    }
  } else {
    std::cout << "mode undeifined in AFopr_Domainwall_5din_dd.\n";
    exit(EXIT_FAILURE);
  }
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::mult_D_sap(AFIELD& v,
                                                  const AFIELD& w,
                                                  const int ieo)
{
#pragma omp barrier

  v.set(0.0);
  m_w1.set(0.0);

  int ith = ThreadManager::get_thread_id();

  if (ith == 0){
    real_t *vp = v.ptr(0);
    real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);
    real_t *yp = m_w1.ptr(0);  // working vector
    real_t *up = m_Ublock.ptr(0);

    int jeo = (ieo + m_Ieo) % 2;
    if(ieo < 0) jeo = ieo;

    //BridgeACC::mult_domainwall_5din_5dir_dirac(
    //                             vp, yp, wp,
    //                             m_mq, m_M0, m_Ns, &m_b[0], &m_c[0], 
    //                             m_Nsize);

    BridgeACC::mult_domainwall_5din_dd_5dir_dirac(
                                vp, yp, wp,
                                m_mq, m_M0, m_Ns, &m_b[0], &m_c[0],
                                m_Nsize, &m_block_size[0], jeo, m_alpha);

    BridgeACC::mult_domainwall_5din_dd_hopb_dirac(
                                 vp, up, yp,
                                 m_Ns, m_bc2, m_Nsize,
                                 &m_block_size[0], jeo, 1);
  }

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::mult_Ddag_sap(AFIELD& v,
                                                     const AFIELD& w,
                                                     const int ieo)
{
#pragma omp barrier

  v.set(0.0);
  m_w1.set(0.0);

  int ith = ThreadManager::get_thread_id();
  if (ith == 0){
    real_t *vp = v.ptr(0);
    real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);
    real_t *yp = m_w1.ptr(0);
    real_t *up = m_Ublock.ptr(0);

    int jeo = (ieo + m_Ieo) % 2;
    if(ieo < 0) jeo = ieo;

    BridgeACC::mult_domainwall_5din_mult_gm5_dirac(
                                         vp, wp, m_Ns, m_Nsize);

    //BridgeACC::mult_domainwall_5din_dd_clear(yp, m_Ns, m_Nsize);

    BridgeACC::mult_domainwall_5din_dd_hopb_dirac(
                                 yp, up, vp,
                                 m_Ns, m_bc2, m_Nsize,
                                 &m_block_size[0], jeo, 0);

    //BridgeACC::mult_domainwall_5din_5dirdag_dirac(
    //                             vp, yp, wp, m_mq, m_M0, m_Ns,
    //                             &m_b[0], &m_c[0], m_Nsize);
 
    //v.set(0.0);

    BridgeACC::mult_domainwall_5din_dd_5dirdag_dirac(
                                vp, yp, wp, m_mq, m_M0, m_Ns,
                                &m_b[0], &m_c[0],
                                m_Nsize, &m_block_size[0], jeo, m_alpha);
  }

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::mult_dd(AFIELD& v,
                                               const AFIELD& w)
{
  if(m_mode == "D") {
    mult_D_dd(v, w);
  }else if(m_mode == "Ddag") {
    mult_Ddag_dd(v, w);
  }else if(m_mode == "DdagD") {
    mult_D_dd(m_v1, w);
    mult_Ddag_dd(v, m_v1);
  }
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::mult_D_dd(AFIELD& v,
                                                 const AFIELD& w)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  if (ith == 0){

    real_t *vp = v.ptr(0);
    real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);
    real_t *yp = m_w1.ptr(0);  // working vector
    real_t *up = m_Ublock.ptr(0);

    BridgeACC::mult_domainwall_5din_5dir_dirac(
                                 vp, yp, wp,
                                 m_mq, m_M0, m_Ns, &m_b[0], &m_c[0], 
                                 m_Nsize);

    BridgeACC::mult_domainwall_5din_hopb_dirac(
                               vp, up, yp,
                               m_Ns, m_bc2, m_Nsize, do_comm, 1);
  }

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::mult_Ddag_dd(AFIELD& v,
                                                    const AFIELD& w)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  if (ith == 0){

    real_t *vp = v.ptr(0);
    real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);
    real_t *yp = m_w1.ptr(0);
    real_t *up = m_Ublock.ptr(0);

    BridgeACC::mult_domainwall_5din_mult_gm5_dirac(
                                           vp, wp, m_Ns, m_Nsize);

    BridgeACC::mult_domainwall_5din_hopb_dirac(
                   yp, up, vp, m_Ns, m_bc2, m_Nsize, do_comm, 0);

    BridgeACC::mult_domainwall_5din_5dirdag_dirac(
                                   vp, yp, wp, m_mq, m_M0, m_Ns,
                                   &m_b[0], &m_c[0], m_Nsize);
  }

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::mult_gm5(AFIELD& v, const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  if (ith == 0)
    BridgeACC::mult_domainwall_5din_mult_gm5_dirac(
                                              vp, wp, m_Ns, m_Nsize);

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::DdagD(AFIELD& v, const AFIELD& w)
{
  assert(w.check_size(m_NinF, m_Nvol, 1));
  assert(v.check_size(m_NinF, m_Nvol, 1));

  //D(m_v1, w);
  //Ddag(v, m_v1);
  D_alt(m_v1, w);
  Ddag_alt(v, m_v1);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::DDdag(AFIELD& v, const AFIELD& w)
{
  assert(w.check_size(m_NinF, m_Nvol, 1));
  assert(v.check_size(m_NinF, m_Nvol, 1));

  Ddag(m_v1, w);
  D(v, m_v1);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::DdagD_prec(AFIELD& v, const AFIELD& w)
{
  assert(w.check_size(m_NinF, m_Nvol, 1));
  assert(v.check_size(m_NinF, m_Nvol, 1));

  D_prec(m_v1, w);
  Ddag_prec(v, m_v1);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::D_prec(AFIELD& v, const AFIELD& w)
{
  //vout.crucial("D_prec is called.\n");
  L_inv(v, w);
  U_inv(m_v2, v);
  D(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::Ddag_prec(AFIELD& v, const AFIELD& w)
{
  //vout.crucial("Ddag_prec is called.\n");
  Ddag(v, w);
  Udag_inv(m_v2, v);
  Ldag_inv(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::Prec(AFIELD& v, const AFIELD& w)
{
  L_inv(m_v2, w);
  U_inv(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::Precdag(AFIELD& v, const AFIELD& w)
{
  Udag_inv(m_v2, w);
  Ldag_inv(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::H(AFIELD& v, const AFIELD& w)
{
  D(m_v2, w);
  mult_gm5R(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::Hdag(AFIELD& v, const AFIELD& w)
{
  mult_gm5R(m_v2, w);
  Ddag(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::mult_gm5R(AFIELD& v, const AFIELD& w)
{
#pragma omp barrier

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

  int ith = ThreadManager::get_thread_id();

  if (ith == 0)
    BridgeACC::mult_domainwall_5din_mult_gm5R_dirac(
                                              vp, wp, m_Ns, m_Nsize);

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::mult_R(AFIELD& v, const AFIELD& w)
{
  assert(w.check_size(m_NinF, m_Nvol, 1));
  assert(v.check_size(m_NinF, m_Nvol, 1));

#pragma omp barrier

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

  int ith = ThreadManager::get_thread_id();

  if (ith == 0)
    BridgeACC::mult_domainwall_5din_mult_R(vp, wp, m_Ns, m_Nsize);

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::D(AFIELD& v, const AFIELD& w)
{
#pragma omp barrier

  //vout.crucial("D is called.\n");

  int ith = ThreadManager::get_thread_id();

  if (ith == 0){

    real_t *vp = v.ptr(0);
    real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);
    real_t *yp = m_w1.ptr(0);  // working vector
    real_t *up = m_U.ptr(0);

    BridgeACC::mult_domainwall_5din_5dir_dirac(
                                   vp, yp, wp, m_mq, m_M0, m_Ns,
                                   &m_b[0], &m_c[0], m_Nsize);

    if (do_comm_any > 0) {

      real_t *buf1_xp = (real_t *)chsend_dn[0].ptr();
      real_t *buf1_xm = (real_t *)chsend_up[0].ptr();
      real_t *buf1_yp = (real_t *)chsend_dn[1].ptr();
      real_t *buf1_ym = (real_t *)chsend_up[1].ptr();
      real_t *buf1_zp = (real_t *)chsend_dn[2].ptr();
      real_t *buf1_zm = (real_t *)chsend_up[2].ptr();
      real_t *buf1_tp = (real_t *)chsend_dn[3].ptr();
      real_t *buf1_tm = (real_t *)chsend_up[3].ptr();

      BridgeACC::mult_domainwall_5din_hop1_dirac(
                          buf1_xp, buf1_xm, buf1_yp, buf1_ym,
                          buf1_zp, buf1_zm, buf1_tp, buf1_tm,
                          up, yp,
                          m_Ns, m_bc, m_Nsize, do_comm);

      chset_recv.start();
      chset_send.start();

    }

    BridgeACC::mult_domainwall_5din_hopb_dirac(
                             vp, up, yp,
                             m_Ns, m_bc2, m_Nsize, do_comm, 1);

    if (do_comm_any > 0) {

      chset_recv.wait();
      chset_send.wait();

      real_t *buf2_xp = (real_t *)chrecv_up[0].ptr();
      real_t *buf2_xm = (real_t *)chrecv_dn[0].ptr();
      real_t *buf2_yp = (real_t *)chrecv_up[1].ptr();
      real_t *buf2_ym = (real_t *)chrecv_dn[1].ptr();
      real_t *buf2_zp = (real_t *)chrecv_up[2].ptr();
      real_t *buf2_zm = (real_t *)chrecv_dn[2].ptr();
      real_t *buf2_tp = (real_t *)chrecv_up[3].ptr();
      real_t *buf2_tm = (real_t *)chrecv_dn[3].ptr();

      BridgeACC::mult_domainwall_5din_hop2_dirac(
                                 vp, up, yp,
                                 buf2_xp, buf2_xm, buf2_yp, buf2_ym,
                                 buf2_zp, buf2_zm, buf2_tp, buf2_tm,
                                 m_Ns, m_bc, m_Nsize, do_comm);
    }

  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::Ddag(AFIELD& v, const AFIELD& w)
{
#pragma omp barrier
  // vout.crucial("Ddag is called.\n");

  int ith = ThreadManager::get_thread_id();
  if (ith == 0){
  
    real_t *vp = v.ptr(0);
    real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);
    real_t *yp = m_w1.ptr(0);
    real_t *up = m_U.ptr(0);

    BridgeACC::mult_domainwall_5din_mult_gm5_dirac(
                                           vp, wp, m_Ns, m_Nsize);

    if (do_comm_any > 0) {

      real_t *buf1_xp = (real_t *)chsend_dn[0].ptr();
      real_t *buf1_xm = (real_t *)chsend_up[0].ptr();
      real_t *buf1_yp = (real_t *)chsend_dn[1].ptr();
      real_t *buf1_ym = (real_t *)chsend_up[1].ptr();
      real_t *buf1_zp = (real_t *)chsend_dn[2].ptr();
      real_t *buf1_zm = (real_t *)chsend_up[2].ptr();
      real_t *buf1_tp = (real_t *)chsend_dn[3].ptr();
      real_t *buf1_tm = (real_t *)chsend_up[3].ptr();

      BridgeACC::mult_domainwall_5din_hop1_dirac(
                     buf1_xp, buf1_xm, buf1_yp, buf1_ym,
                     buf1_zp, buf1_zm, buf1_tp, buf1_tm,
                     up, vp, m_Ns, m_bc, m_Nsize, do_comm);

      chset_recv.start();
      chset_send.start();
    }

    //BridgeACC::mult_domainwall_5din_dd_clear(yp, m_Ns, m_Nsize);

    BridgeACC::mult_domainwall_5din_hopb_dirac(
                   yp, up, vp, m_Ns, m_bc2, m_Nsize, do_comm, 0);

    if (do_comm_any > 0) {

      chset_recv.wait();
      chset_send.wait();

      real_t *buf2_xp = (real_t *)chrecv_up[0].ptr();
      real_t *buf2_xm = (real_t *)chrecv_dn[0].ptr();
      real_t *buf2_yp = (real_t *)chrecv_up[1].ptr();
      real_t *buf2_ym = (real_t *)chrecv_dn[1].ptr();
      real_t *buf2_zp = (real_t *)chrecv_up[2].ptr();
      real_t *buf2_zm = (real_t *)chrecv_dn[2].ptr();
      real_t *buf2_tp = (real_t *)chrecv_up[3].ptr();
      real_t *buf2_tm = (real_t *)chrecv_dn[3].ptr();

      BridgeACC::mult_domainwall_5din_hop2_dirac(
                               yp, up, vp,
                               buf2_xp, buf2_xm, buf2_yp, buf2_ym,
                               buf2_zp, buf2_zm, buf2_tp, buf2_tm,
                               m_Ns, m_bc, m_Nsize, do_comm);
    }

    BridgeACC::mult_domainwall_5din_5dirdag_dirac(
                                   vp, yp, wp, m_mq, m_M0, m_Ns,
                                   &m_b[0], &m_c[0], m_Nsize);

  }
#pragma omp barrier

}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::L_inv(AFIELD& v, const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

  int ith = ThreadManager::get_thread_id();

  if (ith == 0)
    BridgeACC::mult_domainwall_5din_L_inv_dirac(
              vp, wp, m_Ns, m_Nsize, &m_e[0], &m_dpinv[0], &m_dm[0]);

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::U_inv(AFIELD& v, const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

  int ith = ThreadManager::get_thread_id();

  if (ith == 0)
    BridgeACC::mult_domainwall_5din_U_inv_dirac(
              vp, wp, m_Ns, m_Nsize, &m_f[0], &m_dpinv[0], &m_dm[0]);

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::Ldag_inv(AFIELD& v,
                                             const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

  int ith = ThreadManager::get_thread_id();

  if (ith == 0)
    BridgeACC::mult_domainwall_5din_Ldag_inv_dirac(
              vp, wp, m_Ns, m_Nsize, &m_e[0], &m_dpinv[0], &m_dm[0]);
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::Udag_inv(AFIELD& v,
                                             const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

  int ith = ThreadManager::get_thread_id();

  if (ith == 0)
    BridgeACC::mult_domainwall_5din_Udag_inv_dirac(
              vp, wp, m_Ns, m_Nsize, &m_f[0], &m_dpinv[0], &m_dm[0]);
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::D_alt(AFIELD& v, const AFIELD& w)
{
#pragma omp barrier

  //vout.crucial("D_alt is called.\n");

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);
  real_t *yp = m_w1.ptr(0);  // working vector

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){
    BridgeACC::mult_domainwall_5din_5dir_dirac(
                                   vp, yp, wp, m_mq, m_M0, m_Ns,
                                   &m_b[0], &m_c[0], m_Nsize);
  }

  int isap = 0; // no SAP
  int idag = 1; // D
  for(int mu = 0; mu < m_Ndim; ++mu){
    mult_up(v, w, isap, idag, mu);
    mult_dn(v, w, isap, idag, mu);
  }

  //for(int mu = 0; mu < m_Ndim; ++mu){
  //  mult_up4(vp, yp, isap, idag, mu);
  //  mult_dn4(vp, yp, isap, idag, mu);
  //}

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::Ddag_alt(AFIELD& v, const AFIELD& w)
{
#pragma omp barrier

  //vout.crucial("Ddag_alt is called.\n");

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);
  real_t *yp = m_w1.ptr(0);  // working vector

  int ith, nth;
  set_thread(ith, nth);

  m_w1.set(0.0);

  if(ith == 0){

    BridgeACC::mult_domainwall_5din_5dirdag_dirac(
                                   vp, yp, wp, m_mq, m_M0, m_Ns,
                                   &m_b[0], &m_c[0], m_Nsize);
  }

  int isap =  0; // no SAP
  int idag = -1; // Ddag
  for(int mu = 0; mu < m_Ndim; ++mu){
    mult_up(v, w, isap, idag, mu);
    mult_dn(v, w, isap, idag, mu);
  }

  /*
  if(ith == 0){
    BridgeACC::mult_domainwall_5din_mult_gm5_dirac(
                                           vp, wp, m_Ns, m_Nsize);
  }

  m_w1.set(0.0);

  int isap =  0; // no SAP
  int idag = -1; // Ddag
  for(int mu = 0; mu < m_Ndim; ++mu){
    mult_up4(yp, vp, isap, idag, mu);
    mult_dn4(yp, vp, isap, idag, mu);
  }

  if(ith == 0){
    BridgeACC::mult_domainwall_5din_5dirdag_dirac(
                                   vp, yp, wp, m_mq, m_M0, m_Ns,
                                   &m_b[0], &m_c[0], m_Nsize);
  }
  */

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::mult_dup(AFIELD &v,
                                                const AFIELD &w,
                                                const int mu)
{
  if(m_mode == "D") {

    int isap =  1; // SAP
    int idag =  1; // D
    mult_up(v, w, isap, idag, mu);

  } else if (m_mode == "Ddag") {

    int isap =  1; // SAP
    int idag = -1; // D
    mult_up(v, w, isap, idag, mu);

  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::mult_ddn(AFIELD &v,
                                                const AFIELD &w,
                                                const int mu)
{
  if(m_mode == "D") {

    int isap =  1; // SAP
    int idag =  1; // D
    mult_dn(v, w, isap, idag, mu);

  } else if (m_mode == "Ddag") {

    int isap =  1; // SAP
    int idag = -1; // D
    mult_dn(v, w, isap, idag, mu);

  }
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::mult_up(
                                         AFIELD& v, const AFIELD& w,
                                         int isap, int idag, int mu)
{  // isap = 0: no SAP, isap != 0: SAP
   // idag = 1: D, idag = -1: D^dag

#pragma omp barrier

  int ith, nth;
  set_thread(ith, nth);

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);
  real_t *yp = m_w2.ptr(0);  // working vector
  real_t *xp = m_w3.ptr(0);  // working vector
  real_t *zp = m_w4.ptr(0);  // working vector

  if(idag == 1){

    if(ith == 0){
      BridgeACC::mult_domainwall_5din_5dir_dirac(
                                   xp, yp, wp, m_mq, m_M0, m_Ns,
                                   &m_b[0], &m_c[0], m_Nsize);
    }

    mult_up4(vp, yp, isap, idag, mu);
    
  }else{

    if(ith == 0){
      BridgeACC::mult_domainwall_5din_mult_gm5_dirac(
                                           xp, wp, m_Ns, m_Nsize);
    }

    m_w2.set(0.0);
    mult_up4(yp, xp, isap, idag, mu);

    m_w4.set(0.0);
    if(ith == 0){
      BridgeACC::mult_domainwall_5din_5dirdag_dirac(
                                   xp, yp, zp, m_mq, m_M0, m_Ns,
                                   &m_b[0], &m_c[0], m_Nsize);
    }

    axpy(v, real_t(1.0), m_w3);

  }

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::mult_dn(
                                         AFIELD& v, const AFIELD& w,
                                         int isap, int idag, int mu)
{  // isap = 0: no SAP, isap != 0: SAP
   // idag = 1: D, idag = -1: D^dag

#pragma omp barrier

  int ith, nth;
  set_thread(ith, nth);

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);
  real_t *yp = m_w2.ptr(0);  // working vector
  real_t *xp = m_w3.ptr(0);  // working vector
  real_t *zp = m_w4.ptr(0);  // working vector

  if(idag == 1){

    if(ith == 0){
      BridgeACC::mult_domainwall_5din_5dir_dirac(
                                   xp, yp, wp, m_mq, m_M0, m_Ns,
                                   &m_b[0], &m_c[0], m_Nsize);
    }

    mult_dn4(vp, yp, isap, idag, mu);
    
  }else{

    if(ith == 0){
      BridgeACC::mult_domainwall_5din_mult_gm5_dirac(
                                           xp, wp, m_Ns, m_Nsize);
    }

    m_w2.set(0.0);
    mult_dn4(yp, xp, isap, idag, mu);

    m_w4.set(0.0);
    if(ith == 0){
      BridgeACC::mult_domainwall_5din_5dirdag_dirac(
                                   xp, yp, zp, m_mq, m_M0, m_Ns,
                                   &m_b[0], &m_c[0], m_Nsize);
    }

    axpy(v, real_t(1.0), m_w3);

  }

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::mult_up4(
                                         real_t *v2, real_t *v1,
                                         int isap, int idag, int mu)
{  // isap = 0: no SAP, isap != 0: SAP
   // idag = 1: D, idag = -1: D^dag  --- dummy

#pragma omp barrier

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){

    AIndex_lex<real_t,AFIELD::IMPL> index_lex;

    real_t* buf1 = (real_t*)chsend_dn[mu].ptr();
    real_t* buf2 = (real_t*)chrecv_up[mu].ptr();
    real_t* up = m_U.ptr(0);
    if(isap != 0) up = m_Uhop.ptr(0);

    if(do_comm[mu] > 0){
      if(mu == 0){
        BridgeACC::mult_domainwall_5din_xp1(buf1, up, v1,
                                                  m_Ns, m_bc, m_Nsize);
      }else if(mu == 1){
        BridgeACC::mult_domainwall_5din_yp1(buf1, up, v1,
                                                  m_Ns, m_bc, m_Nsize);
      }else if(mu == 2){
        BridgeACC::mult_domainwall_5din_zp1(buf1, up, v1,
                                                  m_Ns, m_bc, m_Nsize);
      }else if(mu == 3){
        BridgeACC::mult_domainwall_5din_tp1_dirac(buf1, up, v1,
                                                  m_Ns, m_bc, m_Nsize);
      }

      BridgeACC::copy_from_device(buf1, m_Nbdsize[mu]);

      chrecv_up[mu].start();
      chsend_dn[mu].start();
    }

    if(mu == 0){
      BridgeACC::mult_domainwall_5din_xpb(v2, up, v1,
                                m_Ns, m_bc2, m_Nsize, do_comm, 1);
    }else if(mu == 1){
      BridgeACC::mult_domainwall_5din_ypb(v2, up, v1,
                                m_Ns, m_bc2, m_Nsize, do_comm, 1);
    }else if(mu == 2){
      BridgeACC::mult_domainwall_5din_zpb(v2, up, v1,
                                m_Ns, m_bc2, m_Nsize, do_comm, 1);
    }else if(mu == 3){
      BridgeACC::mult_domainwall_5din_tpb_dirac(v2, up, v1,
                                m_Ns, m_bc2, m_Nsize, do_comm, 1);
    }

    if(do_comm[mu] > 0){
      chsend_dn[mu].wait();
      chrecv_up[mu].wait();

      BridgeACC::copy_to_device(buf2, m_Nbdsize[mu]);
      if(mu == 0){
        BridgeACC::mult_domainwall_5din_xp2(v2, up, buf2,
                                            m_Ns, m_bc, m_Nsize);
      }else if(mu == 1){
        BridgeACC::mult_domainwall_5din_yp2(v2, up, buf2,
                                            m_Ns, m_bc, m_Nsize);
      }else if(mu == 2){
        BridgeACC::mult_domainwall_5din_zp2(v2, up, buf2,
                                            m_Ns, m_bc, m_Nsize);
      }else if(mu == 3){
        BridgeACC::mult_domainwall_5din_tp2_dirac(v2, up, buf2,
                                            m_Ns, m_bc, m_Nsize);
      }
    }

  }
#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::mult_dn4(
                                         real_t *v2, real_t *v1,
                                         int isap, int idag, int mu)
{  // isap = 0: no SAP, isap != 0: SAP
   // idag = 1: D, idag = -1: D^dag

#pragma omp barrier

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){

    AIndex_lex<real_t,AFIELD::IMPL> index_lex;

    real_t *buf1 = (real_t*)chsend_up[mu].ptr();
    real_t *buf2 = (real_t*)chrecv_dn[mu].ptr();
    real_t *up = m_U.ptr(0);
    if(isap != 0) up = m_Uhop.ptr(0);

    if(do_comm[mu] > 0){
      if(mu == 0){
        BridgeACC::mult_domainwall_5din_xm1(buf1, up, v1,
                                            m_Ns, m_bc, m_Nsize);
      }else if(mu == 1){
        BridgeACC::mult_domainwall_5din_ym1(buf1, up, v1,
                                            m_Ns, m_bc, m_Nsize);
      }else if(mu == 2){
        BridgeACC::mult_domainwall_5din_zm1(buf1, up, v1,
                                            m_Ns, m_bc, m_Nsize);
      }else if(mu == 3){
        BridgeACC::mult_domainwall_5din_tm1_dirac(buf1, up, v1,
                                            m_Ns, m_bc, m_Nsize);
      }

      BridgeACC::copy_from_device(buf1, m_Nbdsize[mu]);

      chrecv_dn[mu].start();
      chsend_up[mu].start();
    }

    if(mu == 0){
      BridgeACC::mult_domainwall_5din_xmb(v2, up, v1,
                                  m_Ns, m_bc2, m_Nsize, do_comm, 1);
    }else if(mu == 1){
      BridgeACC::mult_domainwall_5din_ymb(v2, up, v1,
                                  m_Ns, m_bc2, m_Nsize, do_comm, 1);
    }else if(mu == 2){
      BridgeACC::mult_domainwall_5din_zmb(v2, up, v1,
                                  m_Ns, m_bc2, m_Nsize, do_comm, 1);
    }else if(mu == 3){
      BridgeACC::mult_domainwall_5din_tmb_dirac(v2, up, v1,
                                  m_Ns, m_bc2, m_Nsize, do_comm, 1);
    }

    if(do_comm[mu] > 0){
      chsend_up[mu].wait();
      chrecv_dn[mu].wait();

      BridgeACC::copy_to_device(buf2, m_Nbdsize[mu]);
      if(mu == 0){
          BridgeACC::mult_domainwall_5din_xm2(v2, up, buf2,
                                              m_Ns, m_bc, m_Nsize);
      }else if(mu == 1){
          BridgeACC::mult_domainwall_5din_ym2(v2, up, buf2,
                                              m_Ns, m_bc, m_Nsize);
      }else if(mu == 2){
          BridgeACC::mult_domainwall_5din_zm2(v2, up, buf2,
                                              m_Ns, m_bc, m_Nsize);
      }else if(mu == 3){
          BridgeACC::mult_domainwall_5din_tm2_dirac(v2, up, buf2,
                                              m_Ns, m_bc, m_Nsize);
      }
    }

  }
#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
double AFopr_Domainwall_5din_dd<AFIELD>::flop_count(std::string mode)
{
  int    Lvol  = CommonParameters::Lvol();
  double vsite = static_cast<double>(Lvol);
  double vNs   = static_cast<double>(m_Ns);
  int    Nc    = CommonParameters::Nc();
  int    Nd    = CommonParameters::Nd();

  double flop_Wilson;
  double flop_LU_inv;
  double axpy1 = static_cast<double>(2 * m_NinF);
  double scal1 = static_cast<double>(1 * m_NinF);
  if (m_repr == "Dirac") {
    flop_Wilson = static_cast<double>(
      Nc * Nd * (4 + 6 * (4 * Nc + 2) + 2 * (4 * Nc + 1))) * vsite;
    flop_LU_inv = static_cast<double>(Nc * Nd * (2 + (vNs - 1) * 26)) * vsite;
  } else if (m_repr == "Chiral") {
    flop_Wilson = static_cast<double>(
      Nc * Nd * (4 + 8 * (4 * Nc + 2))) * vsite;
    flop_LU_inv = static_cast<double>(Nc * Nd * (2 + (vNs - 1) * 10)) * vsite;
  }

  // Note that m_NinF := m_Nvcd * m_Ns;
  double flop_DW = vNs * flop_Wilson + vsite * (6 * axpy1 + 2 * scal1);
  // In Ddag case, flop_Wilson + 7 axpy which equals flop_DW.

  //  double flop_LU_inv = 2.0 * vsite *
  //            ((3.0*axpy1 + scal1)*(vNs-1.0) + axpy1 + 2.0*scal1);

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
    vout.crucial(m_vl, "Error at %s: input mode is undefined.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  return flop;
}

//====================================================================
template<typename AFIELD>
double AFopr_Domainwall_5din_dd<AFIELD>::flop_count_sap()
{
  // not yet implemented.
  return 0.0;
}

//============================================================END=====
