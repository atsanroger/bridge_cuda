/*!
      @file    afopr_Domainwall_5din-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2604 $
*/

template<typename AFIELD>
const std::string AFopr_Domainwall_5din<AFIELD>::class_name
  = "AFopr_Domainwall_5din";
//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din<AFIELD>::init(const Parameters& params)
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

  m_Nxv  = m_Nx / VLENX;
  m_Nyv  = m_Ny / VLENY;
  m_Nstv = m_Nvol / VLEN;

  if (VLENX * m_Nxv != m_Nx) {
    vout.crucial(m_vl, "%s: Nx must be multiple of VLENX.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
  if (VLENY * m_Nyv != m_Ny) {
    vout.crucial(m_vl, "%s: Ny must be multiple of VLENY.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  vout.general(m_vl, "  VLENX = %2d  Nxv  = %d\n", VLENX, m_Nxv);
  vout.general(m_vl, "  VLENY = %2d  Nyv  = %d\n", VLENY, m_Nyv);
  vout.general(m_vl, "  VLEN  = %2d  Nstv = %d\n", VLEN, m_Nstv);

  m_Nsize[0] = m_Nxv;
  m_Nsize[1] = m_Nyv;
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
void AFopr_Domainwall_5din<AFIELD>::tidyup()
{
  // need to do nothing.
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din<AFIELD>::setup_channels()
{
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
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din<AFIELD>::set_parameters(
  const Parameters& params)
{
  string           gmset_type;
  double           mq, M0;
  int              Ns;
  std::vector<int> bc;
  double           b, c;

  int err_optional = 0;
  err_optional += params.fetch_string("gamma_matrix_type", gmset_type);
  if (!err_optional) {
    if (gmset_type != m_repr) {
      vout.crucial(m_vl, "%s: unsupported gamma_matrix_type: %s\n",
                   class_name.c_str(), gmset_type.c_str());
      exit(EXIT_FAILURE);
    }
  }

  std::string kernel_type;
  err_optional += params.fetch_string("kernel_type", kernel_type);
  if (!err_optional) {
    if (kernel_type != "Wilson") {
      vout.crucial(m_vl, "%s: unsupported kernel_type: %s\n",
                   class_name.c_str(), kernel_type.c_str());
      exit(EXIT_FAILURE);
    }
  }
  m_kernel_type = kernel_type;

  int err = 0;
  err += params.fetch_double("quark_mass", mq);
  err += params.fetch_double("domain_wall_height", M0);
  err += params.fetch_int("extent_of_5th_dimension", Ns);
  err += params.fetch_int_vector("boundary_condition", bc);
  err += params.fetch_double("coefficient_b", b);
  err += params.fetch_double("coefficient_c", c);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(mq, M0, Ns, bc, b, c);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din<AFIELD>::set_parameters(
  const double mq,
  const double M0,
  const int Ns,
  const vector<int> bc,
  const double b,
  const double c)
{
  assert(bc.size() == m_Ndim);

#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  if (ith == 0) {
    m_M0   = real_t(M0);
    m_mq   = real_t(mq);
    m_Ns   = Ns;
    m_NinF = m_Nvcd * m_Ns;

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

  vout.general(m_vl, "Parameters of %s:\n", class_name.c_str());
  vout.general(m_vl, "  mq   = %8.4f\n", double(m_mq));
  vout.general(m_vl, "  M0   = %8.4f\n", double(m_M0));
  vout.general(m_vl, "  Ns   = %4d\n", m_Ns);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }
  vout.general(m_vl, "  coefficients b = %16.10f  c = %16.10f\n",
               double(m_b[0]), double(m_c[0]));

  set_precond_parameters();

  // 5D working vectors.
  if (m_w1.nin() != m_NinF) {
    if (ith == 0) {
      m_w1.reset(m_NinF, m_Nvol, 1);
      m_v1.reset(m_NinF, m_Nvol, 1);
      m_v2.reset(m_NinF, m_Nvol, 1);
    }
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din<AFIELD>::get_parameters(
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

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din<AFIELD>::set_precond_parameters()
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

    for (int is = 0; is < m_Ns - 1; ++is) {
      m_dpinv[is] = 1.0 / m_dp[is];
    }
    m_dpinv[m_Ns - 1] = 1.0 / (m_dp[m_Ns - 1] + m_g);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din<AFIELD>::set_coefficients(
                                   const std::vector<double> vec_b,
                                   const std::vector<double> vec_c)
{
#pragma omp barrier

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
    }
  }

#pragma omp barrier

  for (int is = 0; is < m_Ns; ++is) {
    vout.general(m_vl, "b[%2d] = %16.10f  c[%2d] = %16.10f\n",
                 is, double(m_b[is]), is, double(m_c[is]));
  }

  set_precond_parameters();
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din<AFIELD>::set_config(Field *u)
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
void AFopr_Domainwall_5din<AFIELD>::set_config_omp(Field *u)
{
  vout.detailed(m_vl, "  set_config_omp is called.\n");

#pragma omp parallel
  {
    set_config_impl(u);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din<AFIELD>::set_config_impl(Field *u)
{
#pragma omp barrier

  AIndex_lex<real_t, AFIELD::IMPL> index_lex;
  convert_gauge(index_lex, m_U, *u);

  QXS_Gauge::set_boundary(m_U, m_boundary);

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din<AFIELD>::convert(AFIELD& v,
                                            const Field& w)
{ // the following implementation only applies Dirac repr.
#pragma omp barrier

  int ith, nth, isite, nsite;
  set_threadtask_mult(ith, nth, isite, nsite, m_Nvol);

  AIndex_lex<real_t, AFIELD::IMPL> index;

  for (int site = isite; site < nsite; ++site) {
    for (int is = 0; is < m_Ns; ++is) {
      for (int ic = 0; ic < NC; ++ic) {
        for (int id = 0; id < ND2; ++id) {
          for (int iri = 0; iri < 2; ++iri) {
            int in_org = iri + 2 * (ic + NC * id);
            int in_alt = iri + 2 * (id + ND * ic) + NVCD * is;
            v.set(index.idx(in_alt, m_NinF, site, 0),
                  real_t(w.cmp(in_org, site, is)));
          }
        }

        for (int id = ND2; id < ND; ++id) {
          for (int iri = 0; iri < 2; ++iri) {
            int in_org = iri + 2 * (ic + NC * id);
            int in_alt = iri + 2 * (id + ND * ic) + NVCD * is;
            v.set(index.idx(in_alt, m_NinF, site, 0),
                  real_t(-w.cmp(in_org, site, is)));
          }
        }
      }
    }
  } // site-llop

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din<AFIELD>::reverse(Field& v,
                                            const AFIELD& w)
{ // the following implementation only applies Dirac repr.
#pragma omp barrier

  int ith, nth, isite, nsite;
  set_threadtask_mult(ith, nth, isite, nsite, m_Nvol);

  AIndex_lex<real_t, AFIELD::IMPL> index;

  for (int site = isite; site < nsite; ++site) {
    for (int is = 0; is < m_Ns; ++is) {
      for (int ic = 0; ic < NC; ++ic) {
        for (int id = 0; id < ND2; ++id) {
          for (int iri = 0; iri < 2; ++iri) {
            int in_org = iri + 2 * (ic + NC * id);
            int in_alt = iri + 2 * (id + ND * ic) + NVCD * is;
            v.set(in_org, site, is,
                  double(  w.cmp(index.idx(in_alt, m_NinF, site, 0))));
          }
        }

        for (int id = ND2; id < ND; ++id) {
          for (int iri = 0; iri < 2; ++iri) {
            int in_org = iri + 2 * (ic + NC * id);
            int in_alt = iri + 2 * (id + ND * ic) + NVCD * is;
            v.set(in_org, site, is,
                  double( -w.cmp(index.idx(in_alt, m_NinF, site, 0))));
          }
        }
      }
    }
  } // site-loopa

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din<AFIELD>::set_mode(std::string mode)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_mode = mode;

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din<AFIELD>::mult(AFIELD& v, const AFIELD& w,
                                         std::string mode)
{
  assert(w.check_size(m_NinF, m_Nvol, 1));
  assert(v.check_size(m_NinF, m_Nvol, 1));

  if (mode == "D") {
    D(v, w);
  } else if (mode == "Ddag") {
    Ddag(v, w);
  } else if (mode == "DdagD") {
    DdagD(v, w);
  } else if (mode == "D_prec") {
    D_prec(v, w);
  } else if (mode == "Ddag_prec") {
    Ddag_prec(v, w);
  } else if (mode == "DdagD_prec") {
    DdagD_prec(v, w);
  } else if (mode == "Prec") {
    Prec(v, w);
  } else if (mode == "Precdag") {
    Precdag(v, w);
  } else {
    vout.crucial(m_vl, "%s: undefined mode = %s\n",
                 class_name.c_str(), mode.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din<AFIELD>::mult_dag(AFIELD& v,
                                             const AFIELD& w,
                                             std::string mode)
{
  assert(w.check_size(m_NinF, m_Nvol, 1));
  assert(v.check_size(m_NinF, m_Nvol, 1));

  if (mode == "D") {
    Ddag(v, w);
  } else if (mode == "Ddag") {
    D(v, w);
  } else if (mode == "DdagD") {
    DdagD(v, w);
  } else if (mode == "D_prec") {
    Ddag_prec(v, w);
  } else if (mode == "Ddag_prec") {
    D_prec(v, w);
  } else if (mode == "DdagD_prec") {
    DdagD_prec(v, w);
  } else if (mode == "Prec") {
    Precdag(v, w);
  } else if (mode == "Precdag") {
    Prec(v, w);
  } else {
    std::cout << "mode undeifined in AFopr_Domainwall_5din.\n";
    abort();
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din<AFIELD>::mult_gm5(AFIELD& v, const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

#pragma omp barrier

  BridgeQXS::mult_domainwall_5din_mult_gm5_dirac(vp, wp, m_Ns, m_Nsize);

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din<AFIELD>::DdagD(AFIELD& v, const AFIELD& w)
{
  D(m_v1, w);
  Ddag(v, m_v1);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din<AFIELD>::DdagD_prec(AFIELD& v,
                                               const AFIELD& w)
{
  D_prec(m_v1, w);
  Ddag_prec(v, m_v1);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din<AFIELD>::D_prec(AFIELD& v,
                                           const AFIELD& w)
{
#pragma omp barrier
  L_inv(v, w);
  U_inv(m_v2, v);
  D(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din<AFIELD>::Ddag_prec(AFIELD& v,
                                              const AFIELD& w)
{
  Ddag(v, w);
  Udag_inv(m_v2, v);
  Ldag_inv(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din<AFIELD>::Prec(AFIELD& v,
                                         const AFIELD& w)
{
#pragma omp barrier
  L_inv(m_v2, w);
  U_inv(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din<AFIELD>::Precdag(AFIELD& v,
                                            const AFIELD& w)
{
#pragma omp barrier
  Udag_inv(m_v2, w);
  Ldag_inv(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din<AFIELD>::D(AFIELD& v, const AFIELD& w)
{
  mult_D(v, w);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din<AFIELD>::Ddag(AFIELD& v, const AFIELD& w)
{
  mult_Ddag(v, w);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din<AFIELD>::mult_D(AFIELD& v, const AFIELD& w)
{
  int Nin4 = VLEN * NVCD;
  int Nin5 = Nin4 * m_Ns;

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);
  real_t *yp = m_w1.ptr(0);  // working vector
  real_t *up = m_U.ptr(0);

  int ith = ThreadManager::get_thread_id();

#pragma omp barrier

  BridgeQXS::mult_domainwall_5din_5dir_dirac(vp, yp, wp,
                                             m_mq, m_M0, m_Ns, &m_boundary[0],
                                             &m_b[0], &m_c[0],
                                             m_Nsize, do_comm);

#pragma omp barrier

  if (do_comm_any > 0) {
    if (ith == 0) chset_recv.start();

    real_t *buf1_xp = (real_t *)chsend_dn[0].ptr();
    real_t *buf1_xm = (real_t *)chsend_up[0].ptr();
    real_t *buf1_yp = (real_t *)chsend_dn[1].ptr();
    real_t *buf1_ym = (real_t *)chsend_up[1].ptr();
    real_t *buf1_zp = (real_t *)chsend_dn[2].ptr();
    real_t *buf1_zm = (real_t *)chsend_up[2].ptr();
    real_t *buf1_tp = (real_t *)chsend_dn[3].ptr();
    real_t *buf1_tm = (real_t *)chsend_up[3].ptr();

    BridgeQXS::mult_domainwall_5din_hop1_dirac(
      buf1_xp, buf1_xm, buf1_yp, buf1_ym,
      buf1_zp, buf1_zm, buf1_tp, buf1_tm,
      up, yp,
      m_mq, m_M0, m_Ns, &m_boundary[0],
      m_Nsize, do_comm);

#pragma omp barrier

    if (ith == 0) chset_send.start();
  }

  BridgeQXS::mult_domainwall_5din_hopb_dirac(vp, up, yp,
                                             m_mq, m_M0, m_Ns, &m_boundary[0],
                                             &m_b[0], &m_c[0],
                                             m_Nsize, do_comm);

  if (do_comm_any > 0) {
    if (ith == 0) chset_recv.wait();

#pragma omp barrier

    real_t *buf2_xp = (real_t *)chrecv_up[0].ptr();
    real_t *buf2_xm = (real_t *)chrecv_dn[0].ptr();
    real_t *buf2_yp = (real_t *)chrecv_up[1].ptr();
    real_t *buf2_ym = (real_t *)chrecv_dn[1].ptr();
    real_t *buf2_zp = (real_t *)chrecv_up[2].ptr();
    real_t *buf2_zm = (real_t *)chrecv_dn[2].ptr();
    real_t *buf2_tp = (real_t *)chrecv_up[3].ptr();
    real_t *buf2_tm = (real_t *)chrecv_dn[3].ptr();

    BridgeQXS::mult_domainwall_5din_hop2_dirac(vp, up, yp,
                                               buf2_xp, buf2_xm, buf2_yp, buf2_ym,
                                               buf2_zp, buf2_zm, buf2_tp, buf2_tm,
                                               m_mq, m_M0, m_Ns, &m_boundary[0],
                                               m_Nsize, do_comm);

    if (ith == 0) chset_send.wait();
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din<AFIELD>::mult_Ddag(AFIELD& v, const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);
  real_t *yp = m_w1.ptr(0);
  real_t *up = m_U.ptr(0);

  int ith = ThreadManager::get_thread_id();

  BridgeQXS::mult_domainwall_5din_mult_gm5_dirac(vp, wp, m_Ns, m_Nsize);

#pragma omp barrier

  if (do_comm_any > 0) {
    if (ith == 0) chset_recv.start();

    real_t *buf1_xp = (real_t *)chsend_dn[0].ptr();
    real_t *buf1_xm = (real_t *)chsend_up[0].ptr();
    real_t *buf1_yp = (real_t *)chsend_dn[1].ptr();
    real_t *buf1_ym = (real_t *)chsend_up[1].ptr();
    real_t *buf1_zp = (real_t *)chsend_dn[2].ptr();
    real_t *buf1_zm = (real_t *)chsend_up[2].ptr();
    real_t *buf1_tp = (real_t *)chsend_dn[3].ptr();
    real_t *buf1_tm = (real_t *)chsend_up[3].ptr();

    BridgeQXS::mult_domainwall_5din_hop1_dirac(
      buf1_xp, buf1_xm, buf1_yp, buf1_ym,
      buf1_zp, buf1_zm, buf1_tp, buf1_tm,
      up, vp,
      m_mq, m_M0, m_Ns, &m_boundary[0],
      m_Nsize, do_comm);

#pragma omp barrier

    if (ith == 0) chset_send.start();
  }

  BridgeQXS::mult_domainwall_5din_clear(yp, m_Ns, m_Nsize);

  BridgeQXS::mult_domainwall_5din_hopb_dirac(yp, up, vp,
                                             m_mq, m_M0, m_Ns, &m_boundary[0],
                                             &m_b[0], &m_c[0],
                                             m_Nsize, do_comm);

  if (do_comm_any > 0) {
    if (ith == 0) chset_recv.wait();

#pragma omp barrier

    real_t *buf2_xp = (real_t *)chrecv_up[0].ptr();
    real_t *buf2_xm = (real_t *)chrecv_dn[0].ptr();
    real_t *buf2_yp = (real_t *)chrecv_up[1].ptr();
    real_t *buf2_ym = (real_t *)chrecv_dn[1].ptr();
    real_t *buf2_zp = (real_t *)chrecv_up[2].ptr();
    real_t *buf2_zm = (real_t *)chrecv_dn[2].ptr();
    real_t *buf2_tp = (real_t *)chrecv_up[3].ptr();
    real_t *buf2_tm = (real_t *)chrecv_dn[3].ptr();

    BridgeQXS::mult_domainwall_5din_hop2_dirac(yp, up, vp,
                                               buf2_xp, buf2_xm, buf2_yp, buf2_ym,
                                               buf2_zp, buf2_zm, buf2_tp, buf2_tm,
                                               m_mq, m_M0, m_Ns, &m_boundary[0],
                                               m_Nsize, do_comm);

    if (ith == 0) chset_send.wait();
  }

#pragma omp barrier

  BridgeQXS::mult_domainwall_5din_5dirdag_dirac(vp, yp, wp,
                                                m_mq, m_M0, m_Ns, &m_boundary[0],
                                                &m_b[0], &m_c[0],
                                                m_Nsize, do_comm);

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din<AFIELD>::L_inv(AFIELD& v,
                                          const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);
  BridgeQXS::mult_domainwall_5din_L_inv_dirac(vp, wp,
                                              m_Ns, m_Nsize,
                                              &m_e[0], &m_dpinv[0], &m_dm[0]);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din<AFIELD>::U_inv(AFIELD& v,
                                          const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);
  BridgeQXS::mult_domainwall_5din_U_inv_dirac(
    vp, wp, m_Ns, m_Nsize,
    &m_f[0], &m_dpinv[0], &m_dm[0]);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din<AFIELD>::Ldag_inv(AFIELD& v,
                                             const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

  BridgeQXS::mult_domainwall_5din_Ldag_inv_dirac(
    vp, wp, m_Ns, m_Nsize,
    &m_e[0], &m_dpinv[0], &m_dm[0]);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din<AFIELD>::Udag_inv(AFIELD& v,
                                             const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

  BridgeQXS::mult_domainwall_5din_Udag_inv_dirac(
    vp, wp, m_Ns, m_Nsize,
    &m_f[0], &m_dpinv[0], &m_dm[0]);
}


//====================================================================
template<typename AFIELD>
double AFopr_Domainwall_5din<AFIELD>::flop_count(std::string mode)
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


//============================================================END=====
