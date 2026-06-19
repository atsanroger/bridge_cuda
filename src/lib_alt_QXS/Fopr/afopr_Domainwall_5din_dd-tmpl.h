/*!
      @file    afopr_Domainwall_5din_dd-tmpl.h
      @brief
      @author  Issaku Kanamori (kanamori)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

template<typename AFIELD>
const std::string AFopr_Domainwall_5din_dd<AFIELD>::class_name
  = "AFopr_Domainwall_5din_dd";
//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::init(const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  int req_comm = 1;  // set 1 if communication forced any time
  //int req_comm = 0;  // set 1 if communication forced any time

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

  // gauge configuration with block condition.
  m_Ublock.reset(NDF, m_Nvol, m_Ndim);

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
  // need to do nothing.
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::setup_channels()
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
void AFopr_Domainwall_5din_dd<AFIELD>::set_parameters(
  const Parameters& params)
{
  string           gmset_type;
  double           mq, M0;
  int              Ns;
  std::vector<int> bc;
  double           b, c;
  std::vector<int> block_size;

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
  err += params.fetch_int_vector("block_size", block_size);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(mq, M0, Ns, bc, b, c, block_size);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::set_parameters(
  const double mq,
  const double M0,
  const int Ns,
  const vector<int> bc,
  const double b,
  const double c,
  const std::vector<int> block_size)
{
  assert(bc.size() == m_Ndim);

#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  if (ith == 0) {
    m_M0   = real_t(M0);
    m_mq   = real_t(mq);
    m_Ns   = Ns;
    m_NinF = m_Nvcd * m_Ns;

    if (m_boundary.size() != m_Ndim) { m_boundary.resize(m_Ndim); }
    for (int mu = 0; mu < m_Ndim; ++mu) {
      m_boundary[mu] = bc[mu];
    }

    if (m_block_size.size() != m_Ndim) { m_block_size.resize(m_Ndim); }
    for (int mu = 0; mu < m_Ndim; ++mu) {
      m_block_size[mu] = block_size[mu];
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
  vout.general(m_vl, "  mq   = %8.4f\n", m_mq);
  vout.general(m_vl, "  M0   = %8.4f\n", m_M0);
  vout.general(m_vl, "  Ns   = %4d\n", m_Ns);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }
  vout.general(m_vl, "  coefficients b = %16.10f  c = %16.10f\n",
               m_b[0], m_c[0]);

  set_precond_parameters();

  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  block_size[%d] = %2d\n", mu, m_block_size[mu]);
  }

  int rem = m_Nx % m_block_size[0]
            + m_Ny % m_block_size[1]
            + m_Nz % m_block_size[2]
            + m_Nt % m_block_size[3];
  if (rem != 0) {
    vout.crucial(m_vl, "%s: block_size is irelevant.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  if (ith == 0) {
    // note that m_block_sizev is simple array
    m_block_sizev[0] = m_block_size[0] / VLENX;
    m_block_sizev[1] = m_block_size[1] / VLENY;
    m_block_sizev[2] = m_block_size[2];
    m_block_sizev[3] = m_block_size[3];
  }
#pragma omp barrier

  if ((m_block_sizev[0] * VLENX != m_block_size[0]) ||
      (m_block_sizev[1] * VLENY != m_block_size[1])) {
    vout.crucial(m_vl, "%s: bad blocksize in XY: must be divided by"
                       " VLENX=%d, VLENY=%d but are %d, %d\n",
                 class_name.c_str(), VLENX, VLENY,
                 m_block_size[0], m_block_size[1]);
    exit(EXIT_FAILURE);
  }

  if (ith == 0) {
    int NBx  = m_Nx / m_block_size[0];
    int NBy  = m_Ny / m_block_size[1];
    int NBz  = m_Nz / m_block_size[2];
    int NBt  = m_Nt / m_block_size[3];
    int ipex = Communicator::ipe(0);
    int ipey = Communicator::ipe(1);
    int ipez = Communicator::ipe(2);
    int ipet = Communicator::ipe(3);
    m_Ieo = (NBx * ipex + NBy * ipey + NBz * ipez + NBt * ipet) % 2;
  }

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
void AFopr_Domainwall_5din_dd<AFIELD>::set_coefficients(
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
                 is, m_b[is], is, m_c[is]);
  }

  set_precond_parameters();
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

  AIndex_lex<real_t, AFIELD::IMPL> index_lex;
  convert_gauge(index_lex, m_U, *u);

  QXS_Gauge::set_boundary(m_U, m_boundary);

  copy(m_Ublock, m_U);
#pragma omp barrier

  set_block_config(m_Ublock);

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::set_block_config(AFIELD& U)
{
#pragma omp barrier

  AIndex_lex<real_t, AFIELD::IMPL> index;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nvol);

  for (int site = is; site < ns; ++site) {
    int ix = site % m_Nx;
    int iy = (site / m_Nx) % m_Ny;
    int iz = (site / (m_Nx * m_Ny)) % m_Nz;
    int it = site / (m_Nx * m_Ny * m_Nz);

    int mu = 0;
    if ((ix + 1) % m_block_size[mu] == 0) {
      for (int in = 0; in < NDF; ++in) {
        int i = index.idx_G(in, site, mu);
        U.set(i, real_t(0.0));
      }
    }

    mu = 1;
    if ((iy + 1) % m_block_size[mu] == 0) {
      for (int in = 0; in < NDF; ++in) {
        int i = index.idx_G(in, site, mu);
        U.set(i, real_t(0.0));
      }
    }

    mu = 2;
    if ((iz + 1) % m_block_size[mu] == 0) {
      for (int in = 0; in < NDF; ++in) {
        int i = index.idx_G(in, site, mu);
        U.set(i, real_t(0.0));
      }
    }

    mu = 3;
    if ((it + 1) % m_block_size[mu] == 0) {
      for (int in = 0; in < NDF; ++in) {
        int i = index.idx_G(in, site, mu);
        U.set(i, real_t(0.0));
      }
    }
  }

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::convert(AFIELD& v,
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
void AFopr_Domainwall_5din_dd<AFIELD>::reverse(Field& v,
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
void AFopr_Domainwall_5din_dd<AFIELD>::set_mode(std::string mode)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_mode = mode;

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::mult(AFIELD& v, const AFIELD& w,
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
void AFopr_Domainwall_5din_dd<AFIELD>::mult_dag(AFIELD& v,
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
    std::cout << "mode undeifined in AFopr_Domainwall_5din_dd.\n";
    abort();
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::mult_gm5(AFIELD& v, const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

#pragma omp barrier

  BridgeQXS::mult_domainwall_5din_mult_gm5_dirac(vp, wp, m_Ns, m_Nsize);

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::scal_local(real_t *vp, real_t a)
{
  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  svbool_t pg = set_predicate();
  for (int site = is; site < ns; ++site) {
    for(int is=0; is < m_Ns; ++is){
      real_t *vp2 = vp + VLEN*site*m_Ns*NVCD + VLEN*NVCD*is;
      for(int in=0; in < NVCD; ++in){
        svreal_t val;
        load_vec(pg, val, vp2 + VLEN*in);
        scal_vec(pg, val, a);
        save_vec(pg, vp2 + VLEN*in, val);
      }
    }
  }
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::DdagD(AFIELD& v, const AFIELD& w)
{
  D(m_v1, w);
  Ddag(v, m_v1);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::DdagD_prec(AFIELD& v,
                                               const AFIELD& w)
{
  D_prec(m_v1, w);
  Ddag_prec(v, m_v1);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::D_prec(AFIELD& v,
                                           const AFIELD& w)
{
#pragma omp barrier
  L_inv(v, w);
  U_inv(m_v2, v);
  D(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::Ddag_prec(AFIELD& v,
                                              const AFIELD& w)
{
  Ddag(v, w);
  Udag_inv(m_v2, v);
  Ldag_inv(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::Prec(AFIELD& v,
                                         const AFIELD& w)
{
#pragma omp barrier
  L_inv(m_v2, w);
  U_inv(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::Precdag(AFIELD& v,
                                            const AFIELD& w)
{
#pragma omp barrier
  Udag_inv(m_v2, w);
  Ldag_inv(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::D(AFIELD& v, const AFIELD& w)
{
  mult_D(v, w);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::Ddag(AFIELD& v, const AFIELD& w)
{
  mult_Ddag(v, w);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::mult_D(AFIELD& v, const AFIELD& w)
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
void AFopr_Domainwall_5din_dd<AFIELD>::mult_Ddag(AFIELD& v, const AFIELD& w)
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
void AFopr_Domainwall_5din_dd<AFIELD>::mult_sap(AFIELD& v, const AFIELD& w,
                                                const int ieo)
{

  if (m_mode == "D") {
    mult_D_sap(v, w, ieo);
  } else if (m_mode == "Ddag") {
    mult_Ddag_sap(v, w, ieo);
  } else if (m_mode == "DdagD") {
    mult_D_sap(m_v1, w, ieo);
    mult_Ddag_sap(v, m_v1, ieo);
  } else {
    std::cout << "mode undeifined for mult_sap in AFopr_Domainwall_5din_dd.\n";
    abort();
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::mult_D_sap(AFIELD& v, const AFIELD& w,
                                                  const int ieo)
{
  int Nin4 = VLEN * NVCD;
  int Nin5 = Nin4 * m_Ns;

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);
  real_t *yp = m_w1.ptr(0);  // working vector
  real_t *up = m_Ublock.ptr(0);

  int jeo = (ieo + m_Ieo) % 2;

#pragma omp barrier
  BridgeQXS::mult_domainwall_5din_dd_5dir_dirac(vp, yp, wp,
                                                m_mq, m_M0, m_Ns, &m_boundary[0],
                                                &m_b[0], &m_c[0],
                                                m_Nsize, m_block_sizev, jeo);
#pragma omp barrier
  BridgeQXS::mult_domainwall_5din_dd_hopb_dirac(vp, up, yp,
                                                m_mq, m_M0, m_Ns, &m_boundary[0],
                                                &m_b[0], &m_c[0],
                                                m_Nsize, m_block_sizev, jeo);
#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::mult_Ddag_sap(AFIELD& v, const AFIELD& w,
                                                  const int ieo)
{
  int Nin4 = VLEN * NVCD;
  int Nin5 = Nin4 * m_Ns;

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);
  real_t *yp = m_w1.ptr(0);  // working vector
  real_t *up = m_Ublock.ptr(0);

  int jeo = (ieo + m_Ieo) % 2;

  BridgeQXS::mult_domainwall_5din_dd_mult_gm5_dirac(vp, wp,
                                                    m_Ns, m_Nsize,
                                                    m_block_sizev, jeo);
  BridgeQXS::mult_domainwall_5din_clear(yp, m_Ns, m_Nsize);

#pragma omp barrier
  BridgeQXS::mult_domainwall_5din_dd_hopb_dirac(yp, up, vp,
                                                m_mq, m_M0, m_Ns, &m_boundary[0],
                                                &m_b[0], &m_c[0],
                                                m_Nsize, m_block_sizev, jeo);
#pragma omp barrier
  BridgeQXS::mult_domainwall_5din_dd_5dirdag_dirac(vp, yp, wp,
                                                   m_mq, m_M0, m_Ns, &m_boundary[0],
                                                   &m_b[0], &m_c[0],
                                                   m_Nsize, m_block_sizev, jeo);
#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::mult_dd(AFIELD& v, const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);
  real_t *yp = m_w1.ptr(0);  // working vector
  real_t *zp = m_v2.ptr(0);  // working vector
  real_t *up = m_Ublock.ptr(0);

#pragma omp barrier

  if (m_mode == "D") {

    BridgeQXS::mult_domainwall_5din_dd_5dir_dirac(
                      vp, yp, wp, m_mq, m_M0, m_Ns, &m_boundary[0],
                      &m_b[0], &m_c[0], m_Nsize, m_block_sizev, -1);
#pragma omp barrier
    BridgeQXS::mult_domainwall_5din_dd_hopb_dirac(
                      vp, up, yp, m_mq, m_M0, m_Ns, &m_boundary[0],
                      &m_b[0], &m_c[0], m_Nsize, m_block_sizev,-1);

  } else if (m_mode == "D_prec") {
    L_inv(v, w);
    U_inv(m_v2, v);

    BridgeQXS::mult_domainwall_5din_dd_5dir_dirac(
                      vp, yp, zp, m_mq, m_M0, m_Ns, &m_boundary[0],
                      &m_b[0], &m_c[0], m_Nsize, m_block_sizev, -1);
#pragma omp barrier
    BridgeQXS::mult_domainwall_5din_dd_hopb_dirac(
                      vp, up, yp, m_mq, m_M0, m_Ns, &m_boundary[0],
                      &m_b[0], &m_c[0], m_Nsize, m_block_sizev,-1);
  }

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::mult_dup(AFIELD& v, const AFIELD& w, const int mu){

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);
  real_t *up = m_U.ptr(0);
  real_t *ublockp = m_Ublock.ptr(0);

  real_t *tmp = m_v1.ptr(0);
#pragma omp barrier

  if (m_mode == "D") {

    BridgeQXS::mult_domainwall_5din_clear(vp, m_Ns, m_Nsize);

    if (mu == 0) {
      mult_xp(vp, wp, ublockp);
      scal_local(vp, real_t(-1.0));
      mult_xp(vp, wp, up);
    } else if (mu == 1) {
      mult_yp(vp, wp, ublockp);
      scal_local(vp, real_t(-1.0));
      mult_yp(vp, wp, up);
    } else if (mu == 2) {
      mult_zp(vp, wp, ublockp);
      scal_local(vp, real_t(-1.0));
      mult_zp(vp, wp, up);
    } else if (mu == 3) {
      mult_tp(vp, wp, ublockp);
      scal_local(vp, real_t(-1.0));
      mult_tp(vp, wp, up);
    }

  } else if (m_mode == "D_prec") {

    L_inv(v, w);
    U_inv(m_v2, v);

    real_t *zp = m_v2.ptr(0);

    BridgeQXS::mult_domainwall_5din_clear(vp, m_Ns, m_Nsize);

    if (mu == 0) {
      mult_xp(vp, zp, ublockp);
      scal_local(vp, real_t(-1.0));
      mult_xp(vp, zp, up);
    } else if (mu == 1) {
      mult_yp(vp, zp, ublockp);
      scal_local(vp, real_t(-1.0));
      mult_yp(vp, zp, up);
    } else if (mu == 2) {
      mult_zp(vp, zp, ublockp);
      scal_local(vp, real_t(-1.0));
      mult_zp(vp, zp, up);
    } else if (mu == 3) {
      mult_tp(vp, zp, ublockp);
      scal_local(vp, real_t(-1.0));
      mult_tp(vp, zp, up);
    }

  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::mult_ddn(AFIELD& v, const AFIELD& w, const int mu){

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);
  real_t *up = m_U.ptr(0);
  real_t *ublockp = m_Ublock.ptr(0);

#pragma omp barrier

  if (m_mode == "D") {

    BridgeQXS::mult_domainwall_5din_clear(vp, m_Ns, m_Nsize);

    if (mu == 0) {
      mult_xm(vp, wp, ublockp);
      scal_local(vp, real_t(-1.0));
      mult_xm(vp, wp, up);
    } else if (mu == 1) {
      mult_ym(vp, wp, ublockp);
      scal_local(vp, real_t(-1.0));
      mult_ym(vp, wp, up);
    } else if (mu == 2) {
      mult_zm(vp, wp, ublockp);
      scal_local(vp, real_t(-1.0));
      mult_zm(vp, wp, up);
    } else if (mu == 3) {
      mult_tm(vp, wp, ublockp);
      scal_local(vp, real_t(-1.0));
      mult_tm(vp, wp, up);
    }

  } else if (m_mode == "D_prec") {

    L_inv(v, w);
    U_inv(m_v2, v);

    real_t *zp = m_v2.ptr(0);

    BridgeQXS::mult_domainwall_5din_clear(vp, m_Ns, m_Nsize);

    if (mu == 0) {
      mult_xm(vp, zp, ublockp);
      scal_local(vp, real_t(-1.0));
      mult_xm(vp, zp, up);
    } else if (mu == 1) {
      mult_ym(vp, zp, ublockp);
      scal_local(vp, real_t(-1.0));
      mult_ym(vp, zp, up);
    } else if (mu == 2) {
      mult_zm(vp, zp, ublockp);
      scal_local(vp, real_t(-1.0));
      mult_zm(vp, zp, up);
    } else if (mu == 3) {
      mult_tm(vp, zp, ublockp);
      scal_local(vp, real_t(-1.0));
      mult_tm(vp, zp, up);
    }

  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::mult_xp(real_t *vp, real_t *wp,
                                               real_t *up){

  constexpr int idir = 0;

  real_t *yp = m_w1.ptr(0);  // working vector
  //  real_t *up = m_U.ptr(0);

  int ith = ThreadManager::get_thread_id();

  constexpr int Nin4 = VLEN * NVCD;
  const int Ns   = m_Ns;
  const int Nin5 = Nin4 * Ns;

  int Nxv = m_Nsize[0];
  int Nyv = m_Nsize[1];
  int Nz  = m_Nsize[2];
  int Nt  = m_Nsize[3];
  int Nstv = Nxv * Nyv * Nz * Nt;
  int Nst  = Nstv * VLEN;
  int NvU  = NDF * Nst;

  svbool_t pg1_xp, pg2_xp, pg1_xm, pg2_xm;
  set_predicate_xp(pg1_xp, pg2_xp);
  set_predicate_xm(pg1_xm, pg2_xm);

  //  svint_t svidx_xp, svidx_xm; // do not work here
  //  set_index_xp(svidx_xp);
  //  set_index_xp(svidx_xm);

#pragma omp barrier
  BridgeQXS::accum_mult_domainwall_5din_5dir_dirac(vp, yp, wp,
                                                   m_mq, m_M0, m_Ns,
                                                   &m_boundary[0],
                                                   &m_b[0], &m_c[0],
                                                   m_Nsize, do_comm);
#pragma omp barrier

  if (do_comm[idir] == 1) {
    if (ith == 0) { chrecv_up[idir].start(); }

    real_t *buf1_xp = (real_t *)chsend_dn[idir].ptr();

    int Nyzt = Nyv * Nz * Nt;

    int ith_dummy, nth, site0, site1;
    set_threadtask(ith_dummy, nth, site0, site1, Nyzt);

    for (int iyzt = site0; iyzt < site1; ++iyzt) {
      int    ix   = 0;
      int    site = ix + Nxv * iyzt;
      real_t *yp2 = &yp[Nin5 * site];
      int    ibf  = VLENY * NVC * ND2 * Ns * iyzt;
      svint_t svidx_xm;
      set_index_xm(svidx_xm);  // need to be here
      for (int is = 0; is < Ns; ++is) {
        mult_wilson_xp1(pg2_xm, svidx_xm,
                        &buf1_xp[ibf + VLENY * NVC * ND2 * is], &yp2[Nin4 * is]);
      }
    }
#pragma omp barrier
    if (ith == 0) {
      chsend_dn[idir].start();
      chrecv_up[idir].wait();
      chsend_dn[idir].wait();
    }
#pragma omp barrier
  }

  int ith_dummy, nth, site0, site1;
  set_threadtask(ith_dummy, nth, site0, site1, Nstv);

  for (int site = site0; site < site1; ++site) {
    int ix   = site % Nxv;
    int iyzt = site / Nxv;

    real_t *yp2 = &yp[Nin5 * site];
    real_t *vp2 = &vp[Nin5 * site];

    Vsimd_t vL[NVCD * Ns];
    load_vec(vL, vp2, NVCD * Ns);

    if ((ix < Nxv - 1) || (do_comm[idir] == 0)) {
      int ix2 = (ix + 1) % Nxv;
      int nei = ix2 + Nxv * iyzt;
      real_t *ypn = &yp[Nin5 * nei];
      real_t *u   = &up[VLEN * NDF * site + NvU * idir];
      for (int is = 0; is < Ns; ++is) {
        mult_wilson_xpb(pg1_xp, pg2_xp, &vL[NVCD * is].v[0], u,
                        &yp2[Nin4 * is], &ypn[Nin4 * is]);
      }
    } else {
      real_t *buf2_xp = (real_t *)chrecv_up[idir].ptr();
      real_t *u2 = &up[VLEN * NDF * site + NvU * idir];
      int    ibf = VLENY * NVC * ND2 * Ns * iyzt;
      for (int is = 0; is < Ns; ++is) {
        svint_t svidx_xp;
        set_index_xp(svidx_xp); // need to be here
        mult_wilson_xp2(pg1_xp, pg2_xp, svidx_xp,
                        &vL[NVCD * is].v[0], u2,
                        &yp2[Nin4 * is], &buf2_xp[ibf + VLENY * NVC * ND2 * is]);
      }
    }

    save_vec(vp2, vL, NVCD*Ns);

  }

#pragma omp barrier

}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::mult_xm(real_t *vp, real_t *wp,
                                               real_t *up){

  constexpr int idir = 0;

  real_t *yp = m_w1.ptr(0);  // working vector
  //  real_t *up = m_U.ptr(0);

  int ith = ThreadManager::get_thread_id();

  constexpr int Nin4 = VLEN * NVCD;
  const int Ns   = m_Ns;
  const int Nin5 = Nin4 * Ns;

  int Nxv = m_Nsize[0];
  int Nyv = m_Nsize[1];
  int Nz  = m_Nsize[2];
  int Nt  = m_Nsize[3];
  int Nstv = Nxv * Nyv * Nz * Nt;
  int Nst  = Nstv * VLEN;
  int NvU  = NDF * Nst;

  svbool_t pg1_xp, pg2_xp, pg1_xm, pg2_xm;
  set_predicate_xp(pg1_xp, pg2_xp);
  set_predicate_xm(pg1_xm, pg2_xm);

  //  svint_t svidx_xm, svidx_xp; // do not work here
  //  set_index_xp(svidx_xp);
  //  set_index_xm(svidx_xm);


#pragma omp barrier
  BridgeQXS::accum_mult_domainwall_5din_5dir_dirac(vp, yp, wp,
                                                   m_mq, m_M0, m_Ns,
                                                   &m_boundary[0],
                                                   &m_b[0], &m_c[0],
                                                   m_Nsize, do_comm);
#pragma omp barrier

  if (do_comm[idir] == 1) {
    if (ith == 0) { chrecv_dn[idir].start(); }

    real_t *buf1_xm = (real_t *)chsend_up[idir].ptr();

    int Nyzt = Nyv * Nz * Nt;

    int ith_dummy, nth, site0, site1;
    set_threadtask(ith_dummy, nth, site0, site1, Nyzt);

    for (int iyzt = site0; iyzt < site1; ++iyzt) {
      int    ix   = Nxv - 1;
      int    site = ix + Nxv * iyzt;
      real_t *yp2 = &yp[Nin5 * site];
      int    ibf  = VLENY * NVC * ND2 * Ns * iyzt;
      real_t *u2  = &up[VLEN * NDF * site + NvU * idir];
      svint_t svidx_xp;
      set_index_xp(svidx_xp); // need to be here
      for (int is = 0; is < Ns; ++is) {
        mult_wilson_xm1(pg2_xp, svidx_xp,
                        &buf1_xm[ibf + VLENY * NVC * ND2 * is], u2, &yp2[Nin4 * is]);
      }
    }
#pragma omp barrier
    if (ith == 0) {
      chsend_up[idir].start();
      chrecv_dn[idir].wait();
      chsend_up[idir].wait();
    }
#pragma omp barrier
  }

  int ith_dummy, nth, site0, site1;
  set_threadtask(ith_dummy, nth, site0, site1, Nstv);


  for (int site = site0; site < site1; ++site) {
    int ix   = site % Nxv;
    int iyzt = site / Nxv;

    real_t *yp2 = &yp[Nin5 * site];
    real_t *vp2 = &vp[Nin5 * site];

    Vsimd_t vL[NVCD * Ns];
    load_vec(vL, vp2, NVCD * Ns);

    if ((ix > 0) || (do_comm[idir] == 0)) {
      int ix2 = (ix - 1 + Nxv) % Nxv;
      int nei = ix2 + Nxv * iyzt;
      real_t *ypn = &yp[Nin5 * nei];
      real_t *ux  = &up[VLEN * NDF * site + NvU * idir];
      real_t *un  = &up[VLEN * NDF * nei  + NvU * idir];
      for (int is = 0; is < Ns; ++is) {
        mult_wilson_xmb(pg1_xm, pg2_xm, &vL[NVCD * is], ux, un,
                        &yp2[Nin4 * is], &ypn[Nin4 * is]);
      }
    } else {
      real_t *buf2_xm = (real_t *)chrecv_dn[idir].ptr();
      real_t *u2 = &up[VLEN * NDF * site + NvU * idir];
      int    ibf = VLENY * NVC * ND2 * Ns * iyzt;
      svint_t svidx_xm;
      set_index_xm(svidx_xm); // need to be here
      for (int is = 0; is < Ns; ++is) {
        mult_wilson_xm2(pg1_xm, pg2_xm, svidx_xm,
                        &vL[NVCD * is].v[0], u2,
                        &yp2[Nin4 * is], &buf2_xm[ibf + VLENY * NVC * ND2 * is]);
      }
    }
    save_vec(vp2, vL, NVCD*Ns);

  }

#pragma omp barrier

}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::mult_yp(real_t *vp, real_t *wp,
                                               real_t *up){
  constexpr int idir = 1;

  real_t *yp = m_w1.ptr(0);  // working vector
  //  real_t *up = m_U.ptr(0);

  int ith = ThreadManager::get_thread_id();

  constexpr int Nin4 = VLEN * NVCD;
  const int Ns   = m_Ns;
  const int Nin5 = Nin4 * Ns;

  int Nxv = m_Nsize[0];
  int Nyv = m_Nsize[1];
  int Nz = m_Nsize[2];
  int Nt = m_Nsize[3];
  int Nxy  = Nxv * Nyv;
  int Nstv = Nxv * Nyv * Nz * Nt;
  int Nst = Nstv * VLEN;
  int NvU   = NDF * Nst;

  svbool_t pg1_yp, pg2_yp, pg1_ym, pg2_ym;
  set_predicate_yp(pg1_yp, pg2_yp);
  set_predicate_ym(pg1_ym, pg2_ym);

#pragma omp barrier
  BridgeQXS::accum_mult_domainwall_5din_5dir_dirac(vp, yp, wp,
                                                   m_mq, m_M0, m_Ns,
                                                   &m_boundary[0],
                                                   &m_b[0], &m_c[0],
                                                   m_Nsize, do_comm);
#pragma omp barrier

  if (do_comm[idir] == 1) {
    if (ith == 0) { chrecv_up[idir].start(); }

    real_t *buf1_yp = (real_t *)chsend_dn[idir].ptr();


    int Nxzt = Nxv * Nz * Nt;

    int ith_dummy, nth, site0, site1;
    set_threadtask(ith_dummy, nth, site0, site1, Nxzt);

    for (int ixzt = site0; ixzt < site1; ++ixzt) {
      int ix  = ixzt % Nxv;
      int izt = ixzt / Nxv;

      int    iy   = 0;
      int    site = ix + Nxv * (iy + Nyv * izt);
      real_t *yp2 = &yp[Nin5 * site];

      int ibf = VLENX * NVC * ND2 * Ns * (ix + Nxv * izt);
      for (int is = 0; is < Ns; ++is) {
        mult_wilson_yp1(pg2_ym,
                        &buf1_yp[ibf + VLENX * NVC * ND2 * is], &yp2[Nin4 * is]);
      }
    }
#pragma omp barrier
    if (ith == 0) {
      chsend_dn[idir].start();
      chrecv_up[idir].wait();
      chsend_dn[idir].wait();
    }
#pragma omp barrier
  }

  int ith_dummy, nth, site0, site1;
  set_threadtask(ith_dummy, nth, site0, site1, Nstv);

  for (int site = site0; site < site1; ++site) {
    int ix   = site % Nxv;
    int iyzt = site / Nxv;
    int iy   = iyzt % Nyv;
    int izt  = site / Nxy;

    real_t *yp2 = &yp[Nin5 * site];
    real_t *vp2 = &vp[Nin5 * site];

    Vsimd_t vL[NVCD * Ns];
    load_vec(vL, vp2, NVCD * Ns);

    if ((iy < Nyv - 1) || (do_comm[idir] == 0)) {
      int iy2 = (iy + 1) % Nyv;
      int nei = ix + Nxv * (iy2 + Nyv * izt);
      real_t *ypn = &yp[Nin5 * nei];
      real_t *u   = &up[VLEN * NDF * site + NvU * idir];
      for (int is = 0; is < Ns; ++is) {
        mult_wilson_ypb(pg1_yp, pg2_yp, &vL[NVCD * is].v[0], u,
                        &yp2[Nin4 * is], &ypn[Nin4 * is]);
      }
    } else {
      real_t *buf2_yp = (real_t *)chrecv_up[idir].ptr();
      int    ibf = VLENX * NVC * ND2 * Ns * (ix + Nxv * izt);
      real_t *u2 = &up[VLEN * NDF * site + NvU * idir];
      for (int is = 0; is < Ns; ++is) {
        mult_wilson_yp2(pg1_yp, pg2_yp,
                        &vL[NVCD * is].v[0], u2,
                        &yp2[Nin4 * is], &buf2_yp[ibf + VLENX * NVC * ND2 * is]);
      }
    }
    save_vec(vp2, vL, NVCD*Ns);

  }

#pragma omp barrier

}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::mult_ym(real_t *vp, real_t *wp,
                                               real_t *up){

  constexpr int idir = 1;

  real_t *yp = m_w1.ptr(0);  // working vector
  //  real_t *up = m_U.ptr(0);

  int ith = ThreadManager::get_thread_id();

  constexpr int Nin4 = VLEN * NVCD;
  const int Ns   = m_Ns;
  const int Nin5 = Nin4 * Ns;

  int Nxv = m_Nsize[0];
  int Nyv = m_Nsize[1];
  int Nz = m_Nsize[2];
  int Nt = m_Nsize[3];
  int Nxy  = Nxv * Nyv;
  int Nstv = Nxv * Nyv * Nz * Nt;
  int Nst = Nstv * VLEN;
  int NvU   = NDF * Nst;

  svbool_t pg1_yp, pg2_yp, pg1_ym, pg2_ym;
  set_predicate_yp(pg1_yp, pg2_yp);
  set_predicate_ym(pg1_ym, pg2_ym);


#pragma omp barrier
  BridgeQXS::accum_mult_domainwall_5din_5dir_dirac(vp, yp, wp,
                                                   m_mq, m_M0, m_Ns,
                                                   &m_boundary[0],
                                                   &m_b[0], &m_c[0],
                                                   m_Nsize, do_comm);
#pragma omp barrier

  if (do_comm[idir] == 1) {
    if (ith == 0) { chrecv_dn[idir].start(); }

    real_t *buf1_ym = (real_t *)chsend_up[idir].ptr();

    int Nxzt = Nxv * Nz * Nt;

    int ith, nth, site0, site1;
    set_threadtask(ith, nth, site0, site1, Nxzt);

    for (int ixzt = site0; ixzt < site1; ++ixzt) {
      int ix  = ixzt % Nxv;
      int izt = ixzt / Nxv;

      int    iy   = Nyv - 1;
      int    site = ix + Nxv * (iy + Nyv * izt);
      real_t *yp2 = &yp[Nin5 * site];

      int    ibf = VLENX * NVC * ND2 * Ns * (ix + Nxv * izt);
      real_t *u2 = &up[VLEN * NDF * site + NvU * idir];
      for (int is = 0; is < Ns; ++is) {
        mult_wilson_ym1(pg2_yp,
                        &buf1_ym[ibf + VLENX * NVC * ND2 * is], u2, &yp2[Nin4 * is]);
      }
    }

#pragma omp barrier
    if (ith == 0) {
      chsend_up[idir].start();
      chrecv_dn[idir].wait();
      chsend_up[idir].wait();
    }
#pragma omp barrier
  }

  int ith_dummy, nth, site0, site1;
  set_threadtask(ith_dummy, nth, site0, site1, Nstv);


  for (int site = site0; site < site1; ++site) {
    int ix   = site % Nxv;
    int iyzt = site / Nxv;
    int iy   = iyzt % Nyv;
    int izt  = site / Nxy;

    real_t *yp2 = &yp[Nin5 * site];
    real_t *vp2 = &vp[Nin5 * site];

    Vsimd_t vL[NVCD * Ns];
    load_vec(vL, vp2, NVCD * Ns);

    if ((iy > 0) || (do_comm[idir] == 0)) {
      int iy2 = (iy - 1 + Nyv) % Nyv;
      int nei = ix + Nxv * (iy2 + Nyv * izt);
      real_t *ypn = &yp[Nin5 * nei];
      real_t *ux  = &up[VLEN * NDF * site + NvU * idir];
      real_t *un  = &up[VLEN * NDF * nei + NvU * idir];
      for (int is = 0; is < Ns; ++is) {
        mult_wilson_ymb(pg1_ym, pg2_ym, &vL[NVCD * is].v[0], ux, un,
                        &yp2[Nin4 * is], &ypn[Nin4 * is]);
      }
    } else {
      real_t *buf2_ym = (real_t *)chrecv_dn[idir].ptr();
      int    ibf = VLENX * NVC * ND2 * Ns * (ix + Nxv * izt);
        real_t *u2 = &up[VLEN * NDF * site + NvU * idir];
        for (int is = 0; is < Ns; ++is) {
          mult_wilson_ym2(pg1_ym, pg2_ym,
                          &vL[NVCD * is].v[0], u2,
                          &yp2[Nin4 * is], &buf2_ym[ibf + VLENX * NVC * ND2 * is]);
        }
    }
    save_vec(vp2, vL, NVCD*Ns);

  }

#pragma omp barrier

}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::mult_zp(real_t *vp, real_t *wp,
                                               real_t *up){
  constexpr int idir = 2;

  real_t *yp = m_w1.ptr(0);  // working vector
  //  real_t *up = m_U.ptr(0);

  int ith = ThreadManager::get_thread_id();

  constexpr int Nin4 = VLEN * NVCD;
  constexpr int Nin4H = VLEN * NVC * ND2;
  const int Ns   = m_Ns;
  const int Nin5 = Nin4 * Ns;
  const int Nin5H = Nin4H * Ns;

  int Nxv = m_Nsize[0];
  int Nyv = m_Nsize[1];
  int Nz = m_Nsize[2];
  int Nt = m_Nsize[3];
  int Nxy  = Nxv * Nyv;
  int Nstv = Nxv * Nyv * Nz * Nt;
  int Nst = Nstv * VLEN;
  int NvU = NDF * Nst;


#pragma omp barrier
  BridgeQXS::accum_mult_domainwall_5din_5dir_dirac(vp, yp, wp,
                                                   m_mq, m_M0, m_Ns,
                                                   &m_boundary[0],
                                                   &m_b[0], &m_c[0],
                                                   m_Nsize, do_comm);
#pragma omp barrier

  if (do_comm[idir] == 1) {
    if (ith == 0) { chrecv_up[idir].start(); }

    real_t *buf1_zp = (real_t *)chsend_dn[idir].ptr();

    int Nxy  = Nxv * Nyv;
    int Nxyt = Nxv * Nyv * Nt;

    int ith_dummy, nth, site0, site1;
    set_threadtask(ith_dummy, nth, site0, site1, Nxyt);

    for (int ixyt = site0; ixyt < site1; ++ixyt) {
      int ixy  = ixyt % Nxy;
      int it   = ixyt / Nxy;
      int iz   = 0;
      int site = ixy + Nxy * (iz + Nz * it);
      real_t *yp2 = &yp[Nin5 * site];
      int    ibf  = Nin5H * (ixy + Nxy * it);
      for (int is = 0; is < Ns; ++is) {
        mult_wilson_zp1(&buf1_zp[ibf + Nin4H * is], &yp2[Nin4 * is]);
      }
    }

#pragma omp barrier
    if (ith == 0) {
      chsend_dn[idir].start();
      chrecv_up[idir].wait();
      chsend_dn[idir].wait();
    }
#pragma omp barrier
  }

  int ith_dummy, nth, site0, site1;
  set_threadtask(ith_dummy, nth, site0, site1, Nstv);

  for (int site = site0; site < site1; ++site) {
    int ixy  = site % Nxy;
    int izt  = site / Nxy;
    int iz   = izt % Nz;
    int it   = izt / Nz;

    real_t *yp2 = &yp[Nin5 * site];
    real_t *vp2 = &vp[Nin5 * site];

    Vsimd_t vL[NVCD * Ns];
    load_vec(vL, vp2, NVCD * Ns);

    if ((iz < Nz - 1) || (do_comm[idir] == 0)) {
      int iz2 = (iz + 1) % Nz;
      int nei = ixy + Nxy * (iz2 + Nz * it);
      real_t *ypn = &yp[Nin5 * nei];
      real_t *u   = &up[VLEN * NDF * site + NvU * idir];
      for (int is = 0; is < Ns; ++is) {
        mult_wilson_zpb(&vL[NVCD * is].v[0], u, &ypn[Nin4 * is]);
      }
    } else {
      real_t *buf2_zp = (real_t *)chrecv_up[idir].ptr();
      int    ibf = Nin5H * (ixy + Nxy * it);
      real_t *u2 = &up[VLEN * NDF * site + NvU * idir];
      for (int is = 0; is < Ns; ++is) {
        mult_wilson_zp2(&vL[NVCD * is].v[0], u2, &buf2_zp[ibf + Nin4H * is]);
        }
    }
    save_vec(vp2, vL, NVCD*Ns);

  }

#pragma omp barrier

}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::mult_zm(real_t *vp, real_t *wp,
                                               real_t *up){

  constexpr int idir = 2;

  real_t *yp = m_w1.ptr(0);  // working vector
  //  real_t *up = m_U.ptr(0);

  int ith = ThreadManager::get_thread_id();

  constexpr int Nin4 = VLEN * NVCD;
  constexpr int Nin4H = VLEN * NVC * ND2;
  const int Ns   = m_Ns;
  const int Nin5 = Nin4 * Ns;
  const int Nin5H = Nin4H * Ns;

  int Nxv = m_Nsize[0];
  int Nyv = m_Nsize[1];
  int Nz = m_Nsize[2];
  int Nt = m_Nsize[3];
  int Nxy  = Nxv * Nyv;
  int Nstv = Nxv * Nyv * Nz * Nt;
  int Nst = Nstv * VLEN;
  int NvU   = NDF * Nst;


#pragma omp barrier
  BridgeQXS::accum_mult_domainwall_5din_5dir_dirac(vp, yp, wp,
                                                   m_mq, m_M0, m_Ns,
                                                   &m_boundary[0],
                                                   &m_b[0], &m_c[0],
                                                   m_Nsize, do_comm);
#pragma omp barrier

  if (do_comm[idir] == 1) {
    if (ith == 0) { chrecv_dn[idir].start(); }

    real_t *buf1_zm = (real_t *)chsend_up[idir].ptr();
    int Nxy  = Nxv * Nyv;
    int Nxyt = Nxv * Nyv * Nt;

    int ith_dummy, nth, site0, site1;
    set_threadtask(ith_dummy, nth, site0, site1, Nxyt);

    for (int ixyt = site0; ixyt < site1; ++ixyt) {
      int ixy = ixyt % Nxy;
      int it  = ixyt / Nxy;
      int    iz   = Nz - 1;
      int    site = ixy + Nxy * (iz + Nz * it);
      real_t *yp2 = &yp[Nin5 * site];
      int    ibf  = Nin5H * (ixy + Nxy * it);
      real_t *u2  = &up[VLEN * NDF * site + NvU * idir];
      for (int is = 0; is < Ns; ++is) {
        mult_wilson_zm1(&buf1_zm[ibf + Nin4H * is], u2, &yp2[Nin4 * is]);
      }
    }

#pragma omp barrier
    if (ith == 0) {
      chsend_up[idir].start();
      chrecv_dn[idir].wait();
      chsend_up[idir].wait();
    }
#pragma omp barrier
  }

  int ith_dummy, nth, site0, site1;
  set_threadtask(ith_dummy, nth, site0, site1, Nstv);


  for (int site = site0; site < site1; ++site) {

    int ixy  = site % Nxy;
    int izt  = site / Nxy;
    int iz   = izt % Nz;
    int it   = izt / Nz;

    real_t *yp2 = &yp[Nin5 * site];
    real_t *vp2 = &vp[Nin5 * site];

    Vsimd_t vL[NVCD * Ns];
    load_vec(vL, vp2, NVCD * Ns);

    if ((iz > 0) || (do_comm[idir] == 0)) {
      int iz2 = (iz - 1 + Nz) % Nz;
      int nei = ixy + Nxy * (iz2 + Nz * it);
      real_t *ypn = &yp[Nin5 * nei];
      real_t *u   = &up[VLEN * NDF * nei + NvU * idir];
      for (int is = 0; is < Ns; ++is) {
        mult_wilson_zmb(&vL[NVCD * is].v[0], u, &ypn[Nin4 * is]);
      }
    } else {
      real_t *buf2_zm = (real_t *)chrecv_dn[idir].ptr();
      int ibf = Nin5H * (ixy + Nxy * it);
      for (int is = 0; is < Ns; ++is) {
        mult_wilson_zm2(&vL[NVCD * is].v[0], &buf2_zm[ibf + Nin4H * is]);
        }
    }
    save_vec(vp2, vL, NVCD*Ns);

  }

#pragma omp barrier

}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::mult_tp(real_t *vp, real_t *wp,
                                               real_t *up){
  constexpr int idir = 3;

  real_t *yp = m_w1.ptr(0);  // working vector
  //  real_t *up = m_U.ptr(0);

  int ith = ThreadManager::get_thread_id();

  constexpr int Nin4 = VLEN * NVCD;
  constexpr int Nin4H = VLEN * NVC * ND2;
  const int Ns   = m_Ns;
  const int Nin5 = Nin4 * Ns;
  const int Nin5H = Nin4H * Ns;

  int Nxv = m_Nsize[0];
  int Nyv = m_Nsize[1];
  int Nz = m_Nsize[2];
  int Nt = m_Nsize[3];
  int Nxyz = Nxv * Nyv * Nz;
  int Nstv = Nxv * Nyv * Nz * Nt;
  int Nst = Nstv * VLEN;
  int NvU = NDF * Nst;


#pragma omp barrier
  BridgeQXS::accum_mult_domainwall_5din_5dir_dirac(vp, yp, wp,
                                             m_mq, m_M0, m_Ns, &m_boundary[0],
                                             &m_b[0], &m_c[0],
                                             m_Nsize, do_comm);
#pragma omp barrier

  if (do_comm[idir] == 1) {
    if (ith == 0) { chrecv_up[idir].start(); }

    real_t *buf1_tp = (real_t *)chsend_dn[idir].ptr();

    int ith_dummy, nth, site0, site1;
    set_threadtask(ith_dummy, nth, site0, site1, Nxyz);

    for (int ixyz = site0; ixyz < site1; ++ixyz) {
      int    it   = 0;
      int    site = ixyz + Nxyz * it;
      real_t *yp2 = &yp[Nin5 * site];
      int    ibf  = Nin5H * ixyz;
      for (int is = 0; is < Ns; ++is) {
        mult_wilson_tp1_dirac(&buf1_tp[ibf + Nin4H * is], &yp2[Nin4 * is]);
      }
    }
#pragma omp barrier
    if (ith == 0) {
      chsend_dn[idir].start();
      chrecv_up[idir].wait();
      chsend_dn[idir].wait();
    }
#pragma omp barrier
  }

  int ith_dummy, nth, site0, site1;
  set_threadtask(ith_dummy, nth, site0, site1, Nstv);

  for (int site = site0; site < site1; ++site) {
    int it   = site / Nxyz;
    int ixyz = site % Nxyz;

    real_t *yp2 = &yp[Nin5 * site];
    real_t *vp2 = &vp[Nin5 * site];

    Vsimd_t vL[NVCD * Ns];
    load_vec(vL, vp2, NVCD * Ns);

    if ((it < Nt - 1) || (do_comm[idir] == 0)) {
      int it2 = (it + 1) % Nt;
      int nei = ixyz + Nxyz * it2;
      real_t *ypn = &yp[Nin5 * nei];
      real_t *u   = &up[VLEN * NDF * site + NvU * idir];
      for (int is = 0; is < Ns; ++is) {
        mult_wilson_tpb_dirac(&vL[NVCD * is].v[0], u, &ypn[Nin4 * is]);
      }
    } else {
      real_t *buf2_tp = (real_t *)chrecv_up[idir].ptr();
      int    ibf = Nin5H * ixyz;
      real_t *u2 = &up[VLEN * NDF * site + NvU * idir];
      for (int is = 0; is < Ns; ++is) {
        mult_wilson_tp2_dirac(&vL[NVCD * is].v[0], u2, &buf2_tp[ibf + Nin4H * is]);
      }
    }
    save_vec(vp2, vL, NVCD*Ns);

  }

#pragma omp barrier

}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::mult_tm(real_t *vp, real_t *wp,
                                               real_t *up){

  constexpr int idir = 3;

  real_t *yp = m_w1.ptr(0);  // working vector
  //  real_t *up = m_U.ptr(0);

  int ith = ThreadManager::get_thread_id();

  constexpr int Nin4 = VLEN * NVCD;
  constexpr int Nin4H = VLEN * NVC * ND2;
  const int Ns   = m_Ns;
  const int Nin5 = Nin4 * Ns;
  const int Nin5H = Nin4H * Ns;

  int Nxv = m_Nsize[0];
  int Nyv = m_Nsize[1];
  int Nz = m_Nsize[2];
  int Nt = m_Nsize[3];
  int Nxyz = Nxv * Nyv * Nz;
  int Nstv = Nxv * Nyv * Nz * Nt;
  int Nst = Nstv * VLEN;
  int NvU   = NDF * Nst;


#pragma omp barrier
  BridgeQXS::accum_mult_domainwall_5din_5dir_dirac(vp, yp, wp,
                                                   m_mq, m_M0, m_Ns,
                                                   &m_boundary[0],
                                                   &m_b[0], &m_c[0],
                                                   m_Nsize, do_comm);
#pragma omp barrier

  if (do_comm[idir] == 1) {
    if (ith == 0) { chrecv_dn[idir].start(); }

    real_t *buf1_tm = (real_t *)chsend_up[idir].ptr();

    int ith_dummy, nth, site0, site1;
    set_threadtask(ith_dummy, nth, site0, site1, Nxyz);

    for (int ixyz = site0; ixyz < site1; ++ixyz) {
      int    it   = Nt - 1;
      int    site = ixyz + Nxyz * it;
      real_t *yp2 = &yp[Nin5 * site];
      int    ibf  = Nin5H * ixyz;
      real_t *u2  = &up[VLEN * NDF * site + NvU * idir];
      for (int is = 0; is < Ns; ++is) {
        mult_wilson_tm1_dirac(&buf1_tm[ibf + Nin4H * is], u2, &yp2[Nin4 * is]);
      }
    }

#pragma omp barrier
    if (ith == 0) {
      chsend_up[idir].start();
      chrecv_dn[idir].wait();
      chsend_up[idir].wait();
    }
#pragma omp barrier
  }

  int ith_dummy, nth, site0, site1;
  set_threadtask(ith_dummy, nth, site0, site1, Nstv);


  for (int site = site0; site < site1; ++site) {

    int it   = site / Nxyz;
    int ixyz = site % Nxyz;

    real_t *yp2 = &yp[Nin5 * site];
    real_t *vp2 = &vp[Nin5 * site];

    Vsimd_t vL[NVCD * Ns];
    load_vec(vL, vp2, NVCD * Ns);

    if ((it > 0) || (do_comm[idir] == 0)) {
      int it2 = (it - 1 + Nt) % Nt;
      int nei = ixyz + Nxyz * it2;
      real_t *ypn = &yp[Nin5 * nei];
      real_t *u   = &up[VLEN * NDF * nei + NvU * idir];
      for (int is = 0; is < Ns; ++is) {
        mult_wilson_tmb_dirac(&vL[NVCD * is].v[0], u, &ypn[Nin4 * is]);
      }
    } else {
      real_t *buf2_tm = (real_t *)chrecv_dn[idir].ptr();
      int ibf = Nin5H * ixyz;
      for (int is = 0; is < Ns; ++is) {
        mult_wilson_tm2_dirac(&vL[NVCD * is].v[0], &buf2_tm[ibf + Nin4H * is]);
      }
    }
    save_vec(vp2, vL, NVCD*Ns);

  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_dd<AFIELD>::L_inv(AFIELD& v,
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
void AFopr_Domainwall_5din_dd<AFIELD>::U_inv(AFIELD& v,
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
void AFopr_Domainwall_5din_dd<AFIELD>::Ldag_inv(AFIELD& v,
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
void AFopr_Domainwall_5din_dd<AFIELD>::Udag_inv(AFIELD& v,
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
double AFopr_Domainwall_5din_dd<AFIELD>::flop_count_sap(std::string mode)
{
  int    Lvol  = CommonParameters::Lvol();
  double vsite = static_cast<double>(Lvol);
  double vNs   = static_cast<double>(m_Ns);
  int    Nc    = CommonParameters::Nc();
  int    Nd    = CommonParameters::Nd();

  return 0;
}
//============================================================END=====
