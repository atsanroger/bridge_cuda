/*!
      @file    afopr_Clover_QWS_dd-tmpl.h
      @brief
      @author  Issaku Kanamori (kanamori)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2023-03-20 10:52:44 +0900 (Mon, 20 Mar 2023) $
      @version $LastChangedRevision: 2499 $
*/

#include "lib_alt_QXS/Fopr/afopr_Clover_QWS_dd.h"

#ifdef USE_QWSLIB
#include "lib_alt_QXS/extra/qws_lib.h"
#endif

template<typename AFIELD>
const std::string AFopr_Clover_QWS_dd<AFIELD>::class_name = "AFopr_Clover_QWS_dd";


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::init()
{
  ThreadManager::assert_single_thread(class_name);

  m_repr = "Dirac";  // now only the Dirac repr is available.

  //int req_comm = 1;  // set 1 if communication forced any time
  int req_comm = 0;  // set 0 if communication called in necessary

  m_Nc  = CommonParameters::Nc();
  m_Nd  = CommonParameters::Nd();
  m_Nvc = 2 * m_Nc;
  m_Ndf = 2 * m_Nc * m_Nc;

  m_Nx   = CommonParameters::Nx();
  m_Ny   = CommonParameters::Ny();
  m_Nz   = CommonParameters::Nz();
  m_Nt   = CommonParameters::Nt();
  m_Nst  = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();

  m_Nxv  = m_Nx / VLENX;
  m_Nyv  = m_Ny / VLENY;
  m_Nstv = m_Nst / VLEN;

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

  m_bdsize.resize(m_Ndim);
  int Nd2 = m_Nd / 2;
  m_bdsize[0] = m_Nvc * Nd2 * m_Ny * m_Nz * m_Nt;
  m_bdsize[1] = m_Nvc * Nd2 * m_Nx * m_Nz * m_Nt;
  m_bdsize[2] = m_Nvc * Nd2 * m_Nx * m_Ny * m_Nt;
  m_bdsize[3] = m_Nvc * Nd2 * m_Nx * m_Ny * m_Nz;

  setup_channels();

  m_fopr_csw        = new Fopr_CloverTerm(m_repr);
  m_fopr_csw_chiral = new Fopr_CloverTerm("Chiral");

  // gauge configuration.
  m_U.reset(NDF, m_Nst, m_Ndim);

  // gauge configuration with block condition.
  m_Ublock.reset(NDF, m_Nst, m_Ndim);

  // working vectors.
  int NinF = 2 * m_Nc * m_Nd;
  m_v1.reset(NinF, m_Nst, 1);
  m_v2.reset(NinF, m_Nst, 1);

  int Ndm2 = m_Nd * m_Nd / 2;
  m_T.reset(NDF, m_Nst, Ndm2);

  m_T_qws.reset(72, m_Nst, 1);  // 72 = 2 * (Nc*Nd/2)^2
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::setup_channels()
{
  chsend_up.resize(m_Ndim);
  chrecv_up.resize(m_Ndim);
  chsend_dn.resize(m_Ndim);
  chrecv_dn.resize(m_Ndim);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    size_t Nvsize = m_bdsize[mu] * sizeof(real_t);

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
void AFopr_Clover_QWS_dd<AFIELD>::tidyup()
{
  ThreadManager::assert_single_thread(class_name);
#ifdef USE_QWSLIB
  //  QWS_lib::tidyup_qws();
#endif
  delete m_fopr_csw;
  delete m_fopr_csw_chiral;
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  double           kappa;
  std::vector<int> bc;
  std::vector<int> block_size;

  int err = 0;
  err += params.fetch_double("hopping_parameter", kappa);
  err += params.fetch_int_vector("boundary_condition", bc);
  err += params.fetch_int_vector("block_size", block_size);
  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(real_t(kappa), bc, block_size);

  Parameters params_csw_chiral = params;
  params_csw_chiral.set_string("gamma_matrix_type", "Chiral");
  m_fopr_csw_chiral->set_parameters(params_csw_chiral);

#ifdef CHIRAL_ROTATION
  // nothing to do for the other rep.
#else
  m_fopr_csw->set_parameters(params);
#endif
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::set_parameters(
  const real_t CKs,
  const std::vector<int> bc,
  const std::vector<int> block_size)
{
  //  ThreadManager::assert_single_thread(class_name);
  int ith = ThreadManager::get_thread_id();
#pragma omp barrier

  if (ith == 0) {
    m_CKs = CKs;

    assert(bc.size() == m_Ndim);
    m_boundary.resize(m_Ndim);
    m_block_size.resize(m_Ndim);
    for (int mu = 0; mu < m_Ndim; ++mu) {
      m_boundary[mu]   = bc[mu];
      m_block_size[mu] = block_size[mu];
    }
  }
#pragma omp barrier

  //- print input parameters
  vout.general(m_vl, "Parameters of %s:\n", class_name.c_str());
  vout.general(m_vl, "  CKs  = %8.4f\n", m_CKs);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }
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
    if ((m_block_sizev[0] * VLENX != m_block_size[0]) ||
        (m_block_sizev[1] * VLENY != m_block_size[1])) {
      vout.crucial(m_vl, "%s: bad blocksize in XY: must be divided by VLENX=%d, VLENY=%d but are %d, %d\n", class_name.c_str(), VLENX, VLENY, m_block_size[0], m_block_size[1]);
      exit(EXIT_FAILURE);
    }

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

#ifdef USE_QWSLIB
  QWS_lib::init_qws(&m_boundary[0]);
#endif

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::set_config(Field *u)
{
  ThreadManager::assert_single_thread(class_name);

  m_conf = u;

  m_fopr_csw->set_config(u);
  m_fopr_csw_chiral->set_config(u);


#pragma omp parallel
  {
    AIndex_lex<real_t, AFIELD::IMPL> index;
    convert_gauge(index, m_U, *u);

#ifdef USE_QWSLIB
#pragma omp barrier
#pragma omp master
    {
      // the configuration with fermion boundary condition must
      // be set only for spatial component: temporal boundary
      // condition is incorporated in the library.
      int         num_reconst = 0;
      std::string str_reconst;

      if (sizeof(real_t) == 4) {
        num_reconst = SU3_RECONSTRUCT_S;
        str_reconst = "SU3_RECONSTRUCT_S";
      } else if (sizeof(real_t) == 8) {
        num_reconst = SU3_RECONSTRUCT_D;
        str_reconst = "SU3_RECONSTRUCT_D";
      }
      for (int mu = 0; mu < m_Ndim - 1; ++mu) {
        if (m_boundary[mu] != 1) {
          if (num_reconst == 18) {
            set_boundary_config(m_U, mu);
          } else {
            vout.crucial("%s: QWS is not available for\n", class_name.c_str());
            vout.crucial("  boundary[%d] = %d\n", mu, m_boundary[mu]);
            vout.crucial("  %s = %d\n", str_reconst.c_str(), num_reconst);
            exit(EXIT_FAILURE);
          }
        }
      }
    } // omp master
#pragma omp barrier

    double kappa = (double)m_CKs;
    if (sizeof(real_t) == 8) {
      qws_loadg_dd_bridge_((double *)m_U.ptr(0), &kappa);
    } else if (sizeof(real_t) == 4) {
      qws_loadgs_dd_bridge_((float *)m_U.ptr(0), &kappa);
    }
#pragma omp barrier
#pragma omp master
    {
      // now temporal boundary condition is set for Bridge++ code.
      int mu = m_Ndim - 1;
      if (m_boundary[mu] != 1) {
        set_boundary_config(m_U, mu);
      }
    }
#else  // USE_QWSLIB
#pragma omp barrier
    set_boundary();
#endif

#pragma omp barrier
    copy(m_Ublock, m_U);
#pragma omp barrier
    set_block_config(m_Ublock);
  } // omp parallel

#ifdef CHIRAL_ROTATION
  set_csw_chrot(); // omp parallelized inside set_csw().
#else
  set_csw();       // omp parallelized inside set_csw().
#endif
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::set_boundary_config(AFIELD& U,
                                                      const int mu)
{
  int ipe[4], Nsize[4], Lsize[4], j[4];
  AIndex_lex<real_t, AFIELD::IMPL> index_lex;

  for (int i = 0; i < m_Ndim; ++i) {
    ipe[i]   = Communicator::ipe(i);
    Nsize[i] = CommonParameters::Nsize(i);
    Lsize[i] = CommonParameters::Lsize(i);
  }

  for (int t = 0; t < m_Nt; ++t) {
    j[3] = ipe[3] * Nsize[3] + t;
    for (int z = 0; z < m_Nz; ++z) {
      j[2] = ipe[2] * Nsize[2] + z;
      for (int y = 0; y < m_Ny; ++y) {
        j[1] = ipe[1] * Nsize[1] + y;
        for (int x = 0; x < m_Nx; ++x) {
          j[0] = ipe[0] * Nsize[0] + x;
          int site = x + m_Nx * (y + m_Ny * (z + m_Nz * t));

          if (j[mu] == Lsize[mu] - 1) {
            double bc = (double)(m_boundary[mu]);
            for (int in = 0; in < m_Ndf; ++in) {
              int    i  = index_lex.idx_G(in, site, mu);
              double uv = bc * U.cmp(i);
              U.set(i, uv);
            }
          }
        }
      }
    }
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::set_boundary()
{
  AIndex_lex<real_t, AFIELD::IMPL> index;

  int mu   = 0;
  int ipex = Communicator::ipe(mu);
  int npex = Communicator::npe(mu);

  if ((m_boundary[mu] != 1) && (ipex == npex - 1)) {
    real_t bc   = real_t(m_boundary[mu]);
    int    Nyzt = m_Ny * m_Nz * m_Nt;

    int ith, nth, is, ns;
    set_threadtask_mult(ith, nth, is, ns, Nyzt);

    for (int iyzt = is; iyzt < ns; ++iyzt) {
      int site = m_Nx - 1 + m_Nx * iyzt;
      for (int in = 0; in < NDF; ++in) {
        int    i  = index.idx_G(in, site, mu);
        real_t uv = m_U.cmp(i);
        uv = uv * bc;
        m_U.set(i, uv);
      }
    }
  }

  mu = 1;
  int ipey = Communicator::ipe(mu);
  int npey = Communicator::npe(mu);

  if ((m_boundary[mu] != 1) && (ipey == npey - 1)) {
    real_t bc   = real_t(m_boundary[mu]);
    int    Nxzt = m_Nx * m_Nz * m_Nt;

    int ith, nth, is, ns;
    set_threadtask_mult(ith, nth, is, ns, Nxzt);

    for (int ixzt = is; ixzt < ns; ++ixzt) {
      int ix   = ixzt % m_Nx;
      int izt  = ixzt / m_Nx;
      int site = ix + m_Nx * (m_Ny - 1 + m_Ny * izt);
      for (int in = 0; in < NDF; ++in) {
        int    i  = index.idx_G(in, site, mu);
        real_t uv = m_U.cmp(i);
        uv = uv * bc;
        m_U.set(i, uv);
      }
    }
  }

  mu = 2;
  int ipez = Communicator::ipe(mu);
  int npez = Communicator::npe(mu);

  if ((m_boundary[mu] != 1) && (ipez == npez - 1)) {
    real_t bc   = real_t(m_boundary[mu]);
    int    Nxy  = m_Nx * m_Ny;
    int    Nxyt = Nxy * m_Nt;

    int ith, nth, is, ns;
    set_threadtask_mult(ith, nth, is, ns, Nxyt);

    for (int ixyt = is; ixyt < ns; ++ixyt) {
      int ixy  = ixyt % Nxy;
      int it   = ixyt / Nxy;
      int site = ixy + Nxy * (m_Nz - 1 + m_Nz * it);
      for (int in = 0; in < NDF; ++in) {
        int    i  = index.idx_G(in, site, mu);
        real_t uv = m_U.cmp(i);
        uv = uv * bc;
        m_U.set(i, uv);
      }
    }
  }

  mu = 3;
  int ipet = Communicator::ipe(mu);
  int npet = Communicator::npe(mu);

  if ((m_boundary[mu] != 1) && (ipet == npet - 1)) {
    real_t bc   = real_t(m_boundary[mu]);
    int    Nxyz = m_Nx * m_Ny * m_Nz;

    int ith, nth, is, ns;
    set_threadtask_mult(ith, nth, is, ns, Nxyz);

    for (int ixyz = is; ixyz < ns; ++ixyz) {
      int site = ixyz + Nxyz * (m_Nt - 1);
      for (int in = 0; in < NDF; ++in) {
        int    i  = index.idx_G(in, site, mu);
        real_t uv = m_U.cmp(i);
        uv = uv * bc;
        m_U.set(i, uv);
      }
    }
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::set_block_config(AFIELD& U)
{
  AIndex_lex<real_t, AFIELD::IMPL> index;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nst);

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
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::set_csw()
{
  // The Dirac representation is assumed.
  // Current implementation makes use of the corelib Fopr_CloverTerm
  // that requires change of the gamma-matrix convention from
  // Bridge++ to QWS (BQCD) by multiplying gamma_4 before and after
  // m_fopr_csw->mult().
  // For numerical efficiency, proper implementation is necessaty.
  //                                        [08 Mar 2021 H.Matsufuru]

  Field_F v(m_Nst, 1), w(m_Nst, 1);

  AIndex_lex<real_t, AFIELD::IMPL> index;

  m_fopr_csw->set_mode("D");

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nst);

#pragma omp parallel
  {
    for (int id = 0; id < m_Nd / 2; ++id) {
      for (int ic = 0; ic < m_Nc; ++ic) {
        w.set(0.0);
        for (int site = is; site < ns; ++site) {
          w.set_r(ic, id, site, 0, 1.0);
        }
#pragma omp barrier

        m_fopr_csw->mult(v, w);
#pragma omp barrier

        for (int site = is; site < ns; ++site) {
          for (int ic2 = 0; ic2 < m_Nc; ++ic2) {
            real_t vt_r  = v.cmp_r(ic2, 0, site, 0);
            real_t vt_i  = v.cmp_i(ic2, 0, site, 0);
            int    idx_r = index.idx_Gr(ic, ic2, site, id + 0);
            int    idx_i = index.idx_Gi(ic, ic2, site, id + 0);
            m_T.set(idx_r, vt_r);
            m_T.set(idx_i, vt_i);
          }

          for (int ic2 = 0; ic2 < m_Nc; ++ic2) {
            real_t vt_r  = v.cmp_r(ic2, 1, site, 0);
            real_t vt_i  = v.cmp_i(ic2, 1, site, 0);
            int    idx_r = index.idx_Gr(ic, ic2, site, id + 4);
            int    idx_i = index.idx_Gi(ic, ic2, site, id + 4);
            m_T.set(idx_r, vt_r);
            m_T.set(idx_i, vt_i);
          }

          for (int ic2 = 0; ic2 < m_Nc; ++ic2) {
            real_t vt_r  = -v.cmp_r(ic2, 2, site, 0);
            real_t vt_i  = -v.cmp_i(ic2, 2, site, 0);
            int    idx_r = index.idx_Gr(ic, ic2, site, id + 2);
            int    idx_i = index.idx_Gi(ic, ic2, site, id + 2);
            m_T.set(idx_r, vt_r);
            m_T.set(idx_i, vt_i);
          }

          for (int ic2 = 0; ic2 < m_Nc; ++ic2) {
            real_t vt_r  = -v.cmp_r(ic2, 3, site, 0);
            real_t vt_i  = -v.cmp_i(ic2, 3, site, 0);
            int    idx_r = index.idx_Gr(ic, ic2, site, id + 6);
            int    idx_i = index.idx_Gi(ic, ic2, site, id + 6);
            m_T.set(idx_r, vt_r);
            m_T.set(idx_i, vt_i);
          }
        } // site loop
      }
    }
#pragma omp barrier

    real_t kappaR = 1.0 / m_CKs;
    scal(m_T, kappaR);
  } // omp parallel

  set_csw_qws();
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::set_csw_chrot()
{
  // This method set the clover term matrix assuming the chiral
  // representation for gamma-matrix in Bridge++ convention.
  // The conversion from the Bridge++ to QWS convention and
  // from Dirac to chiral representations are assumed to be cared
  // in BridgeQXS functions.
  //                                        [17 Jul 2021 H.Matsufuru]

  Field_F v(m_Nst, 1), w(m_Nst, 1);

  AIndex_lex<real_t, AFIELD::IMPL> index;

  const int Nin = NDF * ND * 2;

  m_fopr_csw_chiral->set_mode("D");

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nst);

  //11 22 33 44 55 66   12 34 56  23 45   13 24 35 46  15 26  14 36 25  16
  constexpr int idx[72] = {
    0,  -1,  6,  7, 16, 17, 28, 29, 24, 25, 34, 35,
    -1, -1,  1, -1, 12, 13, 18, 19, 32, 33, 26, 27,
    -1, -1, -1, -5,  2, -1,  8,  9, 20, 21, 30, 31,
    -1, -1, -1, -1, -1, -1,  3, -1, 14, 15, 22, 23,
    -1, -1, -1, -1, -1, -1, -1, -1,  4, -1, 10, 11,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  5, -1,
  };
  vout.general("hoge: chiral rep for Clover term\n");
#pragma omp parallel
  {
    for (int id = 0; id < m_Nd / 2; ++id) {
      for (int ic = 0; ic < m_Nc; ++ic) {
        w.set(0.0);
#pragma omp barrier

        for (int site = is; site < ns; ++site) {
          w.set_r(ic, id, site, 0, 1.0);
          w.set_r(ic, id + 2, site, 0, 1.0);
        }
#pragma omp barrier

        m_fopr_csw->mult(v, w);
#pragma omp barrier

        for (int site = is; site < ns; ++site) {
          for (int id2 = 0; id2 < m_Nd; ++id2) {
            for (int ic2 = 0; ic2 < m_Nc; ++ic2) {
              real_t vt_r = 0.5 * v.cmp_r(ic2, id2, site, 0);
              real_t vt_i = 0.5 * v.cmp_i(ic2, id2, site, 0);
              // int in = ic2 + NC*(ic + NC*(id + 2*id2));
              // int idx_r = index.idx(2*in,   Nin, site, 0);
              // int idx_i = index.idx(2*in+1, Nin, site, 0);
              // m_T.set(idx_r, vt_r);
              // m_T.set(idx_i, vt_i);
              int i    = ic2 + m_Nc * (id2 % 2);
              int j    = ic + m_Nc * id;
              int ij   = m_Nc * 2 * i + j;
              int in_r = idx[2 * ij];
              int in_i = idx[2 * ij + 1];
              if (in_r >= 0) {
                in_r += 36 * (id2 / 2);
                int idx_r = index.idx(in_r, Nin, site, 0);
                m_T.set(idx_r, vt_r);
              }
              if (in_i >= 0) {
                in_i += 36 * (id2 / 2);
                int idx_i = index.idx(in_i, Nin, site, 0);
                m_T.set(idx_i, vt_i);
              }
            }
          }
        } // site loop
#pragma omp barrier
      }
    }
#pragma omp barrier

    real_t kappaR = 1.0 / m_CKs;
    scal(m_T, kappaR);
  } // omp parallel
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::set_csw_qws()
{
  //
  //   T = 1 + clov
  //   T(Dirac,QWS)  = gamma4 T(Dirac,Brdige) gamma4
  //   T(Chiral,QWS) = U T(Dirac,QWS) Udag
  //                 = U gamma4 T(Dirac,Bridge) gamma4 Udag
  // in qws,
  //   T(Dirac,QWS) psi(Dirac,QWS)
  //     = Udag U T(Dirac,QWS) Udag U psi(Dirac,QWS)
  //     = Udag T(Chiral,QWS) U psi(Dirac, QWS)
  //     = Vdag  T(Chiral,QWS)/2 V psi(Dirac,QWS)
  //   with
  //    V = sqrt(2) U
  //      = [ 1  0  1  0 ]
  //        [ 0  1  0  1 ]
  //        [ 1  0 -1  0 ]
  //        [ 0  1  0 -1 ]
  //   and T(chiral,QWS)/2 is stored in the memory
  //
  //  For Bridge++, the conversion matrix is
  //   U'
  //   = sqrt(2)^{-1} [ 1  0  1  0 ]
  //                  [ 0  1  0  1 ]
  //                  [-1  0  1  0 ]
  //                  [ 0 -1  0  1 ]
  //    N.B. since
  //      U gamma_mu(Dirac,Bridge) Udag  = - gamma_mu(Chiral,Bridge)
  //    and thus  U T(Dirac,Bridge) Udag = T(Chiral,Bridge),
  //    one can use U instead for the purpose here
  //
  //  Using U',
  //     T(Chiral,Bridge) = U' T(Dirac,Bridge) U'dag
  //  and
  //   T(Chiral,QWS) = U T(Dirac,QWS) Udag
  //                 = U gamma4 T(Dirac,Bridge) gamma4 Udag
  //                 = U gamma4 U'dag T(Chiral,Bridge) U' gamma4 Udag
  //  and the trasnformation matrix in the chiral rep. is
  //   U gamma4 U'dag =  [ 0  0 -1  0 ]
  //                     [ 0  0  0 -1 ]
  //                     [ 1  0  0  0 ]
  //                     [ 0  1  0  0 ]
  //   in 2x2 block notation,
  //    U gamma4 U'dag [ T1 0  ] U' gamma4 Udag
  //                   [ 0  T2 ]
  //    =  [ T2 0  ]
  //       [ 0  T1 ]
  //   so the upper and lower blocks are interchanged between Bridge and QWS clover term
  //
  //    N.B.
  //     U gamma4 Udag =  [ 0  0  1  0 ]
  //                      [ 0  0  0  1 ]
  //                      [ 1  0  0  0 ]
  //                      [ 0  1  0  0 ]
  //     which also gives
  //     U gamma4 Udag [ T1 0  ] Ugamma4 Udag
  //                   [ 0  T2 ]
  //        =  [ T2 0  ]
  //           [ 0  T1 ]
  //
  // [  20. Jul.2021 I.Kanamori]

#ifdef USE_QWSLIB
  Field_F v(m_Nst, 1), w(m_Nst, 1);
  Field   T(2 * 6 * 6, m_Nst, 1);
  Field   Tinv(2 * 6 * 6, m_Nst, 1);

  AIndex_lex<real_t, AFIELD::IMPL>        index;
  AIndex_lex_QWS_dd<real_t, AFIELD::IMPL> index_qws;

  m_fopr_csw_chiral->set_mode("D");

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nst);

#pragma omp parallel
  {
    int idx[36]
      = { 0,   6,  8, 10, 12, 14,
          -1,  1, 16, 18, 20, 22,
          -1, -1,  2, 24, 26, 28,
          -1, -1, -1,  3, 30, 32,
          -1, -1, -1, -1,  4, 34,
          -1, -1, -1, -1, -1, 5 };
    constexpr double factor = 0.5;

    T.set(0.0);
    for (int ud = 0; ud < 2; ++ud) {
      for (int id = 0; id < m_Nd / 2; ++id) {
        for (int ic = 0; ic < m_Nc; ++ic) {
          w.set(0.0);
          for (int site = is; site < ns; ++site) {
            w.set_r(ic, 2 * ud + id, site, 0, 1.0);
          }
#pragma omp barrier
          m_fopr_csw_chiral->mult(v, w);
#pragma omp barrier
          int j = 3 * id + ic;
          for (int site = is; site < ns; ++site) {
            int ud2    = 1 - ud;
            int offset = ud2 * 36;
            // diagonal: real
            double vt_r = v.cmp_r(ic, 2 * ud + id, site, 0);
            T.set(offset + j, site, 0, factor * (1.0 - vt_r));
            for (int i = 0; i < 6; ++i) {
              int ij = idx[6 * i + j];
              if (ij > 5) {
                int    id2  = i / 3 + 2 * ud;
                int    ic2  = i % 3;
                double vt_r = v.cmp_r(ic2, id2, site, 0);
                double vt_i = v.cmp_i(ic2, id2, site, 0);
                T.set(offset + ij, site, 0, -factor * vt_r);
                T.set(offset + ij + 1, site, 0, -factor * vt_i);
              }
            } // i
          }   // site
        }     // ic
      }       // id
    }         // ud
#pragma omp barrier
    solve_csw_qws_inv(Tinv, T);
    ::convert(index, m_T_qws, Tinv);
#pragma omp barrier
    if (sizeof(real_t) == 8) {
      qws_loadc_dd_bridge_((double *)m_T_qws.ptr(0));
    } else if (sizeof(real_t) == 4) {
      qws_loadcs_dd_bridge_((float *)m_T_qws.ptr(0));
    }
    m_T_qws.set(0.0);
    ::convert(index_qws, m_T_qws, T);
  } // omp parallel
#endif // USE_QWSLIB
}


//====================================================================
template<class AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::solve_csw_qws_inv(Field& Tinv, const Field& T)
{
  // inversion with LU factorization + overall rescaling
  //   in Driac rep., in order to absorb 1/sqrt{2} in the uniary trf,
  //   the rescaling factor must be set to 1/4.  [default]
  //   in Chiral rep., the factor must be set to 1

#ifdef USE_QWSLIB
  constexpr double factor = 0.25;

  typedef dcomplex complex_t;

  constexpr int NS = 4; // 4 sprinor
  constexpr int N  = 6; // matrix size to invert  [ 6x6 in each block]


  // size of Feild
  assert(NCOL == CommonParameters::Nc());
  assert(NS == CommonParameters::Nd());
  assert(N == NS * NCOL / 2);
  assert(T.check_size(2 * 36, T.nvol(), 1));
  assert(Tinv.check_size(2 * 36, T.nvol(), 1));
  const size_t vol = T.nvol();

  const double *pT    = T.ptr(0);
  double       *pTinv = Tinv.ptr(0);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, T.nvol());

  complex_t LU[N * N];

  // working area to store the inversion
  complex_t workT[N * N];


  for (int s = is; s < ns; ++s) {
    size_t offset_clov = 2 * s * N * N;

    const double *clov     = pT + offset_clov;
    double       *clov_inv = pTinv + offset_clov;

    for (int block = 0; block < 2; block++) {
      // fill U part
      // diagonal
      for (int i = 0; i < N; i++) {
        LU[N * i + i] = clov[i];
      }
      // offdiag
      int ii = N;
      for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
          LU[N * i + j] = complex_t(clov[ii], clov[ii + 1]);
          ii           += 2;
        }
      }

      // LU decomposition
      for (int i = 0; i < N; i++) {
        double diag_re  = LU[N * i + i].real(); // (re im) w/ im \simeq 0 --> re
        double diag_inv = 1.0 / diag_re;
        LU[N * i + i] = diag_inv;               // preparaton for inv.
        for (int k = i + 1; k < N; k++) {
          complex_t c = LU[N * i + k] * diag_inv;
          LU[N * k + i] = c;                // L  w/o c.c.
          for (int j = i + 1; j < N; j++) { // U
            LU[N * k + j] -= conj(c) * LU[N * i + j];
          }
        }
      } // i

      for (int i = 0; i < N; i++) { // for i-th vector
        complex_t *x = &workT[N * i];
        for (int j = 0; j < N; j++) {
          x[j] = 0.0;
        }
        x[i] = factor;

        // solving L
        for (int j = 1; j < N; j++) {
          complex_t x_j = x[j];
          for (int k = i; k < j; k++) {
            x_j -= conj(LU[N * j + k]) * x[k];
          }
          x[j] = x_j;
        }
        // solving U
        for (int j = N - 1; j >= i; j--) {
          complex_t x_j = x[j];
          for (int k = j + 1; k < N; k++) {
            x_j -= LU[N * j + k] * x[k];
          }
          x_j *= LU[N * j + j].real();
          x[j] = x_j;
        }
      }

      // pack to the result
      ii = 0;
      for (int i = 0; i < N; i++) {
        clov_inv[ii] = workT[N * i + i].real();
        ii++;
      }
      for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
          complex_t val = conj(workT[i * N + j]);
          clov_inv[ii] = val.real();
          ii++;
          clov_inv[ii] = val.imag();
          ii++;
        }
      }
      clov     += N * N;
      clov_inv += N * N;
    } // block
  }   //s
#endif // USE_QWSLIB
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::convert(AFIELD& v, const AFIELD& w)
{
#ifdef USE_QWSLIB
  if (sizeof(real_t) == 4) {
    qws_loadfs_dd_bridge_((scs_t *)v.ptr(0), (const float *)w.ptr(0));
  } else if (sizeof(real_t) == 8) {
    qws_loadfd_dd_bridge_((scd_t *)v.ptr(0), (const double *)w.ptr(0));
  } else {
    vout.crucial(" %s: cannot happen, unknown sizeof(real_t)=%ul in convert()\n", class_name.c_str(), sizeof(real_t));
    exit(EXIT_FAILURE);
  }
#else
  copy(v, w);
#endif
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::convert(AFIELD& v, const Field& w)
{
  AIndex_lex<real_t, AFIELD::IMPL> index_lex;
  convert_spinor(index_lex, v, w);
  mult_gm4(m_v2, v);
  convert(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::reverse(AFIELD& v, const AFIELD& w)
{
#ifdef USE_QWSLIB
  if (sizeof(real_t) == 4) {
    qws_storefs_dd_bridge_((float *)v.ptr(0), (const scs_t *)w.ptr(0));
  } else if (sizeof(real_t) == 8) {
    qws_storefd_dd_bridge_((double *)v.ptr(0), (const scd_t *)w.ptr(0));
  } else {
    vout.crucial(" %s: cannot happen, unknown sizeof(real_t)=%ul in reverse()\n", class_name.c_str(), sizeof(real_t));
    exit(EXIT_FAILURE);
  }
#else
  copy(v, w);
#endif // USE_QWSLIB
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::reverse(Field& v, const AFIELD& w)
{
  AIndex_lex<real_t, AFIELD::IMPL> index_lex;

  reverse(m_v1, w);
  mult_gm4(m_v2, m_v1);
  reverse_spinor(index_lex, v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::mult_up(int mu, AFIELD& v, const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

  if (mu == 0) {
    mult_xp(vp, wp, 0);
  } else if (mu == 1) {
    mult_yp(vp, wp, 0);
  } else if (mu == 2) {
    mult_zp(vp, wp, 0);
  } else if (mu == 3) {
    mult_tp(vp, wp, 0);
  } else {
    vout.crucial(m_vl, "%s: mult_up for %d direction is undefined.",
                 class_name.c_str(), mu);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::mult_dn(int mu, AFIELD& v, const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

  if (mu == 0) {
    mult_xm(vp, wp, 0);
  } else if (mu == 1) {
    mult_ym(vp, wp, 0);
  } else if (mu == 2) {
    mult_zm(vp, wp, 0);
  } else if (mu == 3) {
    mult_tm(vp, wp, 0);
  } else {
    vout.crucial(m_vl, "%s: mult_dn for %d direction is undefined.",
                 class_name.c_str(), mu);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
std::string AFopr_Clover_QWS_dd<AFIELD>::get_mode() const
{
  return m_mode;
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::set_mode(std::string mode)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_mode = mode;

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::mult(AFIELD& v, const AFIELD& w)
{
  if (m_mode == "D") {
    return D(v, w);
  } else if (m_mode == "DdagD") {
    return DdagD(v, w);
  } else if (m_mode == "Ddag") {
    return Ddag(v, w);
  } else if (m_mode == "H") {
    return H(v, w);
  } else {
    vout.crucial(m_vl, "%s: mode undefined.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::mult_dag(AFIELD& v, const AFIELD& w)
{
  if (m_mode == "D") {
    return Ddag(v, w);
  } else if (m_mode == "DdagD") {
    return DdagD(v, w);
  } else if (m_mode == "Ddag") {
    return D(v, w);
  } else if (m_mode == "H") {
    return H(v, w);
  } else {
    vout.crucial(m_vl, "%s: mode undefined.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::D(AFIELD& v, const AFIELD& w)
{
#ifdef USE_QWSLIB
  mult_D_qws(v, w);
#else
  mult_D(v, w);
  //mult_D_alt(v, w);
#endif
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::DdagD(AFIELD& v, const AFIELD& w)
{
  D(m_v2, w);
  mult_gm5(v, m_v2);
  D(m_v2, v);
  mult_gm5(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::Ddag(AFIELD& v, const AFIELD& w)
{
  mult_gm5(v, w);
  D(m_v2, v);
  mult_gm5(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::project_chiral(AFIELD& v, const AFIELD& w,
                                                 const int ch)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

  mult_gm5(vp, wp);
  if (ch == 1) {
    aypx(real_t(1.0), vp, wp);
  } else {
    aypx(real_t(-1.0), vp, wp);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::mult_gm5(AFIELD& v, const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

#pragma omp barrier

  mult_gm5(vp, wp);
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::mult_gm5(real_t *v, real_t *w)
{ // Dirac representation.
#ifdef USE_QXS_ACLE
  BridgeQXS::mult_wilson_gm5_dirac(v, w, m_Nsize);
#else
  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);
  Vsimd_t wt[2];

  for (int site = is; site < ns; ++site) {
    for (int ic = 0; ic < NC; ++ic) {
      for (int id = 0; id < ND; ++id) {
        int idx1 = 2 * (id + ND * ic) + NVCD * site;
        load_vec(wt, &w[VLEN * idx1], 2);
        scal_vec(wt, (real_t)(-1.0), 2);
        int id2  = (id + 2) % ND;
        int idx2 = 2 * (id2 + ND * ic) + NVCD * site;
        save_vec(&v[VLEN * idx2], wt, 2);
      }
    }
  }
#endif

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::scal_local(real_t *vp, real_t a)
{
  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  for (int site = is; site < ns; ++site) {
    Vsimd_t vt[NVCD];
    load_vec(vt, &vp[VLEN * NVCD * site], NVCD);
    scal_vec(vt, a, NVCD);
    save_vec(&vp[VLEN * NVCD * site], vt, NVCD);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::mult_gm4(AFIELD& v, const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  Vsimd_t wt[2];

  for (int site = is; site < ns; ++site) {
    for (int ic = 0; ic < NC; ++ic) {
      for (int id = 0; id < ND2; ++id) {
        int idx1 = 2 * (id + ND * ic) + NVCD * site;
        load_vec(wt, &wp[VLEN * idx1], 2);
        save_vec(&vp[VLEN * idx1], wt, 2);
      }

      for (int id = ND2; id < ND; ++id) {
        int idx1 = 2 * (id + ND * ic) + NVCD * site;
        load_vec(wt, &wp[VLEN * idx1], 2);
        scal_vec(wt, real_t(-1.0), 2);
        save_vec(&vp[VLEN * idx1], wt, 2);
      }
    }
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::mult_csw(real_t *v2, real_t *v1)
{                               // Dirac representation is assumed.
  real_t *u = m_T.ptr(0);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nstv);

#pragma omp barrier

  for (int site = is; site < ns; ++site) {
    Vsimd_t v2v[NVCD];
    load_vec(v2v, &v2[VLEN * NVCD * site], NVCD);

    Vsimd_t v1v[NVCD];
    load_vec(v1v, &v1[VLEN * NVCD * site], NVCD);

    Vsimd_t ut[NDF], wt1[2], wt2[2];

    for (int jd = 0; jd < ND2; ++jd) {
      for (int id = 0; id < ND; ++id) {
        int ig = VLEN * NDF * (site + m_Nstv * (id + ND * jd));
        load_vec(ut, &u[ig], NDF);
        for (int ic = 0; ic < NC; ++ic) {
          int ic2 = 2 * ic;
          int id2 = (id + ND2) % ND;
          mult_ctv(wt1, &ut[ic2], &v1v[2 * id], NC);
          mult_ctv(wt2, &ut[ic2], &v1v[2 * id2], NC);
          int icd1 = 2 * (jd + ND * ic);
          int icd2 = 2 * (jd + ND2 + ND * ic);
          axpy_vec(&v2v[icd1], real_t(1.0), wt1, 2);
          axpy_vec(&v2v[icd2], real_t(1.0), wt2, 2);
        }
      }
    }

    save_vec(&v2[VLEN * NVCD * site], v2v, NVCD);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::mult_clv_inv(AFIELD& v)
{                               // Dirac representation is assumed.
#ifdef USE_QWSLIB
  // (1+csw)^{-1}
  if (sizeof(real_t) == 4) {
    int deo = 0;
    clv_s_dd_(&deo, (scs_t *)v.ptr(0));
    deo = 1;
    clv_s_dd_(&deo, (scs_t *)v.ptr(0));
  } else if (sizeof(real_t) == 8) {
    int deo = 0;
    clv_d_dd_(&deo, (scd_t *)v.ptr(0));
    deo = 1;
    clv_d_dd_(&deo, (scd_t *)v.ptr(0));
  }
#endif
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::mult_clv(AFIELD& v)
{                               // Dirac representation is assumed.
#ifdef USE_QWSLIB
  // (1+csw)
  if (sizeof(real_t) == 4) {
    int deo = 0;
    clv_matvec_s_dd_(&deo, (scs_t *)v.ptr(0), (clvs_t *)m_T_qws.ptr(0), (scs_t *)v.ptr(0));
    deo = 1;
    clv_matvec_s_dd_(&deo, (scs_t *)v.ptr(0), (clvs_t *)m_T_qws.ptr(0), (scs_t *)v.ptr(0));
  } else if (sizeof(real_t) == 8) {
    int deo = 0;
    clv_matvec_d_dd_(&deo, (scd_t *)v.ptr(0), (clvd_t *)m_T_qws.ptr(0), (scd_t *)v.ptr(0));
    deo = 1;
    clv_matvec_d_dd_(&deo, (scd_t *)v.ptr(0), (clvd_t *)m_T_qws.ptr(0), (scd_t *)v.ptr(0));
  }
#endif
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::mult_D(AFIELD& v, const AFIELD& w)
{
  real_t *v2 = v.ptr(0);
  real_t *v1 = const_cast<AFIELD *>(&w)->ptr(0);
  real_t *up = m_U.ptr(0);
  real_t *ct = m_T.ptr(0);

#pragma omp barrier

  if (do_comm_any > 0) {
    real_t *buf1_xp = (real_t *)chsend_dn[0].ptr();
    real_t *buf1_xm = (real_t *)chsend_up[0].ptr();
    real_t *buf1_yp = (real_t *)chsend_dn[1].ptr();
    real_t *buf1_ym = (real_t *)chsend_up[1].ptr();
    real_t *buf1_zp = (real_t *)chsend_dn[2].ptr();
    real_t *buf1_zm = (real_t *)chsend_up[2].ptr();
    real_t *buf1_tp = (real_t *)chsend_dn[3].ptr();
    real_t *buf1_tm = (real_t *)chsend_up[3].ptr();

    BridgeQXS::mult_wilson_1_dirac(buf1_xp, buf1_xm, buf1_yp, buf1_ym,
                                   buf1_zp, buf1_zm, buf1_tp, buf1_tm,
                                   up, v1, &m_boundary[0], m_Nsize, do_comm);

#pragma omp barrier

#pragma omp master
    {
      chset_recv.start();
      chset_send.start();
    }
  }

  BridgeQXS::mult_clover_bulk_dirac(v2, up, ct, v1,
                                    m_CKs, &m_boundary[0], m_Nsize, do_comm);

  if (do_comm_any > 0) {
#pragma omp master
    {
      chset_recv.wait();
      chset_send.wait();
    }

#pragma omp barrier

    real_t *buf2_xp = (real_t *)chrecv_up[0].ptr();
    real_t *buf2_xm = (real_t *)chrecv_dn[0].ptr();
    real_t *buf2_yp = (real_t *)chrecv_up[1].ptr();
    real_t *buf2_ym = (real_t *)chrecv_dn[1].ptr();
    real_t *buf2_zp = (real_t *)chrecv_up[2].ptr();
    real_t *buf2_zm = (real_t *)chrecv_dn[2].ptr();
    real_t *buf2_tp = (real_t *)chrecv_up[3].ptr();
    real_t *buf2_tm = (real_t *)chrecv_dn[3].ptr();

    BridgeQXS::mult_wilson_2_dirac(v2, up, v1,
                                   buf2_xp, buf2_xm, buf2_yp, buf2_ym,
                                   buf2_zp, buf2_zm, buf2_tp, buf2_tm,
                                   m_CKs, &m_boundary[0], m_Nsize, do_comm);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::mult_sap(AFIELD& v, const AFIELD& w,
                                           const int ieo)
{
  real_t *v2 = v.ptr(0);
  real_t *v1 = const_cast<AFIELD *>(&w)->ptr(0);
  real_t *up = m_Ublock.ptr(0);
  real_t *ct = m_T.ptr(0);

  int jeo = (ieo + m_Ieo) % 2;

#pragma omp barrier

  BridgeQXS::mult_clover_dd_dirac(v2, up, ct, v1,
                                  m_CKs, &m_boundary[0],
                                  m_Nsize, m_block_sizev, jeo);

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::mult_dd(AFIELD& v, const AFIELD& w)
{
  //#ifdef USE_QWSLIB
#if 0
#else
  real_t *v2 = v.ptr(0);
  real_t *v1 = const_cast<AFIELD *>(&w)->ptr(0);
  real_t *up = m_Ublock.ptr(0);
  real_t *ct = m_T.ptr(0);

#pragma omp barrier

  BridgeQXS::mult_clover_dd_dirac(v2, up, ct, v1,
                                  m_CKs, &m_boundary[0],
                                  m_Nsize, m_block_sizev, -1);

  // vector varialbes at the boundary of the domain are not treated here...
  //  BridgeQXS::mult_clover_bulk_dirac(v2, up, ct, v1,
  //                         m_CKs, &m_boundary[0], m_Nsize, do_comm);

#pragma omp barrier
#endif
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::mult_D_qws(AFIELD& v, const AFIELD& w)
{
#ifdef USE_QWSLIB
  if (sizeof(real_t) == 4) {
    //    vout.general("hoge: calling ddd_s_noprl_\n");
    ddd_s_noprl_((scs_t *)v.ptr(0), (scs_t *)const_cast<AFIELD *>(&w)->ptr(0));
  } else if (sizeof(real_t) == 8) {
    //    vout.general("hoge: calling ddd_d_\n");
    int nth = ThreadManager::get_num_threads();
    if (nth > 1) {
      vout.crucial(m_vl, "%s: mult_D_qws cannot be called in thread parallel region\n", class_name.c_str());
      abort();
    }
    ddd_d_((scd_t *)v.ptr(0), (scd_t *)const_cast<AFIELD *>(&w)->ptr(0));
  }
#pragma omp barrier
  mult_clv(v);
#pragma omp barrier
  vout.detailed("mult_D_qws, done.\n");
#ifdef DEBUG
  {
    double w2 = w.norm2();
    double v2 = v.norm2();
    vout.general("  mult_D_qws: (in)  w2=%e\n", w2);
    vout.general("  mult_D_qws: (out) v2=%e\n", v2);
  }
#endif
#endif
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::mult_D_alt(AFIELD& v, const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

  clear(vp);

  mult_xp(vp, wp, 0);
  mult_xm(vp, wp, 0);
  mult_yp(vp, wp, 0);
  mult_ym(vp, wp, 0);
  mult_zp(vp, wp, 0);
  mult_zm(vp, wp, 0);
  mult_tp(vp, wp, 0);
  mult_tm(vp, wp, 0);

  mult_csw(vp, wp);

  aypx(-m_CKs, vp, wp);

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::mult_dup(AFIELD& v, const AFIELD& w,
                                           const int mu)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

#pragma omp barrier

  clear(vp);

  if (mu == 0) {
    mult_xp(vp, wp, 1);
    scal_local(vp, real_t(-1.0));
    mult_xp(vp, wp, 0);
  } else if (mu == 1) {
    mult_yp(vp, wp, 1);
    scal_local(vp, real_t(-1.0));
    mult_yp(vp, wp, 0);
  } else if (mu == 2) {
    mult_zp(vp, wp, 1);
    scal_local(vp, real_t(-1.0));
    mult_zp(vp, wp, 0);
  } else if (mu == 3) {
    mult_tp(vp, wp, 1);
    scal_local(vp, real_t(-1.0));
    mult_tp(vp, wp, 0);
  }

  scal_local(vp, -m_CKs);

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::mult_ddn(AFIELD& v, const AFIELD& w,
                                           const int mu)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

#pragma omp barrier

  clear(vp);

  if (mu == 0) {
    mult_xm(vp, wp, 1);
    scal_local(vp, real_t(-1.0));
    mult_xm(vp, wp, 0);
  } else if (mu == 1) {
    mult_ym(vp, wp, 1);
    scal_local(vp, real_t(-1.0));
    mult_ym(vp, wp, 0);
  } else if (mu == 2) {
    mult_zm(vp, wp, 1);
    scal_local(vp, real_t(-1.0));
    mult_zm(vp, wp, 0);
  } else if (mu == 3) {
    mult_tm(vp, wp, 1);
    scal_local(vp, real_t(-1.0));
    mult_tm(vp, wp, 0);
  }

  scal_local(vp, -m_CKs);

#pragma omp barrier
}


//====================================================================

/*
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::mult_block_hop(AFIELD &v, const AFIELD &w,
                                             int mu)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);

#pragma omp barrier

  clear(vp);

  if(mu == 0){
    mult_xp(vp, wp, 1);
    scal_local(vp, real_t(-1.0));
    mult_xp(vp, wp, 0);
  }else if(mu == 1){
    mult_yp(vp, wp, 1);
    scal_local(vp, real_t(-1.0));
    mult_yp(vp, wp, 0);
  }else if(mu == 2){
    mult_zp(vp, wp, 1);
    scal_local(vp, real_t(-1.0));
    mult_zp(vp, wp, 0);
  }else if(mu == 3){
    mult_tp(vp, wp, 1);
    scal_local(vp, real_t(-1.0));
    mult_tp(vp, wp, 0);
  }

  scal_local(vp, -m_CKs);

#pragma omp barrier

}
*/

//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::H(AFIELD& v, const AFIELD& w)
{
  D(m_v2, w);
  mult_gm5(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::aypx(real_t a, real_t *v, real_t *w)
{
  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  Vsimd_t vt[NVCD], wt[NVCD];

  for (int site = is; site < ns; ++site) {
    load_vec(vt, &v[VLEN * NVCD * site], NVCD);
    load_vec(wt, &w[VLEN * NVCD * site], NVCD);
    aypx_vec(a, vt, wt, NVCD);
    save_vec(&v[VLEN * NVCD * site], vt, NVCD);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::clear(real_t *v)
{
  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  Vsimd_t vt[NVCD];
  clear_vec(vt, NVCD);

  for (int site = is; site < ns; ++site) {
    save_vec(&v[VLEN * NVCD * site], vt, NVCD);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::mult_xp(real_t *v2, real_t *v1, int isap)
{
  int idir = 0;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  Vsimd_t v2v[NVCD];

  real_t *buf1 = (real_t *)chsend_dn[0].ptr();
  real_t *buf2 = (real_t *)chrecv_up[0].ptr();

  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);
  if (isap != 0) u = m_Ublock.ptr(m_Ndf * m_Nst * idir);

#pragma omp barrier

  if (do_comm[0] > 0) {
    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nxv;
      int iyzt = site / m_Nxv;
      if (ix == 0) {
        int ibf = VLENY * NVC * ND2 * iyzt;
        mult_wilson_xp1(&buf1[ibf], &v1[VLEN * NVCD * site]);
      }
    }

#pragma omp barrier

#pragma omp master
    {
      chrecv_up[0].start();
      chsend_dn[0].start();
      chrecv_up[0].wait();
      chsend_dn[0].wait();
    }
#pragma omp barrier
  } // if(do_comm[0] == 1)

  for (int site = is; site < ns; ++site) {
    int ix   = site % m_Nxv;
    int iyzt = site / m_Nxv;

    Vsimd_t v2v[NVCD];
    clear_vec(v2v, NVCD);

    real_t zL[VLEN * NVCD];

    if ((ix < m_Nxv - 1) || (do_comm[0] == 0)) {
      int nei = ix + 1 + m_Nxv * iyzt;
      if (ix == m_Nxv - 1) nei = 0 + m_Nxv * iyzt;
      shift_vec2_xbw(zL, &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei], NVCD);
      mult_wilson_xpb(v2v, &u[VLEN * NDF * site], zL);
    } else {
      int ibf = VLENY * NVC * ND2 * iyzt;
      shift_vec0_xbw(zL, &v1[VLEN * NVCD * site], NVCD);
      mult_wilson_xpb(v2v, &u[VLEN * NDF * site], zL);
      mult_wilson_xp2(v2v, &u[VLEN * NDF * site], &buf2[ibf]);
    }

    add_vec(&v2[VLEN * NVCD * site], v2v, NVCD);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::mult_xm(real_t *v2, real_t *v1, int isap)
{
  int idir = 0;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  Vsimd_t v2v[NVCD];

  real_t *buf1 = (real_t *)chsend_up[0].ptr();
  real_t *buf2 = (real_t *)chrecv_dn[0].ptr();

  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);
  if (isap != 0) u = m_Ublock.ptr(m_Ndf * m_Nst * idir);

#pragma omp barrier

  if (do_comm[0] > 0) {
    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nxv;
      int iyzt = site / m_Nxv;
      if (ix == m_Nxv - 1) {
        int ibf = VLENY * NVC * ND2 * iyzt;
        mult_wilson_xm1(&buf1[ibf], &u[VLEN * NDF * site],
                        &v1[VLEN * NVCD * site]);
      }
    }

#pragma omp barrier
#pragma omp master
    {
      chrecv_dn[0].start();
      chsend_up[0].start();
      chrecv_dn[0].wait();
      chsend_up[0].wait();
    }
#pragma omp barrier
  } // end of if(do_comm[0] > 0)

  for (int site = is; site < ns; ++site) {
    int ix   = site % m_Nxv;
    int iyzt = site / m_Nxv;

    real_t zL[VLEN * NVCD];
    real_t uL[VLEN * NDF];

    clear_vec(v2v, NVCD);

    if ((ix > 0) || (do_comm[0] == 0)) {
      int nei = ix - 1 + m_Nxv * iyzt;
      if (ix == 0) nei = m_Nxv - 1 + m_Nxv * iyzt;
      shift_vec2_xfw(zL, &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei], NVCD);
      shift_vec2_xfw(uL, &u[VLEN * NDF * site], &u[VLEN * NDF * nei], NDF);
      mult_wilson_xmb(v2v, uL, zL);
    } else {
      int ibf = VLENY * NVC * ND2 * iyzt;
      shift_vec0_xfw(zL, &v1[VLEN * NVCD * site], NVCD);
      shift_vec0_xfw(uL, &u[VLEN * NDF * site], NDF);
      mult_wilson_xmb(v2v, uL, zL);
      mult_wilson_xm2(v2v, &buf2[ibf]);
    }

    add_vec(&v2[VLEN * NVCD * site], v2v, NVCD);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::mult_yp(real_t *v2, real_t *v1, int isap)
{
  int idir = 1;
  int Nxy  = m_Nxv * m_Nyv;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  real_t *buf1 = (real_t *)chsend_dn[1].ptr();
  real_t *buf2 = (real_t *)chrecv_up[1].ptr();

  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);
  if (isap != 0) u = m_Ublock.ptr(m_Ndf * m_Nst * idir);

#pragma omp barrier

  if (do_comm[1] > 0) {
    for (int site = is; site < ns; ++site) {
      int ix  = site % m_Nxv;
      int iy  = (site / m_Nxv) % m_Nyv;
      int izt = site / Nxy;
      if (iy == 0) {
        int ibf = VLENX * NVC * ND2 * (ix + m_Nxv * izt);
        mult_wilson_yp1(&buf1[ibf], &v1[VLEN * NVCD * site]);
      }
    }

#pragma omp barrier

#pragma omp master
    {
      chrecv_up[1].start();
      chsend_dn[1].start();
      chrecv_up[1].wait();
      chsend_dn[1].wait();
    }

#pragma omp barrier
  }  // end of if(do_comm[1] > 0)

  for (int site = is; site < ns; ++site) {
    int ix  = site % m_Nxv;
    int iy  = (site / m_Nxv) % m_Nyv;
    int izt = site / Nxy;

    Vsimd_t v2v[NVCD];
    clear_vec(v2v, NVCD);

    real_t zL[VLEN * NVCD];

    if ((iy < m_Nyv - 1) || (do_comm[1] == 0)) {
      int iy2 = (iy + 1) % m_Nyv;
      int nei = ix + m_Nxv * (iy2 + m_Nyv * izt);
      shift_vec2_ybw(zL, &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei], NVCD);
      mult_wilson_ypb(v2v, &u[VLEN * NDF * site], zL);
    } else {
      int ibf = VLENX * NVC * ND2 * (ix + m_Nxv * izt);
      shift_vec0_ybw(zL, &v1[VLEN * NVCD * site], NVCD);
      mult_wilson_ypb(v2v, &u[VLEN * NDF * site], zL);
      mult_wilson_yp2(v2v, &u[VLEN * NDF * site], &buf2[ibf]);
    }

    add_vec(&v2[VLEN * NVCD * site], v2v, NVCD);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::mult_ym(real_t *v2, real_t *v1, int isap)
{
  int idir = 1;
  int Nxy  = m_Nxv * m_Nyv;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  real_t *buf1 = (real_t *)chsend_up[1].ptr();
  real_t *buf2 = (real_t *)chrecv_dn[1].ptr();

  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);
  if (isap != 0) u = m_Ublock.ptr(m_Ndf * m_Nst * idir);

#pragma omp barrier

  if (do_comm[1] > 0) {
    for (int site = is; site < ns; ++site) {
      int ix  = site % m_Nxv;
      int iy  = (site / m_Nxv) % m_Nyv;
      int izt = site / Nxy;
      if (iy == m_Nyv - 1) {
        int ibf = VLENX * NVC * ND2 * (ix + m_Nxv * izt);
        mult_wilson_ym1(&buf1[ibf], &u[VLEN * NDF * site],
                        &v1[VLEN * NVCD * site]);
      }
    }

#pragma omp barrier

#pragma omp master
    {
      chrecv_dn[1].start();
      chsend_up[1].start();
      chrecv_dn[1].wait();
      chsend_up[1].wait();
    }

#pragma omp barrier
  }

  for (int site = is; site < ns; ++site) {
    int ix  = site % m_Nxv;
    int iy  = (site / m_Nxv) % m_Nyv;
    int izt = site / Nxy;

    Vsimd_t v2v[NVCD];
    clear_vec(v2v, NVCD);

    real_t zL[VLEN * NVCD];
    real_t uL[VLEN * NDF];

    if ((iy != 0) || (do_comm[idir] == 0)) {
      int iy2 = (iy - 1 + m_Nyv) % m_Nyv;
      int nei = ix + m_Nxv * (iy2 + m_Nyv * izt);
      shift_vec2_yfw(zL, &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei], NVCD);
      shift_vec2_yfw(uL, &u[VLEN * NDF * site], &u[VLEN * NDF * nei], NDF);
      mult_wilson_ymb(v2v, uL, zL);
    } else {
      int ibf = VLENX * NVC * ND2 * (ix + m_Nxv * izt);
      shift_vec0_yfw(zL, &v1[VLEN * NVCD * site], NVCD);
      shift_vec0_yfw(uL, &u[VLEN * NDF * site], NDF);
      mult_wilson_ymb(v2v, uL, zL);
      mult_wilson_ym2(v2v, &buf2[ibf]);
    }

    add_vec(&v2[VLEN * NVCD * site], v2v, NVCD);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::mult_zp(real_t *v2, real_t *v1, int isap)
{
  int idir = 2;
  int Nxy  = m_Nxv * m_Nyv;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  real_t *buf1 = (real_t *)chsend_dn[2].ptr();
  real_t *buf2 = (real_t *)chrecv_up[2].ptr();

  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);
  if (isap != 0) u = m_Ublock.ptr(m_Ndf * m_Nst * idir);

#pragma omp barrier

  if (do_comm[2] > 0) {
    for (int site = is; site < ns; ++site) {
      int ixy = site % Nxy;
      int iz  = (site / Nxy) % m_Nz;
      int it  = site / (Nxy * m_Nz);
      if (iz == 0) {
        int ibf = VLEN * NVC * ND2 * (ixy + Nxy * it);
        mult_wilson_zp1(&buf1[ibf], &v1[VLEN * NVCD * site]);
      }
    }

#pragma omp barrier

#pragma omp master
    {
      chrecv_up[2].start();
      chsend_dn[2].start();
      chrecv_up[2].wait();
      chsend_dn[2].wait();
    }

#pragma omp barrier
  }

  for (int site = is; site < ns; ++site) {
    int ixy = site % Nxy;
    int iz  = (site / Nxy) % m_Nz;
    int it  = site / (Nxy * m_Nz);

    Vsimd_t v2v[NVCD];
    clear_vec(v2v, NVCD);

    if ((iz != m_Nz - 1) || (do_comm[2] == 0)) {
      int iz2 = (iz + 1) % m_Nz;
      int nei = ixy + Nxy * (iz2 + m_Nz * it);
      mult_wilson_zpb(v2v, &u[VLEN * NDF * site], &v1[VLEN * NVCD * nei]);
    } else {
      int ibf = VLEN * NVC * ND2 * (ixy + Nxy * it);
      mult_wilson_zp2(v2v, &u[VLEN * NDF * site], &buf2[ibf]);
    }

    add_vec(&v2[VLEN * NVCD * site], v2v, NVCD);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::mult_zm(real_t *v2, real_t *v1, int isap)
{
  int idir = 2;
  int Nxy  = m_Nxv * m_Nyv;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  real_t *buf1 = (real_t *)chsend_up[2].ptr();
  real_t *buf2 = (real_t *)chrecv_dn[2].ptr();

  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);
  if (isap != 0) u = m_Ublock.ptr(m_Ndf * m_Nst * idir);

#pragma omp barrier

  if (do_comm[2] > 0) {
    for (int site = is; site < ns; ++site) {
      int ixy = site % Nxy;
      int iz  = (site / Nxy) % m_Nz;
      int it  = site / (Nxy * m_Nz);
      if (iz == m_Nz - 1) {
        int ibf = VLEN * NVC * ND2 * (ixy + Nxy * it);
        mult_wilson_zm1(&buf1[ibf], &u[VLEN * NDF * site],
                        &v1[VLEN * NVCD * site]);
      }
    }

#pragma omp barrier

#pragma omp master
    {
      chrecv_dn[2].start();
      chsend_up[2].start();
      chrecv_dn[2].wait();
      chsend_up[2].wait();
    }

#pragma omp barrier
  }

  for (int site = is; site < ns; ++site) {
    int ixy = site % Nxy;
    int iz  = (site / Nxy) % m_Nz;
    int it  = site / (Nxy * m_Nz);

    Vsimd_t v2v[NVCD];
    clear_vec(v2v, NVCD);

    if ((iz > 0) || (do_comm[2] == 0)) {
      int iz2 = (iz - 1 + m_Nz) % m_Nz;
      int nei = ixy + Nxy * (iz2 + m_Nz * it);
      mult_wilson_zmb(v2v, &u[VLEN * NDF * nei], &v1[VLEN * NVCD * nei]);
    } else {
      int ibf = VLEN * NVC * ND2 * (ixy + Nxy * it);
      mult_wilson_zm2(v2v, &buf2[ibf]);
    }

    add_vec(&v2[VLEN * NVCD * site], v2v, NVCD);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::mult_tp(real_t *v2, real_t *v1, int isap)
{
  int idir = 3;
  int Nxyz = m_Nxv * m_Nyv * m_Nz;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  real_t *buf1 = (real_t *)chsend_dn[3].ptr();
  real_t *buf2 = (real_t *)chrecv_up[3].ptr();

  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);
  if (isap != 0) u = m_Ublock.ptr(m_Ndf * m_Nst * idir);

#pragma omp barrier

  if (do_comm[3] > 0) {
    for (int site = is; site < ns; ++site) {
      int ixyz = site % Nxyz;
      int it   = site / Nxyz;
      if (it == 0) {
        mult_wilson_tp1_dirac(&buf1[VLEN * NVC * ND2 * ixyz],
                              &v1[VLEN * NVCD * site]);
      }
    }

#pragma omp barrier

#pragma omp master
    {
      chrecv_up[3].start();
      chsend_dn[3].start();
      chrecv_up[3].wait();
      chsend_dn[3].wait();
    }

#pragma omp barrier
  }

  for (int site = is; site < ns; ++site) {
    int ixyz = site % Nxyz;
    int it   = site / Nxyz;

    Vsimd_t v2v[NVCD];
    clear_vec(v2v, NVCD);

    if ((it < m_Nt - 1) || (do_comm[3] == 0)) {
      int it2 = (it + 1) % m_Nt;
      int nei = ixyz + Nxyz * it2;
      mult_wilson_tpb_dirac(v2v, &u[VLEN * NDF * site],
                            &v1[VLEN * NVCD * nei]);
    } else {
      mult_wilson_tp2_dirac(v2v, &u[VLEN * NDF * site],
                            &buf2[VLEN * NVC * ND2 * ixyz]);
    }

    add_vec(&v2[VLEN * NVCD * site], v2v, NVCD);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_QWS_dd<AFIELD>::mult_tm(real_t *v2, real_t *v1, int isap)
{
  int idir = 3;
  int Nxyz = m_Nxv * m_Nyv * m_Nz;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  real_t *buf1 = (real_t *)chsend_up[3].ptr();
  real_t *buf2 = (real_t *)chrecv_dn[3].ptr();

  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);
  if (isap != 0) u = m_Ublock.ptr(m_Ndf * m_Nst * idir);

#pragma omp barrier

  if (do_comm[3] > 0) {
    for (int site = is; site < ns; ++site) {
      int ixyz = site % Nxyz;
      int it   = site / Nxyz;
      if (it == m_Nt - 1) {
        mult_wilson_tm1_dirac(&buf1[VLEN * NVC * ND2 * ixyz],
                              &u[VLEN * NDF * site], &v1[VLEN * NVCD * site]);
      }
    }

#pragma omp barrier

#pragma omp master
    {
      chrecv_dn[3].start();
      chsend_up[3].start();
      chrecv_dn[3].wait();
      chsend_up[3].wait();
    }
#pragma omp barrier
  }

  for (int site = is; site < ns; ++site) {
    int ixyz = site % Nxyz;
    int it   = site / Nxyz;

    Vsimd_t v2v[NVCD];
    clear_vec(v2v, NVCD);

    if ((it > 0) || (do_comm[3] == 0)) {
      int it2 = (it - 1 + m_Nt) % m_Nt;
      int nei = ixyz + Nxyz * it2;
      mult_wilson_tmb_dirac(v2v, &u[VLEN * NDF * nei],
                            &v1[VLEN * NVCD * nei]);
    } else {
      mult_wilson_tm2_dirac(v2v, &buf2[VLEN * NVC * ND2 * ixyz]);
    }

    add_vec(&v2[VLEN * NVCD * site], v2v, NVCD);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
double AFopr_Clover_QWS_dd<AFIELD>::flop_count(const std::string mode)
{
  // The following counting explicitly depends on the implementation.
  // It will be recalculated when the code is modified.
  // The present counting is based on rev.1107. [24 Aug 2014 H.Matsufuru]

  int    Lvol = CommonParameters::Lvol();
  double flop_site, flop;

  if (m_repr == "Dirac") {
    flop_site = static_cast<double>(
      m_Nc * m_Nd * (4 + 6 * (4 * m_Nc + 2) + 2 * (4 * m_Nc + 1)));
  } else if (m_repr == "Chiral") {
    flop_site = static_cast<double>(
      m_Nc * m_Nd * (4 + 8 * (4 * m_Nc + 2)));
  } else {
    vout.crucial(m_vl, "%s: input repr is undefined.\n");
    abort();
  }

  flop = flop_site * static_cast<double>(Lvol);
  if ((mode == "DdagD") || (mode == "DDdag")) flop *= 2.0;

  return flop;
}


//====================================================================
template<typename AFIELD>
double AFopr_Clover_QWS_dd<AFIELD>::flop_count_sap()
{
  // The following implentastion was copied from SIMD implementation
  // by Kanamori-san. To be confirmed. [08 May 2021 H.Matsufuru]

  int NPE = CommonParameters::NPE();
  int Nc  = CommonParameters::Nc();
  int Nd  = CommonParameters::Nd();

  int    NBx   = m_Nx / m_block_size[0];
  int    NBy   = m_Ny / m_block_size[1];
  int    NBz   = m_Nz / m_block_size[2];
  int    NBt   = m_Nt / m_block_size[3];
  size_t nvol2 = NBx * NBy * NBz * NBt / 2;

  int block_x = m_block_size[0];
  int block_y = m_block_size[1];
  int block_z = m_block_size[2];
  int block_t = m_block_size[3];

  int clover_site = 4 * Nc * Nc * Nd * Nd;
  int hop_x_site  = Nc * Nd * (4 * Nc + 2);
  int hop_y_site  = Nc * Nd * (4 * Nc + 2);
  int hop_z_site  = Nc * Nd * (4 * Nc + 2);
  int hop_t_site  = Nc * Nd * (4 * Nc + 1);
  int accum_site  = 4 * Nc * Nd;

  // Dirac representation is assumed.
  double flop_x = 2.0 * static_cast<double>(hop_x_site
                                            * (block_x - 1) * block_y * block_z * block_t);
  double flop_y = 2.0 * static_cast<double>(hop_y_site
                                            * block_x * (block_y - 1) * block_z * block_t);
  double flop_z = 2.0 * static_cast<double>(hop_z_site
                                            * block_x * block_y * (block_z - 1) * block_t);
  double flop_t = 2.0 * static_cast<double>(hop_t_site
                                            * block_x * block_y * block_z * (block_t - 1));
  double flop_sap_mult = flop_x + flop_y + flop_z + flop_t
                         + static_cast<double>(clover_site + accum_site);

  return flop_sap_mult * static_cast<double>(nvol2)
         * static_cast<double>(NPE);
}


//============================================================END=====
