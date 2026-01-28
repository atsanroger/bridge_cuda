/*!
      @file    afopr_CloverTerm-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2550 $
*/

#include "lib/ResourceManager/threadManager.h"

// QXS index ordering is assumed:
#define ID1R  0
#define ID1I  1
#define ID2R  2
#define ID2I  3
#define ID3R  4
#define ID3I  5
#define ID4R  6
#define ID4I  7

template<typename AFIELD>
const std::string AFopr_CloverTerm<AFIELD>::class_name = "AFopr_CloverTerm";

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::init(const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  // switch of coomunication
  //int req_comm = 1;  // set 1 if communication forced any time
  int req_comm = 0;  // set 0 if communication called in necessary

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

  m_Nc = CommonParameters::Nc();
  if (m_Nc != 3) {
    vout.crucial("%s: only applicable to Nc = 3\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

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

  m_bdsize.resize(m_Ndim);
  int Nd2 = m_Nd / 2;
  m_bdsize[0] = m_Nvc * Nd2 * m_Ny * m_Nz * m_Nt;
  m_bdsize[1] = m_Nvc * Nd2 * m_Nx * m_Nz * m_Nt;
  m_bdsize[2] = m_Nvc * Nd2 * m_Nx * m_Ny * m_Nt;
  m_bdsize[3] = m_Nvc * Nd2 * m_Nx * m_Ny * m_Nz;

  setup_channels();

  vout.increase_indent();

  Parameters params_ct = params;

#ifdef CHIRAL_ROTATION
  params_ct.set_string("gamma_matrix_type", "Chiral");
  m_fopr_csw = new Fopr_CloverTerm(params_ct);
#else
  m_fopr_csw = new Fopr_CloverTerm(params_ct);
#endif

  set_parameters(params);

  vout.decrease_indent();

  // gauge configuration.
  m_U.reset(NDF, m_Nst, m_Ndim);

  // working vectors.
  int NinF = 2 * m_Nc * m_Nd;
  m_v2.reset(NinF, m_Nst, 1);

  int Ndm2 = m_Nd * m_Nd / 2;
  //  m_T.reset(NDF, m_Nst, Ndm2);
  m_T.reset(NDF * Ndm2, m_Nst, 1);

  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());
}


//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::setup_channels()
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
void AFopr_CloverTerm<AFIELD>::tidyup()
{
  ThreadManager::assert_single_thread(class_name);

  delete m_fopr_csw;
}


//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::set_parameters(const Parameters& params)
{
  double           kappa, csw;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_double("hopping_parameter", kappa);
  err += params.fetch_double("clover_coefficient", csw);
  err += params.fetch_int_vector("boundary_condition", bc);
  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(real_t(kappa), real_t(csw), bc);

  Parameters params_csw = params;
#ifdef CHIRAL_ROTATION
  params_csw.set_string("gamma_matrix_type", "Chiral");
#endif
  m_fopr_csw->set_parameters(params_csw);
}


//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::set_parameters(const real_t CKs,
                                          const real_t csw,
                                          const std::vector<int> bc)
{
  assert(bc.size() == m_Ndim);

#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  if (ith == 0) {
    m_CKs = CKs;
    m_csw = csw;
    m_boundary.resize(m_Ndim);
    for (int mu = 0; mu < m_Ndim; ++mu) {
      m_boundary[mu] = bc[mu];
    }
  }

  vout.general(m_vl, "%s: set parameters\n", class_name.c_str());
  vout.general(m_vl, "  gamma-matrix type = %s\n", m_repr.c_str());
  vout.general(m_vl, "  kappa = %8.4f\n", m_CKs);
  vout.general(m_vl, "  cSW   = %12.8f\n", m_csw);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::get_parameters(Parameters& params) const
{
  params.set_double("hopping_parameter", double(m_CKs));
  params.set_double("clover_coefficient", double(m_csw));
  params.set_int_vector("boundary_condition", m_boundary);
  params.set_string("gamma_matrix_type", m_repr);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::set_config(Field *u)
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
void AFopr_CloverTerm<AFIELD>::set_config_omp(Field *u)
{
  vout.detailed(m_vl, "  set_config_omp is called.\n");

#pragma omp parallel
  {
    set_config_impl(u);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::set_config_impl(Field *u)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  if (ith == 0) m_conf = u;

  AIndex_lex<real_t, AFIELD::IMPL> index;

  convert_gauge(index, m_U, *u);

  QXS_Gauge::set_boundary(m_U, m_boundary);

  m_fopr_csw->set_config(u);

#ifdef CHIRAL_ROTATION
  set_csw_chrot();
#else
  set_csw();
#endif
}


//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::set_csw()
{
  // The Dirac representation is assumed.
  // Current implementation makes use of the corelib Fopr_CloverTerm
  // that requires change of the gamma-matrix convention from
  // Bridge++ to QXS (BQCD) by multiplying gamma_4 before and after
  // m_fopr_csw->mult().
  // For numerical efficiency, proper implementation is necessaty.
  //                                        [08 Mar 2021 H.Matsufuru]

#pragma omp barrier

  AIndex_lex<real_t, AFIELD::IMPL> index;

  const int Nin = NDF * ND * 2;

  m_fopr_csw->set_mode("D");

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nst);

  for (int id = 0; id < m_Nd / 2; ++id) {
    for (int ic = 0; ic < m_Nc; ++ic) {
      m_w1.set(0.0);
#pragma omp barrier

      for (int site = is; site < ns; ++site) {
        m_w1.set_r(ic, id, site, 0, 1.0);
      }
#pragma omp barrier

      m_fopr_csw->mult(m_w2, m_w1);
#pragma omp barrier

      for (int site = is; site < ns; ++site) {
        for (int ic2 = 0; ic2 < m_Nc; ++ic2) {
          real_t vt_r = m_w2.cmp_r(ic2, 0, site, 0);
          real_t vt_i = m_w2.cmp_i(ic2, 0, site, 0);
          int in = ic2 + NC * (ic + NC * (id + 0));
          int idx_r = index.idx(2 * in, Nin, site, 0);
          int idx_i = index.idx(2 * in + 1, Nin, site, 0);
          m_T.set(idx_r, vt_r);
          m_T.set(idx_i, vt_i);
        }

        for (int ic2 = 0; ic2 < m_Nc; ++ic2) {
          real_t vt_r = m_w2.cmp_r(ic2, 1, site, 0);
          real_t vt_i = m_w2.cmp_i(ic2, 1, site, 0);
          int in = ic2 + NC * (ic + NC * (id + 4));
          int idx_r = index.idx(2 * in, Nin, site, 0);
          int idx_i = index.idx(2 * in + 1, Nin, site, 0);
          m_T.set(idx_r, vt_r);
          m_T.set(idx_i, vt_i);
        }

        for (int ic2 = 0; ic2 < m_Nc; ++ic2) {
          real_t vt_r = -m_w2.cmp_r(ic2, 2, site, 0);
          real_t vt_i = -m_w2.cmp_i(ic2, 2, site, 0);
          int in = ic2 + NC * (ic + NC * (id + 2));
          int idx_r = index.idx(2 * in, Nin, site, 0);
          int idx_i = index.idx(2 * in + 1, Nin, site, 0);
          m_T.set(idx_r, vt_r);
          m_T.set(idx_i, vt_i);
        }

        for (int ic2 = 0; ic2 < m_Nc; ++ic2) {
          real_t vt_r = -m_w2.cmp_r(ic2, 3, site, 0);
          real_t vt_i = -m_w2.cmp_i(ic2, 3, site, 0);
          int in = ic2 + NC * (ic + NC * (id + 6));
          int idx_r = index.idx(2 * in, Nin, site, 0);
          int idx_i = index.idx(2 * in + 1, Nin, site, 0);
          m_T.set(idx_r, vt_r);
          m_T.set(idx_i, vt_i);
        }
      } // site loop
#pragma omp barrier
    }
  }
#pragma omp barrier

  real_t kappaR = 1.0 / m_CKs;
  scal(m_T, kappaR);

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::set_csw_chrot()
{
  // This method set the clover term matrix assuming the chiral
  // representation for gamma-matrix in Bridge++ convention.
  // The conversion from the Bridge++ to QWS convention and
  // from Dirac to chiral representations are assumed to be cared
  // in BridgeQXS functions.
  //                                        [17 Jul 2021 H.Matsufuru]

  // The following implementation is by I.Kanamori

  AIndex_lex<real_t, AFIELD::IMPL> index;

  const int Nin = NDF * ND * 2;

  m_fopr_csw->set_mode("D");

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

  for (int id = 0; id < m_Nd / 2; ++id) {
    for (int ic = 0; ic < m_Nc; ++ic) {
      m_w1.set(0.0);
#pragma omp barrier

      for (int site = is; site < ns; ++site) {
        m_w1.set_r(ic, id, site, 0, 1.0);
        m_w1.set_r(ic, id + 2, site, 0, 1.0);
      }
#pragma omp barrier

      m_fopr_csw->mult(m_w2, m_w1);
#pragma omp barrier

      for (int site = is; site < ns; ++site) {
        for (int id2 = 0; id2 < m_Nd; ++id2) {
          for (int ic2 = 0; ic2 < m_Nc; ++ic2) {
            real_t vt_r = 0.5 * m_w2.cmp_r(ic2, id2, site, 0);
            real_t vt_i = 0.5 * m_w2.cmp_i(ic2, id2, site, 0);
            int    i    = ic2 + m_Nc * (id2 % 2);
            int    j    = ic + m_Nc * id;
            int    ij   = m_Nc * 2 * i + j;
            int    in_r = idx[2 * ij];
            int    in_i = idx[2 * ij + 1];
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

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::convert(AFIELD& v, const Field& w)
{
  AIndex_lex<real_t, AFIELD::IMPL> index_lex;
  convert_spinor(index_lex, m_v2, w);

  mult_gm4(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::reverse(Field& v, const AFIELD& w)
{
  mult_gm4(m_v2, w);

  AIndex_lex<real_t, AFIELD::IMPL> index_lex;
  reverse_spinor(index_lex, v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::mult_up(int mu, AFIELD& v, const AFIELD& w)
{
  vout.crucial(m_vl, "%s: mult_up is undefined.", class_name.c_str());
  exit(EXIT_FAILURE);
}


//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::mult_dn(int mu, AFIELD& v, const AFIELD& w)
{
  vout.crucial(m_vl, "%s: mult_dn is undefined.", class_name.c_str());
  exit(EXIT_FAILURE);
}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::set_mode(std::string mode)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_mode = mode;

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
std::string AFopr_CloverTerm<AFIELD>::get_mode() const
{
  return m_mode;
}


//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::mult(AFIELD& v, const AFIELD& w)
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
void AFopr_CloverTerm<AFIELD>::mult_dag(AFIELD& v, const AFIELD& w)
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
void AFopr_CloverTerm<AFIELD>::D(AFIELD& v, const AFIELD& w)
{
  mult_D(v, w);
  //mult_D_alt(v, w);
}


//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::DdagD(AFIELD& v, const AFIELD& w)
{
  D(m_v2, w);
  mult_gm5(v, m_v2);
  D(m_v2, v);
  mult_gm5(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::Ddag(AFIELD& v, const AFIELD& w)
{
  mult_gm5(v, w);
  D(m_v2, v);
  mult_gm5(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::mult_gm5(AFIELD& v, const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

#pragma omp barrier

  BridgeQXS::mult_wilson_gm5_dirac(vp, wp, m_Nsize);

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::mult_gm4(AFIELD& v, const AFIELD& w)
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
void AFopr_CloverTerm<AFIELD>::mult_csw(real_t *v2, real_t *v1)
{                               // Dirac representation is assumed.
  real_t *u = m_T.ptr(0);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nstv);

#pragma omp barrier

  for (int site = is; site < ns; ++site) {
    Vsimd_t v2v[NVCD];
    // clear_vec(v2v, NVCD);
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
void AFopr_CloverTerm<AFIELD>::mult_D(AFIELD& v, const AFIELD& w)
{
  real_t *v2 = v.ptr(0);
  real_t *v1 = const_cast<AFIELD *>(&w)->ptr(0);

  mult_csw(v2, v1);

}


//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::H(AFIELD& v, const AFIELD& w)
{
  D(m_v2, w);
  mult_gm5(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::mult_isigma(AFIELD& v, const AFIELD& w,
                                           const int mu, const int nu)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

  if(mu == 0){
    if(nu == 1) mult_gm12(vp, wp,  1.0);
    if(nu == 2) mult_gm31(vp, wp, -1.0);
    if(nu == 3) mult_gm14_dirac(vp, wp,  1.0);
  }else if(mu == 1){
    if(nu == 0) mult_gm12(vp, wp, -1.0);
    if(nu == 2) mult_gm23(vp, wp,  1.0);
    if(nu == 3) mult_gm24_dirac(vp, wp,  1.0);
  }else if(mu == 2){
    if(nu == 0) mult_gm31(vp, wp,  1.0);
    if(nu == 1) mult_gm23(vp, wp, -1.0);
    if(nu == 3) mult_gm34_dirac(vp, wp,  1.0);
  }else if(mu == 3){
    if(nu == 0) mult_gm14_dirac(vp, wp, -1.0);
    if(nu == 1) mult_gm24_dirac(vp, wp, -1.0);
    if(nu == 2) mult_gm34_dirac(vp, wp, -1.0);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::mult_gm12(real_t* vp, const real_t* wp,
					 const real_t jd)
{
  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  svbool_t pg = set_predicate();
  svreal_t wt1r, wt1i, wt2r, wt2i, wt3r, wt3i, wt4r, wt4i;

  for (int site = is; site < ns; ++site) {
    for (int ic = 0; ic < NC; ++ic) {
      int idx = 2 * ND * ic + NVCD * site;

      load_vec(pg, wt1r, &wp[VLEN*(ID1R + idx)]);
      scal_vec(pg, wt1r,  jd);
      save_vec(pg, &vp[VLEN*(ID1I + idx)], wt1r);
      
      load_vec(pg, wt1i, &wp[VLEN*(ID1I + idx)]);
      scal_vec(pg, wt1i, -jd);
      save_vec(pg, &vp[VLEN*(ID1R + idx)], wt1i);
      
      load_vec(pg, wt2r, &wp[VLEN*(ID2R + idx)]);
      scal_vec(pg, wt2r, -jd);
      save_vec(pg, &vp[VLEN*(ID2I + idx)], wt2r);
      
      load_vec(pg, wt2i, &wp[VLEN*(ID2I + idx)]);
      scal_vec(pg, wt2i,  jd);
      save_vec(pg, &vp[VLEN*(ID2R + idx)], wt2i);
      
      load_vec(pg, wt3r, &wp[VLEN*(ID3R + idx)]);
      scal_vec(pg, wt3r,  jd);
      save_vec(pg, &vp[VLEN*(ID3I + idx)], wt3r);
      
      load_vec(pg, wt3i, &wp[VLEN*(ID3I + idx)]);
      scal_vec(pg, wt3i, -jd);
      save_vec(pg, &vp[VLEN*(ID3R + idx)], wt3i);
      
      load_vec(pg, wt4r, &wp[VLEN*(ID4R + idx)]);
      scal_vec(pg, wt4r, -jd);
      save_vec(pg, &vp[VLEN*(ID4I + idx)], wt4r);
      
      load_vec(pg, wt4i, &wp[VLEN*(ID4I + idx)]);
      scal_vec(pg, wt4i,  jd);
      save_vec(pg, &vp[VLEN*(ID4R + idx)], wt4i);
      
    }
  }

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::mult_gm23(real_t* vp, const real_t* wp,
					 const real_t jd)
{
  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  svbool_t pg = set_predicate();
  svreal_t wt1r, wt1i, wt2r, wt2i, wt3r, wt3i, wt4r, wt4i;

  for (int site = is; site < ns; ++site) {
    for (int ic = 0; ic < NC; ++ic) {
      int idx = 2 * ND * ic + NVCD * site;

      load_vec(pg, wt1r, &wp[VLEN*(ID1R + idx)]);
      scal_vec(pg, wt1r,  jd);
      save_vec(pg, &vp[VLEN*(ID2I + idx)], wt1r);
      
      load_vec(pg, wt1i, &wp[VLEN*(ID1I + idx)]);
      scal_vec(pg, wt1i, -jd);
      save_vec(pg, &vp[VLEN*(ID2R + idx)], wt1i);
      
      load_vec(pg, wt2r, &wp[VLEN*(ID2R + idx)]);
      scal_vec(pg, wt2r,  jd);
      save_vec(pg, &vp[VLEN*(ID1I + idx)], wt2r);
      
      load_vec(pg, wt2i, &wp[VLEN*(ID2I + idx)]);
      scal_vec(pg, wt2i, -jd);
      save_vec(pg, &vp[VLEN*(ID1R + idx)], wt2i);
      
      load_vec(pg, wt3r, &wp[VLEN*(ID3R + idx)]);
      scal_vec(pg, wt3r,  jd);
      save_vec(pg, &vp[VLEN*(ID4I + idx)], wt3r);
      
      load_vec(pg, wt3i, &wp[VLEN*(ID3I + idx)]);
      scal_vec(pg, wt3i, -jd);
      save_vec(pg, &vp[VLEN*(ID4R + idx)], wt3i);
      
      load_vec(pg, wt4r, &wp[VLEN*(ID4R + idx)]);
      scal_vec(pg, wt4r,  jd);
      save_vec(pg, &vp[VLEN*(ID3I + idx)], wt4r);
      
      load_vec(pg, wt4i, &wp[VLEN*(ID4I + idx)]);
      scal_vec(pg, wt4i, -jd);
      save_vec(pg, &vp[VLEN*(ID3R + idx)], wt4i);
      
    }
  }

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::mult_gm31(real_t* vp, const real_t* wp,
					 const real_t jd)
{
  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  svbool_t pg = set_predicate();
  svreal_t wt1r, wt1i, wt2r, wt2i, wt3r, wt3i, wt4r, wt4i;

  for (int site = is; site < ns; ++site) {
    for (int ic = 0; ic < NC; ++ic) {
      int idx = 2 * ND * ic + NVCD * site;

      load_vec(pg, wt1r, &wp[VLEN*(ID1R + idx)]);
      scal_vec(pg, wt1r, -jd);
      save_vec(pg, &vp[VLEN*(ID2R + idx)], wt1r);
      
      load_vec(pg, wt1i, &wp[VLEN*(ID1I + idx)]);
      scal_vec(pg, wt1i, -jd);
      save_vec(pg, &vp[VLEN*(ID2I + idx)], wt1i);
      
      load_vec(pg, wt2r, &wp[VLEN*(ID2R + idx)]);
      scal_vec(pg, wt2r,  jd);
      save_vec(pg, &vp[VLEN*(ID1R + idx)], wt2r);
      
      load_vec(pg, wt2i, &wp[VLEN*(ID2I + idx)]);
      scal_vec(pg, wt2i,  jd);
      save_vec(pg, &vp[VLEN*(ID1I + idx)], wt2i);
      
      load_vec(pg, wt3r, &wp[VLEN*(ID3R + idx)]);
      scal_vec(pg, wt3r, -jd);
      save_vec(pg, &vp[VLEN*(ID4R + idx)], wt3r);
      
      load_vec(pg, wt3i, &wp[VLEN*(ID3I + idx)]);
      scal_vec(pg, wt3i, -jd);
      save_vec(pg, &vp[VLEN*(ID4I + idx)], wt3i);
      
      load_vec(pg, wt4r, &wp[VLEN*(ID4R + idx)]);
      scal_vec(pg, wt4r,  jd);
      save_vec(pg, &vp[VLEN*(ID3R + idx)], wt4r);
      
      load_vec(pg, wt4i, &wp[VLEN*(ID4I + idx)]);
      scal_vec(pg, wt4i,  jd);
      save_vec(pg, &vp[VLEN*(ID3I + idx)], wt4i);
      
    }
  }

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::mult_gm14_dirac(
                                         real_t* vp, const real_t* wp,
					 const real_t jd)
{
  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  svbool_t pg = set_predicate();
  svreal_t wt1r, wt1i, wt2r, wt2i, wt3r, wt3i, wt4r, wt4i;

  for (int site = is; site < ns; ++site) {
    for (int ic = 0; ic < NC; ++ic) {
      int idx = 2 * ND * ic + NVCD * site;

      load_vec(pg, wt1r, &wp[VLEN*(ID1R + idx)]);
      scal_vec(pg, wt1r, -jd);
      save_vec(pg, &vp[VLEN*(ID4I + idx)], wt1r);
      
      load_vec(pg, wt1i, &wp[VLEN*(ID1I + idx)]);
      scal_vec(pg, wt1i,  jd);
      save_vec(pg, &vp[VLEN*(ID4R + idx)], wt1i);
      
      load_vec(pg, wt2r, &wp[VLEN*(ID2R + idx)]);
      scal_vec(pg, wt2r, -jd);
      save_vec(pg, &vp[VLEN*(ID3I + idx)], wt2r);
      
      load_vec(pg, wt2i, &wp[VLEN*(ID2I + idx)]);
      scal_vec(pg, wt2i,  jd);
      save_vec(pg, &vp[VLEN*(ID3R + idx)], wt2i);
      
      load_vec(pg, wt3r, &wp[VLEN*(ID3R + idx)]);
      scal_vec(pg, wt3r, -jd);
      save_vec(pg, &vp[VLEN*(ID2I + idx)], wt3r);
      
      load_vec(pg, wt3i, &wp[VLEN*(ID3I + idx)]);
      scal_vec(pg, wt3i,  jd);
      save_vec(pg, &vp[VLEN*(ID2R + idx)], wt3i);
      
      load_vec(pg, wt4r, &wp[VLEN*(ID4R + idx)]);
      scal_vec(pg, wt4r, -jd);
      save_vec(pg, &vp[VLEN*(ID1I + idx)], wt4r);
      
      load_vec(pg, wt4i, &wp[VLEN*(ID4I + idx)]);
      scal_vec(pg, wt4i,  jd);
      save_vec(pg, &vp[VLEN*(ID1R + idx)], wt4i);
      
    }
  }

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::mult_gm24_dirac(
                                         real_t* vp, const real_t* wp,
					 const real_t jd)
{
  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  svbool_t pg = set_predicate();
  svreal_t wt1r, wt1i, wt2r, wt2i, wt3r, wt3i, wt4r, wt4i;

  for (int site = is; site < ns; ++site) {
    for (int ic = 0; ic < NC; ++ic) {
      int idx = 2 * ND * ic + NVCD * site;

      load_vec(pg, wt1r, &wp[VLEN*(ID1R + idx)]);
      scal_vec(pg, wt1r,  jd);
      save_vec(pg, &vp[VLEN*(ID4R + idx)], wt1r);
      
      load_vec(pg, wt1i, &wp[VLEN*(ID1I + idx)]);
      scal_vec(pg, wt1i,  jd);
      save_vec(pg, &vp[VLEN*(ID4I + idx)], wt1i);
      
      load_vec(pg, wt2r, &wp[VLEN*(ID2R + idx)]);
      scal_vec(pg, wt2r, -jd);
      save_vec(pg, &vp[VLEN*(ID3R + idx)], wt2r);
      
      load_vec(pg, wt2i, &wp[VLEN*(ID2I + idx)]);
      scal_vec(pg, wt2i, -jd);
      save_vec(pg, &vp[VLEN*(ID3I + idx)], wt2i);
      
      load_vec(pg, wt3r, &wp[VLEN*(ID3R + idx)]);
      scal_vec(pg, wt3r,  jd);
      save_vec(pg, &vp[VLEN*(ID2R + idx)], wt3r);
      
      load_vec(pg, wt3i, &wp[VLEN*(ID3I + idx)]);
      scal_vec(pg, wt3i,  jd);
      save_vec(pg, &vp[VLEN*(ID2I + idx)], wt3i);
      
      load_vec(pg, wt4r, &wp[VLEN*(ID4R + idx)]);
      scal_vec(pg, wt4r, -jd);
      save_vec(pg, &vp[VLEN*(ID1R + idx)], wt4r);
      
      load_vec(pg, wt4i, &wp[VLEN*(ID4I + idx)]);
      scal_vec(pg, wt4i, -jd);
      save_vec(pg, &vp[VLEN*(ID1I + idx)], wt4i);
      
    }
  }

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::mult_gm34_dirac(
                                         real_t* vp, const real_t* wp,
					 const real_t jd)
{
  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  svbool_t pg = set_predicate();
  svreal_t wt1r, wt1i, wt2r, wt2i, wt3r, wt3i, wt4r, wt4i;

  for (int site = is; site < ns; ++site) {
    for (int ic = 0; ic < NC; ++ic) {
      int idx = 2 * ND * ic + NVCD * site;

      load_vec(pg, wt1r, &wp[VLEN*(ID1R + idx)]);
      scal_vec(pg, wt1r, -jd);
      save_vec(pg, &vp[VLEN*(ID3I + idx)], wt1r);
      
      load_vec(pg, wt1i, &wp[VLEN*(ID1I + idx)]);
      scal_vec(pg, wt1i,  jd);
      save_vec(pg, &vp[VLEN*(ID3R + idx)], wt1i);
      
      load_vec(pg, wt2r, &wp[VLEN*(ID2R + idx)]);
      scal_vec(pg, wt2r,  jd);
      save_vec(pg, &vp[VLEN*(ID4I + idx)], wt2r);
      
      load_vec(pg, wt2i, &wp[VLEN*(ID2I + idx)]);
      scal_vec(pg, wt2i, -jd);
      save_vec(pg, &vp[VLEN*(ID4R + idx)], wt2i);
      
      load_vec(pg, wt3r, &wp[VLEN*(ID3R + idx)]);
      scal_vec(pg, wt3r, -jd);
      save_vec(pg, &vp[VLEN*(ID1I + idx)], wt3r);
      
      load_vec(pg, wt3i, &wp[VLEN*(ID3I + idx)]);
      scal_vec(pg, wt3i,  jd);
      save_vec(pg, &vp[VLEN*(ID1R + idx)], wt3i);
      
      load_vec(pg, wt4r, &wp[VLEN*(ID4R + idx)]);
      scal_vec(pg, wt4r,  jd);
      save_vec(pg, &vp[VLEN*(ID2I + idx)], wt4r);
      
      load_vec(pg, wt4i, &wp[VLEN*(ID4I + idx)]);
      scal_vec(pg, wt4i, -jd);
      save_vec(pg, &vp[VLEN*(ID2R + idx)], wt4i);
      
    }
  }

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
double AFopr_CloverTerm<AFIELD>::flop_count(const std::string mode)
{
  // The following counting explicitly depends on the implementation.
  // It will be recalculated when the code is modified.
  /*
  int    Lvol = CommonParameters::Lvol();
  double flop_wilson, flop_clover, flop_site, flop;

  if (m_repr == "Dirac") {
    flop_wilson = static_cast<double>(
      m_Nc * m_Nd * (4                        // aypx
                     + 6 * (4 * m_Nc + 2)     // spatial hopping
                     + 2 * (4 * m_Nc + 1)));  // temporal hopping

    // clover term mult assumes rotation to chiral repr.
    flop_clover = static_cast<double>(
      m_Nc * m_Nd * (2                                 // Dirac -> chiral
                     + 2 * (2 * (m_Nc * m_Nd - 1) + 1) // clover term mult
                     + 2                               // chiral -> Dirac
                     + 2));                            // addition to vector
  } else if (m_repr == "Chiral") {
    flop_wilson = static_cast<double>(
      m_Nc * m_Nd * (4 + 8 * (4 * m_Nc + 2)));

    flop_clover = static_cast<double>(
      m_Nc * m_Nd * (2 * (2 * (m_Nc * m_Nd - 1) + 1) // clover term mult
                     + 2));                          // addition to vector
  } else {
    vout.crucial(m_vl, "%s: input repr is undefined.\n");
    exit(EXIT_FAILURE);
  }

  flop_site = flop_wilson + flop_clover;

  flop = flop_site * static_cast<double>(Lvol);
  if ((mode == "DdagD") || (mode == "DDdag")) flop *= 2.0;
  */
  double flop = 0.0;  // not yet implemented.
  return flop;
}


//============================================================END=====
