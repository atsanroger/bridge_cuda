/*!
      @file    afopr_Clover_eo-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#include "lib/Field/index_eo.h"
#include "lib/ResourceManager/threadManager.h"

template<typename AFIELD>
const std::string AFopr_Clover_eo<AFIELD>::class_name = "AFopr_Clover_eo";

//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::init(const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  int req_comm = 1;  // set 1 if communication forced any time
  //int req_comm = 0;

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

  m_Nd = CommonParameters::Nd();
  m_vl = CommonParameters::Vlevel();

  m_Ndf  = 2 * m_Nc * m_Nc;
  m_Nvc  = m_Nc * 2;
  m_Ndim = CommonParameters::Ndim();
  m_Nx   = CommonParameters::Nx();
  m_Ny   = CommonParameters::Ny();
  m_Nz   = CommonParameters::Nz();
  m_Nt   = CommonParameters::Nt();
  m_Nst  = CommonParameters::Nvol();

  m_Nx2   = m_Nx / 2;
  m_Nst2  = m_Nst / 2;
  m_Nx2v  = m_Nx2 / VLENX;
  m_Nyv   = m_Ny / VLENY;
  m_Nst2v = m_Nst2 / VLEN;

  // condition check
  if (m_Nx % 2 != 0) {
    vout.crucial(m_vl, "%s: Nx must be even.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
  if (m_Nx % (2 * VLENX) != 0) {
    vout.crucial(m_vl, "%s: Nx must be mulriple of 2*VLEN.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
  if (m_Ny % (VLENY) != 0) {
    vout.crucial(m_vl, "%s: Ny must be multiple of VLENY.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  vout.general(m_vl, "  VLENX = %2d  Nx2v  = %d\n", VLENX, m_Nx2v);
  vout.general(m_vl, "  VLENY = %2d  Nyv   = %d\n", VLENY, m_Nyv);
  vout.general(m_vl, "  VLEN  = %2d  Nst2v = %d\n", VLEN, m_Nst2v);

  m_Leo.resize(m_Ny * m_Nz * m_Nt);

  int ipe3 = Communicator::ipe(3);
  int ipe2 = Communicator::ipe(2);
  int ipe1 = Communicator::ipe(1);
  for (int t = 0; t < m_Nt; ++t) {
    for (int z = 0; z < m_Nz; ++z) {
      for (int y = 0; y < m_Ny; ++y) {
        int t2 = ipe3 * m_Nt + t;
        int z2 = ipe2 * m_Nz + z;
        int y2 = ipe1 * m_Ny + y;
        m_Leo[y + m_Ny * (z + m_Nz * t)] = (y2 + z2 + t2) % 2;
      }
    }
  }

  m_Nsize[0] = m_Nx2v;
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

  // setup of communication buffers
  m_bdsize.resize(m_Ndim);
  int Nd2 = m_Nd / 2;
  m_bdsize[0] = m_Nvc * Nd2 * ((m_Ny * m_Nz * m_Nt + 1) / 2);
  m_bdsize[1] = m_Nvc * Nd2 * m_Nx2 * m_Nz * m_Nt;
  m_bdsize[2] = m_Nvc * Nd2 * m_Nx2 * m_Ny * m_Nt;
  m_bdsize[3] = m_Nvc * Nd2 * m_Nx2 * m_Ny * m_Nz;

  setup_channels();

  // gauge configuration.
  m_Ueo.reset(m_Ndf, m_Nst, m_Ndim);
  m_Ulex.reset(m_Ndf, m_Nst, m_Ndim);

  int Ndm2 = m_Nd * m_Nd / 2;

  vout.increase_indent();

  Parameters params_ct = params;

#ifdef CHIRAL_ROTATION
  m_Te_inv.reset(NDF * Ndm2 / 2, m_Nst2, 1);
  m_To_inv.reset(NDF * Ndm2 / 2, m_Nst2, 1);
  params_ct.set_string("gamma_matrix_type", "Chiral");
  m_fopr_ct = new Fopr_CloverTerm_eo(params_ct);
#else
  m_Te_inv.reset(NDF * Ndm2, m_Nst2, 1);
  m_To_inv.reset(NDF * Ndm2, m_Nst2, 1);
  m_fopr_ct = new Fopr_CloverTerm_eo(params_ct);
#endif

  set_parameters(params);

  vout.decrease_indent();

  // working vectors.
  int NinF = 2 * m_Nc * m_Nd;
  m_v1.reset(NinF, m_Nst2, 1);
  m_v2.reset(NinF, m_Nst2, 1);

  m_z1.reset(NinF, m_Nst, 1);  // used in convert/revese.

  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::setup_channels()
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
void AFopr_Clover_eo<AFIELD>::tidyup()
{
  ThreadManager::assert_single_thread(class_name);

  delete m_fopr_ct;
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
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

  set_parameters(real_t(kappa), csw, bc);

  Parameters params_ct = params;
#ifdef CHIRAL_ROTATION
  params_ct.set_string("gamma_matrix_type", "Chiral");
#endif
  m_fopr_ct->set_parameters(params_ct);
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::set_parameters(const real_t CKs,
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
void AFopr_Clover_eo<AFIELD>::get_parameters(Parameters& params) const
{
  params.set_double("hopping_parameter", double(m_CKs));
  params.set_double("clover_coefficient", double(m_csw));
  params.set_int_vector("boundary_condition", m_boundary);
  params.set_string("gamma_matrix_type", m_repr);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::set_config(Field *u)
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
void AFopr_Clover_eo<AFIELD>::set_config_omp(Field *u)
{
  vout.detailed(m_vl, "  set_config_omp is called.\n");

#pragma omp parallel
  {
    set_config_impl(u);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::set_config_impl(Field *u)
{
#pragma omp barrier

  int ith, nth, is, ns;
  ith = ThreadManager::get_thread_id();

  if (ith == 0) m_conf = u;

  AIndex_lex<real_t, QXS> index_lex;

  convert_gauge(index_lex, m_Ulex, *u);

  QXS_Gauge::set_boundary(m_Ulex, m_boundary);

  AIndex_lex<real_t, QXS> index_lex2;
  AIndex_eo<real_t, QXS>  index_eo;

  set_threadtask_mult(ith, nth, is, ns, m_Nst);

  for (int ex = 0; ex < m_Ndim; ++ex) {
    for (int site = is; site < ns; ++site) {
      for (int in = 0; in < m_Ndf; ++in) {
        int iv1 = index_lex2.idx(in, m_Ndf, site, ex);
        int iv2 = index_eo.idx(in, m_Ndf, site, ex);
        m_Ueo.e(iv2) = m_Ulex.cmp(iv1);
      }
    }
  }

  // setup of clover term inverse
  m_fopr_ct->set_config(u);

#ifdef CHIRAL_ROTATION
  set_csw_chrot(m_Te_inv, 0);
  set_csw_chrot(m_To_inv, 1);
#else
  set_csw();
#endif

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::set_csw()
{
  vout.crucial(m_vl, "%s: set_csw() has not been implemented.\n",
               class_name.c_str());
  exit(EXIT_FAILURE);
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::set_csw_chrot(AFIELD& m_T_inv, int ieo)
{
  // This method set the clover term matrix assuming the chiral
  // representation for gamma-matrix in Bridge++ convention.
  // The conversion from the Bridge++ to QWS convention and
  // from Dirac to chiral representations are assumed to be cared
  // in BridgeQXS functions.
  //                                       [04 Sep 2021 H.Matsufuru]

  //  Field_F v(m_Nst2, 1), w(m_Nst2, 1);

  AIndex_eo<real_t, AFIELD::IMPL> index;

  //  const int Nin = NDF * ND * 2;
  const int Nin = NDF * ND;

  m_fopr_ct->set_mode("D");

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nst2);

//11 22 33 44 55 66   12 34 56  23 45   13 24 35 46  15 26  14 36 25  16
  constexpr int idx[72] = {
    0,  -1,  6,  7, 16, 17, 28, 29, 24, 25, 34, 35,
    -1, -1,  1, -1, 12, 13, 18, 19, 32, 33, 26, 27,
    -1, -1, -1, -5,  2, -1,  8,  9, 20, 21, 30, 31,
    -1, -1, -1, -1, -1, -1,  3, -1, 14, 15, 22, 23,
    -1, -1, -1, -1, -1, -1, -1, -1,  4, -1, 10, 11,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  5, -1,
  };
//#pragma omp parallel
//  {
  for (int id = 0; id < m_Nd / 2; ++id) {
    for (int ic = 0; ic < m_Nc; ++ic) {
      m_w1.set(0.0);
#pragma omp barrier

      for (int site = is; site < ns; ++site) {
        m_w1.set_r(ic, id, site, 0, 1.0);
        m_w1.set_r(ic, id + 2, site, 0, 1.0);
      }
#pragma omp barrier

      m_fopr_ct->mult_csw_inv(m_w2, m_w1, ieo);
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
              int idx_r = index.idxh(in_r, Nin, site, 0);
              m_T_inv.set(idx_r, vt_r);
            }
            if (in_i >= 0) {
              in_i += 36 * (id2 / 2);
              int idx_i = index.idxh(in_i, Nin, site, 0);
              m_T_inv.set(idx_i, vt_i);
            }
          }
        }
      }   // site loop
#pragma omp barrier
    }
  }
#pragma omp barrier

  //  } // omp parallel
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::convert(AFIELD& v, const Field& w)
{
  AIndex_lex<real_t, AFIELD::IMPL> index_lex;
  convert_spinor(index_lex, m_z1, w);

  mult_gm4(v, m_z1);
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::reverse(Field& v, const AFIELD& w)
{
  mult_gm4(m_z1, w);

  AIndex_lex<real_t, AFIELD::IMPL> index_lex;
  reverse_spinor(index_lex, v, m_z1);
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::set_mode(std::string mode)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_mode = mode;

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::mult(AFIELD& v, const AFIELD& w)
{
  if (m_mode == "D") {
    return D(v, w);
  } else if (m_mode == "DdagD") {
    return DdagD(v, w);
  } else if (m_mode == "Ddag") {
    return Ddag(v, w);
  } else {
    vout.crucial(m_vl, "%s: mode undefined.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::mult_dag(AFIELD& v, const AFIELD& w)
{
  if (m_mode == "D") {
    return Ddag(v, w);
  } else if (m_mode == "DdagD") {
    return DdagD(v, w);
  } else if (m_mode == "Ddag") {
    return D(v, w);
  } else {
    vout.crucial(m_vl, "%s: mode undefined.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::mult_gm5(AFIELD& v, const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

#pragma omp barrier

  BridgeQXS::mult_wilson_gm5_dirac(vp, wp, m_Nsize);

#pragma omp barrier

  //#pragma omp barrier
  //  mult_gm5(vp, wp);
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::mult_gm4(AFIELD& v, const AFIELD& w)
{
  int Nst = w.nvol();
  if (Nst != v.nvol()) {
    vout.crucial(m_vl, "%s: sizes of fields inconsisutent in mult_gm4.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nst);

#pragma omp barrier

  // Dirac representation.
  for (int site = is; site < ns; ++site) {
    for (int ic = 0; ic < NC; ++ic) {
      vp[m_index.idxh_SPr(ic, 0, site, 0)] = wp[m_index.idxh_SPr(ic, 0, site, 0)];
      vp[m_index.idxh_SPi(ic, 0, site, 0)] = wp[m_index.idxh_SPi(ic, 0, site, 0)];
      vp[m_index.idxh_SPr(ic, 1, site, 0)] = wp[m_index.idxh_SPr(ic, 1, site, 0)];
      vp[m_index.idxh_SPi(ic, 1, site, 0)] = wp[m_index.idxh_SPi(ic, 1, site, 0)];
      vp[m_index.idxh_SPr(ic, 2, site, 0)] = -wp[m_index.idxh_SPr(ic, 2, site, 0)];
      vp[m_index.idxh_SPi(ic, 2, site, 0)] = -wp[m_index.idxh_SPi(ic, 2, site, 0)];
      vp[m_index.idxh_SPr(ic, 3, site, 0)] = -wp[m_index.idxh_SPr(ic, 3, site, 0)];
      vp[m_index.idxh_SPi(ic, 3, site, 0)] = -wp[m_index.idxh_SPi(ic, 3, site, 0)];
    }
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::mult(AFIELD& v, const AFIELD& w,
                                   const std::string mode)
{
  if (mode == "Dee_inv") {
    mult_cswinv(v, w, 0);
  } else if (mode == "Doo_inv") {
    mult_cswinv(v, w, 1);
  } else if (mode == "Deo") {
    Meo(v, w, 0);
  } else if (mode == "Doe") {
    Meo(v, w, 1);
  } else {
    vout.crucial(m_vl, "%s: illegal mode is given to mult with mode\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::DdagD(AFIELD& v, const AFIELD& w)
{
  D(m_v2, w);
  Ddag(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::Ddag(AFIELD& v, const AFIELD& w)
{
#pragma omp barrier

  mult_gm5(m_v1, w);
  mult_cswinv(v, m_v1, 0);
  Meo(m_v1, v, 1);
  mult_cswinv(v, m_v1, 1);
  Meo(m_v1, v, 0);
  mult_gm5(v, m_v1);
  aypx(-1.0, v, w);
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::D(AFIELD& v, const AFIELD& w)
{
#pragma omp barrier
  Meo(m_v1, w, 1);
  mult_cswinv(v, m_v1, 1);
  Meo(m_v1, v, 0);
  mult_cswinv(v, m_v1, 0);
  aypx(-1.0, v, w);
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::aypx(real_t a, AFIELD& v, const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

#pragma omp barrier

  aypx(a, vp, wp);
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::Meo(AFIELD& v, const AFIELD& w,
                                  const int ieo)
{
  mult_Meo_qxs(v, w, w, ieo, 0);
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::Meo(AFIELD& v,
                                  const AFIELD& w, const AFIELD& x,
                                  const int ieo, const int iflag)
{
  mult_Meo_qxs(v, w, x, ieo, iflag);
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::mult_Meo_qxs(AFIELD& v,
                                           const AFIELD& w, const AFIELD& x,
                                           const int ieo, const int iflag)
{
  real_t *v2 = v.ptr(0);
  real_t *v1 = const_cast<AFIELD *>(&w)->ptr(0);
  real_t *xp = const_cast<AFIELD *>(&x)->ptr(0);

  real_t *up = m_Ueo.ptr(0);

  int ith = ThreadManager::get_thread_id();

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
    BridgeQXS::mult_wilson_eo_1_dirac(buf1_xp, buf1_xm, buf1_yp, buf1_ym,
                                      buf1_zp, buf1_zm, buf1_tp, buf1_tm,
                                      up, v1, &m_boundary[0],
                                      m_Nsize, do_comm, &m_Leo[0], ieo, iflag);

#pragma omp barrier

    if (ith == 0) chset_send.start();
  }

  BridgeQXS::mult_wilson_eo_bulk_dirac(v2, up, v1, xp, m_CKs, &m_boundary[0],
                                       m_Nsize, do_comm, &m_Leo[0], ieo, iflag);

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

    BridgeQXS::mult_wilson_eo_2_dirac(v2, up, v1, xp,
                                      buf2_xp, buf2_xm, buf2_yp, buf2_ym,
                                      buf2_zp, buf2_zm, buf2_tp, buf2_tm,
                                      m_CKs, &m_boundary[0],
                                      m_Nsize, do_comm, &m_Leo[0], ieo, iflag);

    if (ith == 0) chset_send.wait();
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::mult_cswinv(AFIELD& v, const AFIELD& w,
                                          const int ieo)
{
  real_t *v2 = v.ptr(0);
  real_t *v1 = const_cast<AFIELD *>(&w)->ptr(0);
  real_t *ct;
  if (ieo == 0) {
    ct = m_Te_inv.ptr(0);
  } else {
    ct = m_To_inv.ptr(0);
  }

#pragma omp barrier

#ifdef CHIRAL_ROTATION
  BridgeQXS::mult_clover_cswinv_dirac_chrot(v2, ct, v1, m_Nsize);
#endif

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::aypx(real_t a, real_t *__restrict v,
                                   real_t *__restrict w)
{
  int Ncd = m_Nvc * m_Nd;
#ifdef USE_QXS_ACLE
  svbool_t pg = set_predicate();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nst2v);

  for (int site = is; site < ns; ++site) {
    real_t *__restrict v2 = v + VLEN * Ncd * site;
    real_t *__restrict w2 = w + VLEN * Ncd * site;
    for (int icd = 0; icd < Ncd / 2; ++icd) {
      svreal_t vtr, vti;
      svreal_t wtr, wti;
      load_vec(pg, wtr, &w2[2 * VLEN * icd]);
      load_vec(pg, wti, &w2[2 * VLEN * icd + VLEN]);
      load_vec(pg, vtr, &v2[2 * VLEN * icd]);
      load_vec(pg, vti, &v2[2 * VLEN * icd + VLEN]);
      aypx_vec(pg, a, vtr, wtr);
      aypx_vec(pg, a, vti, wti);
      save_vec(pg, &v2[2 * VLEN * icd], vtr);
      save_vec(pg, &v2[2 * VLEN * icd + VLEN], vti);
    }
  }
#else
  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nst2);
  for (int i = is; i < ns; ++i) {
    for (int icd = 0; icd < Ncd; ++icd) {
      int i2 = icd + Ncd * i;
      v[i2] = a * v[i2] + w[i2];
    }
  }
#endif

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::clear(real_t *v)
{
  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nst2);

  real_t zero = 0.0;
  int    Ncd  = m_Nvc * m_Nd;

  for (int i = is; i < ns; ++i) {
    for (int icd = 0; icd < Ncd; ++icd) {
      v[icd + Ncd * i] = zero;
    }
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::scal(real_t *v, const real_t a)
{
  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nst2);

  int Ncd = m_Nvc * m_Nd;

  for (int i = is; i < ns; ++i) {
    for (int icd = 0; icd < Ncd; ++icd) {
      v[icd + Ncd * i] *= a;
    }
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
double AFopr_Clover_eo<AFIELD>::flop_count(const std::string mode)
{
  // The following counting explicitly depends on the implementation.
  // It will be recalculated when the code is modified.
  // The present counting is copied from AFopr_Clover in QXS impl.

  size_t Lvol = CommonParameters::Lvol();
  double flop_wilson, flop_clover, flop_site, flop;

  if (m_repr == "Dirac") {
    flop_wilson = static_cast<double>(
      m_Nc * m_Nd * (4                          // aypx
                     + 6 * (4 * m_Nc + 2)       // spatial hopping
                     + 2 * (4 * m_Nc + 1)));    // temporal hopping

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
      m_Nc * m_Nd * (2 * (2 * (m_Nc * m_Nd - 1) + 1)   // clover term mult
                     + 2));                            // addition to vector
  } else {
    vout.crucial(m_vl, "%s: input repr is undefined.\n");
    exit(EXIT_FAILURE);
  }

  flop_site = flop_wilson + flop_clover;

  flop = flop_site * static_cast<double>(Lvol);
  if ((mode == "DdagD") || (mode == "DDdag")) flop *= 2.0;

  return flop;
}


//============================================================END=====
