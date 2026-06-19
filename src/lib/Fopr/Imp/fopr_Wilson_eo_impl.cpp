/*!
        @file    fopr_Wilson_eo_impl.cpp

        @brief

        @author  Satoru UEDA (maintailed by H.Matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "fopr_Wilson_eo_impl.h"

#include "Fopr/fopr_thread-inc.h"

#if defined USE_GROUP_SU3
#include "fopr_Wilson_impl_SU3-inc.h"
#elif defined USE_GROUP_SU2
#include "fopr_Wilson_impl_SU2-inc.h"
#elif defined USE_GROUP_SU_N
#include "fopr_Wilson_impl_SU_N-inc.h"
#endif

#include "Fopr/Imp/fopr_Wilson_impl_common-inc.h"

namespace Imp {
#ifdef USE_FACTORY_AUTOREGISTER
  namespace {
    bool init = Fopr_Wilson_eo::register_factory();
  }
#endif

  const std::string Fopr_Wilson_eo::class_name = "Imp::Fopr_Wilson_eo";

//====================================================================
  void Fopr_Wilson_eo::init(const Parameters& params)
  {
    ThreadManager::assert_single_thread(class_name);

    m_vl = CommonParameters::Vlevel();

    vout.general(m_vl, "%s: construction\n", class_name.c_str());
    vout.increase_indent();

    setup();

    std::string repr;
    if (!params.fetch_string("gamma_matrix_type", repr)) {
      m_repr = repr;
    } else {
      m_repr = "Dirac"; // default gamma-matrix type
      vout.general(m_vl, "gamma_matrix_type is not given: defalt = %s\n",
                   m_repr.c_str());
    }

    if ((m_repr != "Dirac") && (m_repr != "Chiral")) {
      vout.crucial("Error at %s: unsupported gamma-matrix type: %s\n",
                   class_name.c_str(), m_repr.c_str());
      exit(EXIT_FAILURE);
    }

    set_parameters(params);

    vout.decrease_indent();
    vout.general(m_vl, "%s: construction finished.\n",
                 class_name.c_str());
  }


//====================================================================
  void Fopr_Wilson_eo::init(const std::string repr)
  {
    ThreadManager::assert_single_thread(class_name);

    m_vl = CommonParameters::Vlevel();

    vout.general(m_vl, "%s: construction (obsolete)\n",
                 class_name.c_str());

    setup();

    m_repr = repr;
    if ((m_repr != "Dirac") && (m_repr != "Chiral")) {
      vout.crucial(m_vl, "Error at %s: input repr is undefined.\n",
                   class_name.c_str());
      exit(EXIT_FAILURE);
    }

    vout.general(m_vl, "  gamma-matrix type is set to %s\n",
                 m_repr.c_str());

    vout.general(m_vl, "%s: construction finished.\n",
                 class_name.c_str());
  }


//====================================================================
  void Fopr_Wilson_eo::setup()
  {
    m_Nc  = CommonParameters::Nc();
    m_Nd  = CommonParameters::Nd();
    m_Nvc = 2 * m_Nc;
    m_Ndf = 2 * m_Nc * m_Nc;
    const int Nvcd = m_Nvc * m_Nd;

    m_Nvol  = CommonParameters::Nvol();
    m_Nvol2 = m_Nvol / 2;
    m_Ndim  = CommonParameters::Ndim();

    m_Nx = CommonParameters::Nx();
    m_Ny = CommonParameters::Ny();
    m_Nz = CommonParameters::Nz();
    m_Nt = CommonParameters::Nt();

    if ((m_Nx % 2) != 0) {
      vout.crucial(m_vl, "Error at %s: Nx=%d must be even.\n",
                   class_name.c_str(), m_Nx);
      exit(EXIT_FAILURE);
    }
    m_Nx2 = m_Nx / 2;

    if ((m_Ny % 2) != 0) {
      vout.crucial(m_vl, "Error at %s: Ny=%d must be even.\n",
                   class_name.c_str(), m_Ny);
      exit(EXIT_FAILURE);
    }

    m_Nvol = CommonParameters::Nvol();
    m_Ndim = CommonParameters::Ndim();
    m_boundary.resize(m_Ndim);
    m_boundary_each_node.resize(m_Ndim);

    m_yzt_eo.resize(m_Ny * m_Nz * m_Nt);

    for (int it = 0; it < m_Nt; ++it) {
      for (int iz = 0; iz < m_Nz; ++iz) {
        for (int iy = 0; iy < m_Ny; ++iy) {
          int it_global = it + Communicator::ipe(3) * m_Nt;
          int iz_global = iz + Communicator::ipe(2) * m_Nz;
          int iy_global = iy + Communicator::ipe(1) * m_Ny;
          m_yzt_eo[iy + m_Ny * (iz + m_Nz * it)]
            = (iy_global + iz_global + it_global) % 2;
        }
      }
    }

    //m_Ueo = new Field_G(m_Nvol, m_Ndim);
    m_Ueo.reset(m_Nvol, m_Ndim);

    const int Nvx = m_Nvc * 2 * ((m_Ny * m_Nz * m_Nt + 1) / 2);
    vcp1_xp = new double[Nvx];
    vcp2_xp = new double[Nvx];
    vcp1_xm = new double[Nvx];
    vcp2_xm = new double[Nvx];

    const int Nvy = m_Nvc * 2 * m_Nx2 * m_Nz * m_Nt;
    vcp1_yp = new double[Nvy];
    vcp2_yp = new double[Nvy];
    vcp1_ym = new double[Nvy];
    vcp2_ym = new double[Nvy];

    const int Nvz = m_Nvc * 2 * m_Nx2 * m_Ny * m_Nt;
    vcp1_zp = new double[Nvz];
    vcp2_zp = new double[Nvz];
    vcp1_zm = new double[Nvz];
    vcp2_zm = new double[Nvz];

    const int Nvt = m_Nvc * 2 * m_Nx2 * m_Ny * m_Nz;
    vcp1_tp = new double[Nvt];
    vcp2_tp = new double[Nvt];
    vcp1_tm = new double[Nvt];
    vcp2_tm = new double[Nvt];

    m_v1.reset(Nvcd, m_Nvol2, 1);
    m_v2.reset(Nvcd, m_Nvol2, 1);
    m_w1.reset(Nvcd, m_Nvol2, 1);
    m_w2.reset(Nvcd, m_Nvol2, 1);
  }


//====================================================================
  void Fopr_Wilson_eo::tidyup()
  {
    //  delete m_Ueo;

    delete[]  vcp1_xp;
    delete[]  vcp2_xp;
    delete[]  vcp1_xm;
    delete[]  vcp2_xm;

    delete[]  vcp1_yp;
    delete[]  vcp2_yp;
    delete[]  vcp1_ym;
    delete[]  vcp2_ym;

    delete[]  vcp1_zp;
    delete[]  vcp2_zp;
    delete[]  vcp1_zm;
    delete[]  vcp2_zm;

    delete[]  vcp1_tp;
    delete[]  vcp2_tp;
    delete[]  vcp1_tm;
    delete[]  vcp2_tm;
  }


//====================================================================
  void Fopr_Wilson_eo::set_parameters(const Parameters& params)
  {
#pragma omp barrier
    int         ith = ThreadManager::get_thread_id();
    std::string vlevel;
    if (!params.fetch_string("verbose_level", vlevel)) {
      if (ith == 0) m_vl = vout.set_verbose_level(vlevel);
    }
#pragma omp barrier

    // fetch and check input parameters
    double           kappa;
    std::vector<int> bc;

    int err = 0;
    err += params.fetch_double("hopping_parameter", kappa);
    err += params.fetch_int_vector("boundary_condition", bc);

    if (err) {
      vout.crucial(m_vl, "Error at %s: parameter fetch error.\n",
                   class_name.c_str());
      exit(EXIT_FAILURE);
    }

    set_parameters(kappa, bc);
  }


//====================================================================
  void Fopr_Wilson_eo::set_parameters(const double kappa,
                                      const std::vector<int> bc)
  {
    assert(bc.size() == m_Ndim);

#pragma omp barrier
    int ith = ThreadManager::get_thread_id();
    if (ith == 0) {
      m_kappa    = kappa;
      m_boundary = bc;

      for (int idir = 0; idir < m_Ndim; ++idir) {
        m_boundary_each_node[idir] = 1.0;
        if (Communicator::ipe(idir) == 0) {
          m_boundary_each_node[idir] = m_boundary[idir];
        }
      }
    }
#pragma omp barrier

    //- print input parameters
    vout.general(m_vl, "%s: input parameters\n", class_name.c_str());
    vout.general(m_vl, "  gamma matrix type = %s\n", m_repr.c_str());
    vout.general(m_vl, "  kappa = %12.8f\n", m_kappa);
    for (int mu = 0; mu < m_Ndim; ++mu) {
      vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
    }
  }


//====================================================================
  void Fopr_Wilson_eo::get_parameters(Parameters& params) const
  {
    params.set_double("hopping_parameter", m_kappa);
    params.set_int_vector("boundary_condition", m_boundary);
    params.set_string("gamma_matrix_type", m_repr);

    params.set_string("verbose_level", vout.get_verbose_level(m_vl));
  }


//====================================================================
  void Fopr_Wilson_eo::set_config(Field *U)
  {
    int nth = ThreadManager::get_num_threads();

    vout.detailed(m_vl, "%s: set_config is called: num_threads = %d\n",
                  class_name.c_str(), nth);

    if (nth > 1) {
      set_config_impl(U);
    } else {
      set_config_omp(U);
    }

    vout.detailed(m_vl, "%s: set_config finished\n", class_name.c_str());
  }


//====================================================================
  void Fopr_Wilson_eo::set_config_omp(Field *U)
  {
    vout.detailed(m_vl, "  set_config_omp is called.\n");

#pragma omp parallel
    {
      set_config_impl(U);
    }
  }


//====================================================================
  void Fopr_Wilson_eo::set_config_impl(Field *U)
  {
    m_index.convertField(m_Ueo, *U);
  }


//====================================================================
  void Fopr_Wilson_eo::set_mode(std::string mode)
  {
#pragma omp barrier
    int ith = ThreadManager::get_thread_id();
    if (ith == 0) m_mode = mode;
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::mult(Field& v, const Field& w)
  {
    if (m_mode == "D") {
      D(v, w);
    } else if (m_mode == "Ddag") {
      Ddag(v, w);
    } else if (m_mode == "DdagD") {
      DdagD(v, w);
    } else if (m_mode == "DDdag") {
      DDdag(v, w);
    } else {
      vout.crucial("Error at %s: irrelevant mult mode = %s.\n",
                   class_name.c_str(), m_mode.c_str());
      exit(EXIT_FAILURE);
    }
  }


//====================================================================
  void Fopr_Wilson_eo::mult_dag(Field& v, const Field& w)
  {
    if (m_mode == "D") {
      Ddag(v, w);
    } else if (m_mode == "Ddag") {
      D(v, w);
    } else if (m_mode == "DdagD") {
      DdagD(v, w);
    } else if (m_mode == "DDdag") {
      DDdag(v, w);
    } else {
      vout.crucial("Error at %s: irrelevant mult mode = %s.\n",
                   class_name.c_str(), m_mode.c_str());
      exit(EXIT_FAILURE);
    }
  }


//====================================================================
  void Fopr_Wilson_eo::mult(Field& v, const Field& w,
                            const std::string mode)
  {
    if (mode == "Deo") {
      Meo(v, w, 0);
    } else if (mode == "Doe") {
      Meo(v, w, 1);
    } else {
      vout.crucial("Error at %s: irrelevant mult mode = %s.\n",
                   class_name.c_str(), mode.c_str());
      exit(EXIT_FAILURE);
    }
  }


//====================================================================
  void Fopr_Wilson_eo::mult_dag(Field& v, const Field& w,
                                const std::string mode)
  {
    if (mode == "Deo") {
      Meo(v, w, 1);               // "Doe" for mult
    } else if (mode == "Doe") {
      Meo(v, w, 0);               // "Deo" for mult
    } else {
      vout.crucial("Error at %s: irrelevant mult mode = %s.\n",
                   class_name.c_str(), mode.c_str());
      exit(EXIT_FAILURE);
    }
  }


//====================================================================
  void Fopr_Wilson_eo::preProp(Field& Be, Field& bo, const Field& b)
  {
    m_index.convertField(Be, b, 0);
    m_index.convertField(bo, b, 1);

    if (m_mode == "D") {
      Meo(m_v1, bo, 0);
    } else {
      Mdageo(m_v1, bo, 0);
    }

    axpy(Be, -1.0, m_v1);

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::postProp(Field& x,
                                const Field& xe, const Field& bo)
  {
    if (m_mode == "D") {
      Meo(m_v1, xe, 1);
    } else {
      Mdageo(m_v1, xe, 1);
    }

    aypx(-1.0, m_v1, bo);
#pragma omp barrier

    m_index.reverseField(x, xe, 0);
    m_index.reverseField(x, m_v1, 1);

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::DdagD(Field& v, const Field& f)
  {
    D(m_w1, f);
    Ddag(v, m_w1);
  }


//====================================================================
  void Fopr_Wilson_eo::DDdag(Field& v, const Field& f)
  {
    Ddag(m_w1, f);
    D(v, m_w1);
  }


//====================================================================
  void Fopr_Wilson_eo::D(Field& v, const Field& w)
  {
    assert(w.nex() == 1);

    Meo(m_v1, w, 1);
    Meo(v, m_v1, 0);

    aypx(-1.0, v, w);
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::Ddag(Field& v, const Field& w)
  {
    Mdageo(m_v1, w, 1);
    Mdageo(v, m_v1, 0);

    aypx(-1.0, v, w);
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::Meo(Field& v, const Field& w, const int ieo)
  {
    if (m_repr == "Dirac") {
      Meo_dirac(v, w, ieo);
    } else {
      Meo_chiral(v, w, ieo);
    }
  }


//====================================================================
  void Fopr_Wilson_eo::Meo_dirac(Field& w, const Field& f, const int ieo)
  {
#pragma omp barrier

    w.set(0.0);
#pragma omp barrier

    mult_xp(w, f, ieo);
    mult_xm(w, f, ieo);

    mult_yp(w, f, ieo);
    mult_ym(w, f, ieo);

    mult_zp(w, f, ieo);
    mult_zm(w, f, ieo);

    mult_tp_dirac(w, f, ieo);
    mult_tm_dirac(w, f, ieo);

    scal(w, -m_kappa);
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::Meo_chiral(Field& w, const Field& f, const int ieo)
  {
#pragma omp barrier

    w.set(0.0);
#pragma omp barrier

    mult_xp(w, f, ieo);
    mult_xm(w, f, ieo);

    mult_yp(w, f, ieo);
    mult_ym(w, f, ieo);

    mult_zp(w, f, ieo);
    mult_zm(w, f, ieo);

    mult_tp_chiral(w, f, ieo);
    mult_tm_chiral(w, f, ieo);

    scal(w, -m_kappa);
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::Mdageo(Field& v, const Field& w, const int ieo)
  {
    mult_gm5(v, w);
    Meo(m_v2, v, ieo);
    mult_gm5(v, m_v2);
  }


//====================================================================
  void Fopr_Wilson_eo::Meo_gm5(Field& v, const Field& w, const int ieo)
  {
    Meo(m_v2, w, ieo);
    mult_gm5(v, m_v2);
  }


//====================================================================
  void Fopr_Wilson_eo::mult_gm5(Field& v, const Field& w)
  {
    if (m_repr == "Dirac") {
      gm5_dirac(v, w);
    } else {
      gm5_chiral(v, w);
    }
  }


//====================================================================
  void Fopr_Wilson_eo::gm5_dirac(Field& v, const Field& w)
  {
    double       *vp = v.ptr(0);
    const double *wp = w.ptr(0);

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol2);

    int Nvcd = m_Nvc * m_Nd;

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      mult_gamma5_dirac(&vp[Nvcd * site], &wp[Nvcd * site], m_Nc);
    }

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::gm5_chiral(Field& v, const Field& w)
  {
    double       *vp = v.ptr(0);
    const double *wp = w.ptr(0);

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol2);

    int Nvcd = m_Nvc * m_Nd;

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      mult_gamma5_chiral(&vp[Nvcd * site], &wp[Nvcd * site], m_Nc);
    }

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::gm5p(const int mu, Field& v, const Field& w)
  {
    // this function is probably not to be used.
    // determines  \gamma_5 * (1 - \gamma_\mu) v(x+\hat{x})
    vout.crucial(m_vl, "Error at %s: gm5p is undefined.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }


//====================================================================
  void Fopr_Wilson_eo::mult_xp(Field& v, const Field& w, const int ieo)
  {
    int idir = 0;

    int Nvc2 = m_Nvc * 2;
    int Nvcd = m_Nvc * m_Nd;

    double bc2 = m_boundary_each_node[idir];

    double       *vp = v.ptr(0);
    const double *wp = w.ptr(0);
    const double *up = m_Ueo.ptr(m_Ndf * m_Nvol2 * (ieo + 2 * idir));

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol2);

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ix2  = site % m_Nx2;
      int iyzt = site / m_Nx2;
      int keo  = (ieo + m_yzt_eo[iyzt]) % 2;
      if ((ix2 == 0) && (keo == 1)) {
        int    iyzt2 = iyzt / 2;
        int    in = Nvcd * site;
        int    ib1 = Nvc2 * iyzt2;
        int    ib2 = ib1 + NVC;
        double vt1[NVC], vt2[NVC];
        set_sp2_xp(vt1, vt2, &wp[in], m_Nc);
        for (int ivc = 0; ivc < NVC; ++ivc) {
          vcp1_xp[ivc + ib1] = bc2 * vt1[ivc];
          vcp1_xp[ivc + ib2] = bc2 * vt2[ivc];
        }
      }
    }

#pragma omp barrier

#pragma omp master
    {
      const int Nv = m_Nvc * 2 * ((m_Ny * m_Nz * m_Nt + 1) / 2);
      Communicator::exchange(Nv, vcp2_xp, vcp1_xp, 0, 1, 1);
    }
#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ix2  = site % m_Nx2;
      int iyzt = site / m_Nx2;
      int keo  = (ieo + m_yzt_eo[iyzt]) % 2;
      int ix2n = ix2 + keo;
      int nei  = ix2n + m_Nx2 * iyzt;
      int iv   = Nvcd * site;
      int ig   = m_Ndf * site;

      if (ix2n < m_Nx2) {
        int    in = Nvcd * nei;
        double vt1[NVC], vt2[NVC];
        set_sp2_xp(vt1, vt2, &wp[in], m_Nc);
        for (int ic = 0; ic < m_Nc; ++ic) {
          int    ic2  = ic * NVC;
          double wt1r = mult_uv_r(&up[ic2 + ig], vt1, m_Nc);
          double wt1i = mult_uv_i(&up[ic2 + ig], vt1, m_Nc);
          double wt2r = mult_uv_r(&up[ic2 + ig], vt2, m_Nc);
          double wt2i = mult_uv_i(&up[ic2 + ig], vt2, m_Nc);
          set_sp4_xp(&vp[2 * ic + iv], wt1r, wt1i, wt2r, wt2i, m_Nc);
        }
      } else {
        int iyzt2 = iyzt / 2;
        int ib1   = Nvc2 * iyzt2;
        int ib2   = ib1 + NVC;
        for (int ic = 0; ic < m_Nc; ++ic) {
          int    ic2  = ic * NVC;
          double wt1r = mult_uv_r(&up[ic2 + ig], &vcp2_xp[ib1], m_Nc);
          double wt1i = mult_uv_i(&up[ic2 + ig], &vcp2_xp[ib1], m_Nc);
          double wt2r = mult_uv_r(&up[ic2 + ig], &vcp2_xp[ib2], m_Nc);
          double wt2i = mult_uv_i(&up[ic2 + ig], &vcp2_xp[ib2], m_Nc);
          set_sp4_xp(&vp[2 * ic + iv], wt1r, wt1i, wt2r, wt2i, m_Nc);
        }
      }
    }

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::mult_xm(Field& v, const Field& w, const int ieo)
  {
    int idir = 0;

    int Nvc2 = m_Nvc * 2;
    int Nvcd = m_Nvc * m_Nd;

    double bc2 = m_boundary_each_node[idir];

    double       *vp = v.ptr(0);
    const double *wp = w.ptr(0);
    const double *up = m_Ueo.ptr(m_Ndf * m_Nvol2 * (1 - ieo + 2 * idir));

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol2);

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ix2  = site % m_Nx2;
      int iyzt = site / m_Nx2;
      int keo  = (ieo + m_yzt_eo[iyzt]) % 2;
      if ((ix2 == m_Nx2 - 1) && (keo == 0)) {
        int    iyzt2 = iyzt / 2;
        int    in = Nvcd * site;
        int    ig = m_Ndf * site;
        int    ib1 = Nvc2 * iyzt2;
        int    ib2 = ib1 + NVC;
        double vt1[NVC], vt2[NVC];
        set_sp2_xm(vt1, vt2, &wp[in], m_Nc);
        for (int ic = 0; ic < m_Nc; ++ic) {
          int ic2 = 2 * ic;
          int icr = 2 * ic;
          int ici = 2 * ic + 1;
          vcp1_xm[icr + ib1] = mult_udagv_r(&up[ic2 + ig], vt1, m_Nc);
          vcp1_xm[ici + ib1] = mult_udagv_i(&up[ic2 + ig], vt1, m_Nc);
          vcp1_xm[icr + ib2] = mult_udagv_r(&up[ic2 + ig], vt2, m_Nc);
          vcp1_xm[ici + ib2] = mult_udagv_i(&up[ic2 + ig], vt2, m_Nc);
        }
      }
    }

#pragma omp barrier


#pragma omp master
    {
      const int Nv = m_Nvc * 2 * ((m_Ny * m_Nz * m_Nt + 1) / 2);
      Communicator::exchange(Nv, vcp2_xm, vcp1_xm, 0, -1, 2);
    }
#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ix2  = site % m_Nx2;
      int iyzt = site / m_Nx2;
      int keo  = (ieo + m_yzt_eo[iyzt]) % 2;
      int ix2n = ix2 - (1 - keo);
      int nei  = ix2n + m_Nx2 * iyzt;
      int iv   = Nvcd * site;

      if (ix2n >= 0) {
        int    ig = m_Ndf * nei;
        int    in = Nvcd * nei;
        double vt1[NVC], vt2[NVC];
        set_sp2_xm(vt1, vt2, &wp[in], m_Nc);
        for (int ic = 0; ic < m_Nc; ++ic) {
          int    ic2  = 2 * ic;
          double wt1r = mult_udagv_r(&up[ic2 + ig], vt1, m_Nc);
          double wt1i = mult_udagv_i(&up[ic2 + ig], vt1, m_Nc);
          double wt2r = mult_udagv_r(&up[ic2 + ig], vt2, m_Nc);
          double wt2i = mult_udagv_i(&up[ic2 + ig], vt2, m_Nc);
          set_sp4_xm(&vp[2 * ic + iv], wt1r, wt1i, wt2r, wt2i, m_Nc);
        }
      } else {
        int iyzt2 = iyzt / 2;
        int ib1   = Nvc2 * iyzt2;
        int ib2   = ib1 + NVC;
        for (int ic = 0; ic < m_Nc; ++ic) {
          double wt1r = bc2 * vcp2_xm[2 * ic + ib1];
          double wt1i = bc2 * vcp2_xm[2 * ic + 1 + ib1];
          double wt2r = bc2 * vcp2_xm[2 * ic + ib2];
          double wt2i = bc2 * vcp2_xm[2 * ic + 1 + ib2];
          set_sp4_xm(&vp[2 * ic + iv], wt1r, wt1i, wt2r, wt2i, m_Nc);
        }
      }
    }

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::mult_yp(Field& v, const Field& w, const int ieo)
  {
    int idir = 1;

    int Nvc2 = m_Nvc * 2;
    int Nvcd = m_Nvc * m_Nd;

    double bc2 = m_boundary_each_node[idir];

    double       *vp = v.ptr(0);
    const double *wp = w.ptr(0);
    const double *up = m_Ueo.ptr(m_Ndf * m_Nvol2 * (ieo + 2 * idir));

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol2);

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nx2;
      int iyzt = site / m_Nx2;
      int iy   = iyzt % m_Ny;
      int izt  = iyzt / m_Ny;
      int ixzt = ix + m_Nx2 * izt;
      if (iy == 0) {
        int    in = Nvcd * site;
        int    ix1 = Nvc2 * ixzt;
        int    ix2 = ix1 + NVC;
        double vt1[NVC], vt2[NVC];
        set_sp2_yp(vt1, vt2, &wp[in], m_Nc);
        for (int ivc = 0; ivc < NVC; ++ivc) {
          vcp1_yp[ivc + ix1] = bc2 * vt1[ivc];
          vcp1_yp[ivc + ix2] = bc2 * vt2[ivc];
        }
      }
    }

#pragma omp barrier

#pragma omp master
    {
      const int Nv = m_Nvc * 2 * m_Nx2 * m_Nz * m_Nt;
      Communicator::exchange(Nv, vcp2_yp, vcp1_yp, 1, 1, 3);
    }
#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nx2;
      int iyzt = site / m_Nx2;
      int iy   = iyzt % m_Ny;
      int izt  = iyzt / m_Ny;
      int ixzt = ix + m_Nx2 * izt;
      int nei  = ix + m_Nx2 * (iy + 1 + m_Ny * izt);
      int iv   = Nvcd * site;
      int ig   = m_Ndf * site;

      if (iy < m_Ny - 1) {
        int    in = Nvcd * nei;
        double vt1[NVC], vt2[NVC];
        set_sp2_yp(vt1, vt2, &wp[in], m_Nc);
        for (int ic = 0; ic < m_Nc; ++ic) {
          int    ic2  = ic * NVC;
          double wt1r = mult_uv_r(&up[ic2 + ig], vt1, m_Nc);
          double wt1i = mult_uv_i(&up[ic2 + ig], vt1, m_Nc);
          double wt2r = mult_uv_r(&up[ic2 + ig], vt2, m_Nc);
          double wt2i = mult_uv_i(&up[ic2 + ig], vt2, m_Nc);
          set_sp4_yp(&vp[2 * ic + iv], wt1r, wt1i, wt2r, wt2i, m_Nc);
        }
      } else {
        int ix1 = Nvc2 * ixzt;
        int ix2 = ix1 + NVC;
        for (int ic = 0; ic < m_Nc; ++ic) {
          int    ic2  = ic * NVC;
          double wt1r = mult_uv_r(&up[ic2 + ig], &vcp2_yp[ix1], m_Nc);
          double wt1i = mult_uv_i(&up[ic2 + ig], &vcp2_yp[ix1], m_Nc);
          double wt2r = mult_uv_r(&up[ic2 + ig], &vcp2_yp[ix2], m_Nc);
          double wt2i = mult_uv_i(&up[ic2 + ig], &vcp2_yp[ix2], m_Nc);
          set_sp4_yp(&vp[2 * ic + iv], wt1r, wt1i, wt2r, wt2i, m_Nc);
        }
      }
    }

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::mult_ym(Field& v, const Field& w, const int ieo)
  {
    int idir = 1;

    int Nvc2 = m_Nvc * 2;
    int Nvcd = m_Nvc * m_Nd;

    double bc2 = m_boundary_each_node[idir];

    double       *vp = v.ptr(0);
    const double *wp = w.ptr(0);
    const double *up = m_Ueo.ptr(m_Ndf * m_Nvol2 * (1 - ieo + 2 * idir));

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol2);

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nx2;
      int iyzt = site / m_Nx2;
      int iy   = iyzt % m_Ny;
      int izt  = iyzt / m_Ny;
      int ixzt = ix + m_Nx2 * izt;
      if (iy == m_Ny - 1) {
        int in  = Nvcd * site;
        int ig  = m_Ndf * site;
        int ix1 = Nvc2 * ixzt;
        int ix2 = ix1 + NVC;

        double vt1[NVC], vt2[NVC];
        set_sp2_ym(vt1, vt2, &wp[in], m_Nc);
        for (int ic = 0; ic < m_Nc; ++ic) {
          int ic2 = 2 * ic;
          int icr = 2 * ic;
          int ici = 2 * ic + 1;
          vcp1_ym[icr + ix1] = mult_udagv_r(&up[ic2 + ig], vt1, m_Nc);
          vcp1_ym[ici + ix1] = mult_udagv_i(&up[ic2 + ig], vt1, m_Nc);
          vcp1_ym[icr + ix2] = mult_udagv_r(&up[ic2 + ig], vt2, m_Nc);
          vcp1_ym[ici + ix2] = mult_udagv_i(&up[ic2 + ig], vt2, m_Nc);
        }
      }
    }

#pragma omp barrier

#pragma omp master
    {
      const int Nv = m_Nvc * 2 * m_Nx2 * m_Nz * m_Nt;
      Communicator::exchange(Nv, vcp2_ym, vcp1_ym, 1, -1, 4);
    }
#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nx2;
      int iyzt = site / m_Nx2;
      int iy   = iyzt % m_Ny;
      int izt  = iyzt / m_Ny;
      int ixzt = ix + m_Nx2 * izt;
      int nei  = ix + m_Nx2 * (iy - 1 + m_Ny * izt);
      int iv   = Nvcd * site;

      if (iy > 0) {
        int    ig = m_Ndf * nei;
        int    in = Nvcd * nei;
        double vt1[NVC], vt2[NVC];
        set_sp2_ym(vt1, vt2, &wp[in], m_Nc);
        for (int ic = 0; ic < m_Nc; ++ic) {
          int    ic2  = 2 * ic;
          double wt1r = mult_udagv_r(&up[ic2 + ig], vt1, m_Nc);
          double wt1i = mult_udagv_i(&up[ic2 + ig], vt1, m_Nc);
          double wt2r = mult_udagv_r(&up[ic2 + ig], vt2, m_Nc);
          double wt2i = mult_udagv_i(&up[ic2 + ig], vt2, m_Nc);
          set_sp4_ym(&vp[2 * ic + iv], wt1r, wt1i, wt2r, wt2i, m_Nc);
        }
      } else {
        int ix1 = Nvc2 * ixzt;
        int ix2 = ix1 + NVC;
        for (int ic = 0; ic < m_Nc; ++ic) {
          double wt1r = bc2 * vcp2_ym[2 * ic + ix1];
          double wt1i = bc2 * vcp2_ym[2 * ic + 1 + ix1];
          double wt2r = bc2 * vcp2_ym[2 * ic + ix2];
          double wt2i = bc2 * vcp2_ym[2 * ic + 1 + ix2];
          set_sp4_ym(&vp[2 * ic + iv], wt1r, wt1i, wt2r, wt2i, m_Nc);
        }
      }
    }

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::mult_zp(Field& v, const Field& w, const int ieo)
  {
    int idir = 2;

    int Nvc2 = m_Nvc * 2;
    int Nvcd = m_Nvc * m_Nd;

    double bc2 = m_boundary_each_node[idir];

    double       *vp = v.ptr(0);
    const double *wp = w.ptr(0);
    const double *up = m_Ueo.ptr(m_Ndf * m_Nvol2 * (ieo + 2 * idir));

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol2);

    int Nxy = m_Nx2 * m_Ny;

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ixy  = site % Nxy;
      int izt  = site / Nxy;
      int iz   = izt % m_Nz;
      int it   = izt / m_Nz;
      int ixyt = ixy + Nxy * it;
      if (iz == 0) {
        int    in = Nvcd * site;
        int    ix1 = Nvc2 * ixyt;
        int    ix2 = ix1 + NVC;
        double vt1[NVC], vt2[NVC];
        set_sp2_zp(vt1, vt2, &wp[in], m_Nc);
        for (int ivc = 0; ivc < NVC; ++ivc) {
          vcp1_zp[ivc + ix1] = bc2 * vt1[ivc];
          vcp1_zp[ivc + ix2] = bc2 * vt2[ivc];
        }
      }
    }

#pragma omp barrier

#pragma omp master
    {
      int Nv = m_Nvc * 2 * m_Nx2 * m_Ny * m_Nt;
      Communicator::exchange(Nv, vcp2_zp, vcp1_zp, 2, 1, 5);
    }
#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ixy  = site % Nxy;
      int izt  = site / Nxy;
      int iz   = izt % m_Nz;
      int it   = izt / m_Nz;
      int ixyt = ixy + Nxy * it;
      int nei  = ixy + Nxy * (iz + 1 + m_Nz * it);
      int iv   = Nvcd * site;
      int ig   = m_Ndf * site;

      if (iz < m_Nz - 1) {
        int    in = Nvcd * nei;
        double vt1[NVC], vt2[NVC];
        set_sp2_zp(vt1, vt2, &wp[in], m_Nc);
        for (int ic = 0; ic < m_Nc; ++ic) {
          int    ic2  = ic * NVC;
          double wt1r = mult_uv_r(&up[ic2 + ig], vt1, m_Nc);
          double wt1i = mult_uv_i(&up[ic2 + ig], vt1, m_Nc);
          double wt2r = mult_uv_r(&up[ic2 + ig], vt2, m_Nc);
          double wt2i = mult_uv_i(&up[ic2 + ig], vt2, m_Nc);
          set_sp4_zp(&vp[2 * ic + iv], wt1r, wt1i, wt2r, wt2i, m_Nc);
        }
      } else {
        int ix1 = Nvc2 * ixyt;
        int ix2 = ix1 + NVC;
        for (int ic = 0; ic < m_Nc; ++ic) {
          int    ic2  = ic * NVC;
          double wt1r = mult_uv_r(&up[ic2 + ig], &vcp2_zp[ix1], m_Nc);
          double wt1i = mult_uv_i(&up[ic2 + ig], &vcp2_zp[ix1], m_Nc);
          double wt2r = mult_uv_r(&up[ic2 + ig], &vcp2_zp[ix2], m_Nc);
          double wt2i = mult_uv_i(&up[ic2 + ig], &vcp2_zp[ix2], m_Nc);
          set_sp4_zp(&vp[2 * ic + iv], wt1r, wt1i, wt2r, wt2i, m_Nc);
        }
      }
    }

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::mult_zm(Field& v, const Field& w, const int ieo)
  {
    int idir = 2;

    int Nvc2 = m_Nvc * 2;
    int Nvcd = m_Nvc * m_Nd;

    double bc2 = m_boundary_each_node[idir];

    double       *vp = v.ptr(0);
    const double *wp = w.ptr(0);
    const double *up = m_Ueo.ptr(m_Ndf * m_Nvol2 * (1 - ieo + 2 * idir));

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol2);

    int Nxy = m_Nx2 * m_Ny;

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ixy  = site % Nxy;
      int izt  = site / Nxy;
      int iz   = izt % m_Nz;
      int it   = izt / m_Nz;
      int ixyt = ixy + Nxy * it;
      if (iz == m_Nz - 1) {
        int in  = Nvcd * site;
        int ig  = m_Ndf * site;
        int ix1 = Nvc2 * ixyt;
        int ix2 = ix1 + NVC;

        double vt1[NVC], vt2[NVC];
        set_sp2_zm(vt1, vt2, &wp[in], m_Nc);
        for (int ic = 0; ic < m_Nc; ++ic) {
          int ic2 = 2 * ic;
          int icr = 2 * ic;
          int ici = 2 * ic + 1;
          vcp1_zm[icr + ix1] = mult_udagv_r(&up[ic2 + ig], vt1, m_Nc);
          vcp1_zm[ici + ix1] = mult_udagv_i(&up[ic2 + ig], vt1, m_Nc);
          vcp1_zm[icr + ix2] = mult_udagv_r(&up[ic2 + ig], vt2, m_Nc);
          vcp1_zm[ici + ix2] = mult_udagv_i(&up[ic2 + ig], vt2, m_Nc);
        }
      }
    }

#pragma omp barrier

#pragma omp master
    {
      const int Nv = m_Nvc * 2 * m_Nx2 * m_Ny * m_Nt;
      Communicator::exchange(Nv, vcp2_zm, vcp1_zm, 2, -1, 6);
    }
#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ixy  = site % Nxy;
      int izt  = site / Nxy;
      int iz   = izt % m_Nz;
      int it   = izt / m_Nz;
      int ixyt = ixy + Nxy * it;
      int nei  = ixy + Nxy * (iz - 1 + m_Nz * it);
      int iv   = Nvcd * site;

      if (iz > 0) {
        int    ig = m_Ndf * nei;
        int    in = Nvcd * nei;
        double vt1[NVC], vt2[NVC];
        set_sp2_zm(vt1, vt2, &wp[in], m_Nc);
        for (int ic = 0; ic < m_Nc; ++ic) {
          int    ic2  = 2 * ic;
          double wt1r = mult_udagv_r(&up[ic2 + ig], vt1, m_Nc);
          double wt1i = mult_udagv_i(&up[ic2 + ig], vt1, m_Nc);
          double wt2r = mult_udagv_r(&up[ic2 + ig], vt2, m_Nc);
          double wt2i = mult_udagv_i(&up[ic2 + ig], vt2, m_Nc);
          set_sp4_zm(&vp[2 * ic + iv], wt1r, wt1i, wt2r, wt2i, m_Nc);
        }
      } else {
        int ix1 = Nvc2 * ixyt;
        int ix2 = ix1 + NVC;
        for (int ic = 0; ic < m_Nc; ++ic) {
          double wt1r = bc2 * vcp2_zm[2 * ic + ix1];
          double wt1i = bc2 * vcp2_zm[2 * ic + 1 + ix1];
          double wt2r = bc2 * vcp2_zm[2 * ic + ix2];
          double wt2i = bc2 * vcp2_zm[2 * ic + 1 + ix2];
          set_sp4_zm(&vp[2 * ic + iv], wt1r, wt1i, wt2r, wt2i, m_Nc);
        }
      }
    }

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::mult_tp_dirac(Field& v,
                                     const Field& w, const int ieo)
  {
    int idir = 3;

    int Nvc2 = m_Nvc * 2;
    int Nvcd = m_Nvc * m_Nd;

    double bc2 = m_boundary_each_node[idir];

    double       *vp = v.ptr(0);
    const double *wp = w.ptr(0);
    const double *up = m_Ueo.ptr(m_Ndf * m_Nvol2 * (ieo + 2 * idir));

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol2);

    int Nxyz = m_Nx2 * m_Ny * m_Nz;

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ixyz = site % Nxyz;
      int it   = site / Nxyz;
      if (it == 0) {
        int    in = Nvcd * site;
        int    ix1 = Nvc2 * ixyz;
        int    ix2 = ix1 + NVC;
        double vt1[NVC], vt2[NVC];
        set_sp2_tp_dirac(vt1, vt2, &wp[in], m_Nc);
        for (int ivc = 0; ivc < NVC; ++ivc) {
          vcp1_tp[ivc + ix1] = bc2 * vt1[ivc];
          vcp1_tp[ivc + ix2] = bc2 * vt2[ivc];
        }
      }
    }

#pragma omp barrier

#pragma omp master
    {
      int Nv = m_Nvc * 2 * m_Nx2 * m_Ny * m_Nz;
      Communicator::exchange(Nv, vcp2_tp, vcp1_tp, 3, 1, 7);
    }
#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ixyz = site % Nxyz;
      int it   = site / Nxyz;
      int nei  = ixyz + Nxyz * (it + 1);
      int iv   = Nvcd * site;
      int ig   = m_Ndf * site;

      if (it < m_Nt - 1) {
        int    in = Nvcd * nei;
        double vt1[NVC], vt2[NVC];
        set_sp2_tp_dirac(vt1, vt2, &wp[in], m_Nc);
        for (int ic = 0; ic < m_Nc; ++ic) {
          int    ic2  = ic * NVC;
          double wt1r = mult_uv_r(&up[ic2 + ig], vt1, m_Nc);
          double wt1i = mult_uv_i(&up[ic2 + ig], vt1, m_Nc);
          double wt2r = mult_uv_r(&up[ic2 + ig], vt2, m_Nc);
          double wt2i = mult_uv_i(&up[ic2 + ig], vt2, m_Nc);
          set_sp4_tp_dirac(&vp[2 * ic + iv], wt1r, wt1i, wt2r, wt2i, m_Nc);
        }
      } else {
        int ix1 = Nvc2 * ixyz;
        int ix2 = ix1 + NVC;
        for (int ic = 0; ic < m_Nc; ++ic) {
          int    ic2  = ic * NVC;
          double wt1r = mult_uv_r(&up[ic2 + ig], &vcp2_tp[ix1], m_Nc);
          double wt1i = mult_uv_i(&up[ic2 + ig], &vcp2_tp[ix1], m_Nc);
          double wt2r = mult_uv_r(&up[ic2 + ig], &vcp2_tp[ix2], m_Nc);
          double wt2i = mult_uv_i(&up[ic2 + ig], &vcp2_tp[ix2], m_Nc);
          set_sp4_tp_dirac(&vp[2 * ic + iv], wt1r, wt1i, wt2r, wt2i, m_Nc);
        }
      }
    }

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::mult_tm_dirac(Field& v,
                                     const Field& w, const int ieo)
  {
    int idir = 3;

    int Nvc2 = m_Nvc * 2;
    int Nvcd = m_Nvc * m_Nd;

    double bc2 = m_boundary_each_node[idir];

    double       *vp = v.ptr(0);
    const double *wp = w.ptr(0);
    const double *up = m_Ueo.ptr(m_Ndf * m_Nvol2 * (1 - ieo + 2 * idir));

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol2);

    int Nxyz = m_Nx2 * m_Ny * m_Nz;

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ixyz = site % Nxyz;
      int it   = site / Nxyz;
      if (it == m_Nt - 1) {
        int in  = Nvcd * site;
        int ig  = m_Ndf * site;
        int ix1 = Nvc2 * ixyz;
        int ix2 = ix1 + NVC;

        double vt1[NVC], vt2[NVC];
        set_sp2_tm_dirac(vt1, vt2, &wp[in], m_Nc);
        for (int ic = 0; ic < m_Nc; ++ic) {
          int ic2 = 2 * ic;
          int icr = 2 * ic;
          int ici = 2 * ic + 1;
          vcp1_tm[icr + ix1] = mult_udagv_r(&up[ic2 + ig], vt1, m_Nc);
          vcp1_tm[ici + ix1] = mult_udagv_i(&up[ic2 + ig], vt1, m_Nc);
          vcp1_tm[icr + ix2] = mult_udagv_r(&up[ic2 + ig], vt2, m_Nc);
          vcp1_tm[ici + ix2] = mult_udagv_i(&up[ic2 + ig], vt2, m_Nc);
        }
      }
    }

#pragma omp barrier

#pragma omp master
    {
      const int Nv = m_Nvc * 2 * m_Nx2 * m_Ny * m_Nz;
      Communicator::exchange(Nv, vcp2_tm, vcp1_tm, 3, -1, 8);
    }
#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ixyz = site % Nxyz;
      int it   = site / Nxyz;
      int nei  = ixyz + Nxyz * (it - 1);
      int iv   = Nvcd * site;

      if (it > 0) {
        int    ig = m_Ndf * nei;
        int    in = Nvcd * nei;
        double vt1[NVC], vt2[NVC];
        set_sp2_tm_dirac(vt1, vt2, &wp[in], m_Nc);
        for (int ic = 0; ic < m_Nc; ++ic) {
          int    ic2  = 2 * ic;
          double wt1r = mult_udagv_r(&up[ic2 + ig], vt1, m_Nc);
          double wt1i = mult_udagv_i(&up[ic2 + ig], vt1, m_Nc);
          double wt2r = mult_udagv_r(&up[ic2 + ig], vt2, m_Nc);
          double wt2i = mult_udagv_i(&up[ic2 + ig], vt2, m_Nc);
          set_sp4_tm_dirac(&vp[2 * ic + iv], wt1r, wt1i, wt2r, wt2i, m_Nc);
        }
      } else {
        int ix1 = Nvc2 * ixyz;
        int ix2 = ix1 + NVC;
        for (int ic = 0; ic < m_Nc; ++ic) {
          int    icr  = 2 * ic;
          int    ici  = 2 * ic + 1;
          double wt1r = bc2 * vcp2_tm[icr + ix1];
          double wt1i = bc2 * vcp2_tm[ici + ix1];
          double wt2r = bc2 * vcp2_tm[icr + ix2];
          double wt2i = bc2 * vcp2_tm[ici + ix2];
          set_sp4_tm_dirac(&vp[2 * ic + iv], wt1r, wt1i, wt2r, wt2i, m_Nc);
        }
      }
    }

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::mult_tp_chiral(Field& v,
                                      const Field& w, const int ieo)
  {
    int idir = 3;

    int Nvc2 = m_Nvc * 2;
    int Nvcd = m_Nvc * m_Nd;

    double bc2 = m_boundary_each_node[idir];

    double       *vp = v.ptr(0);
    const double *wp = w.ptr(0);
    const double *up = m_Ueo.ptr(m_Ndf * m_Nvol2 * (ieo + 2 * idir));

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol2);

    int Nxyz = m_Nx2 * m_Ny * m_Nz;

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ixyz = site % Nxyz;
      int it   = site / Nxyz;
      if (it == 0) {
        int    in = Nvcd * site;
        int    ix1 = Nvc2 * ixyz;
        int    ix2 = ix1 + NVC;
        double vt1[NVC], vt2[NVC];
        set_sp2_tp_chiral(vt1, vt2, &wp[in], m_Nc);
        for (int ivc = 0; ivc < NVC; ++ivc) {
          vcp1_tp[ivc + ix1] = bc2 * vt1[ivc];
          vcp1_tp[ivc + ix2] = bc2 * vt2[ivc];
        }
      }
    }

#pragma omp barrier

#pragma omp master
    {
      int Nv = m_Nvc * 2 * m_Nx2 * m_Ny * m_Nz;
      Communicator::exchange(Nv, vcp2_tp, vcp1_tp, 3, 1, 7);
    }
#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ixyz = site % Nxyz;
      int it   = site / Nxyz;
      int nei  = ixyz + Nxyz * (it + 1);
      int iv   = Nvcd * site;
      int ig   = m_Ndf * site;

      if (it < m_Nt - 1) {
        int    in = Nvcd * nei;
        double vt1[NVC], vt2[NVC];
        set_sp2_tp_chiral(vt1, vt2, &wp[in], m_Nc);
        for (int ic = 0; ic < m_Nc; ++ic) {
          int    ic2  = ic * NVC;
          double wt1r = mult_uv_r(&up[ic2 + ig], vt1, m_Nc);
          double wt1i = mult_uv_i(&up[ic2 + ig], vt1, m_Nc);
          double wt2r = mult_uv_r(&up[ic2 + ig], vt2, m_Nc);
          double wt2i = mult_uv_i(&up[ic2 + ig], vt2, m_Nc);
          set_sp4_tp_chiral(&vp[2 * ic + iv], wt1r, wt1i, wt2r, wt2i, m_Nc);
        }
      } else {
        int ix1 = Nvc2 * ixyz;
        int ix2 = ix1 + NVC;
        for (int ic = 0; ic < m_Nc; ++ic) {
          int    ic2  = ic * NVC;
          double wt1r = mult_uv_r(&up[ic2 + ig], &vcp2_tp[ix1], m_Nc);
          double wt1i = mult_uv_i(&up[ic2 + ig], &vcp2_tp[ix1], m_Nc);
          double wt2r = mult_uv_r(&up[ic2 + ig], &vcp2_tp[ix2], m_Nc);
          double wt2i = mult_uv_i(&up[ic2 + ig], &vcp2_tp[ix2], m_Nc);
          set_sp4_tp_chiral(&vp[2 * ic + iv], wt1r, wt1i, wt2r, wt2i, m_Nc);
        }
      }
    }

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::mult_tm_chiral(Field& v,
                                      const Field& w, const int ieo)
  {
    int idir = 3;

    int Nvc2 = m_Nvc * 2;
    int Nvcd = m_Nvc * m_Nd;

    double bc2 = m_boundary_each_node[idir];

    double       *vp = v.ptr(0);
    const double *wp = w.ptr(0);
    const double *up = m_Ueo.ptr(m_Ndf * m_Nvol2 * (1 - ieo + 2 * idir));

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol2);

    int Nxyz = m_Nx2 * m_Ny * m_Nz;

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ixyz = site % Nxyz;
      int it   = site / Nxyz;
      if (it == m_Nt - 1) {
        int in  = Nvcd * site;
        int ig  = m_Ndf * site;
        int ix1 = Nvc2 * ixyz;
        int ix2 = ix1 + NVC;

        double vt1[NVC], vt2[NVC];
        set_sp2_tm_chiral(vt1, vt2, &wp[in], m_Nc);
        for (int ic = 0; ic < m_Nc; ++ic) {
          int ic2 = 2 * ic;
          int icr = 2 * ic;
          int ici = 2 * ic + 1;
          vcp1_tm[icr + ix1] = mult_udagv_r(&up[ic2 + ig], vt1, m_Nc);
          vcp1_tm[ici + ix1] = mult_udagv_i(&up[ic2 + ig], vt1, m_Nc);
          vcp1_tm[icr + ix2] = mult_udagv_r(&up[ic2 + ig], vt2, m_Nc);
          vcp1_tm[ici + ix2] = mult_udagv_i(&up[ic2 + ig], vt2, m_Nc);
        }
      }
    }

#pragma omp barrier

#pragma omp master
    {
      const int Nv = m_Nvc * 2 * m_Nx2 * m_Ny * m_Nz;
      Communicator::exchange(Nv, vcp2_tm, vcp1_tm, 3, -1, 8);
    }
#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ixyz = site % Nxyz;
      int it   = site / Nxyz;
      int nei  = ixyz + Nxyz * (it - 1);
      int iv   = Nvcd * site;

      if (it > 0) {
        int    ig = m_Ndf * nei;
        int    in = Nvcd * nei;
        double vt1[NVC], vt2[NVC];
        set_sp2_tm_chiral(vt1, vt2, &wp[in], m_Nc);
        for (int ic = 0; ic < m_Nc; ++ic) {
          int    ic2  = 2 * ic;
          double wt1r = mult_udagv_r(&up[ic2 + ig], vt1, m_Nc);
          double wt1i = mult_udagv_i(&up[ic2 + ig], vt1, m_Nc);
          double wt2r = mult_udagv_r(&up[ic2 + ig], vt2, m_Nc);
          double wt2i = mult_udagv_i(&up[ic2 + ig], vt2, m_Nc);
          set_sp4_tm_chiral(&vp[2 * ic + iv], wt1r, wt1i, wt2r, wt2i, m_Nc);
        }
      } else {
        int ix1 = Nvc2 * ixyz;
        int ix2 = ix1 + NVC;
        for (int ic = 0; ic < m_Nc; ++ic) {
          double wt1r = bc2 * vcp2_tm[2 * ic + ix1];
          double wt1i = bc2 * vcp2_tm[2 * ic + 1 + ix1];
          double wt2r = bc2 * vcp2_tm[2 * ic + ix2];
          double wt2i = bc2 * vcp2_tm[2 * ic + 1 + ix2];
          set_sp4_tm_chiral(&vp[2 * ic + iv], wt1r, wt1i, wt2r, wt2i, m_Nc);
        }
      }
    }

#pragma omp barrier
  }


//====================================================================
  double Fopr_Wilson_eo::flop_count()
  {
    // Counting of floating point operations in giga unit.
    // The following counting explicitly depends on the implementation.
    // It will be recalculated when the code is modified.
    // The present counting is based on rev.1107. [24 Aug 2014 H.Matsufuru]

    const int Nvol = CommonParameters::Nvol();
    const int NPE  = CommonParameters::NPE();

    int flop_site;

    if (m_repr == "Dirac") {
      flop_site = m_Nc * m_Nd * (6 * (4 * m_Nc + 2) + 2 * (4 * m_Nc + 1));
    } else if (m_repr == "Chiral") {
      flop_site = m_Nc * m_Nd * 8 * (4 * m_Nc + 2);
    } else {
      vout.crucial(m_vl, "Error at %s: input repr is undefined.\n",
                   class_name.c_str());
      exit(EXIT_FAILURE);
    }

    double gflop = flop_site * ((Nvol / 2) * (NPE / 1.0e+9));

    return gflop;
  }


//====================================================================
}
//============================================================END=====
