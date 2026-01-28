/*!
        @file    fopr_Wilson_impl.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "fopr_Wilson_impl.h"

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
    bool init = Fopr_Wilson::register_factory();
  }
#endif

  const std::string Fopr_Wilson::class_name = "Imp::Fopr_Wilson";

//====================================================================
  void Fopr_Wilson::init(const Parameters& params)
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
  void Fopr_Wilson::init(const std::string repr)
  {
    ThreadManager::assert_single_thread(class_name);

    m_vl = CommonParameters::Vlevel();

    vout.general(m_vl, "%s: construction (obsolete)\n",
                 class_name.c_str());
    vout.increase_indent();

    setup();

    m_repr = repr;
    vout.general(m_vl, "gamma-matrix type is set to %s\n",
                 m_repr.c_str());

    vout.decrease_indent();
    vout.general(m_vl, "%s: construction finished.\n",
                 class_name.c_str());
  }


//====================================================================
  void Fopr_Wilson::setup()
  {
    // check_Nc();
    std::string imple_gauge = imple_Nc();
    vout.general(m_vl, "Gauge group implementation: %s\n",
                 imple_gauge.c_str());

    m_Nc  = CommonParameters::Nc();
    m_Nd  = CommonParameters::Nd();
    m_Nvc = 2 * m_Nc;
    m_Ndf = 2 * m_Nc * m_Nc;

    m_Nx = CommonParameters::Nx();
    m_Ny = CommonParameters::Ny();
    m_Nz = CommonParameters::Nz();
    m_Nt = CommonParameters::Nt();

    m_Nvol = CommonParameters::Nvol();
    m_Ndim = CommonParameters::Ndim();
    m_boundary.resize(m_Ndim);
    m_boundary_each_node.resize(m_Ndim);

    m_U = 0;

    const int Nvx = m_Nvc * 2 * m_Ny * m_Nz * m_Nt;
    vcp1_xp = new double[Nvx];
    vcp2_xp = new double[Nvx];
    vcp1_xm = new double[Nvx];
    vcp2_xm = new double[Nvx];

    const int Nvy = m_Nvc * 2 * m_Nx * m_Nz * m_Nt;
    vcp1_yp = new double[Nvy];
    vcp2_yp = new double[Nvy];
    vcp1_ym = new double[Nvy];
    vcp2_ym = new double[Nvy];

    const int Nvz = m_Nvc * 2 * m_Nx * m_Ny * m_Nt;
    vcp1_zp = new double[Nvz];
    vcp2_zp = new double[Nvz];
    vcp1_zm = new double[Nvz];
    vcp2_zm = new double[Nvz];

    const int Nvt = m_Nvc * 2 * m_Nx * m_Ny * m_Nz;
    vcp1_tp = new double[Nvt];
    vcp2_tp = new double[Nvt];
    vcp1_tm = new double[Nvt];
    vcp2_tm = new double[Nvt];

    m_w1.reset(m_Nvc * m_Nd, m_Nvol, 1);
    m_w2.reset(m_Nvc * m_Nd, m_Nvol, 1);
  }


//====================================================================
  void Fopr_Wilson::tidyup()
  {
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
  void Fopr_Wilson::set_parameters(const Parameters& params)
  {
#pragma omp barrier
    int         ith = ThreadManager::get_thread_id();
    std::string vlevel;
    if (!params.fetch_string("verbose_level", vlevel)) {
      if (ith == 0) m_vl = vout.set_verbose_level(vlevel);
    }
#pragma omp barrier

    //- fetch and check input parameters
    double           kappa;
    std::vector<int> bc;

    int err = 0;
    err += params.fetch_double("hopping_parameter", kappa);
    err += params.fetch_int_vector("boundary_condition", bc);

    if (err) {
      vout.crucial("Error at %s: input parameter not found.\n",
                   class_name.c_str());
      exit(EXIT_FAILURE);
    }

    set_parameters(kappa, bc);
  }


//====================================================================
  void Fopr_Wilson::set_parameters(const double kappa,
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

    vout.general(m_vl, "%s: input parameters\n", class_name.c_str());
    vout.general(m_vl, "  gamma-matrix type = %s\n", m_repr.c_str());
    vout.general(m_vl, "  kappa  = %12.8f\n", m_kappa);
    for (int mu = 0; mu < m_Ndim; ++mu) {
      vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
    }
  }


//====================================================================
  void Fopr_Wilson::get_parameters(Parameters& params) const
  {
    params.set_double("hopping_parameter", m_kappa);
    params.set_int_vector("boundary_condition", m_boundary);
    params.set_string("gamma_matrix_type", m_repr);

    params.set_string("verbose_level", vout.get_verbose_level(m_vl));
  }


//====================================================================
  void Fopr_Wilson::set_config(Field *U)
  {
#pragma omp barrier
    m_U = (Field_G *)U;
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson::set_mode(const std::string mode)
  {
#pragma omp barrier
    int ith = ThreadManager::get_thread_id();
    if (ith == 0) m_mode = mode;
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson::mult(Field& v, const Field& w)
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
    } else {
      vout.crucial(m_vl, "Error at %s: undefined mode: %s\n",
                   class_name.c_str(), m_mode.c_str());
      exit(EXIT_FAILURE);
    }
  }


//====================================================================
  void Fopr_Wilson::mult_dag(Field& v, const Field& w)
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
      H(v, w);
    } else {
      vout.crucial(m_vl, "Error at %s: undefined mode: %s\n",
                   class_name.c_str(), m_mode.c_str());
      exit(EXIT_FAILURE);
    }
  }


//====================================================================
  void Fopr_Wilson::mult(Field& v, const Field& w,
                         const std::string mode)
  {
    if (mode == "D") {
      D(v, w);
    } else if (mode == "Ddag") {
      Ddag(v, w);
    } else if (mode == "DdagD") {
      DdagD(v, w);
    } else if (mode == "DDdag") {
      DDdag(v, w);
    } else if (mode == "H") {
      H(v, w);
    } else {
      vout.crucial(m_vl, "Error at %s: undefined mode: %s\n",
                   class_name.c_str(), mode.c_str());
      exit(EXIT_FAILURE);
    }
  }


//====================================================================
  void Fopr_Wilson::mult_dag(Field& v, const Field& w,
                             const std::string mode)
  {
    if (mode == "D") {
      Ddag(v, w);
    } else if (mode == "Ddag") {
      D(v, w);
    } else if (mode == "DdagD") {
      DdagD(v, w);
    } else if (mode == "DDdag") {
      DDdag(v, w);
    } else if (mode == "H") {
      H(v, w);
    } else {
      vout.crucial(m_vl, "Error at %s: undefined mode: %s\n",
                   class_name.c_str(), mode.c_str());
      exit(EXIT_FAILURE);
    }
  }


//====================================================================
  void Fopr_Wilson::mult_up(const int mu, Field& v, const Field& w)
  {
    if (mu == 0) {
      mult_xp(v, w);
    } else if (mu == 1) {
      mult_yp(v, w);
    } else if (mu == 2) {
      mult_zp(v, w);
    } else if (mu == 3) {
      if (m_repr == "Dirac") {
        mult_tp_dirac(v, w);
      } else {
        mult_tp_chiral(v, w);
      }
    } else {
      vout.crucial(m_vl, "Error at %s::mult_up: illegal mu=%d\n",
                   class_name.c_str(), mu);
      exit(EXIT_FAILURE);
    }
  }


//====================================================================
  void Fopr_Wilson::mult_dn(const int mu, Field& v, const Field& w)
  {
    if (mu == 0) {
      mult_xm(v, w);
    } else if (mu == 1) {
      mult_ym(v, w);
    } else if (mu == 2) {
      mult_zm(v, w);
    } else if (mu == 3) {
      if (m_repr == "Dirac") {
        mult_tm_dirac(v, w);
      } else {
        mult_tm_chiral(v, w);
      }
    } else {
      vout.crucial(m_vl, "Error at %s::mult_dn: illegal mu=%d\n",
                   class_name.c_str(), mu);
      exit(EXIT_FAILURE);
    }
  }


//====================================================================
  void Fopr_Wilson::mult_gm5(Field& v, const Field& w)
  {
    if (m_repr == "Dirac") {
      mult_gm5_dirac(v, w);
    } else if (m_repr == "Chiral") {
      mult_gm5_chiral(v, w);
    }
  }


//====================================================================
  void Fopr_Wilson::D(Field& v, const Field& w)
  {
    if (m_repr == "Dirac") {
      D_ex_dirac(v, 0, w, 0);
      //D_ex_dirac_alt(v, 0, w, 0);
    } else if (m_repr == "Chiral") {
      D_ex_chiral(v, 0, w, 0);
      //D_ex_chiral_alt(v, 0, w, 0);
    }
  }


//====================================================================
  void Fopr_Wilson::D_ex(Field& v, const int ex1,
                         const Field& w, const int ex2)
  {
    if (m_repr == "Dirac") {
      D_ex_dirac(v, ex1, w, ex2);
    } else if (m_repr == "Chiral") {
      D_ex_chiral(v, ex1, w, ex2);
    }
  }


//====================================================================
  void Fopr_Wilson::Ddag(Field& v, const Field& w)
  {
    mult_gm5(v, w);
    D(m_w1, v);
    mult_gm5(v, m_w1);
  }


//====================================================================
  void Fopr_Wilson::DdagD(Field& v, const Field& w)
  {
    D(m_w1, w);
    mult_gm5(v, m_w1);
    D(m_w1, v);
    mult_gm5(v, m_w1);
  }


//====================================================================
  void Fopr_Wilson::DDdag(Field& v, const Field& w)
  {
    mult_gm5(m_w1, w);
    D(v, m_w1);
    mult_gm5(m_w1, v);
    D(v, m_w1);
  }


//====================================================================
  void Fopr_Wilson::H(Field& v, const Field& w)
  {
    D(m_w1, w);
    mult_gm5(v, m_w1);
  }


//====================================================================
  void Fopr_Wilson::D_ex_dirac(Field& v, const int ex1,
                               const Field& w, const int ex2)
  {
    const int    Ninvol = m_Nvc * m_Nd * m_Nvol;
    double       *vp    = v.ptr(Ninvol * ex1);
    const double *wp    = w.ptr(Ninvol * ex2);
    const double *up    = m_U->ptr(0);

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol);

    int Nvc2 = m_Nvc * 2;
    int Nvcd = m_Nvc * m_Nd;

    int Nxy  = m_Nx * m_Ny;
    int Nxyz = m_Nx * m_Ny * m_Nz;

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nx;
      int iyzt = site / m_Nx;
      int iy   = iyzt % m_Ny;
      int izt  = iyzt / m_Ny;
      int iz   = izt % m_Nz;
      int it   = izt / m_Nz;

      int ixy  = ix + m_Nx * iy;
      int ixyz = ixy + Nxy * iz;
      int ixyt = ixy + Nxy * it;
      int ixzt = ix + m_Nx * izt;

      int in = Nvcd * site;

      if (ix == 0) {
        int    ix1 = Nvc2 * iyzt;
        int    ix2 = ix1 + m_Nvc;
        double bc2 = m_boundary_each_node[0];
        double vt1[NVC], vt2[NVC];
        set_sp2_xp(vt1, vt2, &wp[in], m_Nc);
        for (int ivc = 0; ivc < NVC; ++ivc) {
          vcp1_xp[ivc + ix1] = bc2 * vt1[ivc];
          vcp1_xp[ivc + ix2] = bc2 * vt2[ivc];
        }
      }

      if (ix == m_Nx - 1) {
        int    ig = m_Ndf * (site + 0 * m_Nvol);
        int    ix1 = Nvc2 * iyzt;
        int    ix2 = ix1 + NVC;
        double vt1[NVC], vt2[NVC];
        set_sp2_xm(vt1, vt2, &wp[in], m_Nc);
        for (int ic = 0; ic < m_Nc; ++ic) {
          int ic2 = 2 * ic;
          int icr = 2 * ic;
          int ici = 2 * ic + 1;
          vcp1_xm[icr + ix1] = mult_udagv_r(&up[ic2 + ig], vt1, m_Nc);
          vcp1_xm[ici + ix1] = mult_udagv_i(&up[ic2 + ig], vt1, m_Nc);
          vcp1_xm[icr + ix2] = mult_udagv_r(&up[ic2 + ig], vt2, m_Nc);
          vcp1_xm[ici + ix2] = mult_udagv_i(&up[ic2 + ig], vt2, m_Nc);
        }
      }

      if (iy == 0) {
        int    ix1 = Nvc2 * ixzt;
        int    ix2 = ix1 + NVC;
        double bc2 = m_boundary_each_node[1];
        double vt1[NVC], vt2[NVC];
        set_sp2_yp(vt1, vt2, &wp[in], m_Nc);
        for (int ivc = 0; ivc < NVC; ++ivc) {
          vcp1_yp[ivc + ix1] = bc2 * vt1[ivc];
          vcp1_yp[ivc + ix2] = bc2 * vt2[ivc];
        }
      }

      if (iy == m_Ny - 1) {
        int    ig = m_Ndf * (site + 1 * m_Nvol);
        int    ix1 = Nvc2 * ixzt;
        int    ix2 = ix1 + NVC;
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

      if (iz == 0) {
        int    ix1 = Nvc2 * ixyt;
        int    ix2 = ix1 + NVC;
        double bc2 = m_boundary_each_node[2];
        double vt1[NVC], vt2[NVC];
        set_sp2_zp(vt1, vt2, &wp[in], m_Nc);
        for (int ivc = 0; ivc < NVC; ++ivc) {
          vcp1_zp[ivc + ix1] = bc2 * vt1[ivc];
          vcp1_zp[ivc + ix2] = bc2 * vt2[ivc];
        }
      }

      if (iz == m_Nz - 1) {
        int    ig = m_Ndf * (site + 2 * m_Nvol);
        int    ix1 = Nvc2 * ixyt;
        int    ix2 = ix1 + NVC;
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

      if (it == 0) {
        int    ix1 = Nvc2 * ixyz;
        int    ix2 = ix1 + NVC;
        double bc2 = m_boundary_each_node[3];
        double vt1[NVC], vt2[NVC];
        set_sp2_tp_dirac(vt1, vt2, &wp[in], m_Nc);
        for (int ivc = 0; ivc < NVC; ++ivc) {
          vcp1_tp[ivc + ix1] = bc2 * vt1[ivc];
          vcp1_tp[ivc + ix2] = bc2 * vt2[ivc];
        }
      }

      if (it == m_Nt - 1) {
        int    ig = m_Ndf * (site + 3 * m_Nvol);
        int    ix1 = Nvc2 * ixyz;
        int    ix2 = ix1 + NVC;
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
      int Nvx = m_Nvc * 2 * m_Ny * m_Nz * m_Nt;
      Communicator::exchange(Nvx, vcp2_xp, vcp1_xp, 0, 1, 1);
      Communicator::exchange(Nvx, vcp2_xm, vcp1_xm, 0, -1, 2);

      int Nvy = m_Nvc * 2 * m_Nx * m_Nz * m_Nt;
      Communicator::exchange(Nvy, vcp2_yp, vcp1_yp, 1, 1, 3);
      Communicator::exchange(Nvy, vcp2_ym, vcp1_ym, 1, -1, 4);

      int Nvz = m_Nvc * 2 * m_Nx * m_Ny * m_Nt;
      Communicator::exchange(Nvz, vcp2_zp, vcp1_zp, 2, 1, 5);
      Communicator::exchange(Nvz, vcp2_zm, vcp1_zm, 2, -1, 6);

      int Nvt = m_Nvc * 2 * m_Nx * m_Ny * m_Nz;
      Communicator::exchange(Nvt, vcp2_tp, vcp1_tp, 3, 1, 7);
      Communicator::exchange(Nvt, vcp2_tm, vcp1_tm, 3, -1, 8);
    }
#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nx;
      int iyzt = site / m_Nx;
      int iy   = iyzt % m_Ny;
      int izt  = iyzt / m_Ny;
      int iz   = izt % m_Nz;
      int it   = izt / m_Nz;

      int ixy  = site % Nxy;
      int ixyz = ixy + Nxy * iz;
      int ixyt = ixy + Nxy * it;
      int ixzt = ix + m_Nx * izt;

      int iv = Nvcd * site;

      for (int ivcd = 0; ivcd < Nvcd; ++ivcd) {
        vp[ivcd + iv] = 0.0;
      }

      if (ix < m_Nx - 1) {
        int    nei = ix + 1 + m_Nx * iyzt;
        int    in = Nvcd * nei;
        int    ig = m_Ndf * (site + 0 * m_Nvol);
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
        int ix1 = Nvc2 * iyzt;
        int ix2 = ix1 + NVC;
        int ig  = m_Ndf * (site + 0 * m_Nvol);
        for (int ic = 0; ic < m_Nc; ++ic) {
          int    ic2  = ic * NVC;
          double wt1r = mult_uv_r(&up[ic2 + ig], &vcp2_xp[ix1], m_Nc);
          double wt1i = mult_uv_i(&up[ic2 + ig], &vcp2_xp[ix1], m_Nc);
          double wt2r = mult_uv_r(&up[ic2 + ig], &vcp2_xp[ix2], m_Nc);
          double wt2i = mult_uv_i(&up[ic2 + ig], &vcp2_xp[ix2], m_Nc);
          set_sp4_xp(&vp[2 * ic + iv], wt1r, wt1i, wt2r, wt2i, m_Nc);
        }
      }

      if (ix > 0) {
        int    nei = ix - 1 + m_Nx * iyzt;
        int    ig = m_Ndf * (nei + 0 * m_Nvol);
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
        int    ix1 = Nvc2 * iyzt;
        int    ix2 = ix1 + NVC;
        double bc2 = m_boundary_each_node[0];
        for (int ic = 0; ic < m_Nc; ++ic) {
          double wt1r = bc2 * vcp2_xm[2 * ic + ix1];
          double wt1i = bc2 * vcp2_xm[2 * ic + 1 + ix1];
          double wt2r = bc2 * vcp2_xm[2 * ic + ix2];
          double wt2i = bc2 * vcp2_xm[2 * ic + 1 + ix2];
          set_sp4_xm(&vp[2 * ic + iv], wt1r, wt1i, wt2r, wt2i, m_Nc);
        }
      }

      if (iy < m_Ny - 1) {
        int    nei = ix + m_Nx * (iy + 1 + m_Ny * izt);
        int    ig = m_Ndf * (site + 1 * m_Nvol);
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
        int ig  = m_Ndf * (site + 1 * m_Nvol);
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

      if (iy > 0) {
        int    nei = ix + m_Nx * (iy - 1 + m_Ny * izt);
        int    ig = m_Ndf * (nei + 1 * m_Nvol);
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
        int    ix1 = Nvc2 * ixzt;
        int    ix2 = ix1 + NVC;
        double bc2 = m_boundary_each_node[1];
        for (int ic = 0; ic < m_Nc; ++ic) {
          double wt1r = bc2 * vcp2_ym[2 * ic + ix1];
          double wt1i = bc2 * vcp2_ym[2 * ic + 1 + ix1];
          double wt2r = bc2 * vcp2_ym[2 * ic + ix2];
          double wt2i = bc2 * vcp2_ym[2 * ic + 1 + ix2];
          set_sp4_ym(&vp[2 * ic + iv], wt1r, wt1i, wt2r, wt2i, m_Nc);
        }
      }

      if (iz < m_Nz - 1) {
        int    nei = ixy + Nxy * (iz + 1 + m_Nz * it);
        int    ig = m_Ndf * (site + 2 * m_Nvol);
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
        int ig  = m_Ndf * (site + 2 * m_Nvol);
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

      if (iz > 0) {
        int    nei = ixy + Nxy * (iz - 1 + m_Nz * it);
        int    ig = m_Ndf * (nei + 2 * m_Nvol);
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
        int    ix1 = Nvc2 * ixyt;
        int    ix2 = ix1 + NVC;
        double bc2 = m_boundary_each_node[2];
        for (int ic = 0; ic < m_Nc; ++ic) {
          double wt1r = bc2 * vcp2_zm[2 * ic + ix1];
          double wt1i = bc2 * vcp2_zm[2 * ic + 1 + ix1];
          double wt2r = bc2 * vcp2_zm[2 * ic + ix2];
          double wt2i = bc2 * vcp2_zm[2 * ic + 1 + ix2];
          set_sp4_zm(&vp[2 * ic + iv], wt1r, wt1i, wt2r, wt2i, m_Nc);
        }
      }

      if (it < m_Nt - 1) {
        int    nei = ixyz + Nxyz * (it + 1);
        int    ig = m_Ndf * (site + 3 * m_Nvol);
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
        int ig  = m_Ndf * (site + 3 * m_Nvol);
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

      if (it > 0) {
        int    nei = ixyz + Nxyz * (it - 1);
        int    ig = m_Ndf * (nei + 3 * m_Nvol);
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
        int    ix1 = Nvc2 * ixyz;
        int    ix2 = ix1 + NVC;
        double bc2 = m_boundary_each_node[3];
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

      for (int ivcd = 0; ivcd < Nvcd; ++ivcd) {
        vp[ivcd + iv] = -m_kappa * vp[ivcd + iv] + wp[ivcd + iv];
      }
    }

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson::D_ex_chiral(Field& v, const int ex1,
                                const Field& w, const int ex2)
  {
    const int    Ninvol = m_Nvc * m_Nd * m_Nvol;
    double       *vp    = v.ptr(Ninvol * ex1);
    const double *wp    = w.ptr(Ninvol * ex2);
    const double *up    = m_U->ptr(0);

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol);

    int Nvc2 = m_Nvc * 2;
    int Nvcd = m_Nvc * m_Nd;

    int Nxy  = m_Nx * m_Ny;
    int Nxyz = m_Nx * m_Ny * m_Nz;

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nx;
      int iyzt = site / m_Nx;
      int iy   = iyzt % m_Ny;
      int izt  = iyzt / m_Ny;
      int iz   = izt % m_Nz;
      int it   = izt / m_Nz;

      int ixy  = ix + m_Nx * iy;
      int ixyz = ixy + Nxy * iz;
      int ixyt = ixy + Nxy * it;
      int ixzt = ix + m_Nx * izt;

      int in = Nvcd * site;

      if (ix == 0) {
        int    ix1 = Nvc2 * iyzt;
        int    ix2 = ix1 + NVC;
        double bc2 = m_boundary_each_node[0];
        double vt1[NVC], vt2[NVC];
        set_sp2_xp(vt1, vt2, &wp[in], m_Nc);
        for (int ivc = 0; ivc < NVC; ++ivc) {
          vcp1_xp[ivc + ix1] = bc2 * vt1[ivc];
          vcp1_xp[ivc + ix2] = bc2 * vt2[ivc];
        }
      }

      if (ix == m_Nx - 1) {
        int    ig = m_Ndf * (site + 0 * m_Nvol);
        int    ix1 = Nvc2 * iyzt;
        int    ix2 = ix1 + NVC;
        double vt1[NVC], vt2[NVC];
        set_sp2_xm(vt1, vt2, &wp[in], m_Nc);
        for (int ic = 0; ic < m_Nc; ++ic) {
          int ic2 = 2 * ic;
          int icr = 2 * ic;
          int ici = 2 * ic + 1;
          vcp1_xm[icr + ix1] = mult_udagv_r(&up[ic2 + ig], vt1, m_Nc);
          vcp1_xm[ici + ix1] = mult_udagv_i(&up[ic2 + ig], vt1, m_Nc);
          vcp1_xm[icr + ix2] = mult_udagv_r(&up[ic2 + ig], vt2, m_Nc);
          vcp1_xm[ici + ix2] = mult_udagv_i(&up[ic2 + ig], vt2, m_Nc);
        }
      }

      if (iy == 0) {
        int    ix1 = Nvc2 * ixzt;
        int    ix2 = ix1 + NVC;
        double bc2 = m_boundary_each_node[1];
        double vt1[NVC], vt2[NVC];
        set_sp2_yp(vt1, vt2, &wp[in], m_Nc);
        for (int ivc = 0; ivc < NVC; ++ivc) {
          vcp1_yp[ivc + ix1] = bc2 * vt1[ivc];
          vcp1_yp[ivc + ix2] = bc2 * vt2[ivc];
        }
      }

      if (iy == m_Ny - 1) {
        int    ig = m_Ndf * (site + 1 * m_Nvol);
        int    ix1 = Nvc2 * ixzt;
        int    ix2 = ix1 + NVC;
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

      if (iz == 0) {
        int    ix1 = Nvc2 * ixyt;
        int    ix2 = ix1 + NVC;
        double bc2 = m_boundary_each_node[2];
        double vt1[NVC], vt2[NVC];
        set_sp2_zp(vt1, vt2, &wp[in], m_Nc);
        for (int ivc = 0; ivc < NVC; ++ivc) {
          vcp1_zp[ivc + ix1] = bc2 * vt1[ivc];
          vcp1_zp[ivc + ix2] = bc2 * vt2[ivc];
        }
      }

      if (iz == m_Nz - 1) {
        int    ig = m_Ndf * (site + 2 * m_Nvol);
        int    ix1 = Nvc2 * ixyt;
        int    ix2 = ix1 + NVC;
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

      if (it == 0) {
        int    ix1 = Nvc2 * ixyz;
        int    ix2 = ix1 + NVC;
        double bc2 = m_boundary_each_node[3];
        double vt1[NVC], vt2[NVC];
        set_sp2_tp_chiral(vt1, vt2, &wp[in], m_Nc);
        for (int ivc = 0; ivc < NVC; ++ivc) {
          vcp1_tp[ivc + ix1] = bc2 * vt1[ivc];
          vcp1_tp[ivc + ix2] = bc2 * vt2[ivc];
        }
      }

      if (it == m_Nt - 1) {
        int    ig = m_Ndf * (site + 3 * m_Nvol);
        int    ix1 = Nvc2 * ixyz;
        int    ix2 = ix1 + NVC;
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
      int Nvx = m_Nvc * 2 * m_Ny * m_Nz * m_Nt;
      Communicator::exchange(Nvx, vcp2_xp, vcp1_xp, 0, 1, 1);
      Communicator::exchange(Nvx, vcp2_xm, vcp1_xm, 0, -1, 2);

      int Nvy = m_Nvc * 2 * m_Nx * m_Nz * m_Nt;
      Communicator::exchange(Nvy, vcp2_yp, vcp1_yp, 1, 1, 3);
      Communicator::exchange(Nvy, vcp2_ym, vcp1_ym, 1, -1, 4);

      int Nvz = m_Nvc * 2 * m_Nx * m_Ny * m_Nt;
      Communicator::exchange(Nvz, vcp2_zp, vcp1_zp, 2, 1, 5);
      Communicator::exchange(Nvz, vcp2_zm, vcp1_zm, 2, -1, 6);

      int Nvt = m_Nvc * 2 * m_Nx * m_Ny * m_Nz;
      Communicator::exchange(Nvt, vcp2_tp, vcp1_tp, 3, 1, 7);
      Communicator::exchange(Nvt, vcp2_tm, vcp1_tm, 3, -1, 8);
    }
#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nx;
      int iyzt = site / m_Nx;
      int iy   = iyzt % m_Ny;
      int izt  = iyzt / m_Ny;
      int iz   = izt % m_Nz;
      int it   = izt / m_Nz;

      int ixy  = site % Nxy;
      int ixyz = ixy + Nxy * iz;
      int ixyt = ixy + Nxy * it;
      int ixzt = ix + m_Nx * izt;

      int iv = Nvcd * site;

      for (int ivcd = 0; ivcd < Nvcd; ++ivcd) {
        vp[ivcd + iv] = 0.0;
      }

      if (ix < m_Nx - 1) {
        int    nei = ix + 1 + m_Nx * iyzt;
        int    in = Nvcd * nei;
        int    ig = m_Ndf * (site + 0 * m_Nvol);
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
        int ix1 = Nvc2 * iyzt;
        int ix2 = ix1 + NVC;
        int ig  = m_Ndf * (site + 0 * m_Nvol);
        for (int ic = 0; ic < m_Nc; ++ic) {
          int    ic2  = ic * NVC;
          double wt1r = mult_uv_r(&up[ic2 + ig], &vcp2_xp[ix1], m_Nc);
          double wt1i = mult_uv_i(&up[ic2 + ig], &vcp2_xp[ix1], m_Nc);
          double wt2r = mult_uv_r(&up[ic2 + ig], &vcp2_xp[ix2], m_Nc);
          double wt2i = mult_uv_i(&up[ic2 + ig], &vcp2_xp[ix2], m_Nc);
          set_sp4_xp(&vp[2 * ic + iv], wt1r, wt1i, wt2r, wt2i, m_Nc);
        }
      }

      if (ix > 0) {
        int    nei = ix - 1 + m_Nx * iyzt;
        int    ig = m_Ndf * (nei + 0 * m_Nvol);
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
        int    ix1 = Nvc2 * iyzt;
        int    ix2 = ix1 + NVC;
        double bc2 = m_boundary_each_node[0];
        for (int ic = 0; ic < m_Nc; ++ic) {
          double wt1r = bc2 * vcp2_xm[2 * ic + ix1];
          double wt1i = bc2 * vcp2_xm[2 * ic + 1 + ix1];
          double wt2r = bc2 * vcp2_xm[2 * ic + ix2];
          double wt2i = bc2 * vcp2_xm[2 * ic + 1 + ix2];
          set_sp4_xm(&vp[2 * ic + iv], wt1r, wt1i, wt2r, wt2i, m_Nc);
        }
      }

      if (iy < m_Ny - 1) {
        int    nei = ix + m_Nx * (iy + 1 + m_Ny * izt);
        int    ig = m_Ndf * (site + 1 * m_Nvol);
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
        int ig  = m_Ndf * (site + 1 * m_Nvol);
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

      if (iy > 0) {
        int    nei = ix + m_Nx * (iy - 1 + m_Ny * izt);
        int    ig = m_Ndf * (nei + 1 * m_Nvol);
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
        int    ix1 = Nvc2 * ixzt;
        int    ix2 = ix1 + NVC;
        double bc2 = m_boundary_each_node[1];
        for (int ic = 0; ic < m_Nc; ++ic) {
          double wt1r = bc2 * vcp2_ym[2 * ic + ix1];
          double wt1i = bc2 * vcp2_ym[2 * ic + 1 + ix1];
          double wt2r = bc2 * vcp2_ym[2 * ic + ix2];
          double wt2i = bc2 * vcp2_ym[2 * ic + 1 + ix2];
          set_sp4_ym(&vp[2 * ic + iv], wt1r, wt1i, wt2r, wt2i, m_Nc);
        }
      }

      if (iz < m_Nz - 1) {
        int    nei = ixy + Nxy * (iz + 1 + m_Nz * it);
        int    ig = m_Ndf * (site + 2 * m_Nvol);
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
        int ig  = m_Ndf * (site + 2 * m_Nvol);
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

      if (iz > 0) {
        int    nei = ixy + Nxy * (iz - 1 + m_Nz * it);
        int    ig = m_Ndf * (nei + 2 * m_Nvol);
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
        int    ix1 = Nvc2 * ixyt;
        int    ix2 = ix1 + NVC;
        double bc2 = m_boundary_each_node[2];
        for (int ic = 0; ic < m_Nc; ++ic) {
          double wt1r = bc2 * vcp2_zm[2 * ic + ix1];
          double wt1i = bc2 * vcp2_zm[2 * ic + 1 + ix1];
          double wt2r = bc2 * vcp2_zm[2 * ic + ix2];
          double wt2i = bc2 * vcp2_zm[2 * ic + 1 + ix2];
          set_sp4_zm(&vp[2 * ic + iv], wt1r, wt1i, wt2r, wt2i, m_Nc);
        }
      }

      if (it < m_Nt - 1) {
        int    nei = ixyz + Nxyz * (it + 1);
        int    ig = m_Ndf * (site + 3 * m_Nvol);
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
        int ig  = m_Ndf * (site + 3 * m_Nvol);
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

      if (it > 0) {
        int    nei = ixyz + Nxyz * (it - 1);
        int    ig = m_Ndf * (nei + 3 * m_Nvol);
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
        int    ix1 = Nvc2 * ixyz;
        int    ix2 = ix1 + NVC;
        double bc2 = m_boundary_each_node[3];
        for (int ic = 0; ic < m_Nc; ++ic) {
          int    icr  = 2 * ic;
          int    ici  = 2 * ic + 1;
          double wt1r = bc2 * vcp2_tm[icr + ix1];
          double wt1i = bc2 * vcp2_tm[ici + ix1];
          double wt2r = bc2 * vcp2_tm[icr + ix2];
          double wt2i = bc2 * vcp2_tm[ici + ix2];
          set_sp4_tm_chiral(&vp[2 * ic + iv], wt1r, wt1i, wt2r, wt2i, m_Nc);
        }
      }

      for (int ivcd = 0; ivcd < Nvcd; ++ivcd) {
        vp[ivcd + iv] = -m_kappa * vp[ivcd + iv] + wp[ivcd + iv];
      }
    }

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson::D_ex_dirac_alt(Field& w, const int ex1,
                                   const Field& f, const int ex2)
  {
    clear(w);
    mult_xp(w, f);
    mult_xm(w, f);
    mult_yp(w, f);
    mult_ym(w, f);
    mult_zp(w, f);
    mult_zm(w, f);
    mult_tp_dirac(w, f);
    mult_tm_dirac(w, f);
    daypx(w, -m_kappa, f);    // w = -m_kappa * w + f.
  }


//====================================================================
  void Fopr_Wilson::D_ex_chiral_alt(Field& w, const int ex1,
                                    const Field& f, const int ex2)
  {
    clear(w);
    mult_xp(w, f);
    mult_xm(w, f);
    mult_yp(w, f);
    mult_ym(w, f);
    mult_zp(w, f);
    mult_zm(w, f);
    mult_tp_chiral(w, f);
    mult_tm_chiral(w, f);
    daypx(w, -m_kappa, f);    // w = -m_kappa * w + f.
  }


//====================================================================
  void Fopr_Wilson::mult_gm5p(const int mu, Field& v, const Field& w)
  {
    clear(m_w2);
    mult_up(mu, m_w2, w);
    mult_gm5(v, m_w2);
  }


//====================================================================
  void Fopr_Wilson::proj_chiral(Field& w, const int ex1,
                                const Field& v, const int ex2, const int ipm)
  {
    double fpm = 0.0;

    if (ipm == 1) {
      fpm = 1.0;
    } else if (ipm == -1) {
      fpm = -1.0;
    } else {
      vout.crucial(m_vl, "Error at %s: illegal chirality = %d\n", class_name.c_str(), ipm);
      exit(EXIT_FAILURE);
    }

    copy(m_w1, 0, v, ex2);
    mult_gm5(m_w2, m_w1);
    axpy(m_w1, 0, fpm, m_w2, 0);
    copy(w, ex1, m_w1, 0);

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson::clear(Field& w)
  {
    double *wp = w.ptr(0);

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol);

    int Nvcd = m_Nvc * m_Nd;

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      for (int ivcd = 0; ivcd < Nvcd; ++ivcd) {
        wp[ivcd + Nvcd * site] = 0.0;
      }
    }

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson::daypx(Field& v,
                          const double fac, const Field& w)
  {
    double       *vp = v.ptr(0);
    const double *wp = w.ptr(0);

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol);

    int Nvcd = m_Nvc * m_Nd;

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      for (int ivcd = 0; ivcd < Nvcd; ++ivcd) {
        vp[ivcd + Nvcd * site]
          = fac * vp[ivcd + Nvcd * site] + wp[ivcd + Nvcd * site];
      }
    }

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson::mult_gm5_dirac(Field& v, const Field& w)
  {
    double       *vp = v.ptr(0);
    const double *wp = w.ptr(0);

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol);

    int Nvcd = m_Nvc * m_Nd;

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      mult_gamma5_dirac(&vp[Nvcd * site], &wp[Nvcd * site], m_Nc);
    }

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson::mult_gm5_chiral(Field& v, const Field& w)
  {
    double       *vp = v.ptr(0);
    const double *wp = w.ptr(0);

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol);

    int Nvcd = m_Nvc * m_Nd;

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      mult_gamma5_chiral(&vp[Nvcd * site], &wp[Nvcd * site], m_Nc);
    }

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson::mult_xp(Field& v, const Field& w)
  {
    int idir = 0;

    int Nvc2 = m_Nvc * 2;
    int Nvcd = m_Nvc * m_Nd;

    double bc2 = m_boundary_each_node[idir];

    double       *vp = v.ptr(0);
    const double *wp = w.ptr(0);
    const double *up = m_U->ptr(m_Ndf * m_Nvol * idir);

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol);

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nx;
      int iyzt = site / m_Nx;
      if (ix == 0) {
        int    in = Nvcd * site;
        int    ix1 = Nvc2 * iyzt;
        int    ix2 = ix1 + NVC;
        double vt1[NVC], vt2[NVC];
        set_sp2_xp(vt1, vt2, &wp[in], m_Nc);
        for (int ivc = 0; ivc < NVC; ++ivc) {
          vcp1_xp[ivc + ix1] = bc2 * vt1[ivc];
          vcp1_xp[ivc + ix2] = bc2 * vt2[ivc];
        }
      }
    }

#pragma omp barrier

#pragma omp master
    {
      const int Nv = m_Nvc * 2 * m_Ny * m_Nz * m_Nt;
      Communicator::exchange(Nv, vcp2_xp, vcp1_xp, 0, 1, 1);
    }
#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nx;
      int iyzt = site / m_Nx;
      int nei  = ix + 1 + m_Nx * iyzt;
      int iv   = Nvcd * site;
      int ig   = m_Ndf * site;

      if (ix < m_Nx - 1) {
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
        int ix1 = Nvc2 * iyzt;
        int ix2 = ix1 + NVC;
        for (int ic = 0; ic < m_Nc; ++ic) {
          int    ic2  = ic * NVC;
          double wt1r = mult_uv_r(&up[ic2 + ig], &vcp2_xp[ix1], m_Nc);
          double wt1i = mult_uv_i(&up[ic2 + ig], &vcp2_xp[ix1], m_Nc);
          double wt2r = mult_uv_r(&up[ic2 + ig], &vcp2_xp[ix2], m_Nc);
          double wt2i = mult_uv_i(&up[ic2 + ig], &vcp2_xp[ix2], m_Nc);
          set_sp4_xp(&vp[2 * ic + iv], wt1r, wt1i, wt2r, wt2i, m_Nc);
        }
      }
    }

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson::mult_xm(Field& v, const Field& w)
  {
    int idir = 0;

    int Nvc2 = m_Nvc * 2;
    int Nvcd = m_Nvc * m_Nd;

    double bc2 = m_boundary_each_node[idir];

    double       *vp = v.ptr(0);
    const double *wp = w.ptr(0);
    const double *up = m_U->ptr(m_Ndf * m_Nvol * idir);

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol);

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nx;
      int iyzt = site / m_Nx;
      if (ix == m_Nx - 1) {
        int in  = Nvcd * site;
        int ig  = m_Ndf * site;
        int ix1 = Nvc2 * iyzt;
        int ix2 = ix1 + NVC;

        double vt1[NVC], vt2[NVC];
        set_sp2_xm(vt1, vt2, &wp[in], m_Nc);
        for (int ic = 0; ic < m_Nc; ++ic) {
          int ic2 = 2 * ic;
          int icr = 2 * ic;
          int ici = 2 * ic + 1;
          vcp1_xm[icr + ix1] = mult_udagv_r(&up[ic2 + ig], vt1, m_Nc);
          vcp1_xm[ici + ix1] = mult_udagv_i(&up[ic2 + ig], vt1, m_Nc);
          vcp1_xm[icr + ix2] = mult_udagv_r(&up[ic2 + ig], vt2, m_Nc);
          vcp1_xm[ici + ix2] = mult_udagv_i(&up[ic2 + ig], vt2, m_Nc);
        }
      }
    }

#pragma omp barrier

#pragma omp master
    {
      const int Nv = m_Nvc * 2 * m_Ny * m_Nz * m_Nt;
      Communicator::exchange(Nv, vcp2_xm, vcp1_xm, 0, -1, 2);
    }
#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nx;
      int iyzt = site / m_Nx;
      int nei  = ix - 1 + m_Nx * iyzt;
      int iv   = Nvcd * site;

      if (ix > 0) {
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
        int ix1 = Nvc2 * iyzt;
        int ix2 = ix1 + NVC;
        for (int ic = 0; ic < m_Nc; ++ic) {
          double wt1r = bc2 * vcp2_xm[2 * ic + ix1];
          double wt1i = bc2 * vcp2_xm[2 * ic + 1 + ix1];
          double wt2r = bc2 * vcp2_xm[2 * ic + ix2];
          double wt2i = bc2 * vcp2_xm[2 * ic + 1 + ix2];
          set_sp4_xm(&vp[2 * ic + iv], wt1r, wt1i, wt2r, wt2i, m_Nc);
        }
      }
    }

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson::mult_yp(Field& v, const Field& w)
  {
    int idir = 1;

    int Nvc2 = m_Nvc * 2;
    int Nvcd = m_Nvc * m_Nd;

    double bc2 = m_boundary_each_node[idir];

    double       *vp = v.ptr(0);
    const double *wp = w.ptr(0);
    const double *up = m_U->ptr(m_Ndf * m_Nvol * idir);

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol);

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nx;
      int iyzt = site / m_Nx;
      int iy   = iyzt % m_Ny;
      int izt  = iyzt / m_Ny;
      int ixzt = ix + m_Nx * izt;
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
      const int Nv = m_Nvc * 2 * m_Nx * m_Nz * m_Nt;
      Communicator::exchange(Nv, vcp2_yp, vcp1_yp, 1, 1, 3);
    }
#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nx;
      int iyzt = site / m_Nx;
      int iy   = iyzt % m_Ny;
      int izt  = iyzt / m_Ny;
      int ixzt = ix + m_Nx * izt;
      int nei  = ix + m_Nx * (iy + 1 + m_Ny * izt);
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
  void Fopr_Wilson::mult_ym(Field& v, const Field& w)
  {
    int idir = 1;

    int Nvc2 = m_Nvc * 2;
    int Nvcd = m_Nvc * m_Nd;

    double bc2 = m_boundary_each_node[idir];

    double       *vp = v.ptr(0);
    const double *wp = w.ptr(0);
    const double *up = m_U->ptr(m_Ndf * m_Nvol * idir);

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol);

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nx;
      int iyzt = site / m_Nx;
      int iy   = iyzt % m_Ny;
      int izt  = iyzt / m_Ny;
      int ixzt = ix + m_Nx * izt;
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
      const int Nv = m_Nvc * 2 * m_Nx * m_Nz * m_Nt;
      Communicator::exchange(Nv, vcp2_ym, vcp1_ym, 1, -1, 4);
    }
#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nx;
      int iyzt = site / m_Nx;
      int iy   = iyzt % m_Ny;
      int izt  = iyzt / m_Ny;
      int ixzt = ix + m_Nx * izt;
      int nei  = ix + m_Nx * (iy - 1 + m_Ny * izt);
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
  void Fopr_Wilson::mult_zp(Field& v, const Field& w)
  {
    int idir = 2;

    int Nvc2 = m_Nvc * 2;
    int Nvcd = m_Nvc * m_Nd;

    double bc2 = m_boundary_each_node[idir];

    double       *vp = v.ptr(0);
    const double *wp = w.ptr(0);
    const double *up = m_U->ptr(m_Ndf * m_Nvol * idir);

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol);

    int Nxy = m_Nx * m_Ny;

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
      int Nv = m_Nvc * 2 * m_Nx * m_Ny * m_Nt;
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
  void Fopr_Wilson::mult_zm(Field& v, const Field& w)
  {
    int idir = 2;

    int Nvc2 = m_Nvc * 2;
    int Nvcd = m_Nvc * m_Nd;

    double bc2 = m_boundary_each_node[idir];

    double       *vp = v.ptr(0);
    const double *wp = w.ptr(0);
    const double *up = m_U->ptr(m_Ndf * m_Nvol * idir);

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol);

    int Nxy = m_Nx * m_Ny;

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
      const int Nv = m_Nvc * 2 * m_Nx * m_Ny * m_Nt;
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
  void Fopr_Wilson::mult_tp_dirac(Field& v, const Field& w)
  {
    int idir = 3;

    int Nvc2 = m_Nvc * 2;
    int Nvcd = m_Nvc * m_Nd;

    double bc2 = m_boundary_each_node[idir];

    double       *vp = v.ptr(0);
    const double *wp = w.ptr(0);
    const double *up = m_U->ptr(m_Ndf * m_Nvol * idir);

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol);

    int Nxyz = m_Nx * m_Ny * m_Nz;

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
      int Nv = m_Nvc * 2 * m_Nx * m_Ny * m_Nz;
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
  void Fopr_Wilson::mult_tm_dirac(Field& v, const Field& w)
  {
    int idir = 3;

    int Nvc2 = m_Nvc * 2;
    int Nvcd = m_Nvc * m_Nd;

    double bc2 = m_boundary_each_node[idir];

    double       *vp = v.ptr(0);
    const double *wp = w.ptr(0);
    const double *up = m_U->ptr(m_Ndf * m_Nvol * idir);

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol);

    int Nxyz = m_Nx * m_Ny * m_Nz;

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
      const int Nv = m_Nvc * 2 * m_Nx * m_Ny * m_Nz;
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
  void Fopr_Wilson::mult_tp_chiral(Field& v, const Field& w)
  {
    int idir = 3;

    int Nvc2 = m_Nvc * 2;
    int Nvcd = m_Nvc * m_Nd;

    double bc2 = m_boundary_each_node[idir];

    double       *vp = v.ptr(0);
    const double *wp = w.ptr(0);
    const double *up = m_U->ptr(m_Ndf * m_Nvol * idir);

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol);

    int Nxyz = m_Nx * m_Ny * m_Nz;

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
      int Nv = m_Nvc * 2 * m_Nx * m_Ny * m_Nz;
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
  void Fopr_Wilson::mult_tm_chiral(Field& v, const Field& w)
  {
    int idir = 3;

    int Nvc2 = m_Nvc * 2;
    int Nvcd = m_Nvc * m_Nd;

    double bc2 = m_boundary_each_node[idir];

    double       *vp = v.ptr(0);
    const double *wp = w.ptr(0);
    const double *up = m_U->ptr(m_Ndf * m_Nvol * idir);

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol);

    int Nxyz = m_Nx * m_Ny * m_Nz;

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
      const int Nv = m_Nvc * 2 * m_Nx * m_Ny * m_Nz;
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
  double Fopr_Wilson::flop_count()
  {
    // Counting of floating point operations in giga unit.
    // The following counting explicitly depends on the implementation.
    // It will be recalculated when the code is modified.
    // The present counting is based on rev.1107. [24 Aug 2014 H.Matsufuru]

    const int Nvol = CommonParameters::Nvol();
    const int NPE  = CommonParameters::NPE();

    int flop_site;

    if (m_repr == "Dirac") {
      flop_site = m_Nc * m_Nd * (4 + 6 * (4 * m_Nc + 2) + 2 * (4 * m_Nc + 1));
    } else if (m_repr == "Chiral") {
      flop_site = m_Nc * m_Nd * (4 + 8 * (4 * m_Nc + 2));
    } else {
      vout.crucial(m_vl, "Error at %s: input repr is undefined.\n",
                   class_name.c_str());
      exit(EXIT_FAILURE);
    }

    double flop = flop_site * (Nvol * NPE);

    if ((m_mode == "DdagD") || (m_mode == "DDdag")) flop *= 2;

    double gflop = flop * 1.e-9;
    return gflop;
  }


//====================================================================
}
//============================================================END=====
