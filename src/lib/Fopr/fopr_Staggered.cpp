/*!
        @file    fopr_Staggered.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
        @version $LastChangedRevision: 2499 $
*/

#include "Fopr/fopr_Staggered.h"

#include "Field/index_lex.h"

#include "Fopr/fopr_thread-inc.h"

#if defined USE_GROUP_SU3
#include "Fopr/Imp/fopr_Wilson_impl_SU3-inc.h"
#elif defined USE_GROUP_SU2
#include "Fopr/Imp/fopr_Wilson_impl_SU2-inc.h"
#elif defined USE_GROUP_SU_N
#include "Fopr/Imp/fopr_Wilson_impl_SU_N-inc.h"
#endif


#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Fopr_Staggered::register_factory();
}
#endif

const std::string Fopr_Staggered::class_name = "Fopr_Staggered";

//====================================================================
void Fopr_Staggered::init(const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();
  vout.general(m_vl, "%s: construction\n", class_name.c_str());
  vout.increase_indent();

  setup();

  set_parameters(params);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());
}


//====================================================================
void Fopr_Staggered::init()
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();
  vout.general(m_vl, "%s: construction (obsolete)\n",
               class_name.c_str());

  setup();

  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());
}


//====================================================================
void Fopr_Staggered::setup()
{
  m_Nc   = CommonParameters::Nc();
  m_Nvc  = 2 * m_Nc;
  m_Ndf  = 2 * m_Nc * m_Nc;
  m_Nvol = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();

  m_Nx = CommonParameters::Nx();
  m_Ny = CommonParameters::Ny();
  m_Nz = CommonParameters::Nz();
  m_Nt = CommonParameters::Nt();

  m_boundary.resize(m_Ndim);

  m_stg_phase.reset(1, m_Nvol, m_Ndim, Element_type::REAL);
  m_parity.reset(1, m_Nvol, 1, Element_type::REAL);
  set_staggered_phase();

  m_U.reset(m_Nvol, m_Ndim);

  m_v1.reset(m_Nvc, m_Nvol, 1);
  m_v2.reset(m_Nvc, m_Nvol, 1);

  int Nvx = m_Nvc * m_Ny * m_Nz * m_Nt;
  vcp1_xp = new double[Nvx];
  vcp2_xp = new double[Nvx];
  vcp1_xm = new double[Nvx];
  vcp2_xm = new double[Nvx];

  int Nvy = m_Nvc * m_Nx * m_Nz * m_Nt;
  vcp1_yp = new double[Nvy];
  vcp2_yp = new double[Nvy];
  vcp1_ym = new double[Nvy];
  vcp2_ym = new double[Nvy];

  int Nvz = m_Nvc * m_Nx * m_Ny * m_Nt;
  vcp1_zp = new double[Nvz];
  vcp2_zp = new double[Nvz];
  vcp1_zm = new double[Nvz];
  vcp2_zm = new double[Nvz];

  int Nvt = m_Nvc * m_Nx * m_Ny * m_Nz;
  vcp1_tp = new double[Nvt];
  vcp2_tp = new double[Nvt];
  vcp1_tm = new double[Nvt];
  vcp2_tm = new double[Nvt];
}


//====================================================================
void Fopr_Staggered::tidyup()
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
void Fopr_Staggered::set_parameters(const Parameters& params)
{
#pragma omp barrier
  int    ith = ThreadManager::get_thread_id();
  string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    if (ith == 0) m_vl = vout.set_verbose_level(vlevel);
  }
#pragma omp barrier

  double           mq;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_double("quark_mass", mq);
  err += params.fetch_int_vector("boundary_condition", bc);

  if (err) {
    vout.crucial(m_vl, "%s: fetch error, input parameter not found.\n",
                 class_name.c_str());
    abort();
  }

  set_parameters(mq, bc);
}


//====================================================================
void Fopr_Staggered::set_parameters(const double mq,
                                    const std::vector<int> bc)
{
  assert(bc.size() == m_Ndim);

  int err = 0;
  // currently no error check is needed
  if (err) {
    vout.crucial(m_vl, "%s: parameter range check failed.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

#pragma omp barrier
  int ith = ThreadManager::get_thread_id();
  if (ith == 0) {
    m_mq = mq;
    for (int mu = 0; mu < m_Ndim; ++mu) {
      m_boundary[mu] = bc[mu];
    }
  }
#pragma omp barrier

  vout.general(m_vl, "%s: input parameters\n", class_name.c_str());
  vout.general(m_vl, "  mq   = %8.4f\n", m_mq);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }
}


//====================================================================
void Fopr_Staggered::get_parameters(Parameters& params) const
{
  params.set_double("quark_mass", m_mq);
  params.set_int_vector("boundary_condition", m_boundary);
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Fopr_Staggered::set_staggered_phase()
{
  int Nt   = CommonParameters::Nt();
  int Nz   = CommonParameters::Nz();
  int Ny   = CommonParameters::Ny();
  int Nx   = CommonParameters::Nx();
  int Nvol = CommonParameters::Nvol();

  Index_lex idx_lex;

  int ipex = Communicator::ipe(0);
  int ipey = Communicator::ipe(1);
  int ipez = Communicator::ipe(2);
  int ipet = Communicator::ipe(3);

  for (int t = 0; t < Nt; ++t) {
    int t2 = t + ipet * Nt;
    for (int z = 0; z < Nz; ++z) {
      int z2 = z + ipez * Nz;
      for (int y = 0; y < Ny; ++y) {
        int y2 = y + ipey * Ny;
        for (int x = 0; x < Nx; ++x) {
          int x2 = x + ipex * Nx;
          int is = idx_lex.site(x, y, z, t);

          m_stg_phase.set(0, is, 0, 1.0);
          m_stg_phase.set(0, is, 1, 1.0);
          m_stg_phase.set(0, is, 2, 1.0);
          m_stg_phase.set(0, is, 3, 1.0);

          m_parity.set(0, is, 0, 1.0);

          if (((x2 + y2 + z2 + t2) % 2) == 1) {
            m_parity.set(0, is, 0, -1.0);
          }

          if ((x2 % 2) == 1) {
            m_stg_phase.set(0, is, 1, -1.0);
          }
          if (((x2 + y2) % 2) == 1) {
            m_stg_phase.set(0, is, 2, -1.0);
          }
          if (((x2 + y2 + z2) % 2) == 1) {
            m_stg_phase.set(0, is, 3, -1.0);
          }

          /*
          // Fortran style
          if (( x2 + 1) % 2 == 1){
            m_stg_phase.set(0, is, 1, -1.0);
          }
          if (((x2 + y2 + 2) % 2) == 1){
            m_stg_phase.set(0, is, 2, -1.0);
          }
          if (((x2 + y2 + z2 + 3) % 2) == 1){
            m_stg_phase.set(0, is, 3, -1.0);
          }
          */
        }
      }
    }
  }
}


//====================================================================
void Fopr_Staggered::mult_staggered_phase(Field& v, int mu)
{
  int Nvol = CommonParameters::Nvol();

  int Nin = v.nin();
  int Nex = v.nex();

  assert(v.nvol() == Nvol);

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = 0; site < Nvol; ++site) {
      for (int in = 0; in < Nin; ++in) {
        double vt = m_stg_phase.cmp(0, site, mu) * v.cmp(in, site, ex);
        v.set(in, site, ex, vt);
      }
    }
  }
}


//====================================================================
void Fopr_Staggered::set_config(Field *U)
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
void Fopr_Staggered::set_config_omp(Field *U)
{
#pragma omp parallel
  {
    set_config_impl(U);
  }
}


//====================================================================
void Fopr_Staggered::set_config_impl(Field *U)
{
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);

#pragma omp barrier

  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int site = is; site < ns; ++site) {
      double ph = m_stg_phase.cmp(0, site, mu);
      for (int cc = 0; cc < m_Ndf; ++cc) {
        double ut = ph * U->cmp(cc, site, mu);
        m_U.set(cc, site, mu, ut);
      }
    }
  }
#pragma omp barrier
}


//====================================================================
void Fopr_Staggered::set_mode(std::string mode)
{
#pragma omp barrier
  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_mode = mode;
#pragma omp barrier
}


//====================================================================
void Fopr_Staggered::mult(Field& v, const Field& w)
{
  if (m_mode == "D") {
    D(v, w);
  } else if (m_mode == "Ddag") {
    Ddag(v, w);
  } else if (m_mode == "DdagD") {
    DdagD(v, w);
  } else if (m_mode == "H") {
    H(v, w);
  } else {
    vout.crucial(m_vl, "%s: mode undeifined.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void Fopr_Staggered::mult_dag(Field& v, const Field& w)
{
  if (m_mode == "D") {
    Ddag(v, w);
  } else if (m_mode == "Ddag") {
    D(v, w);
  } else if (m_mode == "DdagD") {
    DdagD(v, w);
  } else if (m_mode == "H") {
    H(v, w);
  } else {
    vout.crucial(m_vl, "%s: mode undeifined.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void Fopr_Staggered::H(Field& v, const Field& w)
{
  D(v, w);
  mult_gm5(v);
}


//====================================================================
void Fopr_Staggered::DdagD(Field& v, const Field& w)
{
  // H(m_v1, w);
  // H(v, m_v1);
  D(m_v1, w);
  Ddag(v, m_v1);
}


//====================================================================
void Fopr_Staggered::D(Field& v, const Field& w)
{
  assert(w.check_size(m_Nvc, m_Nvol, 1));
  assert(v.check_size(m_Nvc, m_Nvol, 1));

#pragma omp barrier

  m_v2.set(0.0);
#pragma omp barrier

  mult_xp(m_v2, w);
  mult_xm(m_v2, w);

  mult_yp(m_v2, w);
  mult_ym(m_v2, w);

  mult_zp(m_v2, w);
  mult_zm(m_v2, w);

  mult_tp(m_v2, w);
  mult_tm(m_v2, w);

  // hopping normalization
  //  copy(v, m_v2);
  //  aypx(0.5/m_mq, v, w);

  // mass normalization
  copy(v, m_v2);
  scal(v, 0.5);
  axpy(v, m_mq, w);
#pragma omp barrier
}


//====================================================================
void Fopr_Staggered::Ddag(Field& v, const Field& w)
{
  assert(w.check_size(m_Nvc, m_Nvol, 1));
  assert(v.check_size(m_Nvc, m_Nvol, 1));

#pragma omp barrier

  m_v2.set(0.0);
#pragma omp barrier

  mult_xp(m_v2, w);
  mult_xm(m_v2, w);

  mult_yp(m_v2, w);
  mult_ym(m_v2, w);

  mult_zp(m_v2, w);
  mult_zm(m_v2, w);

  mult_tp(m_v2, w);
  mult_tm(m_v2, w);

  // hopping normalization
  // copy(v, m_v2);
  // aypx(-0.5/m_mq, v, w);

  // mass normalization
  copy(v, m_v2);
  scal(v, -0.5);
  axpy(v, m_mq, w);
#pragma omp barrier
}


//====================================================================
void Fopr_Staggered::mult_gm5(Field& v, const Field& w)
{
  int Nin  = w.nin();
  int Nvol = w.nvol();
  int Nex  = w.nex();
  assert(Nvol == m_Nvol);
  assert(v.check_size(Nin, Nvol, Nex));

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);

#pragma omp barrier

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = is; site < ns; ++site) {
      double ph = m_parity.cmp(0, site, 0);
      for (int in = 0; in < Nin; ++in) {
        double vt = ph * w.cmp(in, site, ex);
        v.set(in, site, ex, vt);
      }
    }
  }

#pragma omp barrier
}


//====================================================================
void Fopr_Staggered::mult_gm5(Field& v)
{
  int Nin  = v.nin();
  int Nvol = v.nvol();
  int Nex  = v.nex();
  assert(Nvol == m_Nvol);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);

#pragma omp barrier

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = is; site < ns; ++site) {
      double ph = m_parity.cmp(0, site, 0);
      for (int in = 0; in < Nin; ++in) {
        double vt = ph * v.cmp(in, site, ex);
        v.set(in, site, ex, vt);
      }
    }
  }

#pragma omp barrier
}


//====================================================================
void Fopr_Staggered::mult_xp(Field& v, const Field& w)
{
  int idir = 0;

  double bc2 = 1.0;
  if (Communicator::ipe(idir) == 0) bc2 = double(m_boundary[idir]);

  double       *vp = v.ptr(0);
  const double *wp = w.ptr(0);
  const double *up = m_U.ptr(m_Ndf * m_Nvol * idir);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);

#pragma omp barrier

  for (int site = is; site < ns; ++site) {
    int ix   = site % m_Nx;
    int iyzt = site / m_Nx;
    if (ix == 0) {
      for (int ivc = 0; ivc < m_Nvc; ++ivc) {
        vcp1_xp[ivc + m_Nvc * iyzt] = bc2 * wp[ivc + m_Nvc * site];
      }
    }
  }

#pragma omp barrier
#pragma omp master
  {
    int Nv = m_Nvc * m_Ny * m_Nz * m_Nt;
    Communicator::exchange(Nv, vcp2_xp, vcp1_xp, 0, 1, 1);
  }
#pragma omp barrier

  for (int site = is; site < ns; ++site) {
    int ix   = site % m_Nx;
    int iyzt = site / m_Nx;
    int nei  = ix + 1 + m_Nx * iyzt;

    if (ix < m_Nx - 1) {
      for (int ic = 0; ic < m_Nc; ++ic) {
        int    ic2 = ic * m_Nvc;
        double wtr, wti;
        wtr = mult_uv_r(&up[ic2 + m_Ndf * site], &wp[m_Nvc * nei], m_Nc);
        wti = mult_uv_i(&up[ic2 + m_Ndf * site], &wp[m_Nvc * nei], m_Nc);
        vp[2 * ic + m_Nvc * site]     += wtr;
        vp[2 * ic + 1 + m_Nvc * site] += wti;
      }
    } else {
      for (int ic = 0; ic < m_Nc; ++ic) {
        int    ic2 = ic * m_Nvc;
        double wtr, wti;
        wtr = mult_uv_r(&up[ic2 + m_Ndf * site], &vcp2_xp[m_Nvc * iyzt], m_Nc);
        wti = mult_uv_i(&up[ic2 + m_Ndf * site], &vcp2_xp[m_Nvc * iyzt], m_Nc);
        vp[2 * ic + m_Nvc * site]     += wtr;
        vp[2 * ic + 1 + m_Nvc * site] += wti;
      }
    }
  }

#pragma omp barrier
}


//====================================================================
void Fopr_Staggered::mult_xm(Field& v, const Field& w)
{
  int idir = 0;

  double bc2 = 1.0;
  if (Communicator::ipe(idir) == 0) bc2 = double(m_boundary[idir]);

  double       *vp = v.ptr(0);
  const double *wp = w.ptr(0);
  const double *up = m_U.ptr(m_Ndf * m_Nvol * idir);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);

#pragma omp barrier

  for (int site = is; site < ns; ++site) {
    int ix   = site % m_Nx;
    int iyzt = site / m_Nx;
    if (ix == m_Nx - 1) {
      for (int ic = 0; ic < m_Nc; ++ic) {
        int    ic2 = 2 * ic;
        double wtr, wti;
        wtr = mult_udagv_r(&up[ic2 + m_Ndf * site], &wp[m_Nvc * site], m_Nc);
        wti = mult_udagv_i(&up[ic2 + m_Ndf * site], &wp[m_Nvc * site], m_Nc);
        vcp1_xm[2 * ic + m_Nvc * iyzt]     = wtr;
        vcp1_xm[2 * ic + 1 + m_Nvc * iyzt] = wti;
      }
    }
  }

#pragma omp barrier
#pragma omp master
  {
    int Nv = m_Nvc * m_Ny * m_Nz * m_Nt;
    Communicator::exchange(Nv, vcp2_xm, vcp1_xm, 0, -1, 2);
  }
#pragma omp barrier

  for (int site = is; site < ns; ++site) {
    int ix   = site % m_Nx;
    int iyzt = site / m_Nx;
    int nei  = ix - 1 + m_Nx * iyzt;

    if (ix > 0) {
      for (int ic = 0; ic < m_Nc; ++ic) {
        int    ic2 = 2 * ic;
        double wtr, wti;
        wtr = mult_udagv_r(&up[ic2 + m_Ndf * nei], &wp[m_Nvc * nei], m_Nc);
        wti = mult_udagv_i(&up[ic2 + m_Ndf * nei], &wp[m_Nvc * nei], m_Nc);
        vp[2 * ic + m_Nvc * site]     -= wtr;
        vp[2 * ic + 1 + m_Nvc * site] -= wti;
      }
    } else {
      for (int ivc = 0; ivc < m_Nvc; ++ivc) {
        vp[ivc + m_Nvc * site] -= bc2 * vcp2_xm[ivc + m_Nvc * iyzt];
      }
    }
  }

#pragma omp barrier
}


//====================================================================
void Fopr_Staggered::mult_yp(Field& v, const Field& w)
{
  int idir = 1;

  double bc2 = 1.0;
  if (Communicator::ipe(idir) == 0) bc2 = double(m_boundary[idir]);

  double       *vp = v.ptr(0);
  const double *wp = w.ptr(0);
  const double *up = m_U.ptr(m_Ndf * m_Nvol * idir);

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
      for (int ivc = 0; ivc < m_Nvc; ++ivc) {
        vcp1_yp[ivc + m_Nvc * ixzt] = bc2 * wp[ivc + m_Nvc * site];
      }
    }
  }

#pragma omp barrier
#pragma omp master
  {
    int Nv = m_Nvc * m_Nx * m_Nz * m_Nt;
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

    if (iy < m_Ny - 1) {
      for (int ic = 0; ic < m_Nc; ++ic) {
        int    ic2 = ic * m_Nvc;
        double wtr, wti;
        wtr = mult_uv_r(&up[ic2 + m_Ndf * site], &wp[m_Nvc * nei], m_Nc);
        wti = mult_uv_i(&up[ic2 + m_Ndf * site], &wp[m_Nvc * nei], m_Nc);
        vp[2 * ic + m_Nvc * site]     += wtr;
        vp[2 * ic + 1 + m_Nvc * site] += wti;
      }
    } else {
      for (int ic = 0; ic < m_Nc; ++ic) {
        int    ic2 = ic * m_Nvc;
        double wtr, wti;
        wtr = mult_uv_r(&up[ic2 + m_Ndf * site], &vcp2_yp[m_Nvc * ixzt], m_Nc);
        wti = mult_uv_i(&up[ic2 + m_Ndf * site], &vcp2_yp[m_Nvc * ixzt], m_Nc);
        vp[2 * ic + m_Nvc * site]     += wtr;
        vp[2 * ic + 1 + m_Nvc * site] += wti;
      }
    }
  }

#pragma omp barrier
}


//====================================================================
void Fopr_Staggered::mult_ym(Field& v, const Field& w)
{
  int idir = 1;

  double bc2 = 1.0;
  if (Communicator::ipe(idir) == 0) bc2 = double(m_boundary[idir]);

  double       *vp = v.ptr(0);
  const double *wp = w.ptr(0);
  const double *up = m_U.ptr(m_Ndf * m_Nvol * idir);

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
      for (int ic = 0; ic < m_Nc; ++ic) {
        int    ic2 = 2 * ic;
        double wtr, wti;
        wtr = mult_udagv_r(&up[ic2 + m_Ndf * site], &wp[m_Nvc * site], m_Nc);
        wti = mult_udagv_i(&up[ic2 + m_Ndf * site], &wp[m_Nvc * site], m_Nc);
        vcp1_ym[2 * ic + m_Nvc * ixzt]     = wtr;
        vcp1_ym[2 * ic + 1 + m_Nvc * ixzt] = wti;
      }
    }
  }

#pragma omp barrier
#pragma omp master
  {
    int Nv = m_Nvc * m_Nx * m_Nz * m_Nt;
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

    if (iy > 0) {
      for (int ic = 0; ic < m_Nc; ++ic) {
        int    ic2 = 2 * ic;
        double wtr, wti;
        wtr = mult_udagv_r(&up[ic2 + m_Ndf * nei], &wp[m_Nvc * nei], m_Nc);
        wti = mult_udagv_i(&up[ic2 + m_Ndf * nei], &wp[m_Nvc * nei], m_Nc);
        vp[2 * ic + m_Nvc * site]     -= wtr;
        vp[2 * ic + 1 + m_Nvc * site] -= wti;
      }
    } else {
      for (int ivc = 0; ivc < m_Nvc; ++ivc) {
        vp[ivc + m_Nvc * site] -= bc2 * vcp2_ym[ivc + m_Nvc * ixzt];
      }
    }
  }

#pragma omp barrier
}


//====================================================================
void Fopr_Staggered::mult_zp(Field& v, const Field& w)
{
  int idir = 2;

  double bc2 = 1.0;
  if (Communicator::ipe(idir) == 0) bc2 = double(m_boundary[idir]);

  double       *vp = v.ptr(0);
  const double *wp = w.ptr(0);
  const double *up = m_U.ptr(m_Ndf * m_Nvol * idir);

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
      for (int ivc = 0; ivc < m_Nvc; ++ivc) {
        vcp1_zp[ivc + m_Nvc * ixyt] = bc2 * wp[ivc + m_Nvc * site];
      }
    }
  }

#pragma omp barrier
#pragma omp master
  {
    int Nv = m_Nvc * m_Nx * m_Ny * m_Nt;
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

    if (iz < m_Nz - 1) {
      for (int ic = 0; ic < m_Nc; ++ic) {
        int    ic2 = ic * m_Nvc;
        double wtr, wti;
        wtr = mult_uv_r(&up[ic2 + m_Ndf * site], &wp[m_Nvc * nei], m_Nc);
        wti = mult_uv_i(&up[ic2 + m_Ndf * site], &wp[m_Nvc * nei], m_Nc);
        vp[2 * ic + m_Nvc * site]     += wtr;
        vp[2 * ic + 1 + m_Nvc * site] += wti;
      }
    } else {
      for (int ic = 0; ic < m_Nc; ++ic) {
        int    ic2 = ic * m_Nvc;
        double wtr, wti;
        wtr = mult_uv_r(&up[ic2 + m_Ndf * site], &vcp2_zp[m_Nvc * ixyt], m_Nc);
        wti = mult_uv_i(&up[ic2 + m_Ndf * site], &vcp2_zp[m_Nvc * ixyt], m_Nc);
        vp[2 * ic + m_Nvc * site]     += wtr;
        vp[2 * ic + 1 + m_Nvc * site] += wti;
      }
    }
  }

#pragma omp barrier
}


//====================================================================
void Fopr_Staggered::mult_zm(Field& v, const Field& w)
{
  int idir = 2;

  double bc2 = 1.0;
  if (Communicator::ipe(idir) == 0) bc2 = double(m_boundary[idir]);

  double       *vp = v.ptr(0);
  const double *wp = w.ptr(0);
  const double *up = m_U.ptr(m_Ndf * m_Nvol * idir);

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
      for (int ic = 0; ic < m_Nc; ++ic) {
        int    ic2 = 2 * ic;
        double wtr, wti;
        wtr = mult_udagv_r(&up[ic2 + m_Ndf * site], &wp[m_Nvc * site], m_Nc);
        wti = mult_udagv_i(&up[ic2 + m_Ndf * site], &wp[m_Nvc * site], m_Nc);
        vcp1_zm[2 * ic + m_Nvc * ixyt]     = wtr;
        vcp1_zm[2 * ic + 1 + m_Nvc * ixyt] = wti;
      }
    }
  }

#pragma omp barrier
#pragma omp master
  {
    int Nv = m_Nvc * m_Nx * m_Ny * m_Nt;
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

    if (iz > 0) {
      for (int ic = 0; ic < m_Nc; ++ic) {
        int    ic2 = 2 * ic;
        double wtr, wti;
        wtr = mult_udagv_r(&up[ic2 + m_Ndf * nei], &wp[m_Nvc * nei], m_Nc);
        wti = mult_udagv_i(&up[ic2 + m_Ndf * nei], &wp[m_Nvc * nei], m_Nc);
        vp[2 * ic + m_Nvc * site]     -= wtr;
        vp[2 * ic + 1 + m_Nvc * site] -= wti;
      }
    } else {
      for (int ivc = 0; ivc < m_Nvc; ++ivc) {
        vp[ivc + m_Nvc * site] -= bc2 * vcp2_zm[ivc + m_Nvc * ixyt];
      }
    }
  }

#pragma omp barrier
}


//====================================================================
void Fopr_Staggered::mult_tp(Field& v, const Field& w)
{
  int idir = 3;

  double bc2 = 1.0;
  if (Communicator::ipe(idir) == 0) bc2 = double(m_boundary[idir]);

  double       *vp = v.ptr(0);
  const double *wp = w.ptr(0);
  const double *up = m_U.ptr(m_Ndf * m_Nvol * idir);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);

  int Nxyz = m_Nx * m_Ny * m_Nz;

#pragma omp barrier

  for (int site = is; site < ns; ++site) {
    int ixyz = site % Nxyz;
    int it   = site / Nxyz;
    if (it == 0) {
      for (int ivc = 0; ivc < m_Nvc; ++ivc) {
        vcp1_tp[ivc + m_Nvc * ixyz] = bc2 * wp[ivc + m_Nvc * site];
      }
    }
  }

#pragma omp barrier
#pragma omp master
  {
    int Nv = m_Nvc * m_Nx * m_Ny * m_Nz;
    Communicator::exchange(Nv, vcp2_tp, vcp1_tp, 3, 1, 7);
  }
#pragma omp barrier

  for (int site = is; site < ns; ++site) {
    int ixyz = site % Nxyz;
    int it   = site / Nxyz;
    int nei  = ixyz + Nxyz * (it + 1);

    if (it < m_Nt - 1) {
      for (int ic = 0; ic < m_Nc; ++ic) {
        int    ic2 = ic * m_Nvc;
        double wtr, wti;
        wtr = mult_uv_r(&up[ic2 + m_Ndf * site], &wp[m_Nvc * nei], m_Nc);
        wti = mult_uv_i(&up[ic2 + m_Ndf * site], &wp[m_Nvc * nei], m_Nc);
        vp[2 * ic + m_Nvc * site]     += wtr;
        vp[2 * ic + 1 + m_Nvc * site] += wti;
      }
    } else {
      for (int ic = 0; ic < m_Nc; ++ic) {
        int    ic2 = ic * m_Nvc;
        double wtr, wti;
        wtr = mult_uv_r(&up[ic2 + m_Ndf * site], &vcp2_tp[m_Nvc * ixyz], m_Nc);
        wti = mult_uv_i(&up[ic2 + m_Ndf * site], &vcp2_tp[m_Nvc * ixyz], m_Nc);
        vp[2 * ic + m_Nvc * site]     += wtr;
        vp[2 * ic + 1 + m_Nvc * site] += wti;
      }
    }
  }

#pragma omp barrier
}


//====================================================================
void Fopr_Staggered::mult_tm(Field& v, const Field& w)
{
  int idir = 3;

  double bc2 = 1.0;
  if (Communicator::ipe(idir) == 0) bc2 = double(m_boundary[idir]);

  double       *vp = v.ptr(0);
  const double *wp = w.ptr(0);
  const double *up = m_U.ptr(m_Ndf * m_Nvol * idir);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);

  int Nxyz = m_Nx * m_Ny * m_Nz;

#pragma omp barrier

  for (int site = is; site < ns; ++site) {
    int ixyz = site % Nxyz;
    int it   = site / Nxyz;
    if (it == m_Nt - 1) {
      for (int ic = 0; ic < m_Nc; ++ic) {
        int    ic2 = 2 * ic;
        double wtr, wti;
        wtr = mult_udagv_r(&up[ic2 + m_Ndf * site], &wp[m_Nvc * site], m_Nc);
        wti = mult_udagv_i(&up[ic2 + m_Ndf * site], &wp[m_Nvc * site], m_Nc);
        vcp1_tm[2 * ic + m_Nvc * ixyz]     = wtr;
        vcp1_tm[2 * ic + 1 + m_Nvc * ixyz] = wti;
      }
    }
  }

#pragma omp barrier
#pragma omp master
  {
    int Nv = m_Nvc * m_Nx * m_Ny * m_Nz;
    Communicator::exchange(Nv, vcp2_tm, vcp1_tm, 3, -1, 8);
  }
#pragma omp barrier

  for (int site = is; site < ns; ++site) {
    int ixyz = site % Nxyz;
    int it   = site / Nxyz;
    int nei  = ixyz + Nxyz * (it - 1);

    if (it > 0) {
      for (int ic = 0; ic < m_Nc; ++ic) {
        int    ic2 = 2 * ic;
        double wtr, wti;
        wtr = mult_udagv_r(&up[ic2 + m_Ndf * nei], &wp[m_Nvc * nei], m_Nc);
        wti = mult_udagv_i(&up[ic2 + m_Ndf * nei], &wp[m_Nvc * nei], m_Nc);
        vp[2 * ic + m_Nvc * site]     -= wtr;
        vp[2 * ic + 1 + m_Nvc * site] -= wti;
      }
    } else {
      for (int ivc = 0; ivc < m_Nvc; ++ivc) {
        vp[ivc + m_Nvc * site] -= bc2 * vcp2_tm[ivc + m_Nvc * ixyz];
      }
    }
  }

#pragma omp barrier
}


//====================================================================
double Fopr_Staggered::flop_count(const std::string mode)
{
  // the following flop counting assumes mass normalization.

  int Nvol = m_Nvol;
  int NPE  = CommonParameters::NPE();

  int flop_site = m_Nvc * (3 + 2 * m_Ndim * 2 * m_Nvc);

  double flop_vol = double(Nvol) * double(NPE);

  double gflop = double(flop_site) * flop_vol * 1.0e-9;

  if ((mode == "DdagD") || (mode == "DDdag")) {
    gflop *= 2.0;
  }

  return gflop;
}


//============================================================END=====
