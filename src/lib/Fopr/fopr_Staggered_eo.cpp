/*!
        @file    fopr_Staggered_eo.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate:: 2023-04-04 15:24:44 #$

        @version $LastChangedRevision: 2506 $
*/

#include "fopr_Staggered_eo.h"

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
  bool init = Fopr_Staggered_eo::register_factory();
}
#endif

const std::string Fopr_Staggered_eo::class_name = "Fopr_Staggered_eo";

//====================================================================
void Fopr_Staggered_eo::init(const Parameters& params)
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
void Fopr_Staggered_eo::init()
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();
  vout.general(m_vl, "%s: construction (obsolete)\n",
               class_name.c_str());
  vout.increase_indent();

  setup();

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());
}


//====================================================================
void Fopr_Staggered_eo::setup()
{
  m_Nc    = CommonParameters::Nc();
  m_Nvc   = 2 * m_Nc;
  m_Nvol  = CommonParameters::Nvol();
  m_Nvol2 = m_Nvol / 2;
  m_Ndim  = CommonParameters::Ndim();

  m_boundary.resize(m_Ndim);

  m_Ueo.reset(m_Nvol, m_Ndim);

  m_staggered_phase.reset(1, m_Nvol, m_Ndim);
  set_staggered_phase();

  m_shift_eo = new ShiftField_eo(m_Nvc);

  m_t1.reset(m_Nvc, m_Nvol2, 1);
  m_t2.reset(m_Nvc, m_Nvol2, 1);

  m_v1.reset(m_Nvc, m_Nvol2, 1);
  m_v2.reset(m_Nvc, m_Nvol2, 1);
}


//====================================================================
void Fopr_Staggered_eo::tidyup()
{
  delete m_shift_eo;
}


//====================================================================
void Fopr_Staggered_eo::set_parameters(const Parameters& params)
{
#pragma omp barrier
  int         ith = ThreadManager::get_thread_id();
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    if (ith == 0) m_vl = vout.set_verbose_level(vlevel);
  }
#pragma omp barrier

  //- fetch and check input parameters
  double           mq;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_double("quark_mass", mq);
  err += params.fetch_int_vector("boundary_condition", bc);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(mq, bc);
}


//====================================================================
void Fopr_Staggered_eo::set_parameters(const double mq,
                                       const std::vector<int> bc)
{
  assert(bc.size() == m_Ndim);

  int err = 0;
  // currently no error check is needed.
  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

#pragma omp barrier
  int ith = ThreadManager::get_thread_id();
  if (ith == 0) {
    m_mq       = mq;
    m_boundary = bc;
  }
#pragma omp barrier

  vout.general(m_vl, "%s: input parameters\n", class_name.c_str());
  vout.general(m_vl, "  mq   = %12.8f\n", m_mq);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }
}


//====================================================================
void Fopr_Staggered_eo::get_parameters(Parameters& params) const
{
  params.set_double("quark_mass", m_mq);
  params.set_int_vector("boundary_condition", m_boundary);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Fopr_Staggered_eo::set_staggered_phase()
{
  const int Nt   = CommonParameters::Nt();
  const int Nz   = CommonParameters::Nz();
  const int Ny   = CommonParameters::Ny();
  const int Nx   = CommonParameters::Nx();
  const int Nvol = CommonParameters::Nvol();

  Field     phase_lex(1, Nvol, m_Ndim, Element_type::REAL);
  Index_lex idx_lex;

  const int ipex = Communicator::ipe(0);
  const int ipey = Communicator::ipe(1);
  const int ipez = Communicator::ipe(2);
  const int ipet = Communicator::ipe(3);

  for (int it = 0; it < Nt; ++it) {
    int kt = it + ipet * Nt;
    for (int iz = 0; iz < Nz; ++iz) {
      int kz = iz + ipez * Nz;
      for (int iy = 0; iy < Ny; ++iy) {
        int ky = iy + ipey * Ny;
        for (int ix = 0; ix < Nx; ++ix) {
          int kx = ix + ipex * Nx;

          int is = idx_lex.site(ix, iy, iz, it);

          phase_lex.set(0, is, 0, 1.0);
          phase_lex.set(0, is, 1, 1.0);
          phase_lex.set(0, is, 2, 1.0);
          phase_lex.set(0, is, 3, 1.0);

          if (kx % 2 == 1) phase_lex.set(0, is, 1, -1.0);
          if ((kx + ky) % 2 == 1) phase_lex.set(0, is, 2, -1.0);
          if ((kx + ky + kz) % 2 == 1) phase_lex.set(0, is, 3, -1.0);

          /*
          // the following setting is for check with Fortran
          if (((kx + 1) % 2) == 1) phase_lex.set(0, is, 1, -1.0);
          if (((kx + ky + 2) % 2) == 1) phase_lex.set(0, is, 2, -1.0);
          if (((kx + ky + kz + 3) % 2) == 1) phase_lex.set(0, is, 3, -1.0);
          */
        }
      }
    }
  }

  m_index_eo.convertField(m_staggered_phase, phase_lex);
}


//====================================================================
void Fopr_Staggered_eo::mult_staggered_phase(Field& v,
                                             const int mu, const int ieo)
{
  const int Nvol  = CommonParameters::Nvol();
  const int Nvol2 = Nvol / 2;

  const int Nin = v.nin();
  const int Nex = v.nex();

  assert(v.nvol() == Nvol2);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol2);

#pragma omp barrier
  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = is; site < ns; ++site) {
      int    site2 = m_index_eo.site(site, ieo);
      double ph    = m_staggered_phase.cmp(0, site2, mu);
      for (int in = 0; in < Nin; ++in) {
        double vt = ph * v.cmp(in, site, ex);
        v.set(in, site, ex, vt);
      }
    }
  }
#pragma omp barrier
}


//====================================================================
void Fopr_Staggered_eo::mult(Field& v, const Field& w)
{
  if (m_mode == "D") {
    MeoMoe(v, w);
  } else if (m_mode == "Deo") {
    Meo(v, w, 0, 1);
  } else if (m_mode == "Doe") {
    Meo(v, w, 1, 1);
  } else if (m_mode == "Ddag") {
    MeoMoe(v, w);
  } else if (m_mode == "DdagD") {
    MeoMoe(m_v1, w);
    MeoMoe(v, m_v1);
  } else {
    vout.crucial(m_vl, "Error at %s: mode undeifined.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void Fopr_Staggered_eo::mult_dag(Field& v, const Field& w)
{
  if (m_mode == "D") {
    MeoMoe(v, w);
  } else if (m_mode == "Deo") {
    Meo(v, w, 0, -1);
  } else if (m_mode == "Doe") {
    Meo(v, w, 1, -1);
  } else if (m_mode == "Ddag") {
    MeoMoe(v, w);
  } else if (m_mode == "DdagD") {
    MeoMoe(m_v1, w);
    MeoMoe(v, m_v1);
  } else {
    vout.crucial(m_vl, "Error at %s: mode undeifined.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void Fopr_Staggered_eo::set_mode(std::string mode)
{
#pragma omp barrier
  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_mode = mode;
#pragma omp barrier
}


//====================================================================
void Fopr_Staggered_eo::mult(Field& v, const Field& w,
                             const std::string mode)
{
  if (mode == "D") {
    MeoMoe(v, w);
  } else if (mode == "Deo") {
    Meo(v, w, 0, 1);
  } else if (mode == "Doe") {
    Meo(v, w, 1, 1);
  } else if ((mode == "Dee") || (mode == "Doo")) {
#pragma omp barrier
    copy(v, w);
    scal(v, m_mq);
#pragma omp barrier
  } else {
    vout.crucial(m_vl, "Error at %s: mode undeifined.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void Fopr_Staggered_eo::mult_dag(Field& v, const Field& w,
                                 const std::string mode)
{
  if (mode == "D") {
    MeoMoe(v, w);
  } else if (mode == "Deo") {
    Meo(v, w, 0, -1);
  } else if (mode == "Doe") {
    Meo(v, w, 1, -1);
  } else if ((mode == "Dee") || (mode == "Doo")) {
#pragma omp barrier
    copy(v, w);
    scal(v, m_mq);
#pragma omp barrier
  } else {
    vout.crucial(m_vl, "Error at %s: mode undeifined.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void Fopr_Staggered_eo::set_config(Field *U)
{
  int nth = ThreadManager::get_num_threads();
  int ith = ThreadManager::get_thread_id();

  vout.detailed(m_vl, "%s: set_config is called: num_threads = %d\n",
                class_name.c_str(), nth);

  if (nth > 1) {
    set_config_impl(U);
  } else {
    set_config_omp(U);
  }

  vout.detailed(m_vl, "%s: set_config finished.\n", class_name.c_str());
}


//====================================================================
void Fopr_Staggered_eo::set_config_omp(Field *U)
{
  vout.detailed(m_vl, "  set_config_omp is called.\n");

#pragma omp parallel
  {
    set_config_impl(U);
  }
}


//====================================================================
void Fopr_Staggered_eo::set_config_impl(Field *U)
{
  int Ncc = m_Nc * m_Nc;

  m_index_eo.convertField(m_Ueo, *U);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);
  vout.detailed(m_vl, "  set_config_impl is called: num_threads = %d\n",
                nth);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int site = is; site < ns; ++site) {
      double ph = m_staggered_phase.cmp(0, site, mu);
      for (int cc = 0; cc < Ncc; ++cc) {
        double ur = ph * m_Ueo.cmp_r(cc, site, mu);
        double ui = ph * m_Ueo.cmp_i(cc, site, mu);
        m_Ueo.set_ri(cc, site, mu, ur, ui);
      }
    }
  }
#pragma omp barrier
}


//====================================================================
void Fopr_Staggered_eo::preProp(Field& Be, Field& Bo, const Field& b)
{
  // assuming Be(Nin, Nvol/2, 1), bo(Nin, Nvol/2, 1), b(Nin, Nvol, 1).

  vout.detailed(m_vl, "%s: preProp starts.\n", class_name.c_str());

  m_index_eo.convertField(Be, b, 0);  // even component field
  m_index_eo.convertField(Bo, b, 1);  //  odd component field

  Meo(m_v1, Bo, 0, 1);

  // axpy(m_xe, -1.0, vt);    // hopping normalization
  scal(Be, m_mq);          // mass normalization
  axpy(Be, -1.0, m_v1);    // mass normalization

#pragma omp barrier

  vout.detailed(m_vl, "%s: preProp finished.\n", class_name.c_str());
}


//====================================================================
void Fopr_Staggered_eo::postProp(Field& xq,
                                 const Field& xe, const Field& bo)
{
  vout.detailed(m_vl, "%s: postProp starts.\n", class_name.c_str());

  copy(m_v1, bo);

#pragma omp barrier

  Meo(m_v2, xe, 1, 1);

  axpy(m_v1, -1.0, m_v2);

  scal(m_v1, 1.0 / m_mq);  // needed for mass normalization

#pragma omp barrier

  m_index_eo.reverseField(xq, xe, 0);    // even component
  m_index_eo.reverseField(xq, m_v1, 1);  //  odd component

  vout.detailed(m_vl, "%s: postProp finished.\n", class_name.c_str());
}


//====================================================================
void Fopr_Staggered_eo::MeoMoe(Field& v, const Field& w)
{
  assert(v.check_size(m_Nvc, m_Nvol2, 1));
  assert(w.check_size(m_Nvc, m_Nvol2, 1));

  Meo(m_v2, w, 1, 1);
  Meo(v, m_v2, 0, -1);

  axpy(v, m_mq * m_mq, w); // mass normalization

#pragma omp barrier
}


//====================================================================
void Fopr_Staggered_eo::Meo(Field& v, const Field& w,
                            const int ieo, const int jd)
{
  assert(w.check_size(m_Nvc, m_Nvol2, 1));
  assert(v.check_size(m_Nvc, m_Nvol2, 1));

#pragma omp barrier

  v.set(0.0);

  mult_up(v, w, 0, ieo);
  mult_dn(v, w, 0, ieo);
  mult_up(v, w, 1, ieo);
  mult_dn(v, w, 1, ieo);
  mult_up(v, w, 2, ieo);
  mult_dn(v, w, 2, ieo);
  mult_up(v, w, 3, ieo);
  mult_dn(v, w, 3, ieo);

  double fac = double(jd) * 0.5; // mass normalization
  // double fac = double(jd) * 0.5/m_mq; // hopping normalization
  scal(v, fac);

#pragma omp barrier
}


//====================================================================
void Fopr_Staggered_eo::mult_up(Field& v, const Field& w,
                                const int mu, const int ieo)
{
  m_shift_eo->backward_h(m_t1, w, m_boundary[mu], mu, ieo);

  int          Ndf = m_Nvc * m_Nc;
  double       *vp = v.ptr(0);
  const double *wp = m_t1.ptr(0);
  const double *up = m_Ueo.ptr(Ndf * m_Nvol * mu);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol2);

  for (int site = is; site < ns; ++site) {
    int site2 = m_index_eo.site(site, ieo);
    for (int ic = 0; ic < m_Nc; ++ic) {
      int    ic2 = ic * m_Nvc;
      double wtr, wti;
      wtr = mult_uv_r(&up[ic2 + Ndf * site2], &wp[m_Nvc * site], m_Nc);
      wti = mult_uv_i(&up[ic2 + Ndf * site2], &wp[m_Nvc * site], m_Nc);
      vp[2 * ic + m_Nvc * site]     += wtr;
      vp[2 * ic + 1 + m_Nvc * site] += wti;
    }
  }

#pragma omp barrier
}


//====================================================================
void Fopr_Staggered_eo::mult_dn(Field& v, const Field& w,
                                const int mu, const int ieo)
{
  int          Ndf = m_Nvc * m_Nc;
  double       *vp = m_t1.ptr(0);
  const double *wp = w.ptr(0);
  const double *up = m_Ueo.ptr(Ndf * m_Nvol * mu);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol2);

#pragma omp barrier

  for (int site = is; site < ns; ++site) {
    int site2 = m_index_eo.site(site, 1 - ieo);
    for (int ic = 0; ic < m_Nc; ++ic) {
      int    ic2 = ic * 2;
      double wtr, wti;
      wtr = mult_udagv_r(&up[ic2 + Ndf * site2], &wp[m_Nvc * site], m_Nc);
      wti = mult_udagv_i(&up[ic2 + Ndf * site2], &wp[m_Nvc * site], m_Nc);
      vp[2 * ic + m_Nvc * site]     = wtr;
      vp[2 * ic + 1 + m_Nvc * site] = wti;
    }
  }

#pragma omp barrier

  m_shift_eo->forward_h(m_t2, m_t1, m_boundary[mu], mu, ieo);

  axpy(v, -1.0, m_t2);

#pragma omp barrier
}


//====================================================================
double Fopr_Staggered_eo::flop_count(const std::string mode)
{
  // the following flop counting assumes mass normalization.

  int Nvol2 = m_Nvol2;
  int NPE   = CommonParameters::NPE();

  int flop_site_Meo = m_Nvc * (1 + 2 * m_Ndim * 2 * m_Nvc);

  int flop_site;

  if ((mode == "Deo") || (mode == "Doe")) {
    flop_site = flop_site_Meo;
  } else if (mode == "D") {
    flop_site = m_Nvc * 2 + 2 * flop_site_Meo;
  } else {
    vout.crucial(m_vl, "%s: unsupported mode in flop_count.\n",
                 class_name.c_str());
    flop_site = 0;
  }

  double flop_vol = double(Nvol2) * double(NPE);

  double gflop = double(flop_site) * flop_vol * 1.0e-9;

  return gflop;
}


//============================================================END=====
