/*!
      @file    fopr_NonRelativistic.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
      @version $LastChangedRevision: 2499 $
*/

#include "lib/Fopr/fopr_NonRelativistic.h"

#include "lib/ResourceManager/threadManager.h"
#include "lib/Measurements/Gauge/fieldStrength.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Fopr_NonRelativistic::register_factory();
}
#endif

const std::string Fopr_NonRelativistic::class_name
  = "Fopr_NonRelativistic";

//====================================================================
namespace {  // inline functions
  inline void set_threadtask(int& ith, int& nth, int& is, int& ns,
                             const int size)
  {
    nth = ThreadManager_OpenMP::get_num_threads();
    ith = ThreadManager_OpenMP::get_thread_id();
    is  = size * ith / nth;
    ns  = size * (ith + 1) / nth;
  }


  double mult_Gn_r(const double *g, const double *w, int Nc)
  {
    return g[0] * w[0] - g[1] * w[1]
           + g[2] * w[2] - g[3] * w[3]
           + g[4] * w[4] - g[5] * w[5];
  }


  double mult_Gn_i(const double *g, const double *w, int Nc)
  {
    return g[0] * w[1] + g[1] * w[0]
           + g[2] * w[3] + g[3] * w[2]
           + g[4] * w[5] + g[5] * w[4];
  }


  double mult_Gd_r(const double *g, const double *w, int Nc)
  {
    return g[0] * w[0] + g[1] * w[1]
           + g[6] * w[2] + g[7] * w[3]
           + g[12] * w[4] + g[13] * w[5];
  }


  double mult_Gd_i(const double *g, const double *w, int Nc)
  {
    return g[0] * w[1] - g[1] * w[0]
           + g[6] * w[3] - g[7] * w[2]
           + g[12] * w[5] - g[13] * w[4];
  }
} // end of nameless namespace

//====================================================================
void Fopr_NonRelativistic::init(const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();
  vout.general(m_vl, "%s: construction\n", class_name.c_str());
  vout.increase_indent();

  m_Nc    = CommonParameters::Nc();
  m_Nvc   = 2 * m_Nc;
  m_Nd    = CommonParameters::Nd();
  m_Nd2   = m_Nd / 2;
  m_Nx    = CommonParameters::Nx();
  m_Ny    = CommonParameters::Ny();
  m_Nz    = CommonParameters::Nz();
  m_Nt    = CommonParameters::Nt();
  m_Ltime = CommonParameters::Lt();
  m_Nvol  = CommonParameters::Nvol();
  m_Ndim  = CommonParameters::Ndim();
  m_Nspc  = m_Nx * m_Ny * m_Nz;

  m_boundary.resize(m_Ndim);

  m_U.reset(m_Nvol, m_Ndim);

  m_repr = "Dirac";  // only Dirac representation is available.

  // number of correction terms in this implementation
  m_num_correct = 6;

  set_parameters(params);

  m_Fstr.resize(6);
  for (int j = 0; j < 6; ++j) {
    m_Fstr[j].reset(m_Nvol, 1);
  }

  int NinF2 = 2 * m_Nc * m_Nd2;

  m_Xt.resize(m_Nt);
  for (int it = 0; it < m_Nt; ++it) {
    m_Xt[it].reset(NinF2, m_Nspc, 1);
    m_Xt[it].set(0.0);
  }

  m_Zt1.reset(NinF2, m_Nspc, 1);
  m_Zt2.reset(NinF2, m_Nspc, 1);
  m_Yt1.reset(NinF2, m_Nspc, 1);
  m_Yt2.reset(NinF2, m_Nspc, 1);
  m_Yt3.reset(NinF2, m_Nspc, 1);
  m_Wt1.reset(NinF2, m_Nspc, 1);
  m_Wt2.reset(NinF2, m_Nspc, 1);
  m_Wt3.reset(NinF2, m_Nspc, 1);

  int Nbuf_x = NinF2 * m_Nspc / m_Nx;
  buf1_x.resize(Nbuf_x);
  buf2_x.resize(Nbuf_x);

  int Nbuf_y = NinF2 * m_Nspc / m_Ny;
  buf1_y.resize(Nbuf_y);
  buf2_y.resize(Nbuf_y);

  int Nbuf_z = NinF2 * m_Nspc / m_Nz;
  buf1_z.resize(Nbuf_z);
  buf2_z.resize(Nbuf_z);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());
}


//====================================================================
void Fopr_NonRelativistic::tidyup()
{
  // do nothing.
}


//====================================================================
void Fopr_NonRelativistic::set_parameters(const Parameters& params)
{
  const std::string str_vlevel = params.get_string("verbose_level");

  vout.general(m_vl, "%s: parameter setup:\n", class_name.c_str());

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  double              MQ;
  int                 nstab;
  double              u0;
  std::vector<int>    bc;
  std::vector<double> coeff;
  std::string         evolution_type;
  std::string         correction_terms;

  int err = 0;
  err += params.fetch_double("quark_mass", MQ);
  err += params.fetch_int("stabilization_parameter", nstab);
  err += params.fetch_double("mean_field_value", u0);
  err += params.fetch_int_vector("boundary_condition", bc);
  err += params.fetch_double_vector("coefficients", coeff);
  err += params.fetch_string("evolution_type", evolution_type);
  err += params.fetch_string("correction_terms", correction_terms);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(MQ, nstab, u0, coeff, bc,
                 evolution_type, correction_terms);

  vout.general(m_vl, "%s: parameter setup finished.\n",
               class_name.c_str());
}


//====================================================================
void Fopr_NonRelativistic::set_parameters(
  const double MQ,
  const int nstab,
  const double u0,
  const std::vector<double> coeff,
  const std::vector<int> bc,
  const std::string evolution_type,
  const std::string correction_terms)
{
  assert(bc.size() == m_Ndim);

  m_MQ               = MQ;
  m_u0               = u0;
  m_nstab            = nstab;
  m_coeff            = coeff;
  m_boundary         = bc;
  m_evolution_type   = evolution_type;
  m_correction_terms = correction_terms;

  if (m_num_correct != m_coeff.size()) {
    vout.crucial(m_vl, "Error at %s: inconsistent number of coefficients.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  vout.general(m_vl, "%s: parameters:\n", class_name.c_str());
  vout.general(m_vl, "  quark mass = %12.8f\n", m_MQ);
  vout.general(m_vl, "  stabilization parameter = %2d\n", m_nstab);
  vout.general(m_vl, "  mean-field value = %12.8f\n", m_u0);
  vout.general(m_vl, "  evolution_type   = %s\n",
               m_evolution_type.c_str());
  vout.general(m_vl, "  correction_terms = %s\n",
               m_correction_terms.c_str());
  for (int i = 0; i < m_num_correct; ++i) {
    vout.general(m_vl, "  coeff[%d] = %12.8f\n", i, m_coeff[i]);
  }
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }
}


//====================================================================
void Fopr_NonRelativistic::get_parameters(Parameters& params) const
{
  params.set_double("quark_mass", m_MQ);
  params.set_int("stabilization_parameter", m_nstab);
  params.set_double("mean_field_value", m_u0);
  params.set_int_vector("boundary_condition", m_boundary);
  params.set_double_vector("coefficients", m_coeff);
  params.set_string("evolution_type", m_evolution_type);
  params.set_string("correction_terms", m_correction_terms);
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Fopr_NonRelativistic::set_config(Field *U)
{
  copy(m_U, *U);
  scal(m_U, 1.0 / m_u0);  // applying mean-field improvement.

  // Note: the following definition of field strength is not traceless.
  FieldStrength fstr;
  fstr.construct_Fmunu_1x1(m_Fstr[0], 2, 1, m_U);
  fstr.construct_Fmunu_1x1(m_Fstr[1], 0, 2, m_U);
  fstr.construct_Fmunu_1x1(m_Fstr[2], 1, 0, m_U);
  fstr.construct_Fmunu_1x1(m_Fstr[3], 0, 3, m_U);
  fstr.construct_Fmunu_1x1(m_Fstr[4], 1, 3, m_U);
  fstr.construct_Fmunu_1x1(m_Fstr[5], 2, 3, m_U);
}


//====================================================================
void Fopr_NonRelativistic::set_mode(std::string mode)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_mode = mode;

#pragma omp barrier
}


//====================================================================
void Fopr_NonRelativistic::mult(Field& v, const Field& w)
{
  if (m_mode == "Evolve") {
    evolve(v, w);
  } else if (m_mode == "Rotation") {
    rotation(v, w, 1);
  } else {
    vout.crucial(m_vl, "Error at %s: undefined mode = %s.\n",
                 class_name.c_str(), m_mode.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void Fopr_NonRelativistic::mult_dag(Field& v, const Field& w)
{
  if (m_mode == "Rotation") {
    rotation(v, w, -1);
  } else {
    vout.crucial(m_vl, "Error at %s: undefined mode = %s.\n",
                 class_name.c_str(), m_mode.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void Fopr_NonRelativistic::mult(Field& v, const Field& w,
                                const std::string mode)
{
  if (mode == "Evolve") {
    evolve(v, w);
  } else if (mode == "Rotation") {
    rotation(v, w, 1);
  } else {
    vout.crucial(m_vl, "Error at %s: undefined mode = %s.\n",
                 class_name.c_str(), m_mode.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void Fopr_NonRelativistic::mult_dag(Field& v, const Field& w,
                                    const std::string mode)
{
  if (mode == "Rotation") {
    rotation(v, w, -1);
  } else {
    vout.crucial(m_vl, "Error at %s: undefined mode = %s.\n",
                 class_name.c_str(), m_mode.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void Fopr_NonRelativistic::mult_gm5(Field& v, const Field& w)
{
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);

#pragma omp barrier

  for (int site = is; site < ns; ++site) {
    for (int id = 0; id < m_Nd2; ++id) {
      for (int ivc = 0; ivc < m_Nvc; ++ivc) {
        int    ivcd  = ivc + m_Nvc * id;
        int    ivcd2 = ivc + m_Nvc * (2 + id);
        double w1    = w.cmp(ivcd, site, 0);
        double w2    = w.cmp(ivcd2, site, 0);
        v.set(ivcd2, site, 0, w1);
        v.set(ivcd, site, 0, w2);
      }
    }
  }
#pragma omp barrier
}


//====================================================================
void Fopr_NonRelativistic::mult_up(const int mu, Field& v, const Field& w)
{
  vout.crucial(m_vl, "Error at %s: mult_up undefined.\n",
               class_name.c_str());
  exit(EXIT_FAILURE);
}


//====================================================================
void Fopr_NonRelativistic::mult_dn(const int mu, Field& v, const Field& w)
{
  vout.crucial(m_vl, "Error at %s: mult_dn undefined.\n",
               class_name.c_str());
  exit(EXIT_FAILURE);
}


//====================================================================
void Fopr_NonRelativistic::evolve(Field& v, const Field& w)
{
  vout.detailed(m_vl, "%s: evolution of heavy quark field start.\n",
                class_name.c_str());

  int IPEt = Communicator::ipe(3);

  // setting source time slice and source vector
  set_source(m_Yt1, w);

  int it_src   = m_itime_src % m_Nt;
  int ipet_src = m_itime_src / m_Nt;
  if (ipet_src == IPEt) copy(m_Xt[it_src], m_Yt1);
#pragma omp barrier

  // time evolution
  for (int itime = 1; itime < m_Ltime - 1; ++itime) {
    int itime1 = (itime - 1 + m_itime_src) % m_Ltime;
    int it1    = itime1 % m_Nt;
    int itime2 = (itime + m_itime_src) % m_Ltime;
    int it2    = itime2 % m_Nt;

    evolve_impl(m_Xt[it2], m_Xt[it1], itime2);

    if (itime2 / m_Nt == Communicator::ipe(3)) {
      copy(m_Zt1, m_Xt[it2]);
    } else {
      m_Zt1.set(0.0);
    }
#pragma omp barrier

    double xnorm2 = m_Zt1.norm2();
    vout.detailed(m_vl, "  itime = %d  norm2 = %18.14e\n",
                  itime2, xnorm2);
  }

  // set 4D spinor field

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nspc);

  for (int it = 0; it < m_Nt; ++it) {
    for (int ispc = is; ispc < ns; ++ispc) {
      for (int id = 0; id < m_Nd2; ++id) {
        for (int ivc = 0; ivc < m_Nvc; ++ivc) {
          int    ivcd  = ivc + m_Nvc * id;
          int    ivcd2 = ivc + m_Nvc * (2 + id);
          double Xt1   = m_Xt[it].cmp(ivcd, ispc, 0);
          v.set(ivcd, ispc + m_Nspc * it, 0, Xt1);
          v.set(ivcd2, ispc + m_Nspc * it, 0, 0.0);
        }
      }
    }
  }
#pragma omp barrier
}


//====================================================================
void Fopr_NonRelativistic::rotation(Field& v, const Field& w,
                                    const int jd)
{
#pragma omp barrier

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nspc);

  for (int it = 0; it < m_Nt; ++it) {
    for (int ispc = is; ispc < ns; ++ispc) {
      for (int id = 0; id < m_Nd2; ++id) {
        for (int ivc = 0; ivc < m_Nvc; ++ivc) {
          int    ivcd  = ivc + m_Nvc * id;
          int    ivcd2 = ivc + m_Nvc * (2 + id);
          double w1    = w.cmp(ivcd, ispc + m_Nspc * it, 0);
          double w2    = w.cmp(ivcd2, ispc + m_Nspc * it, 0);
          m_Zt1.set(ivcd, ispc, 0, w1);
          m_Zt2.set(ivcd, ispc, 0, w2);
        }
      }
    }
#pragma omp barrier

    int IPEt  = Communicator::ipe(3);
    int itime = it + m_Nt * IPEt;

    // upper components of v:
    copy(m_Xt[0], m_Zt1);
#pragma omp barrier

    add_R2(m_Xt[0], m_Zt2, -jd, itime);

    if (m_correction_terms == "next-leading") {
      add_R3(m_Xt[0], m_Zt1, itime);
      add_R4(m_Xt[0], m_Zt1, itime);
      add_R5(m_Xt[0], m_Zt2, itime);
    }

    for (int ispc = is; ispc < ns; ++ispc) {
      for (int id = 0; id < m_Nd2; ++id) {
        for (int ivc = 0; ivc < m_Nvc; ++ivc) {
          int    ivcd = ivc + m_Nvc * id;
          double v1   = m_Xt[0].cmp(ivcd, ispc, 0);
          v.set(ivcd, ispc + m_Nspc * it, 0, v1);
        }
      }
    }
#pragma omp barrier

    // lower components of v:
    copy(m_Xt[0], m_Zt2);
#pragma omp barrier

    add_R2(m_Xt[0], m_Zt1, jd, itime);

    if (m_correction_terms == "next-leading") {
      add_R3(m_Xt[0], m_Zt2, itime);
      add_R4(m_Xt[0], m_Zt2, itime);
      add_R5(m_Xt[0], m_Zt1, itime);
    }

    for (int ispc = is; ispc < ns; ++ispc) {
      for (int id = 0; id < m_Nd2; ++id) {
        for (int ivc = 0; ivc < m_Nvc; ++ivc) {
          int    ivcd  = ivc + m_Nvc * id;
          int    ivcd2 = ivc + m_Nvc * (2 + id);
          double v2    = m_Xt[0].cmp(ivcd, ispc, 0);
          v.set(ivcd2, ispc + m_Nspc * it, 0, v2);
        }
      }
    }
#pragma omp barrier
  }
}


//====================================================================
void Fopr_NonRelativistic::add_R2(Field& Xt, const Field& Yt,
                                  const int jd, const int itime)
{
  double fac = -1.0 / (2.0 * m_MQ);

  m_Yt3.set(0.0);

  calc_Delta1(m_Yt1, Yt, 0, itime);
  mult_sigma1(m_Yt2, m_Yt1);
  axpy(m_Yt3, fac, m_Yt2);
#pragma omp barrier

  calc_Delta1(m_Yt1, Yt, 1, itime);
  mult_sigma2(m_Yt2, m_Yt1);
  axpy(m_Yt3, fac, m_Yt2);
#pragma omp barrier

  calc_Delta1(m_Yt1, Yt, 2, itime);
  mult_sigma3(m_Yt2, m_Yt1);
  axpy(m_Yt3, fac, m_Yt2);
#pragma omp barrier

  m_Yt3.xI();
  axpy(Xt, double(jd), m_Yt3);
}


//====================================================================
void Fopr_NonRelativistic::add_R3(Field& Xt, const Field& Yt,
                                  const int itime)
{
  double fac = 1.0 / (8.0 * m_MQ * m_MQ);

  calc_Delta2(m_Yt1, Yt, itime);
  axpy(Xt, fac, m_Yt1);
#pragma omp barrier
}


//====================================================================
void Fopr_NonRelativistic::add_R4(Field& Xt, const Field& Yt,
                                  const int itime)
{
  double fac = 1.0 / (8.0 * m_MQ * m_MQ);

  m_Yt3.set(0.0);
#pragma omp barrier

  mult_F(m_Yt1, Yt, 0, itime);
  mult_sigma1(m_Yt2, m_Yt1);
  axpy(m_Yt3, fac, m_Yt2);
#pragma omp barrier

  mult_F(m_Yt1, Yt, 1, itime);
  mult_sigma2(m_Yt2, m_Yt1);
  axpy(m_Yt3, fac, m_Yt2);
#pragma omp barrier

  mult_F(m_Yt1, Yt, 2, itime);
  mult_sigma3(m_Yt2, m_Yt1);
  axpy(m_Yt3, fac, m_Yt2);
#pragma omp barrier

  axpy(Xt, 1.0, m_Yt3);
#pragma omp barrier
}


//====================================================================
void Fopr_NonRelativistic::add_R5(Field& Xt, const Field& Yt,
                                  const int itime)
{
  double fac = 1.0 / (4.0 * m_MQ * m_MQ);

  m_Yt3.set(0.0);
#pragma omp barrier

  mult_F(m_Yt1, Yt, 3, itime);
  mult_sigma1(m_Yt2, m_Yt1);
  axpy(m_Yt3, fac, m_Yt2);
#pragma omp barrier

  mult_F(m_Yt1, Yt, 4, itime);
  mult_sigma2(m_Yt2, m_Yt1);
  axpy(m_Yt3, fac, m_Yt2);
#pragma omp barrier

  mult_F(m_Yt1, Yt, 5, itime);
  mult_sigma3(m_Yt2, m_Yt1);
  axpy(m_Yt3, fac, m_Yt2);
#pragma omp barrier

  axpy(Xt, -1.0, m_Yt3);
#pragma omp barrier
}


//====================================================================
void Fopr_NonRelativistic::evolve_impl(Field& Xt, const Field& Yt,
                                       const int itime)
{
  if (m_evolution_type == "evolve-A") {
    evolve_typeA(Xt, Yt, itime);
  } else if (m_evolution_type == "evolve-B") {
    evolve_typeB(Xt, Yt, itime);
  } else {
    vout.crucial(m_vl, "Error in %s: unsupported evolution_type: %s\n",
                 class_name.c_str(), m_evolution_type.c_str());
    exit(EXIT_FAILURE);
  }

  ThreadManager_OpenMP::sync_barrier_all();
}


//====================================================================
void Fopr_NonRelativistic::evolve_typeA(Field& Xt, const Field& Yt,
                                        const int itime)
{
#pragma omp barrier

  vout.paranoiac(m_vl, "  evolve-A: itime = %d\n", itime);

  int itime1 = (itime - 1 + m_Ltime) % m_Ltime;

  copy(m_Zt1, Yt);
#pragma omp barrier

  evolve_H0(m_Zt2, m_Zt1, itime1);

  evolve_deltaH(m_Zt1, m_Zt2, itime1);

  evolve_U4(m_Zt2, m_Zt1, itime1);

  evolve_deltaH(m_Zt1, m_Zt2, itime);

  evolve_H0(m_Zt2, m_Zt1, itime);

  int ipet = itime / m_Nt;
  int IPEt = Communicator::ipe(3);

  if (ipet == IPEt) copy(Xt, m_Zt2);
#pragma omp barrier
}


//====================================================================
void Fopr_NonRelativistic::evolve_typeB(Field& Xt, const Field& Yt,
                                        const int itime)
{
#pragma omp barrier

  vout.paranoiac(m_vl, "  evolve-B: itime = %d\n", itime);

  int itime1 = (itime - 1 + m_Ltime) % m_Ltime;

  copy(m_Zt1, Yt);
#pragma omp barrier

  if (itime1 != m_itime_src) {
    evolve_deltaH(m_Zt2, m_Zt1, itime1);
    m_Zt1.set(0.0);
    copy(m_Zt1, m_Zt2);
  }
#pragma omp barrier

  evolve_H0(m_Zt2, m_Zt1, itime1);

  evolve_U4(m_Zt1, m_Zt2, itime1);

  evolve_H0(m_Zt2, m_Zt1, itime);

  int ipet = itime / m_Nt;
  int IPEt = Communicator::ipe(3);

  if (ipet == IPEt) copy(Xt, m_Zt2);
#pragma omp barrier
}


//====================================================================
void Fopr_NonRelativistic::evolve_H0(Field& Xt, const Field& Yt,
                                     const int itime)
{
#pragma omp barrier

  double fac = 1.0 / (4.0 * m_MQ * double(m_nstab));

  copy(m_Yt1, Yt);
#pragma omp barrier

  for (int i = 0; i < m_nstab; ++i) {
    calc_Delta2(Xt, m_Yt1, itime);
    aypx(fac, Xt, m_Yt1);
#pragma omp barrier

    if (i < m_nstab - 1) {
      copy(m_Yt1, Xt);
#pragma omp barrier
    }
  }
}


//====================================================================
void Fopr_NonRelativistic::evolve_U4(Field& Xt, const Field& Yt,
                                     const int itime)
{
  int NPEt = Communicator::npe(3);
  int IPEt = Communicator::ipe(3);

  int it   = itime % m_Nt;
  int ipet = itime / m_Nt;

  ThreadManager_OpenMP::sync_barrier_all();

  if (ipet == IPEt) mult_Gd(m_Yt1, Yt, 3, itime);

  ThreadManager_OpenMP::sync_barrier_all();

  if (it < m_Nt - 1) {
    if (ipet == IPEt) {
      copy(Xt, m_Yt1);
#pragma omp barrier
    }
  } else {
#pragma omp master
    {
      int    size = m_Nvc * m_Nd2 * m_Nspc;
      double *yp1 = m_Yt1.ptr(0);
      double *yp2 = m_Yt2.ptr(0);
      Communicator::exchange(size, yp2, yp1, 3, -1, 7);
    }
#pragma omp barrier

    int ipet_up = (ipet + 1) % NPEt;
    if (ipet_up == IPEt) copy(Xt, m_Yt2);
  }

  ThreadManager_OpenMP::sync_barrier_all();
}


//====================================================================
void Fopr_NonRelativistic::evolve_deltaH(Field& Xt, const Field& Yt,
                                         const int itime)
{
  // Be sure that the number of correction terms implemented here
  // coincides with the variable m_num_correct (currently 6).

#pragma omp barrier

  Xt.set(0.0);
  add_deltaH1(Xt, Yt, itime);

  if (m_correction_terms == "next-leading") {
    add_deltaH2(Xt, Yt, itime);
    add_deltaH3(Xt, Yt, itime);
    add_deltaH5(Xt, Yt, itime);
    add_deltaH46(Xt, Yt, itime);
  }

  double fac;
  if (m_evolution_type == "evolve-A") {
    fac = -0.5;
  } else if (m_evolution_type == "evolve-B") {
    fac = -1.0;
  } else {
    vout.crucial(m_vl, "Error in %s: unsupported evolution_type: %s\n",
                 class_name.c_str(), m_evolution_type.c_str());
    exit(EXIT_FAILURE);
  }

  aypx(fac, Xt, Yt);
#pragma omp barrier
}


//====================================================================
void Fopr_NonRelativistic::add_deltaH1(Field& Xt, const Field& Yt,
                                       const int itime)
{
  double fac = -m_coeff[0] / (2.0 * m_MQ);

  mult_F(m_Yt1, Yt, 0, itime);
  mult_sigma1(m_Yt2, m_Yt1);
  axpy(Xt, fac, m_Yt2);
#pragma omp barrier

  mult_F(m_Yt1, Yt, 1, itime);
  mult_sigma2(m_Yt2, m_Yt1);
  axpy(Xt, fac, m_Yt2);
#pragma omp barrier

  mult_F(m_Yt1, Yt, 2, itime);
  mult_sigma3(m_Yt2, m_Yt1);
  axpy(Xt, fac, m_Yt2);
#pragma omp barrier
}


//====================================================================
void Fopr_NonRelativistic::add_deltaH2(Field& Xt, const Field& Yt,
                                       const int itime)
{
  double fac = m_coeff[1] / (8.0 * m_MQ * m_MQ);

  m_Yt3.set(0.0);

  mult_F(m_Yt1, Yt, 3, itime);
  calc_Delta1(m_Yt2, m_Yt1, 0, itime);
  axpy(m_Yt3, 1.0, m_Yt2);
#pragma omp barrier

  mult_F(m_Yt1, Yt, 4, itime);
  calc_Delta1(m_Yt2, m_Yt1, 1, itime);
  axpy(m_Yt3, 1.0, m_Yt2);
#pragma omp barrier

  mult_F(m_Yt1, Yt, 5, itime);
  calc_Delta1(m_Yt2, m_Yt1, 2, itime);
  axpy(m_Yt3, 1.0, m_Yt2);
#pragma omp barrier

  calc_Delta1(m_Yt1, Yt, 0, itime);
  mult_F(m_Yt2, m_Yt1, 3, itime);
  axpy(m_Yt3, -1.0, m_Yt2);
#pragma omp barrier

  calc_Delta1(m_Yt1, Yt, 1, itime);
  mult_F(m_Yt2, m_Yt1, 4, itime);
  axpy(m_Yt3, -1.0, m_Yt2);
#pragma omp barrier

  calc_Delta1(m_Yt1, Yt, 2, itime);
  mult_F(m_Yt2, m_Yt1, 5, itime);
  axpy(m_Yt3, -1.0, m_Yt2);
#pragma omp barrier

  m_Yt3.xI();
  axpy(Xt, fac, m_Yt3);
#pragma omp barrier
}


//====================================================================
void Fopr_NonRelativistic::add_deltaH3(Field& Xt, const Field& Yt,
                                       const int itime)
{
  double fac = -m_coeff[2] / (8.0 * m_MQ * m_MQ);

  // x-component
  m_Yt3.set(0.0);
#pragma omp barrier

  mult_F(m_Yt1, Yt, 5, itime);          // E_z
  calc_Delta1(m_Yt2, m_Yt1, 1, itime);  // D_y
  axpy(m_Yt3, 1.0, m_Yt2);              // Yt3 += D_y * E_z
#pragma omp barrier

  mult_F(m_Yt1, Yt, 4, itime);          // E_y
  calc_Delta1(m_Yt2, m_Yt1, 2, itime);  // D_z
  axpy(m_Yt3, -1.0, m_Yt2);             // Yt3 -= D_z * E_y
#pragma omp barrier

  calc_Delta1(m_Yt1, Yt, 2, itime);     // D_z
  mult_F(m_Yt2, m_Yt1, 4, itime);       // E_y
  axpy(m_Yt3, -1.0, m_Yt2);             // Yt3 += E_y * D_z
#pragma omp barrier

  calc_Delta1(m_Yt1, Yt, 1, itime);     // D_y
  mult_F(m_Yt2, m_Yt1, 5, itime);       // E_z
  axpy(m_Yt3, 1.0, m_Yt2);              // Yt3 -= E_z * D_y
#pragma omp barrier

  mult_sigma1(m_Yt1, m_Yt3);
  axpy(Xt, fac, m_Yt1);
#pragma omp barrier

  // y-component
  m_Yt3.set(0.0);
#pragma omp barrier

  mult_F(m_Yt1, Yt, 3, itime);          // E_x
  calc_Delta1(m_Yt2, m_Yt1, 2, itime);  // D_z
  axpy(m_Yt3, 1.0, m_Yt2);              // Yt3 += D_z * E_x
#pragma omp barrier

  mult_F(m_Yt1, Yt, 5, itime);          // E_z
  calc_Delta1(m_Yt2, m_Yt1, 0, itime);  // D_x
  axpy(m_Yt3, -1.0, m_Yt2);             // Yt3 -= D_x * E_z
#pragma omp barrier

  calc_Delta1(m_Yt1, Yt, 0, itime);     // D_x
  mult_F(m_Yt2, m_Yt1, 5, itime);       // E_z
  axpy(m_Yt3, -1.0, m_Yt2);             // Yt3 += E_z * D_x
#pragma omp barrier

  calc_Delta1(m_Yt1, Yt, 2, itime);     // D_z
  mult_F(m_Yt2, m_Yt1, 3, itime);       // E_x
  axpy(m_Yt3, 1.0, m_Yt2);              // Yt3 -= E_x * D_z
#pragma omp barrier

  mult_sigma2(m_Yt1, m_Yt3);
  axpy(Xt, fac, m_Yt1);
#pragma omp barrier

  // z-component
  m_Yt3.set(0.0);
#pragma omp barrier

  mult_F(m_Yt1, Yt, 4, itime);          // E_y
  calc_Delta1(m_Yt2, m_Yt1, 0, itime);  // D_x
  axpy(m_Yt3, 1.0, m_Yt2);              // Yt3 += D_x * E_y
#pragma omp barrier

  mult_F(m_Yt1, Yt, 3, itime);          // E_x
  calc_Delta1(m_Yt2, m_Yt1, 1, itime);  // D_y
  axpy(m_Yt3, -1.0, m_Yt2);             // Yt3 -= D_y * E_x
#pragma omp barrier

  calc_Delta1(m_Yt1, Yt, 1, itime);     // D_y
  mult_F(m_Yt2, m_Yt1, 3, itime);       // E_x
  axpy(m_Yt3, -1.0, m_Yt2);             // Yt3 += E_x * D_y
#pragma omp barrier

  calc_Delta1(m_Yt1, Yt, 0, itime);     // D_x
  mult_F(m_Yt2, m_Yt1, 4, itime);       // E_y
  axpy(m_Yt3, 1.0, m_Yt2);              // Yt3 -= E_y * D_x
#pragma omp barrier

  mult_sigma3(m_Yt1, m_Yt3);
  axpy(Xt, fac, m_Yt1);
#pragma omp barrier
}


//====================================================================
void Fopr_NonRelativistic::add_deltaH5(Field& Xt, const Field& Yt,
                                       const int itime)
{
  double fac = m_coeff[4] / (24.0 * m_MQ);

  calc_Delta4(m_Yt1, Yt, itime);
  axpy(Xt, fac, m_Yt1);
#pragma omp barrier
}


//====================================================================
void Fopr_NonRelativistic::add_deltaH46(Field& Xt, const Field& Yt,
                                        const int itime)
{
  double fac4 = -m_coeff[3] / (8.0 * m_MQ * m_MQ * m_MQ);
  double fac6 = -m_coeff[5] / (16.0 * double(m_nstab) * m_MQ * m_MQ);
  double fac  = fac4 + fac6;

  calc_Delta2(m_Yt1, Yt, itime);
  calc_Delta2(m_Yt2, m_Yt1, itime);
  axpy(Xt, fac, m_Yt2);
#pragma omp barrier
}


//====================================================================
void Fopr_NonRelativistic::calc_Delta2(Field& Xt, const Field& Yt,
                                       const int itime)
{
  int ipet = itime / m_Nt;
  int IPEt = Communicator::ipe(3);
  if (ipet != IPEt) return;

#pragma omp barrier

  Xt.set(0.0);
  axpy(Xt, -6.0, Yt);
#pragma omp barrier

  for (int idir = 0; idir < m_Ndim - 1; ++idir) {
    shift_backward(m_Wt1, Yt, idir, itime);
    mult_Gn(m_Wt2, m_Wt1, idir, itime);
    axpy(Xt, 1.0, m_Wt2);
#pragma omp barrier

    mult_Gd(m_Wt1, Yt, idir, itime);
    shift_forward(m_Wt2, m_Wt1, idir, itime);
    axpy(Xt, 1.0, m_Wt2);
#pragma omp barrier
  }
}


//====================================================================
void Fopr_NonRelativistic::calc_Delta4(Field& Xt, const Field& Yt,
                                       const int itime)
{
  int ipet = itime / m_Nt;
  int IPEt = Communicator::ipe(3);
  if (ipet != IPEt) return;

#pragma omp barrier

  Xt.set(0.0);
#pragma omp barrier

  for (int idir = 0; idir < m_Ndim - 1; ++idir) {
    m_Wt3.set(0.0);

    shift_backward(m_Wt1, Yt, idir, itime);
    mult_Gn(m_Wt2, m_Wt1, idir, itime);
    axpy(m_Wt3, 1.0, m_Wt2);
#pragma omp barrier

    mult_Gd(m_Wt1, Yt, idir, itime);
    shift_forward(m_Wt2, m_Wt1, idir, itime);
    axpy(m_Wt3, 1.0, m_Wt2);

    axpy(m_Wt3, -2.0, Yt);
#pragma omp barrier

    shift_backward(m_Wt1, m_Wt3, idir, itime);
    mult_Gn(m_Wt2, m_Wt1, idir, itime);
    axpy(Xt, 1.0, m_Wt2);
#pragma omp barrier

    mult_Gd(m_Wt1, m_Wt3, idir, itime);
    shift_forward(m_Wt2, m_Wt1, idir, itime);
    axpy(Xt, 1.0, m_Wt2);

    axpy(Xt, -2.0, m_Wt3);
#pragma omp barrier
  }
}


//====================================================================
void Fopr_NonRelativistic::calc_Delta1(Field& Xt, const Field& Yt,
                                       const int idir, const int itime)
{
  int ipet = itime / m_Nt;
  int IPEt = Communicator::ipe(3);
  if (ipet != IPEt) return;

#pragma omp barrier

  Xt.set(0.0);
#pragma omp barrier

  shift_backward(m_Wt1, Yt, idir, itime);
  mult_Gn(m_Wt2, m_Wt1, idir, itime);
  axpy(Xt, 0.5, m_Wt2);
#pragma omp barrier

  mult_Gd(m_Wt1, Yt, idir, itime);
  shift_forward(m_Wt2, m_Wt1, idir, itime);
  axpy(Xt, -0.5, m_Wt2);
#pragma omp barrier
}


//====================================================================
void Fopr_NonRelativistic::set_source(Field& Xt, const Field& w)
{
  // here check of Field size is necessary.

  vout.detailed(m_vl, " source norm2 = %e.\n", w.norm2());
  // vout.paranoiac(m_vl, " source norm2 = %e.\n", w.norm2());

#pragma omp master
  {
    int    IPEt    = Communicator::ipe(m_Ndim - 1);
    int    Nvcd2   = m_Nvc * m_Nd2; // only upper spinor components matter
    double epsilon = 1.0e-10;

    // extracting source time slice
    int ndef_src = 0;
    for (int itime = 0; itime < m_Ltime; ++itime) {
      int it   = itime % m_Nt;
      int ipet = itime / m_Nt;

      double norm2t = 0.0;
      if (ipet == IPEt) {
        for (int ispc = 0; ispc < m_Nspc; ++ispc) {
          for (int ivcd = 0; ivcd < Nvcd2; ++ivcd) {
            double wt = w.cmp(ivcd, ispc + m_Nspc * it, 0);
            norm2t += wt * wt;
          }
        }
      }
      norm2t = Communicator::reduce_sum(norm2t);

      if (norm2t > epsilon) {
        m_itime_src = itime;
        ndef_src   += 1;
      }
    }

    if (ndef_src == 0) {
      vout.crucial(m_vl, "Error in %s: source is not set.\n",
                   class_name.c_str());
      exit(EXIT_FAILURE);
    }

    if (ndef_src > 1) {
      vout.crucial(m_vl, "Error in %s: source is multiplly set.\n",
                   class_name.c_str());
      exit(EXIT_FAILURE);
    }

    // setting source vector on single time slice
    int it   = m_itime_src % m_Nt;
    int ipet = m_itime_src / m_Nt;
    if (ipet == IPEt) {
      for (int ispc = 0; ispc < m_Nspc; ++ispc) {
        for (int id = 0; id < m_Nd2; ++id) {
          for (int ivc = 0; ivc < m_Nvc; ++ivc) {
            int    ivcd = ivc + m_Nvc * id;
            double wt   = w.cmp(ivcd, ispc + m_Nspc * it, 0);
            Xt.set(ivcd, ispc, 0, wt);
          }
        }
      }
    } else {
      Xt.set(0.0);
    }
  } // #pragma omp master

  vout.detailed(m_vl, "  source is set on itime = %d.\n", m_itime_src);
  vout.detailed(m_vl, "  source norm2 = %e.\n", Xt.norm2());

#pragma omp barrier
}


//====================================================================
void Fopr_NonRelativistic::shift_forward(Field& Xt, const Field& Yt,
                                         const int idir, const int itime)
{  // Xt(x) = Yt(x-\hat{idir})
  if (idir == 0) {
    shift_xdn(Xt, Yt, itime);
  } else if (idir == 1) {
    shift_ydn(Xt, Yt, itime);
  } else if (idir == 2) {
    shift_zdn(Xt, Yt, itime);
  } else {
    vout.crucial(m_vl, "%s: irrelevant direction = %d\n",
                 class_name.c_str(), idir);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void Fopr_NonRelativistic::shift_backward(Field& Xt, const Field& Yt,
                                          const int idir, const int itime)
{  // Xt(x) = Yt(x+\hat{idir})
  if (idir == 0) {
    shift_xup(Xt, Yt, itime);
  } else if (idir == 1) {
    shift_yup(Xt, Yt, itime);
  } else if (idir == 2) {
    shift_zup(Xt, Yt, itime);
  } else {
    vout.crucial(m_vl, "%s: irrelevant direction = %d\n",
                 class_name.c_str(), idir);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void Fopr_NonRelativistic::shift_xup(Field& Xt, const Field& Yt,
                                     const int itime)
{
  int idir = 0;

  int ipet = itime / m_Nt;
  int IPEt = Communicator::ipe(3);
  if (ipet != IPEt) return;

  double       *xp = Xt.ptr(0);
  const double *yp = Yt.ptr(0);

  int NinF2 = 2 * m_Nc * m_Nd2;
  int Nyz   = m_Ny * m_Nz;

  double bc2 = 1.0;
  if (Communicator::ipe(idir) == 0) bc2 = m_boundary[idir];

#pragma omp barrier

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nyz);

  // boundary
  int ix = 0;
  for (int iyz = is; iyz < ns; ++iyz) {
    int nei = ix + m_Nx * iyz;
    for (int in = 0; in < NinF2; ++in) {
      buf1_x[in + NinF2 * iyz] = bc2 * yp[in + NinF2 * nei];
    }
  }
#pragma omp barrier

#pragma omp master
  {
    const int size = NinF2 * Nyz;
    Communicator::exchange(size, &buf2_x[0], &buf1_x[0], 0, 1, 0);
  }
#pragma omp barrier

  set_threadtask(ith, nth, is, ns, m_Nspc);

  for (int site = is; site < ns; ++site) {
    int ix  = site % m_Nx;
    int iyz = site / m_Nx;
    if (ix < m_Nx - 1) { // bulk
      int nei = ix + 1 + m_Nx * iyz;
      for (int in = 0; in < NinF2; ++in) {
        xp[in + NinF2 * site] = yp[in + NinF2 * nei];
      }
    } else {           // boundary
      for (int in = 0; in < NinF2; ++in) {
        xp[in + NinF2 * site] = buf2_x[in + NinF2 * iyz];
      }
    }
  }

#pragma omp barrier
}


//====================================================================
void Fopr_NonRelativistic::shift_xdn(Field& Xt, const Field& Yt,
                                     const int itime)
{
  int idir = 0;

  int ipet = itime / m_Nt;
  int IPEt = Communicator::ipe(3);
  if (ipet != IPEt) return;

  double       *xp = Xt.ptr(0);
  const double *yp = Yt.ptr(0);

  int NinF2 = 2 * m_Nc * m_Nd2;
  int Nyz   = m_Ny * m_Nz;

  double bc2 = 1.0;
  if (Communicator::ipe(idir) == 0) bc2 = m_boundary[idir];

#pragma omp barrier

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nyz);

  // boundary
  int ix = m_Nx - 1;
  for (int iyz = is; iyz < ns; ++iyz) {
    int nei = ix + m_Nx * iyz;
    for (int in = 0; in < NinF2; ++in) {
      buf1_x[in + NinF2 * iyz] = yp[in + NinF2 * nei];
    }
  }
#pragma omp barrier

#pragma omp master
  {
    const int size = NinF2 * Nyz;
    Communicator::exchange(size, &buf2_x[0], &buf1_x[0], 0, -1, 4);
  }
#pragma omp barrier

  set_threadtask(ith, nth, is, ns, m_Nspc);

  for (int site = is; site < ns; ++site) {
    int ix  = site % m_Nx;
    int iyz = site / m_Nx;
    if (ix > 0) {  // bulk
      int nei = ix - 1 + m_Nx * iyz;
      for (int in = 0; in < NinF2; ++in) {
        xp[in + NinF2 * site] = yp[in + NinF2 * nei];
      }
    } else {       // boundary
      for (int in = 0; in < NinF2; ++in) {
        xp[in + NinF2 * site] = bc2 * buf2_x[in + NinF2 * iyz];
      }
    }
  }
#pragma omp barrier
}


//====================================================================
void Fopr_NonRelativistic::shift_yup(Field& Xt, const Field& Yt,
                                     const int itime)
{
  int idir = 1;

  int ipet = itime / m_Nt;
  int IPEt = Communicator::ipe(3);
  if (ipet != IPEt) return;

  double       *xp = Xt.ptr(0);
  const double *yp = Yt.ptr(0);

  int NinF2 = 2 * m_Nc * m_Nd2;
  int Nxz   = m_Nx * m_Nz;

  double bc2 = 1.0;
  if (Communicator::ipe(idir) == 0) bc2 = m_boundary[idir];

#pragma omp barrier

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nxz);

  // boundary
  int iy = 0;
  for (int ixz = is; ixz < ns; ++ixz) {
    int ix  = ixz % m_Nx;
    int iz  = ixz / m_Nx;
    int nei = ix + m_Nx * (iy + m_Ny * iz);
    for (int in = 0; in < NinF2; ++in) {
      buf1_y[in + NinF2 * ixz] = bc2 * yp[in + NinF2 * nei];
    }
  }
#pragma omp barrier

#pragma omp master
  {
    const int size = NinF2 * Nxz;
    Communicator::exchange(size, &buf2_y[0], &buf1_y[0], 1, 1, 1);
  }
#pragma omp barrier

  set_threadtask(ith, nth, is, ns, m_Nspc);

  for (int site = is; site < ns; ++site) {
    int ix = site % m_Nx;
    int iy = (site / m_Nx) % m_Ny;
    int iz = site / (m_Nx * m_Ny);
    if (iy < m_Ny - 1) { // bulk
      int nei = ix + m_Nx * (iy + 1 + m_Ny * iz);
      for (int in = 0; in < NinF2; ++in) {
        xp[in + NinF2 * site] = yp[in + NinF2 * nei];
      }
    } else {           // boundary
      int ixz = ix + m_Nx * iz;
      for (int in = 0; in < NinF2; ++in) {
        xp[in + NinF2 * site] = buf2_y[in + NinF2 * ixz];
      }
    }
  }
#pragma omp barrier
}


//====================================================================
void Fopr_NonRelativistic::shift_ydn(Field& Xt, const Field& Yt,
                                     const int itime)
{
  int idir = 1;

  int ipet = itime / m_Nt;
  int IPEt = Communicator::ipe(3);
  if (ipet != IPEt) return;

  double       *xp = Xt.ptr(0);
  const double *yp = Yt.ptr(0);

  int NinF2 = 2 * m_Nc * m_Nd2;
  int Nxz   = m_Nx * m_Nz;

  double bc2 = 1.0;
  if (Communicator::ipe(idir) == 0) bc2 = m_boundary[idir];

#pragma omp barrier

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nxz);

  // boundary
  int iy = m_Ny - 1;
  for (int ixz = is; ixz < ns; ++ixz) {
    int ix  = ixz % m_Nx;
    int iz  = ixz / m_Nx;
    int nei = ix + m_Nx * (iy + m_Ny * iz);
    for (int in = 0; in < NinF2; ++in) {
      buf1_y[in + NinF2 * ixz] = yp[in + NinF2 * nei];
    }
  }
#pragma omp barrier

#pragma omp master
  {
    const int size = NinF2 * Nxz;
    Communicator::exchange(size, &buf2_y[0], &buf1_y[0], 1, -1, 5);
  }
#pragma omp barrier

  set_threadtask(ith, nth, is, ns, m_Nspc);

  for (int site = is; site < ns; ++site) {
    int ix = site % m_Nx;
    int iy = (site / m_Nx) % m_Ny;
    int iz = site / (m_Nx * m_Ny);
    if (iy > 0) { // bulk
      int nei = ix + m_Nx * (iy - 1 + m_Ny * iz);
      for (int in = 0; in < NinF2; ++in) {
        xp[in + NinF2 * site] = yp[in + NinF2 * nei];
      }
    } else {      // boundary
      int ixz = ix + m_Nx * iz;
      for (int in = 0; in < NinF2; ++in) {
        xp[in + NinF2 * site] = bc2 * buf2_y[in + NinF2 * ixz];
      }
    }
  }

#pragma omp barrier
}


//====================================================================
void Fopr_NonRelativistic::shift_zup(Field& Xt, const Field& Yt,
                                     const int itime)
{
  int idir = 2;

  int ipet = itime / m_Nt;
  int IPEt = Communicator::ipe(3);
  if (ipet != IPEt) return;

  double       *xp = Xt.ptr(0);
  const double *yp = Yt.ptr(0);

  int    NinF2 = 2 * m_Nc * m_Nd2;
  int    Nxy   = m_Nx * m_Ny;
  double bc2   = 1.0;
  if (Communicator::ipe(idir) == 0) bc2 = m_boundary[idir];

#pragma omp barrier

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nxy);

  // boundary
  int iz = 0;
  for (int ixy = is; ixy < ns; ++ixy) {
    int nei = ixy + Nxy * iz;
    for (int in = 0; in < NinF2; ++in) {
      buf1_z[in + NinF2 * ixy] = bc2 * yp[in + NinF2 * nei];
    }
  }
#pragma omp barrier

#pragma omp master
  {
    const int size = NinF2 * Nxy;
    Communicator::exchange(size, &buf2_z[0], &buf1_z[0], 2, 1, 2);
  }
#pragma omp barrier

  set_threadtask(ith, nth, is, ns, m_Nspc);

  for (int site = is; site < ns; ++site) {
    int ixy = site % Nxy;
    int iz  = site / Nxy;
    if (iz < m_Nz - 1) {  // bulk
      int nei = ixy + Nxy * (iz + 1);
      for (int in = 0; in < NinF2; ++in) {
        xp[in + NinF2 * site] = yp[in + NinF2 * nei];
      }
    } else {            // boundary
      for (int in = 0; in < NinF2; ++in) {
        xp[in + NinF2 * site] = buf2_z[in + NinF2 * ixy];
      }
    }
  }
#pragma omp barrier
}


//====================================================================
void Fopr_NonRelativistic::shift_zdn(Field& Xt, const Field& Yt,
                                     const int itime)
{
  int idir = 2;

  int ipet = itime / m_Nt;
  int IPEt = Communicator::ipe(3);
  if (ipet != IPEt) return;

  double       *xp = Xt.ptr(0);
  const double *yp = Yt.ptr(0);

  int    NinF2 = 2 * m_Nc * m_Nd2;
  int    Nxy   = m_Nx * m_Ny;
  double bc2   = 1.0;
  if (Communicator::ipe(idir) == 0) bc2 = m_boundary[idir];

#pragma omp barrier

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nxy);

  // boundary
  int iz = m_Nz - 1;
  for (int ixy = is; ixy < ns; ++ixy) {
    int nei = ixy + Nxy * iz;
    for (int in = 0; in < NinF2; ++in) {
      buf1_z[in + NinF2 * ixy] = yp[in + NinF2 * nei];
    }
  }
#pragma omp barrier

#pragma omp master
  {
    const int size = NinF2 * Nxy;
    Communicator::exchange(size, &buf2_z[0], &buf1_z[0], 2, -1, 6);
  }
#pragma omp barrier

  set_threadtask(ith, nth, is, ns, m_Nspc);

  // bulk
  for (int site = is; site < ns; ++site) {
    int ixy = site % Nxy;
    int iz  = site / Nxy;
    if (iz > 0) {
      int nei = ixy + Nxy * (iz - 1);
      for (int in = 0; in < NinF2; ++in) {
        xp[in + NinF2 * site] = yp[in + NinF2 * nei];
      }
    } else {  // boundary
      for (int in = 0; in < NinF2; ++in) {
        xp[in + NinF2 * site] = bc2 * buf2_z[in + NinF2 * ixy];
      }
    }
  }
#pragma omp barrier
}


//====================================================================
void Fopr_NonRelativistic::mult_sigma1(Field& Xt, const Field& Yt)
{
  int Nc2 = 2 * m_Nc;
  int Ncd = Nc2 * m_Nd2;

  double       *xp = Xt.ptr(0);
  const double *yp = Yt.ptr(0);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nspc);

#pragma omp barrier

  for (int site = is; site < ns; ++site) {
    int iv = Ncd * site;
    for (int ic = 0; ic < m_Nc; ++ic) {
      xp[2 * ic + Nc2 * 0 + iv]     = yp[2 * ic + Nc2 * 1 + iv];
      xp[2 * ic + 1 + Nc2 * 0 + iv] = yp[2 * ic + 1 + Nc2 * 1 + iv];
      xp[2 * ic + Nc2 * 1 + iv]     = yp[2 * ic + Nc2 * 0 + iv];
      xp[2 * ic + 1 + Nc2 * 1 + iv] = yp[2 * ic + 1 + Nc2 * 0 + iv];
    }
  }
#pragma omp barrier
}


//====================================================================
void Fopr_NonRelativistic::mult_sigma2(Field& Xt, const Field& Yt)
{
  int Nc2 = 2 * m_Nc;
  int Ncd = Nc2 * m_Nd2;

  double       *xp = Xt.ptr(0);
  const double *yp = Yt.ptr(0);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nspc);

#pragma omp barrier

  for (int site = is; site < ns; ++site) {
    int iv = Ncd * site;
    for (int ic = 0; ic < m_Nc; ++ic) {
      xp[2 * ic + Nc2 * 0 + iv]     = yp[2 * ic + 1 + Nc2 * 1 + iv];
      xp[2 * ic + 1 + Nc2 * 0 + iv] = -yp[2 * ic + Nc2 * 1 + iv];
      xp[2 * ic + Nc2 * 1 + iv]     = -yp[2 * ic + 1 + Nc2 * 0 + iv];
      xp[2 * ic + 1 + Nc2 * 1 + iv] = yp[2 * ic + Nc2 * 0 + iv];
    }
  }
#pragma omp barrier
}


//====================================================================
void Fopr_NonRelativistic::mult_sigma3(Field& Xt, const Field& Yt)
{
  int Nc2 = 2 * m_Nc;
  int Ncd = Nc2 * m_Nd2;

  double       *xp = Xt.ptr(0);
  const double *yp = Yt.ptr(0);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nspc);

#pragma omp barrier

  for (int site = is; site < ns; ++site) {
    int iv = Ncd * site;
    for (int ic = 0; ic < m_Nc; ++ic) {
      xp[2 * ic + Nc2 * 0 + iv]     = yp[2 * ic + Nc2 * 0 + iv];
      xp[2 * ic + 1 + Nc2 * 0 + iv] = yp[2 * ic + 1 + Nc2 * 0 + iv];
      xp[2 * ic + Nc2 * 1 + iv]     = -yp[2 * ic + Nc2 * 1 + iv];
      xp[2 * ic + 1 + Nc2 * 1 + iv] = -yp[2 * ic + 1 + Nc2 * 1 + iv];
    }
  }
#pragma omp barrier
}


//====================================================================
void Fopr_NonRelativistic::mult_Gn(Field& Xt, const Field& Yt,
                                   const int idir, const int itime)
{
  int ipet = itime / m_Nt;
  int IPEt = Communicator::ipe(3);
  if (ipet != IPEt) return;

  int          it  = itime % m_Nt;
  const double *up = m_U.ptr(0, 0, idir);

  double       *xp = Xt.ptr(0);
  const double *yp = Yt.ptr(0);

  int Nc2 = 2 * m_Nc;
  int Ndf = 2 * m_Nc * m_Nc;
  int Ncd = 2 * m_Nc * m_Nd2;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nspc);

#pragma omp barrier

  for (int site = is; site < ns; ++site) {
    int ig = Ndf * (site + m_Nspc * it);
    int iv = Ncd * site;
    for (int id = 0; id < m_Nd2; ++id) {
      for (int ic = 0; ic < m_Nc; ++ic) {
        int ig2 = Nc2 * ic + ig;
        int iv2 = Nc2 * id + iv;
        xp[2 * ic + iv2]     = mult_Gn_r(&up[ig2], &yp[iv2], m_Nc);
        xp[2 * ic + 1 + iv2] = mult_Gn_i(&up[ig2], &yp[iv2], m_Nc);
      }
    }
  }
#pragma omp barrier
}


//====================================================================
void Fopr_NonRelativistic::mult_F(Field& Xt, const Field& Yt,
                                  const int icomp, const int itime)
{
  int ipet = itime / m_Nt;
  int IPEt = Communicator::ipe(3);
  if (ipet != IPEt) return;

  int          it  = itime % m_Nt;
  const double *up = m_Fstr[icomp].ptr(0);

  double       *xp = Xt.ptr(0);
  const double *yp = Yt.ptr(0);

  int Nc2 = 2 * m_Nc;
  int Ndf = 2 * m_Nc * m_Nc;
  int Ncd = 2 * m_Nc * m_Nd2;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nspc);

#pragma omp barrier

  for (int site = is; site < ns; ++site) {
    int ig = Ndf * (site + m_Nspc * it);
    int iv = Ncd * site;
    for (int id = 0; id < m_Nd2; ++id) {
      for (int ic = 0; ic < m_Nc; ++ic) {
        int ig2 = Nc2 * ic + ig;
        int iv2 = Nc2 * id + iv;
        xp[2 * ic + iv2]     = mult_Gn_r(&up[ig2], &yp[iv2], m_Nc);
        xp[2 * ic + 1 + iv2] = mult_Gn_i(&up[ig2], &yp[iv2], m_Nc);
      }
    }
  }
#pragma omp barrier
}


//====================================================================
void Fopr_NonRelativistic::mult_Gd(Field& Xt, const Field& Yt,
                                   const int idir, const int itime)
{
  int ipet = itime / m_Nt;
  int IPEt = Communicator::ipe(3);
  if (ipet != IPEt) return;

  int          it  = itime % m_Nt;
  const double *up = m_U.ptr(0, 0, idir);

  double       *xp = Xt.ptr(0);
  const double *yp = Yt.ptr(0);

  int Nc2 = 2 * m_Nc;
  int Ndf = 2 * m_Nc * m_Nc;
  int Ncd = 2 * m_Nc * m_Nd2;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nspc);

#pragma omp barrier

  for (int site = is; site < ns; ++site) {
    int ig = Ndf * (site + m_Nspc * it);
    int iv = Ncd * site;
    for (int id = 0; id < m_Nd2; ++id) {
      for (int ic = 0; ic < m_Nc; ++ic) {
        int ig2 = 2 * ic + ig;
        int iv2 = Nc2 * id + iv;
        xp[2 * ic + iv2]     = mult_Gd_r(&up[ig2], &yp[iv2], m_Nc);
        xp[2 * ic + 1 + iv2] = mult_Gd_i(&up[ig2], &yp[iv2], m_Nc);
      }
    }
  }

#pragma omp barrier
}


//====================================================================
double Fopr_NonRelativistic::flop_count()
{
  return flop_count("Dirac");
}


//====================================================================
double Fopr_NonRelativistic::flop_count(const std::string repr)
{
  // not implemented yet

  const double gflop = 0.0;

  return gflop;
}


//============================================================END=====
