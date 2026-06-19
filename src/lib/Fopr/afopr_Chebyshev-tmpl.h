/*!
        @file    afopr_Chebyshev-tmpl.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
        @version $LastChangedRevision: 2499 $
*/

#include "lib/Fopr/afopr_Chebyshev.h"

#include "ResourceManager/threadManager.h"


template<typename AFIELD>
const std::string AFopr_Chebyshev<AFIELD>::class_name
  = "AFopr_Chebyshev";
//====================================================================
template<typename AFIELD>
void AFopr_Chebyshev<AFIELD>::init(const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);
  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: construction\n", class_name.c_str());
  vout.increase_indent();

  m_NinF = m_fopr->field_nin();
  m_Nvol = m_fopr->field_nvol();
  m_NexF = m_fopr->field_nex();

  vout.general(m_vl, "kernel: %s\n", m_fopr->class_name.c_str());

  set_parameters(params);

  m_dj.resize(3);
  for (int k = 0; k < 3; ++k) {
    m_dj[k].reset(m_NinF, m_Nvol, m_NexF);
  }

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());
}


//====================================================================
template<typename AFIELD>
void AFopr_Chebyshev<AFIELD>::init()
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: construction (obsoletre)\n",
               class_name.c_str());
  vout.increase_indent();

  m_NinF = m_fopr->field_nin();
  m_Nvol = m_fopr->field_nvol();
  m_NexF = m_fopr->field_nex();

  m_dj.resize(3);
  for (int k = 0; k < 3; ++k) {
    m_dj[k].reset(m_NinF, m_Nvol, m_NexF);
  }

  vout.general(m_vl, "kernel: %s\n", m_fopr->class_name.c_str());

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());
}


//====================================================================
template<typename AFIELD>
void AFopr_Chebyshev<AFIELD>::tidyup()
{
  // do nothing.
}


//====================================================================
template<typename AFIELD>
void AFopr_Chebyshev<AFIELD>::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  int    Np;
  double v_threshold, v_max;

  int err = 0;
  err += params.fetch_int("degree_of_polynomial", Np);
  err += params.fetch_double("threshold_value", v_threshold);
  err += params.fetch_double("upper_bound", v_max);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(Np, real_t(v_threshold), real_t(v_max));
}


//====================================================================
template<typename AFIELD>
void AFopr_Chebyshev<AFIELD>::set_parameters(const int Np,
                                             const real_t Vthrs,
                                             const real_t Vmax)
{
  //- range check
  int err = 0;
  err += ParameterCheck::non_negative(Np);
  // NB. Vthrs, Vmax == 0 is allowed.

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  int ith = ThreadManager::get_thread_id();

  if (ith == 0) {
    m_Npcb = Np;

    m_Vthrs = real_t(Vthrs);
    m_Vmax  = real_t(Vmax);

    real_t b_max = Vmax / Vthrs;
    real_t r     = 2.0 / (b_max * b_max - 1.0);
    real_t s     = Vthrs / sqrt(0.5 * r);

    m_Fcb1 = 2.0 / (s * s);
    m_Fcb2 = -(1.0 + r);
  }

  //- print input parameters
  vout.general(m_vl, "%s: input parameters\n", class_name.c_str());
  vout.general(m_vl, "  Np    = %d\n", Np);
  vout.general(m_vl, "  Vthrs = %16.8e\n", Vthrs);
  vout.general(m_vl, "  Vmax  = %16.8e\n", Vmax);
  vout.general(m_vl, "  Fcb1  = %16.8e\n", m_Fcb1);
  vout.general(m_vl, "  Fcb2  = %16.8e\n", m_Fcb2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Chebyshev<AFIELD>::get_parameters(Parameters& params) const
{
  params.set_int("degree_of_polynomial", m_Npcb);
  params.set_double("threshold_value", m_Vthrs);
  params.set_double("upper_bound", m_Vmax);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
template<typename AFIELD>
void AFopr_Chebyshev<AFIELD>::set_config(Field *U)
{
  m_fopr->set_config(U);
}


//====================================================================
template<typename AFIELD>
void AFopr_Chebyshev<AFIELD>::set_mode(std::string mode)
{
  m_mode = mode;
  m_fopr->set_mode(mode);
}


//====================================================================
template<typename AFIELD>
std::string AFopr_Chebyshev<AFIELD>::get_mode() const
{
  return m_mode;
}


//====================================================================
template<typename AFIELD>
void AFopr_Chebyshev<AFIELD>::mult(AFIELD& v, const AFIELD& w)
{
  v.check_size(m_NinF, m_Nvol, m_NexF);
  w.check_size(m_NinF, m_Nvol, m_NexF);

  copy(m_dj[0], w);
  scal(m_dj[0], real_t(-1.0));
  m_dj[1].set(real_t(0.0));

  int jn  = 2;
  int jp1 = 1;
  int jp2 = 0;

  for (int j = m_Npcb; j >= 2; --j) {
    m_fopr->mult(m_dj[jn], m_dj[jp1]);
    scal(m_dj[jn], m_Fcb1);
    axpy(m_dj[jn], m_Fcb2, m_dj[jp1]);

    scal(m_dj[jn], real_t(2.0));
    axpy(m_dj[jn], real_t(-1.0), m_dj[jp2]);

    jn  = (jn + 1) % 3;
    jp1 = (jp1 + 1) % 3;
    jp2 = (jp2 + 1) % 3;
  }

  m_fopr->mult(v, m_dj[jp1]);
  scal(v, m_Fcb1);
  axpy(v, m_Fcb2, m_dj[jp1]);
  axpy(v, real_t(-1.0), m_dj[jp2]);
}


//====================================================================
template<typename AFIELD>
void AFopr_Chebyshev<AFIELD>::mult(real_t& v, const real_t x)
{
  std::vector<real_t> dj(3);

  dj[0] = -1.0;
  dj[1] = 0.0;

  int jn  = 2;
  int jp1 = 1;
  int jp2 = 0;

  for (int j = m_Npcb; j >= 2; --j) {
    dj[jn]  = x * dj[jp1];
    dj[jn] *= m_Fcb1;
    dj[jn] += m_Fcb2 * dj[jp1];

    dj[jn] *= 2.0;
    dj[jn] -= 1.0 * dj[jp2];

    jn  = (jn + 1) % 3;
    jp1 = (jp1 + 1) % 3;
    jp2 = (jp2 + 1) % 3;
  }

  v  = x * dj[jp1];
  v *= m_Fcb1;
  v += m_Fcb2 * dj[jp1];
  v -= dj[jp2];
}


//====================================================================
template<typename AFIELD>
typename AFIELD::real_t AFopr_Chebyshev<AFIELD>::mult(const real_t x)
{
  real_t v;
  mult(v, x);
  return v;
}


//====================================================================
template<typename AFIELD>
double AFopr_Chebyshev<AFIELD>::flop_count()
{
  return 0.0;
}


//============================================================END=====
