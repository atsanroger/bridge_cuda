/*!
        @file    staple_eo.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "staple_eo.h"

#include "Field/field_thread-inc.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Staple_eo::register_factory();
}
#endif

const std::string Staple_eo::class_name = "Staple_eo";

//====================================================================
void Staple_eo::init()
{
  m_vl = CommonParameters::Vlevel();

  m_filename_output = "stdout";

  m_Nc   = CommonParameters::Nc();
  m_Nvol = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();

  int Ndf = 2 * m_Nc * m_Nc;
  m_shift = new ShiftField_eo(Ndf);
}


//====================================================================
void Staple_eo::init(const Parameters& params)
{
  m_vl = CommonParameters::Vlevel();

  m_filename_output = "stdout";

  m_Nc   = CommonParameters::Nc();
  m_Nvol = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();

  int Ndf = 2 * m_Nc * m_Nc;
  m_shift = new ShiftField_eo(Ndf);

  set_parameters(params);
}


//====================================================================
void Staple_eo::tidyup()
{
  delete m_shift;
}


//====================================================================
void Staple_eo::set_parameters(const Parameters& params)
{
  m_filename_output = params.get_string("filename_output");
  if (m_filename_output.empty()) {
    m_filename_output = "stdout";
  }

  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }
}


//====================================================================
void Staple_eo::get_parameters(Parameters& params) const
{
  params.set_string("filename_output", m_filename_output);
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
double Staple_eo::plaquette(const Field_G& U)
{
  Field_G  Ueo(U);
  Index_eo index_eo;

  index_eo.convertField(Ueo, U);

  return (plaq_s(Ueo) + plaq_t(Ueo)) / 2;
}


//====================================================================
double Staple_eo::plaq_s(const Field_G& U)
{
  int nth = ThreadManager::get_num_threads();

  double plaq;

  if (nth > 1) {
    plaq = plaq_s_impl(U);
  } else {
    plaq = plaq_s_omp(U);
  }

  return plaq;
}


//====================================================================
double Staple_eo::plaq_t(const Field_G& U)
{
  int nth = ThreadManager::get_num_threads();

  double plaq;

  if (nth > 1) {
    plaq = plaq_t_impl(U);
  } else {
    plaq = plaq_t_omp(U);
  }

  return plaq;
}


//====================================================================
double Staple_eo::plaq_s_omp(const Field_G& U)
{
  double plaq;

#pragma omp parallel
  {
    double plaq1 = plaq_s_impl(U);
#pragma omp master
    plaq = plaq1;
  }

  return plaq;
}


//====================================================================
double Staple_eo::plaq_t_omp(const Field_G& U)
{
  double plaq;

#pragma omp parallel
  {
    double plaq1 = plaq_t_impl(U);
#pragma omp master
    plaq = plaq1;
  }

  return plaq;
}


//====================================================================
double Staple_eo::plaq_s_impl(const Field_G& Ueo)
{
  const int Ndim_spc = m_Ndim - 1;

  const int NPE = CommonParameters::NPE();

  //  Field_G staple;
  double plaq = 0.0;

  for (int mu = 0; mu < Ndim_spc; ++mu) {
    int nu = (mu + 1) % Ndim_spc;

    copy(m_v1, 0, Ueo, mu);

    upper(m_v2, Ueo, mu, nu);

    double plaq1 = real(dotc(m_v2, m_v1));

    plaq += plaq1;
  }

  /*
  upper(m_staple, Ueo, 0, 1);
  for (int site = 0; site < Nvol; site++) {
    plaq += ReTr(Ueo.mat(site, 0) * m_staple.mat_dag(site));   // P_xy
  }

  upper(m_staple, Ueo, 1, 2);
  for (int site = 0; site < Nvol; site++) {
    plaq += ReTr(Ueo.mat(site, 1) * m_staple.mat_dag(site));   // P_yz
  }

  upper(m_staple, Ueo, 2, 0);
  for (int site = 0; site < Nvol; site++) {
    plaq += ReTr(Ueo.mat(site, 2) * m_staple.mat_dag(site));   // P_zx
  }
  */

  //  plaq = Communicator::reduce_sum(plaq);

  plaq *= 1.0 / (m_Nvol * m_Nc * Ndim_spc);
  plaq *= 1.0 / NPE;
  return plaq;
}


//====================================================================
double Staple_eo::plaq_t_impl(const Field_G& Ueo)
{
  const int Ndim_spc = m_Ndim - 1;

  const int NPE = CommonParameters::NPE();

  //  Field_G staple;
  int    mu   = m_Ndim - 1;
  double plaq = 0.0;

#pragma omp barrier

  copy(m_v1, 0, Ueo, mu);
#pragma omp barrier

  for (int nu = 0; nu < Ndim_spc; ++nu) {
    upper(m_v2, Ueo, mu, nu);
    double plaq1 = real(dotc(m_v2, m_v1));
    plaq += plaq1;
  }

  plaq *= 1.0 / (m_Nvol * m_Nc * Ndim_spc);
  plaq *= 1.0 / NPE;

  return plaq;

  /*
  for (int nu = 0; nu < Ndim - 1; nu++) {
    lower(m_staple, Ueo, 3, nu);
    for (int site = 0; site < Nvol; site++) {
      plaq += ReTr(Ueo.mat(site, 3) * m_staple.mat_dag(site));   // P_zx
    }
  }

  plaq = Communicator::reduce_sum(plaq);

  return plaq / Nvol / NPE / Nc / Ndim_spc;
  */
}


//====================================================================
void Staple_eo::staple(Field_G& W, const Field_G& Ueo, const int mu)
{
#pragma omp barrier

  W.set(0.0);
#pragma omp barrier

  for (int nu = 0; nu < m_Ndim; nu++) {
    if (nu != mu) {
      Field_G u_tmp(W.nvol(), 1);

      upper(u_tmp, Ueo, mu, nu);
      axpy(W, 1.0, u_tmp);
#pragma omp barrier

      lower(u_tmp, Ueo, mu, nu);
      axpy(W, 1.0, u_tmp);
#pragma omp barrier
    }
  }
}


//====================================================================
void Staple_eo::upper(Field_G& c, const Field_G& Ueo,
                      const int mu, const int nu)
{
  // (1)  mu (2)
  //    +-->--+
  // nu |     |
  //   i+     +

#pragma omp barrier

  copy(m_Umu, 0, Ueo, mu);
  copy(m_Unu, 0, Ueo, nu);
#pragma omp barrier

  m_shift->backward(m_w2, m_Unu, mu);
  m_shift->backward(c, m_Umu, nu);

  mult_Field_Gnd(m_w1, 0, c, 0, m_w2, 0);
  mult_Field_Gnn(c, 0, m_Unu, 0, m_w1, 0);

#pragma omp barrier
}


//====================================================================
void Staple_eo::lower(Field_G& c, const Field_G& Ueo,
                      const int mu, const int nu)
{
  //    +     +
  // nu |     |
  //   i+-->--+
  //  (1)  mu (2)

#pragma omp barrier

  copy(m_Umu, 0, Ueo, mu);
  copy(m_Unu, 0, Ueo, nu);
#pragma omp barrier

  m_shift->backward(m_w1, m_Unu, mu);
  mult_Field_Gnn(m_w2, 0, m_Umu, 0, m_w1, 0);
  mult_Field_Gdn(m_w1, 0, m_Unu, 0, m_w2, 0);
  m_shift->forward(c, m_w1, nu);
}


//============================================================END=====
