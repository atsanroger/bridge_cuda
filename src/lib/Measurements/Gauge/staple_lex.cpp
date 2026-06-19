/*!
        @file    staple_lex.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-10-17 13:51:32 #$

        @version $LastChangedRevision: 2550 $
*/

#include "staple_lex.h"

#include "Field/field_thread-inc.h"


#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Staple_lex::register_factory();
}
#endif

const std::string Staple_lex::class_name = "Staple_lex";

//====================================================================
void Staple_lex::init(const Parameters& params)
{
  setup();
  set_parameters(params);
}


//====================================================================
void Staple_lex::init()
{
  setup();
}


//====================================================================
void Staple_lex::setup()
{
  m_vl = CommonParameters::Vlevel();

  m_Nvol = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();
  m_Nc   = CommonParameters::Nc();
  int Ndf = 2 * m_Nc * m_Nc;

  m_shift = new ShiftField_lex(Ndf);

  m_filename_output = "stdout";
}


//====================================================================
void Staple_lex::tidyup()
{
  delete m_shift;
}


//====================================================================
void Staple_lex::set_parameters(const Parameters& params)
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
void Staple_lex::get_parameters(Parameters& params) const
{
  params.set_string("filename_output", m_filename_output);
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
double Staple_lex::plaquette(const Field_G& U)
{
  ThreadManager::assert_single_thread(class_name);

  const double result = (plaq_s(U) + plaq_t(U)) / 2;

  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  vout.detailed(m_vl, "plaq = %20.16e\n", result);

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }

  return result;
}


//====================================================================
double Staple_lex::plaq_s(const Field_G& U)
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
double Staple_lex::plaq_t(const Field_G& U)
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
double Staple_lex::plaq_s_omp(const Field_G& U)
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
double Staple_lex::plaq_t_omp(const Field_G& U)
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
double Staple_lex::plaq_s_impl(const Field_G& U)
{
  const int Ndim_spc = m_Ndim - 1;

  const int NPE = CommonParameters::NPE();

  double plaq = 0.0;

#pragma omp barrier

  for (int mu = 0; mu < Ndim_spc; ++mu) {
    int nu = (mu + 1) % Ndim_spc;

    copy(m_v1, 0, U, mu);

    upper(m_v2, U, mu, nu);

    double plaq1 = real(dotc(m_v2, m_v1));

    plaq += plaq1;
  }

  plaq *= 1.0 / (m_Nvol * m_Nc * Ndim_spc);
  plaq *= 1.0 / NPE;

  return plaq;
}


//====================================================================
double Staple_lex::plaq_t_impl(const Field_G& U)
{
  const int NPE      = CommonParameters::NPE();
  const int Ndim_spc = m_Ndim - 1;

  const int mu = Ndim_spc;

#pragma omp barrier

  copy(m_v1, 0, U, mu);
#pragma omp barrier

  double plaq = 0.0;

  for (int nu = 0; nu < Ndim_spc; ++nu) {
    upper(m_v2, U, mu, nu);
    double plaq1 = real(dotc(m_v2, m_v1));
    plaq += plaq1;
  }

  plaq *= 1.0 / (m_Nvol * m_Nc * Ndim_spc);
  plaq *= 1.0 / NPE;

  return plaq;
}


//====================================================================
void Staple_lex::staple(Field_G& W, const Field_G& U, const int mu)
{
  const int Ndim = CommonParameters::Ndim();

#pragma omp barrier

  W.set(0.0);
#pragma omp barrier

  for (int nu = 0; nu < Ndim; ++nu) {
    if (nu != mu) {
      upper(m_v1, U, mu, nu);
      axpy(W, 1.0, m_v1);

      lower(m_v1, U, mu, nu);
      axpy(W, 1.0, m_v1);
    }
  }

#pragma omp barrier
}


//====================================================================
void Staple_lex::upper(Field_G& Uup, const Field_G& U,
                       const int mu, const int nu)
{
  // (1)  mu (2)
  //    +-->--+
  // nu |     |
  //   i+     +

#pragma omp barrier

  copy(m_Umu, 0, U, mu);
  copy(m_Unu, 0, U, nu);
#pragma omp barrier

  m_shift->backward(m_w2, m_Unu, mu);
  m_shift->backward(Uup, m_Umu, nu);

  mult_Field_Gnd(m_w1, 0, Uup, 0, m_w2, 0);
  mult_Field_Gnn(Uup, 0, U, nu, m_w1, 0);
}


//====================================================================
void Staple_lex::lower(Field_G& Udn, const Field_G& U,
                       const int mu, const int nu)
{
  //    +     +
  // nu |     |
  //   i+-->--+
  //  (1)  mu (2)

#pragma omp barrier

  copy(m_Unu, 0, U, nu);

  m_shift->backward(m_w1, m_Unu, mu);

  mult_Field_Gnn(m_w2, 0, U, mu, m_w1, 0);
  mult_Field_Gdn(m_w1, 0, U, nu, m_w2, 0);

  m_shift->forward(Udn, m_w1, nu);
}


//============================================================END=====
