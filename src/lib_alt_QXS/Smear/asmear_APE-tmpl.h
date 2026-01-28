/*!
        @file    asmear_APE-tmpl.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-07-03 23:32:09 #$

        @version $LastChangedRevision: 2529 $
*/

#include "lib_alt_QXS/Smear/asmear_APE.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Smear_APE::register_factory();
}
#endif

template<typename AFIELD>
const std::string ASmear_APE<AFIELD>::class_name = "ASmear_APE<AFIELD>";

//====================================================================
template<typename AFIELD>
void ASmear_APE<AFIELD>::init(AProjection<AFIELD>* proj,
                              const Parameters& params)
{
  vout.general(m_vl, "%s: construction\n", class_name.c_str());
  vout.increase_indent();

  m_vl = CommonParameters::Vlevel();

  int Nc  = CommonParameters::Nc();
  int Ndf = 2 * Nc * Nc;
  m_Nvol = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();

  m_rho.resize(m_Ndim * m_Ndim);
  for(int i = 0; i < m_Ndim * m_Ndim; ++i){
    m_rho[i] = 0.0;
  }

  m_proj = proj;

  set_parameters(params);

  m_v1.reset(Ndf, m_Nvol, 1);
  m_v2.reset(Ndf, m_Nvol, 1);
  m_v3.reset(Ndf, m_Nvol, 1);

  m_w1.reset(Ndf, m_Nvol, m_Ndim);
  m_w2.reset(Ndf, m_Nvol, m_Ndim);

  m_staple = new AStaple_lex<AFIELD>();

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished\n",
               class_name.c_str());
}

//====================================================================
template<typename AFIELD>
void ASmear_APE<AFIELD>::tidyup()
{
  delete m_staple;
}

//====================================================================
template<typename AFIELD>
void ASmear_APE<AFIELD>::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  double rho1;

  int err = 0;
  err += params.fetch_double("rho_uniform", rho1);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(rho1);
}


//====================================================================
template<typename AFIELD>
void ASmear_APE<AFIELD>::get_parameters(Parameters& params) const
{
  std::vector<double> rho(m_rho.size());

  for (size_t i = 0; i < m_rho.size(); ++i) {
    rho[i] = m_rho[i];
  }
  params.set_double_vector("rho", rho);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));

  return;
}


//====================================================================
template<typename AFIELD>
void ASmear_APE<AFIELD>::set_parameters(const double rho1)
{
  //- print input parameters
  vout.general(m_vl, "%s: parameters\n", class_name.c_str());
  vout.general(m_vl, "  rho = %8.4f\n", rho1);

  //- range check
  // NB. rho == 0 is allowed.

  //- store values
  // m_rho.resize(m_Ndim * m_Ndim);  // already resized in init.
  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      m_rho[mu + nu * m_Ndim] = rho1;
    }
  }
}


//====================================================================
template<typename AFIELD>
void ASmear_APE<AFIELD>::set_parameters(const std::vector<double>& rho)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  rho[%d] = %8.4f\n", mu, rho[mu]);
  }

  // range check
  // NB. rho == 0 is allowed.
  assert(rho.size() == m_Ndim * m_Ndim);

  //- store values
  // m_rho.resize(m_Ndim * m_Ndim);  // already resized in init.
  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      m_rho[mu + nu * m_Ndim] = rho[mu + nu * m_Ndim];
    }
  }
}


//====================================================================
template<typename AFIELD>
void ASmear_APE<AFIELD>::smear(Field_G& Usmear, const Field_G& U)
{
#pragma omp barrier
  // ThreadManager::assert_single_thread(class_name);

  //int Nin = U.nin();
  //AFIELD w1(Nin, m_Nvol, m_Ndim);
  //AFIELD w2(Nin, m_Nvol, m_Ndim);

  AIndex_lex<real_t,AFIELD::IMPL> index;

  convert_gauge(index, m_w1, U);

  smear(m_w2, m_w1);

  reverse_gauge(index, Usmear, m_w2);

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void ASmear_APE<AFIELD>::smear(AFIELD& Usmear, const AFIELD& U)
{
  const int Nvol = CommonParameters::Nvol();

  assert(U.nvol() == Nvol);
  assert(U.nex() == m_Ndim);
  assert(Usmear.nvol() == Nvol);
  assert(Usmear.nex() == m_Ndim);

#pragma omp barrier

  Usmear.set(0.0);
#pragma omp barrier

  for (int mu = 0; mu < m_Ndim; ++mu) {

    m_v1.set(0.0);
    copy(m_v2, 0, U, mu);
#pragma omp barrier

    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu != mu) {
        double rho = m_rho[mu + m_Ndim * nu];
        m_staple->upper(m_v3, U, mu, nu);
        axpy(m_v1, 0, rho, m_v3, 0);
#pragma omp barrier

        m_staple->lower(m_v3, U, mu, nu);
        axpy(m_v1, 0, rho, m_v3, 0);
#pragma omp barrier
      }
    }

    double rho0 = m_rho[mu + m_Ndim * mu];
    m_proj->project(m_v3, rho0, m_v1, m_v2);
    copy(Usmear, mu, m_v3, 0);
#pragma omp barrier

  }

}

//============================================================END=====
