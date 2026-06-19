/*!                                                                             
        @file    afopr_OptimalDomainwall-tmpl.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2022-12-16 15:57:38 #$
        @version $LastChangedRevision: 2422 $
*/

#include "afopr_OptimalDomainwall.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
using namespace std;

#include "lib/Parameters/commonParameters.h"
#include "lib/Communicator/communicator.h"
#include "lib/Tools/math_Sign_Zolotarev_Omega.h"

//#include "Field/afield.h"
//#include "Field/afield-inc.h"


template<typename AFIELD>
const std::string AFopr_OptimalDomainwall<AFIELD>::class_name
                                      = "AFopr_OptimalDomainwall";
//====================================================================
template<typename AFIELD>
void AFopr_OptimalDomainwall<AFIELD>::init(const Parameters& params)
{
  vout.general(m_vl, "Initialization of %s:\n", class_name.c_str());

  int Nc   = CommonParameters::Nc();
  int Nd   = CommonParameters::Nd();
  m_NinF = 2 * Nc * Nd;

  m_Nvol = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();

  m_foprdw = new AFopr_Domainwall_General<AFIELD>(params);

  set_parameters(params);

}

//====================================================================
template<typename AFIELD>
void AFopr_OptimalDomainwall<AFIELD>::tidyup()
{
  delete m_foprdw;
}

//====================================================================
template<typename AFIELD>
void AFopr_OptimalDomainwall<AFIELD>::set_parameters(
                                         const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  string str_gmset_type;
  double mq, M0;
  int    Ns;
  std::vector<int> bc;
  double b, c;
  double lambda_min, lambda_max;

  int err_optional = 0;
  err_optional += params.fetch_string("gamma_matrix_type", str_gmset_type);
  int err = 0;
  err += params.fetch_double("quark_mass", mq);
  err += params.fetch_double("domain_wall_height", M0);
  err += params.fetch_int("extent_of_5th_dimension", Ns);
  err += params.fetch_int_vector("boundary_condition", bc);
  err += params.fetch_double("coefficient_b", b);
  err += params.fetch_double("coefficient_c", c);
  err += params.fetch_double("lower_bound", lambda_min);
  err += params.fetch_double("upper_bound", lambda_max);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(mq, M0, Ns, bc, b, c, lambda_min, lambda_max);

}

//====================================================================
template<typename AFIELD>
void AFopr_OptimalDomainwall<AFIELD>::set_parameters(
                                const double mq,
                                const double M0,
                                const int Ns,
                                const vector<int> bc,
                                const double b,
                                const double c,
                                const double lambda_min,
                                const double lambda_max)
{

  m_M0 = M0;
  m_mq = mq;
  m_Ns = Ns;

  m_boundary.resize(m_Ndim);
  assert(bc.size() == m_Ndim);
  for(int mu = 0; mu < m_Ndim; ++mu){
    m_boundary[mu] = bc[mu];
  }

  vout.general(m_vl, "Parameters specific to %s:\n", class_name.c_str());
  vout.general(m_vl, "  coefficient_b = %8.4f\n", b);
  vout.general(m_vl, "  coefficient_c = %8.4f\n", c);
  vout.general(m_vl, "  lower_bound   = %12.8f\n", lambda_min);
  vout.general(m_vl, "  upper_bound   = %12.8f\n", lambda_max);

  m_b.resize(m_Ns);
  m_c.resize(m_Ns);
  set_optimalDomainwall(b, c, lambda_min, lambda_max);

  m_foprdw->set_coefficients(m_b, m_c);

}

//====================================================================
template<typename AFIELD>
void AFopr_OptimalDomainwall<AFIELD>::set_optimalDomainwall(
                                       const double b,
                                       const double c,
                                       const double lambda_min,
                                       const double lambda_max)
{
  vector<double> omega(m_Ns);
  double bmax = lambda_max/lambda_min;
  double delta;

  Math_Sign_Zolotarev_Omega sign_func(m_Ns, bmax);
  sign_func.get_sign_prms(omega, delta);

  for(int is = 0; is < m_Ns; ++is){
    omega[is] = omega[is]/lambda_min;
  }

  for(int k = 0; k < m_Ns; ++k){
    vout.general(m_vl, "    %2d  %20.14f\n", k, omega[k]);
  }
  vout.general("  max deviation from 1 = %20.14f\n",delta);

  for(int is = 0; is < m_Ns; ++is){
    m_b[is] = 0.5 * (b * omega[is] + c);
    m_c[is] = 0.5 * (b * omega[is] - c);
  }

}

void AFopr_OptimalDomainwall<AFIELD>::set_mode(std::string mode)
{
#pragma omp barrier
  int ith = ThreadManager::get_thread_id();
  if(ith == 0) m_mode = mode;
  m_foprdw->set_mode(mode);
#pragma omp barrier   // redundant as the set_mode above has barriers 

}

//============================================================END=====
