/*!
      @file    afopr_Precond_dd-tmpl.h
      @brief
      @author  Issaku Kanamori (kanamori)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2604 $
*/

#include "afopr_Precond_dd.h"


template<typename AFIELD>
const std::string AFopr_Precond_dd<AFIELD>::class_name
                                                = "AFopr_Precond_dd";

//====================================================================
template<typename AFIELD>
void AFopr_Precond_dd<AFIELD>::init(AFopr_dd<AFIELD>* afopr,
				    ASolver<AFIELD>* asolver,
				    const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  } else {
    m_vl = CommonParameters::Vlevel();
  }

  vout.general(m_vl, "%s: construction\n", class_name.c_str());
  vout.increase_indent();

  m_afopr = afopr;
  m_asolver = asolver;

  set_parameters(params);

  int Nin  = m_afopr->field_nin();
  int Nvol = m_afopr->field_nvol();
  int Nex  = m_afopr->field_nex();

  m_v1.reset(Nin, Nvol, Nex);
  m_w1.reset(Nin, Nvol, Nex);
  m_w2.reset(Nin, Nvol, Nex);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());
}


//====================================================================
template<typename AFIELD>
void AFopr_Precond_dd<AFIELD>::tidyup()
{
  // do nothing
}


//====================================================================
template<typename AFIELD>
void AFopr_Precond_dd<AFIELD>::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);
}


//====================================================================
template<typename AFIELD>
void AFopr_Precond_dd<AFIELD>::set_config(Field *u)
{
  vout.crucial(m_vl, "%s: set_config is called\n", class_name.c_str());
  exit(EXIT_FAILURE);
}



//====================================================================
template<typename AFIELD>
void AFopr_Precond_dd<AFIELD>::mult_gm5(AFIELD& v, const AFIELD& w)
{
  m_afopr->mult_gm5(v, w);
}


//====================================================================
template<typename AFIELD>
void AFopr_Precond_dd<AFIELD>::mult_dup(AFIELD& v, const AFIELD& w, int mu)
{
  int nconv;
  real_t diff;
  m_asolver->solve(m_w1, w, nconv, diff);

  m_afopr->mult_dup(v, m_w1, mu);
}

//====================================================================
template<typename AFIELD>
void AFopr_Precond_dd<AFIELD>::mult_ddn(AFIELD& v, const AFIELD& w, int mu)
{
  int nconv;
  real_t diff;
  m_asolver->solve(m_w1, w, nconv, diff);

  m_afopr->mult_ddn(v, m_w1, mu);
}


//====================================================================
template<typename AFIELD>
void AFopr_Precond_dd<AFIELD>::set_mode(std::string mode)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_mode = mode;

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
std::string AFopr_Precond_dd<AFIELD>::get_mode() const
{
  return m_mode;
}


//====================================================================
template<typename AFIELD>
void AFopr_Precond_dd<AFIELD>::mult(AFIELD& v, const AFIELD& w)
{
  if (m_mode == "D") {
    return D(v, w);
  } else if (m_mode == "DdagD") {
    D(m_v1, w);
    Ddag(v, m_v1);
    return;
  } else if (m_mode == "Ddag") {
    return Ddag(v, w);
  } else {
    vout.crucial(m_vl, "%s: mode undefined.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Precond_dd<AFIELD>::mult(AFIELD& v, const AFIELD& w,
                                                    const string mode)
{
  if (mode == "D") {
    return D(v, w);
  } else if (mode == "DdagD") {
    D(m_v1, w);
    Ddag(v, m_v1);
    return;
  } else if (mode == "Ddag") {
    return Ddag(v, w);
  } else {
    vout.crucial(m_vl, "%s: mode undefined.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Precond_dd<AFIELD>::mult_dag(AFIELD& v, const AFIELD& w)
{
  if (m_mode == "D") {
    return Ddag(v, w);
  } else if (m_mode == "DdagD") {
    D(m_v1, w);
    Ddag(v, m_v1);
    return;
  } else if (m_mode == "Ddag") {
    return D(v, w);
  } else {
    vout.crucial(m_vl, "%s: mode undefined.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
}

//====================================================================
template<typename AFIELD>
void AFopr_Precond_dd<AFIELD>::mult_sap(AFIELD& v, const AFIELD& w, const int ieo){
  if (m_mode == "D") {
    mult_D_sap(v, w, ieo);
  } else if (m_mode == "Ddag") {
    mult_Ddag_sap(v, w, ieo);
  } else if (m_mode == "DdagD") {
    mult_D_sap(m_v1, w, ieo);
    mult_Ddag_sap(v, m_v1, ieo);
  } else if (m_mode == "DDdag") {
    mult_Ddag_sap(m_v1, w, ieo);
    mult_D_sap(v, m_v1, ieo);
  } else {
    std::cout << "mode undeifined for mult_sap in AFopr_Domainwall_5din_dd.\n";
    abort();
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Precond_dd<AFIELD>::mult_dd(AFIELD& v, const AFIELD& w){
  int ieo = -1;
  if (m_mode == "D") {
    mult_D_sap(v, w, ieo);
  } else if (m_mode == "Ddag") {
    mult_Ddag_sap(v, w, ieo);
  } else if (m_mode == "DdagD") {
    mult_D_sap(m_v1, w, ieo);
    mult_Ddag_sap(v, m_v1, ieo);
  } else if (m_mode == "DDdag") {
    mult_Ddag_sap(m_v1, w, ieo);
    mult_D_sap(v, m_v1, ieo);
  } else {
    std::cout << "mode undeifined for mult_sap in AFopr_Domainwall_5din_dd.\n";
    abort();
  }

}


//====================================================================
template<typename AFIELD>
void AFopr_Precond_dd<AFIELD>::D(AFIELD& v, const AFIELD& w)
{
  /*
  m_afopr->set_mode("D");
  m_afopr->mult(m_w1, w);

  int nconv;
  real_t diff;
  m_asolver->get_fopr()->set_mode("D");
  m_asolver->solve(v, m_w1, nconv, diff);
  */

  int nconv;
  real_t diff;
  m_asolver->get_fopr()->set_mode("D");
  m_asolver->solve(m_w1, w, nconv, diff);

  m_afopr->set_mode("D");
  m_afopr->mult(v, m_w1);

}


//====================================================================
template<typename AFIELD>
void AFopr_Precond_dd<AFIELD>::Ddag(AFIELD& v, const AFIELD& w)
{

  int nconv;
  real_t diff;
  m_asolver->get_fopr()->set_mode("Ddag");
  m_asolver->solve(m_w1, w, nconv, diff);

  m_afopr->set_mode("Ddag");
  m_afopr->mult(v, m_w1);

}

//====================================================================
template<typename AFIELD>
void AFopr_Precond_dd<AFIELD>::mult_D_sap(AFIELD& v, const AFIELD& w,
                                                          int ieo)
{
  m_afopr->set_mode("D");
  if(ieo < 0){
    m_afopr->mult_dd(m_w1, w);
  } else {
    m_afopr->mult_sap(m_w1, w, ieo);
  }
  int nconv;
  real_t diff;
  m_asolver->get_fopr()->set_mode("D");
  m_asolver->solve(v, m_w1, nconv, diff);
}


//====================================================================
template<typename AFIELD>
void AFopr_Precond_dd<AFIELD>::mult_Ddag_sap(AFIELD& v, const AFIELD& w,
                                                          int ieo)
{
  int nconv;
  real_t diff;
  m_asolver->get_fopr()->set_mode("D");
  m_asolver->solve(m_w1, w, nconv, diff);

  m_afopr->set_mode("D");
  if(ieo < 0){
    m_afopr->mult_dd(v, m_w1);
  } else {
    m_afopr->mult_sap(v, m_w1, ieo);
  }
}

//====================================================================
template<typename AFIELD>
double AFopr_Precond_dd<AFIELD>::flop_count(const std::string mode)
{

  double flop = 0.0;

  vout.general("Warning: %s: flop_count is not ready.\n",
               class_name.c_str());
  return flop;
}

//====================================================================
template<typename AFIELD>
double AFopr_Precond_dd<AFIELD>::flop_count_sap()
{

  double flop = 0.0;

  vout.general("Warning: %s: flop_count_sap is not ready.\n",
               class_name.c_str());
  return flop;
}


//============================================================END=====
