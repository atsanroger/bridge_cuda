/*!
        @file    afopr_Smeared_alt-tmpl.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-04-17 11:32:37 #$

        @version $LastChangedRevision: 2513 $
*/

#include "lib/Fopr/afopr_Smeared_alt.h"

#include "lib/ResourceManager/threadManager.h"

//#ifdef USE_FACTORY_AUTOREGISTER
//namespace {
//  bool init = AFopr_Smeared_alt::register_factory();
//}
//#endif

template<typename AFIELD>
const std::string AFopr_Smeared_alt<AFIELD>::class_name = "AFopr_Smeared_alt";

//====================================================================
template<typename AFIELD>
void AFopr_Smeared_alt<AFIELD>::init()
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: construction\n", class_name.c_str());

  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}

//====================================================================
template<typename AFIELD>
void AFopr_Smeared_alt<AFIELD>::init(const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: construction\n", class_name.c_str());
  vout.increase_indent();

  set_parameters(params);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}


//====================================================================
template<typename AFIELD>
void AFopr_Smeared_alt<AFIELD>::tidyup()
{
  // do nothing.
}


//====================================================================
template<typename AFIELD>
void AFopr_Smeared_alt<AFIELD>::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Smeared_alt<AFIELD>::get_parameters(Parameters& params) const
{
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
template<typename AFIELD>
void AFopr_Smeared_alt<AFIELD>::set_config(Field *U)
{
  m_dr_smear->set_config(U);

  int Nsmear  = m_dr_smear->get_Nsmear();
  Field* Uptr = m_dr_smear->getptr_smearedConfig(Nsmear);

  m_fopr->set_config(Uptr);
}


//====================================================================
template<typename AFIELD>
double AFopr_Smeared_alt<AFIELD>::flop_count()
{
  double flop_fopr = m_fopr->flop_count();

  return flop_fopr;
}

//============================================================END=====
