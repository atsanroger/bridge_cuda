/*!
        @file    afopr_Smeared-tmpl.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-04-17 11:32:37 #$

        @version $LastChangedRevision: 2513 $
*/

#include "Fopr/afopr_Smeared.h"

//#ifdef USE_FACTORY_AUTOREGISTER
//namespace {
//  bool init = AFopr_Smeared::register_factory();
//}
//#endif

template<typename AFIELD>
const std::string AFopr_Smeared<AFIELD>::class_name = "AFopr_Smeared";

//====================================================================
template<typename AFIELD>
void AFopr_Smeared<AFIELD>::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Smeared<AFIELD>::get_parameters(Parameters& params) const
{
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
template<typename AFIELD>
void AFopr_Smeared<AFIELD>::set_config(Field *U)
{
  m_dr_smear->set_config(U);

  const int Nsmear = m_dr_smear->get_Nsmear();
  Field     *Uptr  = m_dr_smear->getptr_smearedConfig(Nsmear);

  m_fopr->set_config(Uptr);
}


//====================================================================
template<typename AFIELD>
double AFopr_Smeared<AFIELD>::flop_count()
{
  double flop_fopr = m_fopr->flop_count();

  return flop_fopr;
}


//============================================================END=====
