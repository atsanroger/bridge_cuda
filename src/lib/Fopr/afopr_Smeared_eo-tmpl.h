/*!
        @file    afopr_Smeared_eo-tmpl.h

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "Fopr/afopr_Smeared_eo.h"

//#ifdef USE_FACTORY_AUTOREGISTER
//namespace {
//  bool init = Fopr_Smeared_eo::register_factory();
//}
//#endif

template<typename AFIELD>
const std::string AFopr_Smeared_eo<AFIELD>::class_name = "Fopr_Smeared_eo";

//====================================================================
template<typename AFIELD>
void AFopr_Smeared_eo<AFIELD>::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Smeared_eo<AFIELD>::get_parameters(Parameters& params) const
{
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
template<typename AFIELD>
void AFopr_Smeared_eo<AFIELD>::set_config(Field *U)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  m_dr_smear->set_config(U);

  const int Nsmear = m_dr_smear->get_Nsmear();
  Field     *Uptr  = m_dr_smear->getptr_smearedConfig(Nsmear);

  m_fopr_eo->set_config(Uptr);
}


//====================================================================
//============================================================END=====
