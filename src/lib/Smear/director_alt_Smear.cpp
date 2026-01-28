/*!
      @file    director_alt_Smear.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-04-17 11:32:37 #$
      @version $LastChangedRevision: 2513 $
*/

#ifdef USE_ALT_CODE

#include "lib_alt/Smear/director_alt_Smear.h"

#include "lib/ResourceManager/threadManager.h"
#include "lib/Fopr/afopr_Smeared.h"

#include "lib/Field/field.h"
//#include "lib/Field/index_lex.h"


template<>
const std::string Director_alt_Smear<Field>::class_name
                         = "Director_alt_Smear<Field>";

//====================================================================
template<>
void Director_alt_Smear<Field>::get_config(Field& Usmr)
{
  if(m_Nsmear == 0){
    copy(Usmr, *m_U);
  }else{
    copy(Usmr, m_Usmear[m_Nsmear-1]);
  }

}

//====================================================================
template<>
void Director_alt_Smear<Field>::get_config(Field& Usmr, const int ismr)
{
  if(m_Nsmear == 0){
    copy(Usmr, *m_U);
  }else{
    copy(Usmr, m_Usmear[ismr-1]);
  }

}

//====================================================================
// class instanciation for Field

#define CORELIB_INSTANCE
#include "lib_alt/Smear/director_alt_Smear-tmpl.h"
#undef  CORELIB_INSTANCE

template class Director_alt_Smear<Field>;

#endif
//============================================================END=====
