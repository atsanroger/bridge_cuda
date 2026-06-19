/*!
      @file    director_alt_Smear.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-12-24 19:34:22 #$
      @version $LastChangedRevision: 2562 $
*/

#include "lib_alt/Smear/director_alt_Smear.h"

#include "lib/ResourceManager/threadManager.h"
#include "lib/Fopr/afopr_Smeared.h"

#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Field/aindex_lex.h"

#include "lib_alt/Smear/director_alt_Smear-tmpl.h"

//====================================================================
// explicit instanciation for AField<double,QXS>.
template<>
const std::string Director_alt_Smear<AField<double,QXS> >::class_name
                         = "Director_alt_Smear<AField<double,QXS> >";

template class Director_alt_Smear<AField<double,QXS> >;

template<>
const std::string Director_alt_Smear<AField<float,QXS> >::class_name
                         = "Director_alt_Smear<AField<float,QXS> >";

template class Director_alt_Smear<AField<float,QXS> >;

#ifdef USE_QXS_FP16

template<>
const std::string Director_alt_Smear<AField<half,QXS> >::class_name
                         = "Director_alt_Smear<AField<half,QXS> >";

template class Director_alt_Smear<AField<half,QXS> >;

#endif
//============================================================END=====
