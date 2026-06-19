/*!
        @file    action_F_Standard_lex_alt.cpp
        @brief
        @author  Matsufuru, Hideo (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
        @version $LastChangedRevision: 2160 $
*/

#include "lib_alt/Action/Fermion/action_F_Standard_lex_alt.h"

#include "lib/Tools/timer.h"

#include "lib_alt_SIMD2/Field/index_lex_alt.h"
#include "lib_alt_SIMD2/Field/afield.h"
#include "lib_alt_SIMD2/Field/afield-inc.h"


#include "lib_alt/Action/Fermion/action_F_Standard_lex_alt.h"

template<>
const std::string Action_F_Standard_lex_alt<AField<double,SIMD2> >::class_name
                      = "Action_F_Standard_lex_alt<AField<double,SIMD2> >";

// class instanciation.
template class Action_F_Standard_lex_alt<AField<double,SIMD2> >;

//============================================================END=====
