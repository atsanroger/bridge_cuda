/*!
        @file    action_F_alt_Ratio_lex.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-04-17 11:32:37 #$
        @version $LastChangedRevision: 2513 $
*/

#include "lib_alt/Action/Fermion/action_F_alt_Ratio_lex.h"

#include<cassert>

// include files in core library
#include "lib/ResourceManager/threadManager.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/aindex_lex.h"
#include "lib_alt_QXS/Field/afield-inc.h"

#include "lib_alt/Action/Fermion/action_F_alt_Ratio_lex-tmpl.h"

// explicit instanciation for AField<double,QXS>.
template<>
const std::string Action_F_alt_Ratio_lex<AField<double,QXS> >::
       class_name = "Action_F_alt_Ratio_lex<Afield<double,QXS> >";


template class Action_F_alt_Ratio_lex<AField<double,QXS> >;

//============================================================END=====
