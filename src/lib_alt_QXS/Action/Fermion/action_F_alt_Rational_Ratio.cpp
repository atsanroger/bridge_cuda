/*!
        @file    action_F_alt_Rational_Ratio.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-06-29 22:43:06 #$
        @version $LastChangedRevision: 2527 $
*/

#include "lib_alt/Action/Fermion/action_F_alt_Rational_Ratio.h"

#include<cassert>

// include files in core library
#include "lib/ResourceManager/threadManager.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/aindex_lex.h"
#include "lib_alt_QXS/Field/afield-inc.h"

#include "lib_alt/Action/Fermion/action_F_alt_Rational_Ratio-tmpl.h"

// explicit instanciation for AField<double,QXS>.
template<>
const std::string Action_F_alt_Rational_Ratio<AField<double,QXS> >::
     class_name = "Action_F_alt_Rational_Ratio<AField<double,QXS> >";


template class Action_F_alt_Rational_Ratio<AField<double,QXS> >;

//============================================================END=====
