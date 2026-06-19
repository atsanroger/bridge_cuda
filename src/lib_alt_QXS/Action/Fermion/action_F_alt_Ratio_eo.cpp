/*!
        @file    action_F_alt_Ratio_eo.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-05-28 16:32:14 #$
        @version $LastChangedRevision: 2520 $
*/

#include "lib_alt/Action/Fermion/action_F_alt_Ratio_eo.h"

#include<cassert>

// include files in core library
#include "lib/ResourceManager/threadManager.h"
#include "lib/Field/index_eo.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/aindex_lex.h"
#include "lib_alt_QXS/Field/aindex_eo.h"
#include "lib_alt_QXS/Field/aindex_eo-inc.h"
#include "lib_alt_QXS/Field/afield-inc.h"

#include "lib_alt/Action/Fermion/action_F_alt_Ratio_eo-tmpl.h"

// explicit instanciation for AField<double,QXS>.
template<>
const std::string Action_F_alt_Ratio_eo<AField<double,QXS> >::
       class_name = "Action_F_alt_Ratio_eo<Afield<double,QXS> >";


template class Action_F_alt_Ratio_eo<AField<double,QXS> >;

//============================================================END=====
