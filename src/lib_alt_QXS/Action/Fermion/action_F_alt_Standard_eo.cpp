/*!
        @file    action_F_alt_Standard_eo.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-05-23 16:09:22 #$
        @version $LastChangedRevision: 2517 $
*/

#include "lib_alt/Action/Fermion/action_F_alt_Standard_eo.h"

#include<cassert>

// include files in core library
#include "lib/ResourceManager/threadManager.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/aindex_lex.h"
#include "lib_alt_QXS/Field/aindex_eo.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Field/aindex_eo-inc.h"

#include "lib_alt/Action/Fermion/action_F_alt_Standard_eo-tmpl.h"

// explicit instanciation for AField<double,QXS>.
template<>
const std::string Action_F_alt_Standard_eo<AField<double,QXS> >::
       class_name = "Action_F_alt_Standard_eo<Afield<double,QXS> >";


template class Action_F_alt_Standard_eo<AField<double,QXS> >;

//============================================================END=====
