// Time-stamp: <2021-03-24 16:14:58 kanamori>

/*!

        @file    $Id: apreccond_MG.cpp #$

        @brief   multi grid peconditionor with afield (QXS version)

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision: 2499 $

 */

//====================================================================

#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Field/aindex_lex.h"
#include "lib_alt_QXS/Fopr/afopr_Clover_coarse.h"

#include "lib_alt/Solver/aprecond_MG-tmpl.h"

// explicit instanciation.


template<>
const std::string APrecond_MG<AField<double, QXS>, AField<float, QXS> >::class_name
  = "APrecond_MG<AField<double,QXS>, AField<float,QXS> >";
template class APrecond_MG<AField<double, QXS>, AField<float, QXS> >;


//============================================================END=====
