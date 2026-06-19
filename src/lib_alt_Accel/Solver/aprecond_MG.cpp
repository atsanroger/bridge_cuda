/*!
        @file    apreccond_MG.cpp
        @brief   multi grid peconditionor with afield (Accel version)
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate::  $
        @version $LastChangedRevision: 2422 $
*/

#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/afield-inc.h"
#include "lib_alt_Accel/Field/aindex_lex.h"
#include "lib_alt_Accel/Fopr/afopr_Clover_coarse.h"

//#define DEBUG

#include "lib_alt/Solver/aprecond_MG-tmpl.h"

// explicit instanciation.

template<>
const std::string APrecond_MG<AField<double,ACCEL>, AField<float,ACCEL> >
::class_name = "APrecond_MG<AField<double,ACCEL>, AField<float,ACCEL> >";
template class APrecond_MG<AField<double,ACCEL>, AField<float,ACCEL> >;

template<>
const std::string APrecond_MG<AField<double,ACCEL>, AField<double,ACCEL> >
::class_name = "APrecond_MG<AField<double,ACCEL>, AField<double,ACCEL> >";
template class APrecond_MG<AField<double,ACCEL>, AField<double,ACCEL> >;

//============================================================END=====
