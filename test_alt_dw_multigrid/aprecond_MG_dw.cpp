/*!

        @file    $Id: apreccond_MG_dw.cpp #$

        @brief   multi grid peconditionor with afield (QXS version)

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision: 2578 $

 */

//====================================================================

#ifdef USE_ALT_QXS
#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Field/aindex_lex.h"
//#include "lib_alt_QXS/Fopr/afopr_Clover_coarse.h"
typedef AField<float, QXS>    AField_f;
typedef AField<double, QXS>   AField_d;
#endif

#ifdef USE_ALT_ACCEL
#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/afield-inc.h"
#include "lib_alt_Accel/Field/aindex_lex.h"
//#include "lib_alt_Accel/Fopr/afopr_Clover_coarse.h"
typedef AField<float,  ACCEL>  AField_f;
typedef AField<double, ACCEL>  AField_d;
#endif

//#include "lib_alt/Solver/aprecond_MG-tmpl.h"
#include "aprecond_MG_dw-tmpl.h"

// explicit instanciation.


#ifdef USE_ALT_QXS
template<>
const std::string APrecond_MG_dw<AField_d, AField_f>::class_name
    = "APrecond_MG_dw<AField<double,QXS>, AField<float,QXS> >";
#endif

#ifdef USE_ALT_ACCEL
template<>
const std::string APrecond_MG_dw<AField_d, AField_f>::class_name
    = "APrecond_MG_dw<AField<double,ACCEL>, AField<float,ACCEL> >";
#endif

template class APrecond_MG_dw<AField_d, AField_f>;


//============================================================END=====
