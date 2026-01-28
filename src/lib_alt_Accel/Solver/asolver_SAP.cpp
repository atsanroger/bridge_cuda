/*!
        @file    asolver_SAP.cpp
        @brief   SAP solver (Accel version)
        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate: 2023-08-20 14:25:12 +0900 (2023年08月20日 (日)) $
        @version $LastChangedRevision: 2535 $
*/

//====================================================================
#include "lib_alt/Solver/asolver_SAP.h"

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/afield-inc.h"

#define DEBUG

#include "lib_alt/Solver/asolver_SAP-tmpl.h"


//====================================================================
// explicit instanciation

template<>
const std::string ASolver_SAP<AField<float,ACCEL> >::class_name
                              = "ASolver_SAP<AField<float,ACCEL> >";


template class ASolver_SAP<AField<float,ACCEL> >;

//====================================================================
// explicit instanciation

template<>
const std::string ASolver_SAP<AField<double,ACCEL> >::class_name
                              = "ASolver_SAP<AField<double,ACCEL> >";


template class ASolver_SAP<AField<double,ACCEL> >;

//============================================================END=====
