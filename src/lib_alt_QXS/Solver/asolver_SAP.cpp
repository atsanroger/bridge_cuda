/*!
        @file    asolver_SAP.cpp
        @brief   SAP solver (QXS version)
        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate: 2023-03-20 10:52:44 +0900 (Mon, 20 Mar 2023) $
        @version $LastChangedRevision: 2499 $
*/
//====================================================================
#include "lib_alt/Solver/asolver_SAP.h"

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"

//#define DEBUG

#include "lib_alt/Solver/asolver_SAP-tmpl.h"


//====================================================================
// explicit instanciation
template<>
const std::string ASolver_SAP<AField<double, QXS> >::class_name
  = "ASolver_SAP<AField<double,QXS> >";

template class ASolver_SAP<AField<double, QXS> >;

//====================================================================
// explicit instanciation
template<>
const std::string ASolver_SAP<AField<float, QXS> >::class_name
  = "ASolver_SAP<AField<float,QXS> >";

template class ASolver_SAP<AField<float, QXS> >;

//============================================================END=====
