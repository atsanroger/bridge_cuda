/*!
        @file    asolver_SAP_dw2.cpp
        @brief   SAP solver (QXS version)
        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate: 2024-03-17 22:33:12 +0900 (2024年03月17日 (日)) $
        @version $LastChangedRevision: 2586 $
*/
//====================================================================
#include "asolver_SAP_dw2.h"

#include "lib/ResourceManager/threadManager.h"

#ifdef USE_ALT_QXS
#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"
typedef AField<double, QXS>  AField_d;
typedef AField<float, QXS>   AField_f;
#endif

#ifdef USE_ALT_ACCEL
#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/afield-inc.h"
typedef AField<double, ACCEL>  AField_d;
typedef AField<float, ACCEL>   AField_f;
#endif

//#define DEBUG

#include "asolver_SAP_MINRES_dw2.h"

#include "asolver_SAP_dw2-tmpl.h"

//====================================================================
// explicit instanciation


#ifdef USE_ALT_QXS
template<>
const std::string ASolver_SAP_dw2<AField_d>::class_name
                         = "ASolver_SAP_dw2<AField<double,QXS> >";
#endif

#ifdef USE_ALT_ACCEL
template<>
const std::string ASolver_SAP_dw2<AField_d>::class_name
                       = "ASolver_SAP_dw2<AField<double,ACCEL> >";
#endif

template class ASolver_SAP_dw2<AField_d>;

//====================================================================
// explicit instanciation
#ifdef USE_ALT_QXS
template<>
const std::string ASolver_SAP_dw2<AField_f>::class_name
                           = "ASolver_SAP_dw2<AField<float,QXS> >";
#endif

#ifdef USE_ALT_ACCEL
template<>
const std::string ASolver_SAP_dw2<AField_f>::class_name
                         = "ASolver_SAP_dw2<AField<float,ACCEL> >";
#endif

template class ASolver_SAP_dw2<AField_f>;

//============================================================END=====
