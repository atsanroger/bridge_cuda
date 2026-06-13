/*!
        @file    asolver_SAP.cpp
        @brief   SAP solver (QXS version)
        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate: 2024-04-01 12:13:09 +0900 (2024年04月01日 (月)) $
        @version $LastChangedRevision: 2595 $
*/
//====================================================================
#include "asolver_SAP_dw.h"


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

#include "asolver_SAP_dw-tmpl.h"


//====================================================================
// explicit instanciation

#ifdef USE_ALT_QXS
template<>
const std::string ASolver_SAP_dw<AField_d>::class_name
                         = "ASolver_SAP_dw<AField<double,QXS> >";
#endif

#ifdef USE_ALT_ACCEL
template<>
const std::string ASolver_SAP_dw<AField_d>::class_name
                       = "ASolver_SAP_dw<AField<double,ACCEL> >";
#endif

template class ASolver_SAP_dw<AField_d>;

//====================================================================
// explicit instanciation

#ifdef USE_ALT_QXS
template<>
const std::string ASolver_SAP_dw<AField_f>::class_name
                           = "ASolver_SAP_dw<AField<float,QXS> >";
#endif

#ifdef USE_ALT_ACCEL
template<>
const std::string ASolver_SAP_dw<AField_f>::class_name
                         = "ASolver_SAP_dw<AField<float,ACCEL> >";
#endif

template class ASolver_SAP_dw<AField_f>;

//============================================================END=====
