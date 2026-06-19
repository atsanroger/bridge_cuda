/*!
        @file    aindex_eo_double.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-07-23 18:49:38 #$
        @version $LastChangedRevision: 2534 $
*/

#include "lib_alt_Accel/Field/aindex_eo.h"

#include <assert.h>

#include "lib_alt_Accel/Field/afield.h"

#include "lib_alt_Accel/BridgeACC/bridgeACC_AField.h"
#include "lib_alt_Accel/BridgeACC/bridgeACC_Index_eo_alt.h"

typedef double real_t;


// class template
#include "lib_alt_Accel/Field/aindex_eo-tmpl.h"

template<>
const std::string AIndex_eo<double, ACCEL>::class_name
                                 = "AIndex_eo<double, ACCEL>";

// explicit instanciation.
template class AIndex_eo<double, ACCEL>;

//============================================================END=====
