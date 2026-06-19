/*!
        @file    aindex_eo_float.cpp
        @brief
        @author  <Hideo Matsufuru> hideo.matsufuru@kek.jp(matsufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-10-02 13:50:31 #$
        @version $LastChangedRevision: 2543 $
*/

#include "lib_alt_Vector/Field/aindex_eo.h"

#include <assert.h>

#include "lib_alt_Vector/Field/afield.h"

typedef float real_t;

// kernel functions
#include "lib_alt_Vector/BridgeVec/bridgeVec_aindex_eo.h"

// class template
#include "lib_alt_Vector/Field/aindex_eo-tmpl.h"

template<>
const std::string AIndex_eo<float, VECTOR>::class_name
                                 = "AIndex_eo<float, VECTOR>";

// explicit instanciation.
template class AIndex_eo<float, VECTOR>;

//============================================================END=====
