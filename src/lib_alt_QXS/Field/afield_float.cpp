/*!
        @file    afield_float.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2024-01-13 22:19:26 #$
        @version $LastChangedRevision: 2569 $
*/

#include "lib_alt_QXS/Field/afield.h"

#include <cassert>

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_QXS/inline/define_vlen.h"

#define  VLEN     VLENS
#define  VLENX    VLENXS
#define  VLENY    VLENYS

typedef float real_t;

#include "lib_alt_QXS/inline/vsimd_float-inc.h"
#include "lib_alt_QXS/inline/vsimd_common_float-inc.h"

// template definition
#include "lib_alt_QXS/Field/afield-tmpl.h"

template<>
const std::string AField<float, QXS>::class_name = "AField<float, QXS>";

// explicit instanciation.
template class AField<float, QXS>;

//============================================================END=====
