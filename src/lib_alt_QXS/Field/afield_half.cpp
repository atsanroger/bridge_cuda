/*!
        @file    afield_double.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-12-24 19:34:22 #$
        @version $LastChangedRevision: 2562 $
*/
#ifdef USE_QXS_FP16

#include "lib_alt_QXS/Field/afield.h"


#include <cassert>

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_QXS/inline/define_vlen.h"

#define  VLEN   VLENH
#define  VLENX  VLENXH
#define  VLENY  VLENYH

typedef half real_t;

#include "lib_alt_QXS/inline/vsimd_half-inc.h"
#include "lib_alt_QXS/inline/vsimd_common_half-inc.h"

// template definition
#include "lib_alt_QXS/Field/afield-tmpl.h"

template<>
const std::string AField<real_t, QXS>::class_name = "AField<half, QXS>";

// specialized template for FP16
#ifdef USE_QXS_ACLE
#include "lib_alt_QXS/Field/afield_fp16_acle-tmpl.h"
#endif


// explicit instanciation.
template class AField<real_t, QXS>;

#endif
//============================================================END=====
