/*!
      @file    afield_double.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
      @version $LastChangedRevision: 2160 $
*/

#include "lib_alt_SIMD2/Field/afield.h"

#include <cassert>

#include "lib/ResourceManager/threadManager_OpenMP.h"

#include "lib_alt_SIMD2/Field/index_lex_alt.h"

#define  VLEN    VLEND
#define  VLEN2   VLEND2

#include "lib_alt_SIMD2/inline/afield_th-inc.h"

#include "lib_alt_SIMD2/inline/vsimd_double-inc.h"
#include "lib_alt_SIMD2/inline/vsimd_common_double-inc.h"

#include "lib_alt_SIMD2/Field/afield-tmpl.h"

template<>
const std::string AField<double,SIMD2>::class_name = "AField<double,SIMD2>";

// explicit instanciation.
template class AField<double,SIMD2>;

//============================================================END=====
