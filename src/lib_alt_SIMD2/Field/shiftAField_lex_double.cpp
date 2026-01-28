/*!
      @file    shiftAField_lex_double.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
      @version $LastChangedRevision: 2160 $
*/


#include "lib_alt_SIMD2/Field/shiftAField_lex.h"

#include <string>

#include "lib/ResourceManager/threadManager_OpenMP.h"

#include "lib_alt_SIMD2/inline/define_params.h"

#define  VLEN    VLEND
#define  VLEN2   VLEND2

#include "lib_alt_SIMD2/Field/index_lex_alt.h"
#include "lib_alt_SIMD2/Field/afield.h"

#include "lib_alt_SIMD2/inline/define_vlen.h"

// inline function files
#include "lib_alt_SIMD2/inline/vsimd_double-inc.h"
#include "lib_alt_SIMD2/inline/vsimd_common_double-inc.h"
#include "lib_alt_SIMD2/inline/afield_th-inc.h"

// function template files
#include "lib_alt_SIMD2/Field/afield-inc.h"

#include "lib_alt_SIMD2/Field/shiftAField_lex-tmpl.h"

template<>
const std::string ShiftAField_lex<AField<double,SIMD2> >::class_name
                              = "ShiftAField_lex<AField<double,SIMD2> >";

// explicit instanciation.
template class ShiftAField_lex<AField<double,SIMD2> >;


//============================================================END=====
