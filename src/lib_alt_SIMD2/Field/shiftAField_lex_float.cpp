/*!
      @file    shiftAField_lex_float.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
      @version $LastChangedRevision: 2160 $
*/


#include "Field/shiftAField_lex.h"

#include <string>

#include "lib/ResourceManager/threadManager_OpenMP.h"

#include "Field/index_lex_alt.h"
#include "Field/afield.h"

#include "inline/define_params.h"

#define  VLEN    VLENS
#define  VLEN2   VLENS2

// inline function files
#include "inline/afield_th-inc.h"

#include "inline/vsimd_float-inc.h"
#include "inline/vsimd_common_float-inc.h"

// function template files
#include "Field/afield-inc.h"

#include "Field/shiftAField_lex-tmpl.h"

template<>
const std::string ShiftAField_lex<AField<float,SIMD2> >::class_name
                              = "ShiftAField_lex<AField<float,SIMD2> >";

// explicit instanciation.
template class ShiftAField_lex<AField<float,SIMD2> >;

//============================================================END=====
