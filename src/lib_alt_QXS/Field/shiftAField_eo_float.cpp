/*!
      @file    shiftAField_eo_float.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
      @version $LastChangedRevision: 2499 $
*/


#include "lib_alt_QXS/Field/shiftAField_eo.h"

#include <string>

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt_QXS/inline/define_params.h"

#define  VLEN     VLENS
#define  VLENX    VLENXS
#define  VLENY    VLENYS

typedef float real_t;

#include "lib_alt_QXS/Field/aindex_eo.h"
#include "lib_alt_QXS/Field/afield.h"

// inline function files
#include "lib_alt_QXS/inline/vsimd_float-inc.h"
#include "lib_alt_QXS/inline/vsimd_common_float-inc.h"
#include "lib_alt_QXS/inline/afield_th-inc.h"

// function template files
#include "lib_alt_QXS/Field/afield-inc.h"

#include "lib_alt_QXS/Field/shiftAField_eo-tmpl.h"

template<>
const std::string ShiftAField_eo<AField<float, QXS> >::class_name
  = "ShiftAField_eo<AField<float,QXS> >";

// explicit instanciation.
template class ShiftAField_eo<AField<float, QXS> >;

//============================================================END=====
