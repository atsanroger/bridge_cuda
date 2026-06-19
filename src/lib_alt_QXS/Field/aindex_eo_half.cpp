/*!
      @file    aindex_eo_alt_half.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-12-24 19:34:22 #$
      @version $LastChangedRevision: 2562 $
*/
#ifdef USE_QXS_FP16

#include <assert.h>

#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt_QXS/Field/aindex_eo.h"

#define  VLEN     VLENH
#define  VLENX    VLENXH
#define  VLENY    VLENYH

#include "lib_alt_QXS/Field/afield.h"

#include "lib_alt_QXS/Field/aindex_eo-tmpl.h"

// explicit instanciation.
template class AIndex_eo<half, QXS>;

#endif
//============================================================END=====
