/*!
      @file    aindex_eo_alt_double.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-12-24 19:34:22 #$
      @version $LastChangedRevision: 2562 $
*/

#include <assert.h>

#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt_QXS/Field/aindex_eo.h"

#define  VLEN     VLEND
#define  VLENX    VLENXD
#define  VLENY    VLENYD

#include "lib_alt_QXS/Field/afield.h"

#include "lib_alt_QXS/Field/aindex_eo-tmpl.h"

// explicit instanciation.
template class AIndex_eo<double, QXS>;
//============================================================END=====
