/*!
        @file    aforce_F_Clover_Nf2.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-10-11 19:12:24 #$
        @version $LastChangedRevision: 2545 $
*/

#include "lib_alt/Force/Fermion/aforce_F_Clover_Nf2.h"

#include<cassert>

// include files in core library
#include "lib/ResourceManager/threadManager.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt_QXS/inline/define_params.h"

#define  VLEN     VLEND
#define  VLENX    VLENXD
#define  VLENY    VLENYD

typedef double real_t;

#include "lib_alt_QXS/inline/vsimd_double-inc.h"
#include "lib_alt_QXS/inline/vsimd_common_double-inc.h"
#include "lib_alt_QXS/inline/vsimd_Wilson_SU3_double-inc.h"

// include files in alt-code dorectories
#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/aindex_lex.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Field/afield_Spinor-inc.h"
namespace Alt_Spinor = QXS_Spinor;
#include "lib_alt_QXS/Fopr/afopr_Clover.h"


#include "lib_alt/Force/Fermion/aforce_F_Clover_Nf2-tmpl.h"

// explicit instanciation for AField<double,QXS>.
template<>
const std::string AForce_F_Clover_Nf2<AField<double,QXS> >::class_name
                         = "AForce_F_Clover_Nf2<AField<double,QXS> >";


template class AForce_F_Clover_Nf2<AField<double,QXS> >;

//============================================================END=====
