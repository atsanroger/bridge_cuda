/*!
        @file    aforce_F_Domainwall.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-06-29 22:43:06 #$
        @version $LastChangedRevision: 2527 $
*/

#include "lib_alt/Force/Fermion/aforce_F_Domainwall.h"

#include<cassert>

// include files in core library
#include "lib/ResourceManager/threadManager.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

typedef double real_t;

#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/aindex_lex.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Fopr/afopr_Wilson.h"


#include "lib_alt/Force/Fermion/aforce_F_Domainwall-tmpl.h"

// explicit instanciation for AField<double,QXS>.
template<>
const std::string AForce_F_Domainwall<AField<double,QXS> >::class_name
                         = "AForce_F_Domainwall<AField<double,QXS> >";


template class AForce_F_Domainwall<AField<double,QXS> >;

//============================================================END=====
