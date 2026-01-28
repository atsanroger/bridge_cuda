/*!                                                                             
        @file    afopr_Domainwall.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
        @version $LastChangedRevision: 2160 $
*/

#include "lib_alt/Fopr/afopr_Domainwall.h"

#include<cassert>

#include "lib/Field/field.h"

#include "lib_alt_SIMD2/Field/afield.h"
#include "lib_alt_SIMD2/Field/afield-inc.h"


#ifdef USE_FACTORY
namespace {
  AFopr<AField<double,SIMD2> > *create_object_with_params1(
                                      const Parameters& params)
  {
    return new AFopr_Domainwall<AField<double,SIMD2> >(params);
  }
  bool init1 = AFopr<AField<double,SIMD2> >::Factory_params::Register(
                        "Domainwall", create_object_with_params1);

  AFopr<AField<float,SIMD2> > *create_object_with_params2(
                                      const Parameters& params)
  {
    return new AFopr_Domainwall<AField<float,SIMD2> >(params);
  }
  bool init2 = AFopr<AField<float,SIMD2> >::Factory_params::Register(
                        "Domainwall", create_object_with_params2);

}
#endif

#include "lib_alt/Fopr/afopr_Domainwall-tmpl.h"

// explicit instanciation.
template class AFopr_Domainwall<AField<float,SIMD2> >;
template class AFopr_Domainwall<AField<double,SIMD2> >;

//============================================================END=====
