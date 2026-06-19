/*!                                                                             
        @file    afopr_Chebyshev.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
        @version $LastChangedRevision: 2160 $
*/

#include "lib_alt/Fopr/afopr_Chebyshev.h"

#include <cassert>

#include "lib/Field/field.h"

#include "lib_alt_SIMD2/Field/afield.h"
#include "lib_alt_SIMD2/Field/afield-inc.h"


#ifdef USE_FACTORY
namespace {
  AFopr<AField<double,SIMD2> > *create_object_1(
                                      AFopr<AField<double,SIMD2> > *fopr,
                                      const Parameters& params)
  {
    return new AFopr_Chebyshev<AField<double,SIMD2> >(fopr, params);
  }
  bool init1 = AFopr<AField<double,SIMD2> >::Factory_fopr_params::Register(
                        "Chebyshev", create_object_1);

  AFopr<AField<float,SIMD2> > *create_object_2(
                                      AFopr<AField<float,SIMD2> > *fopr,
                                      const Parameters& params)
  {
    return new AFopr_Chebyshev<AField<float,SIMD2> >(fopr, params);
  }
  bool init2 = AFopr<AField<float,SIMD2> >::Factory_fopr_params::Register(
                        "Chebyshev", create_object_2);

}
#endif


#include "lib_alt/Fopr/afopr_Chebyshev-tmpl.h"


// class instanciation.
template class AFopr_Chebyshev<AField<float,SIMD2> >;
template class AFopr_Chebyshev<AField<double,SIMD2> >;

//============================================================END=====
