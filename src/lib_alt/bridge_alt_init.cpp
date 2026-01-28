/*!
        @file    bridge_init_factory_alt.cpp
        @brief
        @author  Hideo Matsufuru  (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-04-17 11:32:37 #$
        @version $LastChangedRevision: 2513 $
*/

#include "lib/Parameters/parameters.h"

// alt-code
#ifdef USE_ALT_QXS
#include "lib_alt_QXS/init_alt_QXS.h"
#endif

#ifdef USE_ALT_ACCEL
#include "lib_alt_Accel/init_alt_Accel.h"
#endif

#ifdef USE_ALT_SIMD2
#include "lib_alt_SIMD2/init_alt_SIMD2.h"
#endif

#ifdef USE_ALT_SIMD
#include "lib_alt_SIMD/init_alt_SIMD.h"
#endif

#ifdef USE_ALT_VECTOR
#include "lib_alt_Vector/init_alt_Vector.h"
#endif


bool bridge_alt_init(Parameters& params)
{
  bool result = true;

#ifdef USE_ALT_VECTOR
  result &= init_alt_Vector();
#endif

#ifdef USE_ALT_ACCEL
  result &= init_alt_Accel(params);
#endif

#ifdef USE_ALT_SIMD2
  result &= init_alt_SIMD2();
#endif

#ifdef USE_ALT_SIMD
  result &= init_alt_SIMD();
#endif

#ifdef USE_ALT_QXS
  result &= init_alt_QXS();
#endif

  return result;
}


bool bridge_alt_fin()
{
  bool result = true;

#ifdef USE_ALT_VECTOR
  result &= fin_alt_Vector();
#endif

#ifdef USE_ALT_ACCEL
  result &= fin_alt_Accel();
#endif

#ifdef USE_ALT_SIMD2
  result &= fin_alt_SIMD2();
#endif

#ifdef USE_ALT_SIMD
  result &= fin_alt_SIMD();
#endif

#ifdef USE_ALT_QXS
  result &= fin_alt_QXS();
#endif

  return result;
}
