/*!
        @file    projection.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:01:33 #$

        @version $LastChangedRevision: 1928 $
*/

#include "lib/Smear/projection.h"

#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else
// setup factories for all subclasses

#include "lib/Smear/projection_Maximum_SU_N.h"
#include "lib/Smear/projection_Stout_SU3.h"

template<>
bool Projection::init_factory()
{
  bool result = true;

  result &= Projection_Maximum_SU_N::register_factory();
  result &= Projection_Stout_SU3::register_factory();

  return result;
}


#endif /* USE_FACTORY_AUTOREGISTER */

#endif /* USE_FACTORY */
