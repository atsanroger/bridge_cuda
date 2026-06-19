/*!
        @file    source.cpp

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2022-12-16 15:57:38 #$

        @version $LastChangedRevision: 2422 $
*/

#include "source.h"

#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else
// setup factories for all subclasses

#include "source_Local.h"
#include "source_Wall.h"
#include "source_Exponential.h"
#include "source_MomentumWall.h"
#include "source_Random.h"

#include "source_Staggered_Wall_Evenodd.h"
#include "source_Staggered_Wall_Cube.h"

bool Source::init_factory()
{
  bool result = true;

  result &= Source_Local::register_factory();
  result &= Source_Wall::register_factory();
  result &= Source_Exponential::register_factory();
  result &= Source_MomentumWall::register_factory();
  result &= Source_Random::register_factory();

  result &= Source_Staggered_Wall_Evenodd::register_factory();
  result &= Source_Staggered_Wall_Cube::register_factory();

  return result;
}


#endif /* USE_FACTORY_AUTOREGISTER */

#endif /* USE_FACTORY */
