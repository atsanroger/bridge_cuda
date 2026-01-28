/*!
        @file    force_G.cpp

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:01:33 #$

        @version $LastChangedRevision: 1928 $
*/

#include "force_G.h"

#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else
// setup factories for all subclasses

#include "force_G_Plaq.h"
#include "force_G_Rectangle.h"
#include "force_G_Plaq_SF.h"
#include "force_G_Rectangle_SF.h"

bool Force_G::init_factory()
{
  bool result = true;

  result &= Force_G_Plaq::register_factory();
  result &= Force_G_Rectangle::register_factory();
  result &= Force_G_Plaq_SF::register_factory();
  result &= Force_G_Rectangle_SF::register_factory();

  return result;
}


#endif /* USE_FACTORY_AUTOREGISTER */

#endif /* USE_FACTORY */
