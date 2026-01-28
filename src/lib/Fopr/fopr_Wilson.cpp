/*!
        @file    fopr_Wilson.cpp

        @brief

        @author  Hideo Matsufuru  (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate: 2013-03-21 15:28:34 #$

        @version $LastChangedRevision: 2499 $
*/

#include "fopr_Wilson.h"

#ifdef USE_FACTORY
namespace Selector_Fopr_Wilson
{
  namespace {
    Fopr *create_object()
    {
      return new Fopr_Wilson();
    }


    Fopr *create_object_with_repr(const std::string& repr)
    {
      return new Fopr_Wilson(repr);
    }


    Fopr *create_object_with_params(const Parameters& params)
    {
      return new Fopr_Wilson(params);
    }
  }

  bool register_factory()
  {
    bool init = true;
    init &= Fopr::Factory_noarg::Register("Wilson", create_object);
    init &= Fopr::Factory_string::Register("Wilson", create_object_with_repr);
    init &= Fopr::Factory_params::Register("Wilson", create_object_with_params);
    return init;
  }


#ifdef USE_FACTORY_AUTOREGISTER
  namespace {
    bool init = register_factory();
  }
#endif
}
#endif
