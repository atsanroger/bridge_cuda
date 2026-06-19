/*!
        @file    fopr_Wilson_eo.cpp

        @brief

        @author  UEDA, Satoru  (sueda)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "fopr_Wilson_eo.h"

#ifdef USE_FACTORY

namespace Selector_Fopr_Wilson_eo
{
  namespace {
    Fopr *create_object()
    {
      return new Fopr_Wilson_eo();
    }


    Fopr *create_object_with_repr(const std::string& repr)
    {
      return new Fopr_Wilson_eo(repr);
    }


    Fopr *create_object_with_params(const Parameters& params)
    {
      return new Fopr_Wilson_eo(params);
    }
  }

  bool register_factory()
  {
    bool init = true;
    init &= Fopr::Factory_noarg::Register("Wilson_eo", create_object);
    init &= Fopr::Factory_string::Register("Wilson_eo", create_object_with_repr);
    init &= Fopr::Factory_params::Register("Wilson_eo", create_object_with_params);
    return init;
  }


#ifdef USE_FACTORY_AUTOREGISTER
  namespace {
    bool init = register_factory();
  }
#endif
}
#endif
