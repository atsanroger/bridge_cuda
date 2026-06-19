/*!
        @file    fopr_WilsonGeneral.cpp

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate: 2013-03-21 15:28:34 #$

        @version $LastChangedRevision: 2499 $
*/

#include "fopr_WilsonGeneral.h"

#ifdef USE_FACTORY

namespace Selector_Fopr_WilsonGeneral
{
  namespace {
    Fopr *create_object()
    {
      return new Fopr_WilsonGeneral();
    }


    Fopr *create_object_with_repr(const std::string& repr)
    {
      return new Fopr_WilsonGeneral(repr);
    }


    Fopr *create_object_with_params(const Parameters& params)
    {
      return new Fopr_WilsonGeneral(params);
    }
  }

  bool register_factory()
  {
    bool init = true;
    init &= Fopr::Factory_noarg::Register("WilsonGeneral", create_object);
    init &= Fopr::Factory_string::Register("WilsonGeneral", create_object_with_repr);
    init &= Fopr::Factory_params::Register("WilsonGeneral", create_object_with_params);
    return init;
  }


#ifdef USE_FACTORY_AUTOREGISTER
  namespace {
    bool init = register_factory();
  }
#endif
}
#endif
