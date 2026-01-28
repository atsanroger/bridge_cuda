/*!
        @file    gaugeFixing.h

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef GAUGEFIXING_INCLUDED
#define GAUGEFIXING_INCLUDED

#include "Field/field_G.h"
#include "Parameters/parameters.h"
#include "IO/bridgeIO.h"

#ifdef USE_FACTORY
#include "Tools/factory.h"
#endif

//! gauge fixing.

/*
    This class fixes a gauge of the configuration
                                        [10 Oct 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                        [21 Mar 2015 Y.Namekawa]
    Staple and RandomNumbers are moved into gaugeFixing
                                        [30 Mar 2016 Y.Namekawa]
 */

class GaugeFixing
{
 public:
  GaugeFixing() {}

  virtual ~GaugeFixing() {}

 private:
  // non-copyable
  GaugeFixing(const GaugeFixing&);
  GaugeFixing& operator=(const GaugeFixing&);

 public:
  virtual void set_parameters(const Parameters& params) = 0;
  virtual void get_parameters(Parameters& params) const = 0;

  virtual void fix(Field_G& Ufix, const Field_G& Uorg) = 0;


#ifdef USE_FACTORY
 public:
  typedef GaugeFixing *(*ProductCreator)();
  typedef GaugeFixing *(*ProductCreator_params)(const Parameters&);

  typedef FactoryTemplate<GaugeFixing, ProductCreator>          Factory;
  typedef FactoryTemplate<GaugeFixing, ProductCreator_params>   Factory_params;

  static GaugeFixing *New(const IdentifierType& subtype)
  {
    ProductCreator p = Factory::Find(subtype);
    return p ? (*p)() : 0;
  }

  static GaugeFixing *New(const IdentifierType& subtype, const Parameters& params)
  {
    ProductCreator_params p = Factory_params::Find(subtype);
    return p ? (*p)(params) : 0;
  }

#ifdef USE_FACTORY_AUTOREGISTER
#else
  static bool init_factory();
#endif
#endif
};
#endif
