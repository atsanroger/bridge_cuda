/*!
        @file    source.h

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef SOURCE_INCLUDED
#define SOURCE_INCLUDED

#include "Parameters/parameters.h"
#include "Field/field.h"

#include "IO/bridgeIO.h"

#ifdef USE_FACTORY
#include "Tools/factory.h"
#endif

//! Base class of source for a linear solver.

/*!
    (Coding history will be recovered from trac.)
    Parameters_Source_All is implemented for source_selector.
                                      [ 2 Feb 2013 Y.Namekawa]
    Bug fix: initialization of verbose level was added.
                                      [23 May 2016 H.Matsufuru]
    Add set_all_color_spin,etc        [ 4 Apr 2017 Y.Namekawa]
 */

class Source
{
 public:
  Source() {}
  virtual ~Source() {}

 private:
  // non-copyable
  Source(const Source&);
  Source& operator=(const Source&);

 public:
  virtual void set_parameters(const Parameters&) = 0;

  virtual void get_parameters(Parameters&) const = 0;

  virtual void set(Field&, const int)            = 0;
  virtual void set(Field&, const int, const int) = 0;
  virtual void set_all_color(Field&, const int)  = 0;
  virtual void set_all_color_spin(Field&)        = 0;

#ifdef USE_FACTORY
 public:
  typedef Source *(*ProductCreator)();
  typedef Source *(*ProductCreator_params)(const Parameters&);

  typedef FactoryTemplate<Source, ProductCreator>          Factory;
  typedef FactoryTemplate<Source, ProductCreator_params>   Factory_params;

  static Source *New(const IdentifierType& subtype)
  {
    ProductCreator p = Factory::Find(subtype);
    return p ? (*p)() : 0;
  }

  static Source *New(const IdentifierType& subtype, const Parameters& params)
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
#endif /* SOURCE_INCLUDED */
