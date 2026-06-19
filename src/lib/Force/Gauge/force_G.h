/*!
        @file    force_G.h

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef FORCE_G_INCLUDED
#define FORCE_G_INCLUDED

#include "Field/field_G.h"
#include "Parameters/parameters.h"

#ifdef USE_FACTORY
#include "Tools/factory.h"
#endif

//! Base class of gauge force calculation.

/*!
    This class defines the interface of gauge force calculation.
                                       [3 Mar 2016 Y.Namekawa]
*/

class Force_G
{
 protected:
  Field_G *m_U;

 public:
  Force_G()
    : m_U(0) {}

  virtual ~Force_G() {}

 private:
  // non-copyable
  Force_G(const Force_G&);
  Force_G& operator=(const Force_G&);

 public:
  virtual void set_parameters(const Parameters&) = 0;

  virtual void get_parameters(Parameters&) const = 0;

  void set_config(Field *U)
  {
    m_U = (Field_G *)U;
  }

  void set_config(Field_G *U)
  {
    m_U = U;
  }

  virtual void force_core(Field&) = 0;

  virtual void force_core(Field& v, Field *U)
  {
    set_config(U);
    force_core(v);
  }

  virtual void force_core(Field& v, Field_G *U)
  {
    set_config(U);
    force_core(v);
  }

#ifdef USE_FACTORY
 public:
  typedef Force_G *(*ProductCreator)();
  typedef Force_G *(*ProductCreator_params)(const Parameters&);

  typedef FactoryTemplate<Force_G, ProductCreator>          Factory;
  typedef FactoryTemplate<Force_G, ProductCreator_params>   Factory_params;

  static Force_G *New(const IdentifierType& subtype)
  {
    ProductCreator p = Factory::Find(subtype);

    return p ? (*p)() : 0;
  }

  static Force_G *New(const IdentifierType& subtype, const Parameters& params)
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
