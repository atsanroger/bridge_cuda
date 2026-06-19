/*!
        @file    staple.h

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-09-08 06:55:16 #$

        @version $LastChangedRevision: 2314 $
*/

#ifndef STAPLE_INCLUDED
#define STAPLE_INCLUDED

#include "Parameters/parameters.h"
#include "Field/field_G.h"

#include "IO/bridgeIO.h"

#ifdef USE_FACTORY
#include "Tools/factory.h"
#endif

//! Base class for Staple construction.

/*!
    This class defines interface of Staple-type classes.
                                    [24 Jan 2017 Y.Namekawa]
 */

class Staple
{
 public:
  Staple() {}
  virtual ~Staple() {}

 private:
  // non-copyable
  Staple(const Staple&);
  Staple& operator=(const Staple&);

 public:
  //! setting parameters.
  virtual void set_parameters(const Parameters& params) = 0;

  //! getting parameters.
  virtual void get_parameters(Parameters& params) const = 0;

  //! constructs upper staple in mu-nu plane.
  virtual void upper(Field_G&, const Field_G&, const int mu, const int nu) = 0;

  //! constructs lower staple in mu-nu plane.
  virtual void lower(Field_G&, const Field_G&, const int mu, const int nu) = 0;

  //! constructs staple in mu-direction (summing up nu-direction).
  virtual void staple(Field_G&, const Field_G&, const int mu) = 0;

  //! calculates plaquette value.
  virtual double plaquette(const Field_G&) = 0;

  //! calculates spatial plaquette value.
  virtual double plaq_s(const Field_G&) = 0;

  //! calculates temporal plaquette value.
  virtual double plaq_t(const Field_G&) = 0;

#ifdef USE_FACTORY
 public:
  typedef Staple *(*ProductCreator)();
  typedef FactoryTemplate<Staple, ProductCreator> Factory;

  static Staple *New(const IdentifierType& subtype)
  {
    ProductCreator p = Factory::Find(subtype);

    return p ? (*p)() : 0;
  }

#ifdef USE_FACTORY_AUTOREGISTER
#else
  static bool init_factory();
#endif
#endif
};
#endif
