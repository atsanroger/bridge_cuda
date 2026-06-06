/*!
        @file    afopr.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef AFOPR_INCLUDED
#define AFOPR_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <string>

#include  "Parameters/commonParameters.h"
#include  "Parameters/parameters.h"
#include  "Field/field.h"
#include  "Field/field_G.h"
#include  "IO/bridgeIO.h"
using Bridge::vout;

#ifdef USE_FACTORY
#include "Tools/factory.h"
#include "Tools/director.h"
#endif

class Field;

/*!
    The base template class of fermion operators.
    AFopr<FIELD> is a fermion operator that acts on a FIELD
    vector.  This template class was first introduced to
    incorporate the alternative fermion classes, but later
    realized as the basis of all the fermion operator by
    incorporating to the Bridge++ core library.
    The implementation basically follows the original Fopr
    class in the core libirary.
                                   [17 Sep 2018 H.Matsufuru]
*/

template<typename AFIELD>
class AFopr
{
 protected:
  static const std::string class_name;

 public:
  AFopr() {}
  virtual ~AFopr() {}

 private:
  //! non-copyable
  AFopr(const AFopr&);
  AFopr& operator=(const AFopr&);

 public:

  //! sets parameters by a Parameter object: to be implemented in a subclass.
  virtual void set_parameters(const Parameters& params)
  {
    vout.crucial("AFopr: set_parameters not implemented.\n");
  }

  //! gets parameters by a Parameter object: to be implemented in a subclass.
  virtual void get_parameters(Parameters& params) const
  {
    vout.crucial("AFopr: get_parameters not implemented.\n");
  }

  //! sets the gauge configuration.
  virtual void set_config(Field *) = 0;

  //! \brief setting the mode of multiplication if necessary.
  //!  Default implementation here is just to avoid irrelevant call.
  virtual void set_mode(std::string mode)
  {
    vout.crucial("Fopr: set_mode not implemented.\n");
    exit(EXIT_FAILURE);
  }

  //! returns the current mult mode.
  virtual std::string get_mode() const
  {
    vout.general("Fopr: get_mode not implemented.\n");
    return std::string();
  }

  //! multiword arithmetic mode: FP = standard, DW = double-word, TW = triple-word.
  enum class MWMode : int { FP = 0, DW = 1, TW = 2 };

  virtual void set_mw_mode(MWMode mode) {}
  virtual MWMode get_mw_mode() const { return MWMode::FP; }

  //! multiplies fermion operator to a given field.
  virtual void mult(AFIELD&, const AFIELD&)
  {
    vout.crucial("AFopr: mult not implemented.\n");
    exit(EXIT_FAILURE);
  }

  //! hermitian conjugate of mult.
  virtual void mult_dag(AFIELD&, const AFIELD&)
  {
    vout.crucial("AFopr: mult_dag not implemented.\n");
    exit(EXIT_FAILURE);
  }

  //! executes mult with specified mode (unchanging internal mode).
  virtual void mult(AFIELD&, const AFIELD&, const std::string mode)
  {
    vout.crucial("AFopr: mult with mode not implemented.\n");
    exit(EXIT_FAILURE);
  }

  //! executes mult_dag with specified mode (unchanging internal mode).
  virtual void mult_dag(AFIELD&, const AFIELD&, const std::string mode)
  {
    vout.crucial("AFopr: mult with mode not implemented.\n");
    exit(EXIT_FAILURE);
  }

  //! multiplies gamma_5 matrix.
  virtual void mult_gm5(AFIELD&, const AFIELD&)
  {
    vout.crucial("AFopr: mult_gm5 not implemented.\n");
    exit(EXIT_FAILURE);
  }

  //! upward nearest neighbor hopping term.
  virtual void mult_up(int mu, AFIELD&, const AFIELD&)
  {
    vout.crucial("AFopr: mult_up not implemented.\n");
    exit(EXIT_FAILURE);
  }

  //! downward nearest neighbor hopping term.
  virtual void mult_dn(int mu, AFIELD&, const AFIELD&)
  {
    vout.crucial("AFopr: mult_dn not implemented.\n");
    exit(EXIT_FAILURE);
  }

  //! normalize propagator if necessary (default: do nothing)
  virtual void normalize_fprop(AFIELD&) {         }

  //! normalize propagator if necessary (default: do nothing)
  virtual void normalize_fopr(AFIELD&) {         }


  //! returns the on-site degree of freedom of the fermion field.
  virtual int field_nin() = 0;

  //! returns the volume of the fermion field.
  virtual int field_nvol() = 0;

  //! returns the external degree of freedom of the fermion field.
  virtual int field_nex() = 0;

  //! returns the number of floating point operations.
  virtual double flop_count()
  {
    vout.crucial("AFopr: mult_dn not implemented.\n");
    return 0.0;
  }

  //! returns the flops per site for specified mode.
  virtual double flop_count(const std::string mode)
  {
    vout.crucial("AFopr: mult_dn not implemented.\n");
    return 0.0;
  }

  //! returns true if additional field conversion is needed.
  virtual bool needs_convert() { return false; }

  //! converts a Field object into other format if necessary.
  virtual void convert(AFIELD&, const Field&)
  {
    vout.crucial("AFopr: convert is not implemented.\n");
    exit(EXIT_FAILURE);
  }

  //! converts an alternative field to a Field object.
  virtual void reverse(Field&, const AFIELD&)
  {
    vout.crucial("AFopr: reverse is not implemented.\n");
    exit(EXIT_FAILURE);
  }

#ifdef USE_FACTORY
 public:
  typedef AFopr *(*ProductCreator_noarg)();
  typedef AFopr *(*ProductCreator_fopr)(AFopr *fopr);

  typedef AFopr *(*ProductCreator_fopr_director_params)(AFopr *fopr, Director *director, const Parameters& params);
  typedef AFopr *(*ProductCreator_fopr_director)(AFopr *fopr, Director *director);
  typedef AFopr *(*ProductCreator_string)(const std::string& arg);
  typedef AFopr *(*ProductCreator_params)(const Parameters& params);
  typedef AFopr *(*ProductCreator_fopr_params)(AFopr *fopr, const Parameters& params);

  typedef FactoryTemplate<AFopr, ProductCreator_noarg>    Factory_noarg;
  typedef FactoryTemplate<AFopr, ProductCreator_fopr>     Factory_fopr;
  typedef FactoryTemplate<AFopr, ProductCreator_fopr_director>
    Factory_fopr_director;
  typedef FactoryTemplate<AFopr, ProductCreator_fopr_director_params>
    Factory_fopr_director_params;
  typedef FactoryTemplate<AFopr, ProductCreator_string>   Factory_string;
  typedef FactoryTemplate<AFopr, ProductCreator_params>   Factory_params;
  typedef FactoryTemplate<AFopr, ProductCreator_fopr_params>
    Factory_fopr_params;

  static AFopr *New(const IdentifierType& subtype)
  {
    ProductCreator_noarg p = Factory_noarg::Find(subtype);
    return p ? (*p)() : 0;
  }

  static AFopr *New(const IdentifierType& subtype, AFopr *fopr)
  {
    ProductCreator_fopr p = Factory_fopr::Find(subtype);
    return p ? (*p)(fopr) : 0;
  }

  static AFopr *New(const IdentifierType& subtype, AFopr *fopr,
                    Director *director)
  {
    ProductCreator_fopr_director p = Factory_fopr_director::Find(subtype);
    return p ? (*p)(fopr, director) : 0;
  }

  static AFopr *New(const IdentifierType& subtype, AFopr *fopr,
                    Director *director, const Parameters& params)
  {
    ProductCreator_fopr_director_params p =
                    Factory_fopr_director_params::Find(subtype);
    return p ? (*p)(fopr, director, params) : 0;
  }

  static AFopr *New(const IdentifierType& subtype, const std::string& arg)
  {
    ProductCreator_string p = Factory_string::Find(subtype);
    return p ? (*p)(arg) : 0;
  }

  static AFopr *New(const IdentifierType& subtype, const Parameters& params)
  {
    ProductCreator_params p = Factory_params::Find(subtype);
    return p ? (*p)(params) : 0;
  }

  static AFopr *New(const IdentifierType& subtype,
                    AFopr *fopr, const Parameters& params)
  {
    ProductCreator_fopr_params p = Factory_fopr_params::Find(subtype);
    return p ? (*p)(fopr, params) : 0;
  }

#ifdef USE_FACTORY_AUTOREGISTER
#else
  static bool init_factory();
#endif
#endif  // USE_FACTORY
};

#endif  // AFOPR_H
