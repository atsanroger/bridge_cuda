#ifndef ASOLVER_H
#define ASOLVER_H

#include <cstdio>
#include <cstdlib>
#include <vector>
using std::vector;
#include <string>
using std::string;

#include  "lib/Parameters/commonParameters.h"
#include  "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib/Fopr/afopr.h"

#ifdef USE_FACTORY
#include "lib/Tools/factory.h"
#include "lib/Tools/director.h"
#endif

template<typename AFIELD>
class ASolver
{
 protected:
  Bridge::VerboseLevel m_vl;

 public:
  typedef typename AFIELD::real_t real_t;

  enum InitialGuess { RHS, GIVEN, ZERO };

  ASolver()
    : m_vl(CommonParameters::Vlevel()) {}

  virtual ~ASolver() {}

  virtual void set_parameters(const Parameters& params) = 0;

  void set_parameter_verboselevel(const Bridge::VerboseLevel vl)
  { m_vl = vl; }

  virtual void set_init_mode(const InitialGuess init_guess)
  {
    vout.crucial("set_init_mode is called\n");
    abort();
  }

  virtual AFopr<AFIELD> *get_fopr() { return 0; }

  virtual void solve(AFIELD& x, const AFIELD& b, int& nconv, real_t& diff) { }

  // added for multi-grid algorithm
  virtual void solve(AFIELD& x, const AFIELD& b, int& nconv, real_t& diff, int ieo) { }

  virtual double flop_count() { return 0.0; }


#ifdef USE_FACTORY
 public:
  typedef ASolver *(*ProductCreator)(AFopr<AFIELD> *);
  typedef FactoryTemplate<ASolver, ProductCreator> Factory_fopr;

  static ASolver *New(const IdentifierType& subtype, AFopr<AFIELD> *afopr)
  {
    ProductCreator p = Factory_fopr::Find(subtype);
    return p ? (*p)(afopr) : 0;
  }

  static ASolver *New(const IdentifierType& subtype,
                      unique_ptr<AFopr<AFIELD> >& afopr)
  {
    ProductCreator p = Factory_fopr::Find(subtype);
    return p ? (*p)(afopr.get()) : 0;
  }

#ifdef USE_FACTORY_AUTOREGISTER
#else
  static bool init_factory();
#endif
#endif  // USE_FACTORY
};

#endif // ASOLVER_H
