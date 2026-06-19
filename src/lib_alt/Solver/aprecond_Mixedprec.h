#ifndef APRECOND_MIXEDPREC_H
#define APRECOND_MIXEDPREC_H

#include <cstdio>
#include <cstdlib>

#include <string>
using std::string;
#include <vector>
using std::vector;

#include "lib/Parameters/commonParameters.h"
#include "lib/Fopr/afopr.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib_alt/Solver/aprecond.h"

#include "lib_alt/Field/afield_base.h"

// temporary
#include "lib_alt/Solver/asolver_BiCGStab.h"


class Field;

template<typename AFIELD, typename AFIELD2>
class APrecond_Mixedprec : public APrecond<AFIELD>
{
 public:
  typedef typename AFIELD::real_t    real_t;
  typedef typename AFIELD2::real_t   real_t2;
  static const std::string class_name;

  enum field_type { LEXICAL, EVEN_ODD };

 private:
  ASolver<AFIELD2> *m_solver;
  AFIELD2 m_v2, m_w2;

  field_type m_field_type;

  double m_accum_flop;

 public:

  APrecond_Mixedprec(ASolver<AFIELD2> *solver)
  {
    m_solver = solver;
    init();
  }

  //!< destructor
  ~APrecond_Mixedprec()
  {
    tidyup();
  }

  void mult(AFIELD&, const AFIELD&);

  void reset_flop_count();

  double flop_count();

 private:
  void init();
  void tidyup();
};

#endif  // APRECOND_MIXEDPREC_H
