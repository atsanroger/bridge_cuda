/*!
      @file    aprecond.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#ifndef APRECOND_H
#define APRECOND_H

#include <cstdio>
#include <cstdlib>

#include <string>
using std::string;

#include <vector>
using std::vector;

#include  "lib/Parameters/commonParameters.h"
#include  "lib/IO/bridgeIO.h"
using Bridge::vout;

class Field;

template<typename AFIELD>
class APrecond
{
 protected:
  Bridge::VerboseLevel m_vl;

 public:

  APrecond()
    : m_vl(CommonParameters::Vlevel()) {}

  virtual ~APrecond() {}

  virtual void mult(AFIELD&, const AFIELD&)
  {
    vout.crucial(m_vl, "this is mult of base class APrecond!\n");
  }

  virtual void reset_flop_count() {  }

  virtual double flop_count() { return 0; }
};

#endif  // AFOPR_H
