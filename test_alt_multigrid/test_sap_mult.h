/*!

        @file    test_sap_mult.cpp

        @brief   test for SAP operator

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: #$

        @version $LastChangedRevision: 2500 $
*/
//====================================================================
#ifndef TEST_SAP_MULT_SIMD_H
#define TEST_SAP_MULT_SIMD_H

#include <vector>
#include <string>

#include "lib/Parameters/parameters.h"
#include "lib/Field/field_G.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

//class Fopr;
class Source;

//====================================================================
class Test_SapMult
{
 public:
  static const std::string class_name;

 private:
  Parameters params_all;
  Bridge::VerboseLevel m_vl;
  unique_ptr<Field_G> U;

 public:

  //! constructor
  Test_SapMult()
    : m_vl(CommonParameters::Vlevel())
  { init(); }

  //! destructor
  ~Test_SapMult() {}

  //! a simple way to test to call something in myGrid
  int test();

 private:

  //! initial setup
  void init();
};
#endif // TEST_SAP_MULT_H

//============================================================END=====
