// Time-stamp: <2019-03-11 15:26:55 kanamori>

/*!

        @file    test_afopr_coarse.h

        @brief   test for coarse operator (alt_Simd version)

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: #$

        @version $LastChangedRevision: 2500 $
*/


#ifndef TEST_MG_SOLVER_H
#define TEST_MG_SOLVER_H

#include <vector>
#include <string>

#include "lib/Parameters/parameters.h"
#include "lib/Field/field_G.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

//class Fopr;
//class Source;

class Test_afopr_coarse
{
 public:
  static const std::string class_name;

 private:
  Parameters params_all;
  Bridge::VerboseLevel m_vl;
  unique_ptr<Field_G> U;

 public:

  //! constructor
  Test_afopr_coarse()
    : m_vl(CommonParameters::Vlevel())
  { init(); }

  //! destructor
  ~Test_afopr_coarse() {}

  int test();
  int test_parallel();

 private:

  //! initial setup
  void init();
};
#endif //
