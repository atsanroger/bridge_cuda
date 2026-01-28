/*!

        @file    test_coarse_operator.h

        @brief

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: #$

        @version $LastChangedRevision: 2500 $
*/


#ifndef TEST_COARSE_OPERATOR_H
#define TEST_COARSE_OPERATOR_H

#include <vector>
#include <string>

#include "lib/Parameters/parameters.h"
#include "lib/Field/field_G.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

class Source;

class Test_Coarse
{
 public:
  static const std::string class_name;

 private:
  Parameters params_all;
  Bridge::VerboseLevel m_vl;
  unique_ptr<Field_G> U;

 public:

  //! constructor
  Test_Coarse()
    : m_vl(CommonParameters::Vlevel())
  { init(); }

  //! destructor
  ~Test_Coarse() {}

  int test();

 private:

  //! initial setup
  void init();
};
#endif //
