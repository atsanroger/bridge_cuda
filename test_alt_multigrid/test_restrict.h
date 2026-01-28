/*!

        @file    test_coarse_operator.h

        @brief

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: #$

        @version $LastChangedRevision: 2500 $
*/


#ifndef TEST_RESTRICT_H
#define TEST_RESTRICT_H

#include <vector>
#include <string>

#include "lib/Parameters/parameters.h"
#include "lib/Field/field_G.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

class Source;

class Test_Restrict
{
 public:
  static const std::string class_name;

 private:
  Parameters params_all;
  Bridge::VerboseLevel m_vl;
  unique_ptr<Field_G> U;

 public:

  //! constructor
  Test_Restrict()
    : m_vl(CommonParameters::Vlevel())
  { init(); }

  //! destructor
  ~Test_Restrict() {}

  int test(std::string);

 private:

  //! initial setup
  void init();
};
#endif //
