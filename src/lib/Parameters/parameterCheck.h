/*!
        @file    parameterCheck.h

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:01:33 #$

        @version $LastChangedRevision: 1928 $
*/

//! checker of input parameters.

/*!
    This checker examines if a parameter satisfies the specified
    condition such as non-zero. The checker returns EXIT_SUCCESS
    or EXIT_FAILURE, instead of bool in the original Aoyama-san's code.
                                      [16 Jun 2013 Y.Namekawa]
*/

#ifndef PARAMETERCHECK_INCLUDED
#define PARAMETERCHECK_INCLUDED

#include "commonParameters.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

namespace ParameterCheck
{
  typedef bool (*valid_double)(const double);
  typedef bool (*valid_int)(const int);
  typedef bool (*valid_double_vector)(const std::vector<double>&);
  typedef bool (*valid_int_vector)(const std::vector<int>&);
  typedef bool (*valid_string)(const std::string&);

  int non_negative(const int v);
  int non_zero(const int v);
  int non_zero(const double v);
  int square_non_zero(const double v);
  int non_NULL(const std::string v);

  int is_satisfied(const bool cond);
}
#endif
