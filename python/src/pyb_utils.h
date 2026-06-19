/*!
        @file    pyb_utils.h

        @brief   pybridge utilities

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#ifndef PYB_UTILS_INCLUDED
#define PYB_UTILS_INCLUDED

#include <iostream>

template <typename T>
std::ostream& operator<<(std::ostream& os, const vector<T>& v)
{
  if (v.size() == 0) return os << "[]";
  os << "[ ";
  int i = 0;
  for (auto item : v) {
    if (i++ > 0) os << ", ";
    os << item;
  }
  os << " ]";
  return os;
}

#endif /* PYB_UTILS_INCLUDED */
