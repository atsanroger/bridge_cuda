/*!
        @file    evalexpr_symbol.h

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:01:33 #$

        @version $LastChangedRevision: 1928 $
*/

#ifndef EVALEXPR_SYMBOL_INCLUDED
#define EVALEXPR_SYMBOL_INCLUDED

#include <map>
#include <string>
#include "IO/bridgeIO.h"

typedef double (*function_t)(double);

class SymbolTable
{
 private:

  enum ValueType { VARIABLE, FUNCTION, };

  struct SymbolRecord
  {
    ValueType type;
    union
    {
      double     val;
      function_t fptr;
    }
              value;
  };

  typedef std::map<std::string, SymbolRecord> SymbolMap_t;

  SymbolMap_t table;

 public:

  bool find_symbol(const std::string& name);

  bool put_symbol(const std::string& name, const double value);
  bool put_symbol(const std::string& name, const function_t tptr);

  double get_symbol_value(const std::string& name) const;
  function_t get_symbol_function(const std::string& name) const;

  void print() const;
};
#endif
