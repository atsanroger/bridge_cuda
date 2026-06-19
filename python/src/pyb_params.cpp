/*!
        @file    pyb_params.cpp

        @brief   translating python dictionary to Parameters class

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#include <vector>
#include "pyb_defs.h"
#include "pyb_params.h"
#include "bridge.h"

Parameters make_parameters(py::dict kwargs)
{
  Parameters params;

  for (auto item : kwargs) {
    string key = item.first.cast<string>();
    // py::print("key: {}, value={}"_s.format(item.first, item.second));
    // std::cout << item.second.ptr()->ob_type->tp_name << std::endl;

    const char *tp = item.second.ptr()->ob_type->tp_name;

    if (tp == nullptr) {
      warn("value type unknown ");

    } else if (strcmp(tp, "int") == 0) {
      params.set_int(key, item.second.cast<int>());

    } else if (strcmp(tp, "float") == 0) {
      params.set_double(key, item.second.cast<double>());

    } else if (strcmp(tp, "str") == 0) {
      params.set_string(key, item.second.cast<string>());

    } else if (strcmp(tp, "list") == 0) {

      auto vv = item.second.cast<py::list>();
      
      if (vv.size() == 0) {
        params.set_string_vector(key, {});
        warn("empty list");
      } else {

        const char *item_type = nullptr;
        for (auto item : vv) {
          item_type = item.ptr()->ob_type->tp_name;
        }
        // warn("item_type = %s", item_type);

        if (item_type == nullptr) {
          warn("unknown item type.");

        } else if (strcmp(item_type, "int") == 0) {
          vector<int> v;
          for (auto item : vv) {
            v.push_back(item.cast<int>());
          }
          params.set_int_vector(key, v);

        } else if (strcmp(item_type, "float") == 0) {
          vector<double> v;
          for (auto item : vv) {
            v.push_back(item.cast<double>());
          }
          params.set_double_vector(key, v);

        } else if (strcmp(item_type, "str") == 0) {
          vector<string> v;
          for (auto item : vv) {
            v.push_back(item.cast<string>());
          }
          params.set_string_vector(key, v);

        } else {
          warn("unknown list type %s", item_type);
        }
      }

    } else if (strcmp(tp, "dict") == 0) {
      warn("dict: unsupported.");

    } else {
      warn("unsupported type %s", tp);
    }
  }

  return params;
}

void pybind_params(py::module& m)
{
  py::class_<Parameters>(m, "Parameters")
    .def(py::init([](py::kwargs kw) {
          return new Parameters(make_parameters(kw));
        }), REDIRECT)

    .def("get_int", &Parameters::get_int, REDIRECT)
    .def("get_double", &Parameters::get_double, REDIRECT)
    .def("get_string", &Parameters::get_string, REDIRECT)
    .def("get_bool", &Parameters::get_bool, REDIRECT)
    .def("get_int_vector", &Parameters::get_int_vector, REDIRECT)
    .def("get_double_vector", &Parameters::get_double_vector, REDIRECT)
    .def("get_string_vector", &Parameters::get_string_vector, REDIRECT)
    ;

#undef PROP
}

