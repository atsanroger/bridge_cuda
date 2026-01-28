/*!
        @file    pyb_rng.cpp

        @brief   python binding of random number generator classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#include "pyb_defs.h"
#include "pyb_params.h"
#include "pyb_rng.h"
#include "bridge.h"
#include "Tools/randomNumbers.h"
#include "Parameters/parameters.h"

void pybind_rng(py::module m)
{
  //! RandomNumbers
  //
  py::class_<RandomNumbers>(m, "RandomNumbers", R"pybdoc(
Random number generator classes

See Also
--------
RandomNumbers_Mseries
RandomNumbers_MT19937
)pybdoc"
#ifdef USE_SFMTLIB
R"pybdoc(
RandomNumbers_SFMT
)pybdoc"
#endif
)
    .def(py::init([](const std::string& subtype, py::kwargs kwargs){

          RandomNumbers *obj = nullptr;

          if (kwargs.contains("state_file")) {
            std::string state_file = kwargs["state_file"].cast<std::string>();
            obj = RandomNumbers::New(subtype, state_file);
          } else if (kwargs.contains("seed")) {
            int seed = kwargs["seed"].cast<int>();
            obj = RandomNumbers::New(subtype, seed);
          }

          // if (obj) obj->set_parameters(make_parameters(kwargs));
          return obj;
        }), py::arg("type"),
      REDIRECT)

    .def("get", &RandomNumbers::get)
    .def("gauss", [](RandomNumbers& self) -> std::tuple<double,double>
         {
           double r1, r2;
           self.gauss(r1, r2);
           return std::tuple<double,double>(r1, r2);
         })

    .def("lex_global", [](RandomNumbers& self,
                          const std::string& type, Field& v) -> void
         {
           return self.lex_global(type, v);
         }, py::arg("type"), py::arg("target"),
         REDIRECT
    )

    .def("gauss_lex_global", &RandomNumbers::gauss_lex_global, REDIRECT)
    .def("uniform_lex_global", &RandomNumbers::uniform_lex_global, REDIRECT)
    .def("U1_lex_global", &RandomNumbers::U1_lex_global, REDIRECT)
    .def("Z2_lex_global", &RandomNumbers::Z2_lex_global, REDIRECT)

    .def("read_file", &RandomNumbers::read_file, REDIRECT)
    .def("write_file", &RandomNumbers::read_file, REDIRECT)

    .def("reset", &RandomNumbers::reset, REDIRECT)
    ;

          
  //! RNG algorithms
  //

  py::class_<RandomNumbers_MT19937, RandomNumbers>(m, "RandomNumbers_MT19937")
    .def(py::init<int>(), py::arg("seed"), REDIRECT)
    // .def(py::init<std::vector<unsigned long> >(), py::arg("seed"), REDIRECT)
    .def(py::init([](std::vector<unsigned long>& seed) {
          return new RandomNumbers_MT19937(seed);
        }), py::arg("seed"),
      REDIRECT)
    .def(py::init<std::string>(), py::arg("state_file"), REDIRECT)
    ;

  py::class_<RandomNumbers_Mseries, RandomNumbers>(m, "RandomNumbers_Mseries")
    .def(py::init<int>(), py::arg("seed"), REDIRECT)
    .def(py::init<std::string>(), py::arg("state_file"), REDIRECT)
    ;

#ifdef USE_SFMTLIB
  py::class_<RandomNumbers_SFMT, RandomNumbers>(m, "RandomNumbers_SFMT")
    .def(py::init<int>(), py::arg("seed"), REDIRECT)
    .def(py::init<std::string>(), py::arg("state_file"), REDIRECT)
    ;
#endif

  //! RandomNumberManager
  //
  py::class_<RandomNumberManager>(m, "RandomNumberManager", R"pybdoc(
Manages a global instance of RandomNumbers class
)pybdoc")
    .def(py::init([](){
          warn("invalid initialization of RandomNumberManager class");
          return nullptr;
        }), REDIRECT)

    // .def_static("initialize", &RandomNumberManager::initialize, REDIRECT)
    .def_static("reset", &RandomNumberManager::reset, REDIRECT)

    .def_static("save_state", &RandomNumberManager::save_state, REDIRECT)
    .def_static("restore_state", &RandomNumberManager::restore_state, REDIRECT)
    ;
  
}
