/*!
        @file    pyb_index.cpp

        @brief   python binding of site index classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#include "pyb_defs.h"
#include "pyb_params.h"
#include "pyb_index.h"
#include "bridge.h"
#include "Field/index_lex.h"
#include "Field/index_eo.h"
#include "Parameters/parameters.h"

void pybind_index(py::module m)
{
  //! Index_lex
  //
  py::class_<Index_lex>(m, "Index_lex", R"pybdoc(
Lexicographical site index
)pybdoc")
    .def(py::init<>(), REDIRECT)
    .def(py::init<int,int,int,int>(), REDIRECT)

    .def("site", &Index_lex::site, REDIRECT)
    ;

  //! Index_eo
  //
  py::class_<Index_eo>(m, "Index_eo", R"pybdoc(
Even-odd site index
)pybdoc")
    .def(py::init<>(), REDIRECT)

    .def("leo", &Index_eo::leo,
         py::arg("y"), py::arg("z"), py::arg("t"))

    .def("site", [](Index_eo& self, int x2, int y, int z, int t, int ieo)
         {
           return self.site(x2, y, z, t, ieo);
         },
         py::arg("x2"), py::arg("y"), py::arg("z"), py::arg("t"), py::arg("ieo"))
    .def("site", [](Index_eo& self, int is, int ieo)
         {
           return self.site(is, ieo);
         },
         py::arg("is"), py::arg("ieo"))

    .def("site_up", &Index_eo::site_up,
         py::arg("x2"), py::arg("y"), py::arg("z"), py::arg("t"), py::arg("ieo"))
    .def("site_dn", &Index_eo::site_dn,
         py::arg("x2"), py::arg("y"), py::arg("z"), py::arg("t"), py::arg("ieo"))

    .def("siteh", &Index_eo::siteh,
         py::arg("x2"), py::arg("y"), py::arg("z"), py::arg("t"))
    .def("siteh_up", &Index_eo::siteh_up,
         py::arg("x2"), py::arg("y"), py::arg("z"), py::arg("t"), py::arg("ieo"))
    .def("siteh_dn", &Index_eo::siteh_dn,
         py::arg("x2"), py::arg("y"), py::arg("z"), py::arg("t"), py::arg("ieo"))

    .def("convert_field", [](Index_eo& self, Field& f_eo, const Field& f_lex)
         {
           return self.convertField(f_eo, f_lex);
         }, py::arg("eo"), py::arg("lex"),
         REDIRECT)

    .def("convert_field", [](Index_eo& self, Field& f_eo, const Field& f_lex, const int ieo)
         {
           return self.convertField(f_eo, f_lex, ieo);
         }, py::arg("eo"), py::arg("lex"), py::arg("ieo"),
         REDIRECT)

    .def("reverse_field", [](Index_eo& self, Field& f_lex, const Field& f_eo)
         {
           return self.reverseField(f_lex, f_eo);
         }, py::arg("lex"), py::arg("eo"),
         REDIRECT)

    .def("reverse_field", [](Index_eo& self, Field& f_lex, const Field& f_eo, const int ieo)
         {
           return self.reverseField(f_lex, f_eo, ieo);
         }, py::arg("lex"), py::arg("eo"), py::arg("ieo"),
         REDIRECT)


    .def("split_field", [](Index_eo& self, Field& f_even, Field& f_odd, const Field& f_eo)
         {
           return self.splitField(f_even, f_odd, f_eo);
         }, py::arg("even"), py::arg("odd"), py::arg("eo"),
         REDIRECT)

    .def("merge_field", [](Index_eo& self, Field& f_eo, const Field& f_even,  const Field& f_odd)
         {
           return self.mergeField(f_eo, f_even, f_odd);
         }, py::arg("eo"), py::arg("even"), py::arg("odd"),
         REDIRECT)
    ;

}
