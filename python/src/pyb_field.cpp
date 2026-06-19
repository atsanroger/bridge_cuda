/*!
        @file    pyb_field.cpp

        @brief   python binding of Field classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#include <pybind11/buffer_info.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "pyb_defs.h"
#include "pyb_field.h"
#include "pyb_tools_su_n.h"
#include "bridge.h"
#include "Field/field.h"
#include "Field/field_G.h"
#include "Field/field_F.h"
#ifdef TRUNK
#include "Field/field_F_1spinor.h"
#endif

void pybind_field(py::module m)
{
  using SU_N::Vec_SU_N;
  using SU_N::Mat_SU_N;
  
  py::enum_<typename Element_type::type>(m, "ElementType")
    .value("REAL", Element_type::REAL)
    .value("COMPLEX", Element_type::COMPLEX)
    .export_values();

  py::class_<Field>(m, "Field", py::buffer_protocol(), R"pybdoc(
Container of Field-type object
)pybdoc")
    .def(py::init<>(), REDIRECT)
    .def(py::init<int, int, int, typename Element_type::type>(),
         py::arg("nin"),
         py::arg("nvol"),
         py::arg("nex"),
         py::arg("element_type") = Element_type::COMPLEX)

    .def("reset", &Field::reset,
         py::arg("nin"),
         py::arg("nvol"),
         py::arg("nex"),
         py::arg("element_type") = Element_type::COMPLEX)

    .def("clone", &Field::clone)
    
    .def_property_readonly("nin", &Field::nin)
    .def_property_readonly("nvol", &Field::nvol)
    .def_property_readonly("nex", &Field::nex)
    .def_property_readonly("element_type", &Field::field_element_type)

    .def_property_readonly("ntot", &Field::ntot)
    .def_property_readonly("size", &Field::size)

    .def("norm2", &Field::norm2)
    .def("norm", &Field::norm)

    .def("is_compat", [](const Field& self, const Field& rhs)->bool {
        return self.check_size(rhs.nin(), rhs.nvol(), rhs.nex());
      })

    .def("cmp", (double (Field::*)(int i) const) &Field::cmp)
    .def("cmp", (double (Field::*)(int iin, int isite, int iex) const) &Field::cmp)
    .def("set", (void (Field::*)(int i, double v)) &Field::set)
    .def("set", (void (Field::*)(int iin, int isite, int iex, double v)) &Field::set)

    .def("fill", (void (Field::*)(double)) &Field::set, REDIRECT)

    .def("setc", (void (Field::*)(double)) &Field::setc, REDIRECT)
    .def("setc", (void (Field::*)(dcomplex)) &Field::setc, REDIRECT)
    
    // buffer_protocol
    .def_buffer([](Field& f) -> py::buffer_info
                {
                  if (f.field_element_type() == Element_type::COMPLEX) {

                    std::vector<ssize_t> shape = {
                      f.nin()/2,
                      f.nvol(),
                      f.nex()
                    };
                    std::vector<ssize_t> strides = {
                      ssize_t(sizeof(dcomplex)),
                      ssize_t(sizeof(dcomplex)) * f.nin()/2,
                      ssize_t(sizeof(dcomplex)) * f.nin()/2 * f.nvol()
                    };

                    return py::buffer_info(
                      (dcomplex *)f.ptr(0),
                      shape,
                      strides,
                      false // readonly
                    );

                  }
                  else if (f.field_element_type() == Element_type::REAL) {

                    std::vector<ssize_t> shape = {
                      f.nin(),
                      f.nvol(),
                      f.nex()
                    };
                    std::vector<ssize_t> strides = {
                      ssize_t(sizeof(double)),
                      ssize_t(sizeof(double)) * f.nin(),
                      ssize_t(sizeof(double)) * f.nin() * f.nvol()
                    };

                    return py::buffer_info(
                      (double *)f.ptr(0),
                      shape,
                      strides,
                      false // readonly
                    );

                  }
                  else {
                    throw std::runtime_error("bad type");
                  }
                })

    // import array
    .def(py::init([](py::buffer b) {
          py::buffer_info info = b.request();

          if (info.ndim != 3) {
            throw std::runtime_error("incompatible buffer dimension. (Nin,Nvol,Nex) expected.");
          }

          if (info.format == py::format_descriptor<dcomplex>::format()) {
            // complex array
            int nin  = info.shape[0] * 2;
            int nvol = info.shape[1];
            int nex  = info.shape[2];

            Field f(nin, nvol, nex, Element_type::COMPLEX);

            for (int k = 0; k < info.shape[2]; ++k) {
              for (int j = 0; j < info.shape[1]; ++j) {
                for (int i = 0; i < info.shape[0]; ++i) {
                  int idx = (
                      i * info.strides[0]
                    + j * info.strides[1]
                    + k * info.strides[2] ) / sizeof(dcomplex);

                  dcomplex v = ((dcomplex *) info.ptr)[idx];

                  f.set(i*2,   j, k, real(v));
                  f.set(i*2+1, j, k, imag(v));
                }
              }
            }

            return f;
          }
          else if (info.format == py::format_descriptor<double>::format()) {
            // real array
            int nin  = info.shape[0];
            int nvol = info.shape[1];
            int nex  = info.shape[2];

            Field f(nin, nvol, nex, Element_type::REAL);

            for (int k = 0; k < info.shape[2]; ++k) {
              for (int j = 0; j < info.shape[1]; ++j) {
                for (int i = 0; i < info.shape[0]; ++i) {
                  int idx = (
                      i * info.strides[0]
                    + j * info.strides[1]
                    + k * info.strides[2] ) / sizeof(double);

                  double v = ((double *) info.ptr)[idx];

                  f.set(i, j, k, v);
                }
              }
            }

            return f;
          }
          else {
            throw std::runtime_error("incompatible format. expected a double or complex array.");
          }
        }), REDIRECT)
    ;

  // blas-like routines
  m
    .def("dot", [](const Field& y, const Field& x) {
        return dot(y, x);
      })
    .def("dotc", [](const Field& y, const Field& x) {
        return dotc(y, x);
      })
    .def("axpy", [](Field& y, const double a, const Field& x) {
        return axpy(y, a, x);
      })
    .def("axpy", [](Field& y, const dcomplex a, const Field& x) {
        return axpy(y, a, x);
      })
    .def("scal", [](Field& y, const double a) {
        return scal(y, a);
      })
    .def("scal", [](Field& y, const dcomplex a) {
        return scal(y, a);
      })
    .def("copy", [](Field& y, const Field& x) {
        return copy(y, x);
      })
    .def("aypx", [](const double a, Field& y, const Field& x) {
        return aypx(a, y, x);
      })
    .def("aypx", [](const dcomplex a, Field& y, const Field& x) {
        return aypx(a, y, x);
      })
    ;

  py::class_<Field_F, Field>(m, "Field_F", R"pybdoc(
Wilson-type fermion field
)pybdoc")
    .def(py::init<>(), REDIRECT)
    .def(py::init<int,int>(),
         py::arg("nvol"),
         py::arg("nex") =1)

    .def("reset", &Field_F::reset)
    .def("clone", &Field_F::clone)

    .def_property_readonly("nc", &Field_F::nc)
    .def_property_readonly("nc2", &Field_F::nc2)
    .def_property_readonly("nd", &Field_F::nd)

    .def("vec", &Field_F::vec)
    .def("set_vec", &Field_F::set_vec)
    .def("add_vec", &Field_F::add_vec)
    .def("clear_vec", &Field_F::clear_vec)

    .def("xI", &Field_F::xI)
    .def("Ix", &Field_F::Ix)
    ;

  typedef void (*field_f_mult)(Field_F& y, const int ex,
                               const Field_G& u, int ex1,
                               const Field_F& x, int ex2);

  typedef void (*field_f_multadd)(Field_F& y, const int ex,
                                  const Field_G& u, int ex1,
                                  const Field_F& x, int ex2,
                                  const double a);

  typedef void (*field_f_multgm)(Field_F& y,
                                 const GammaMatrix& gm,
                                 const Field_F& x);
  typedef void (*field_f_multgmproj)(Field_F& y,
                                     const int pm,
                                     const GammaMatrix& gm,
                                     const Field_F& x);
  typedef void (*field_f_multgmproj2)(Field_F& y,
                                      const double nu_s,
                                      const int pm,
                                      const double r_s,
                                      const GammaMatrix& gm,
                                      const Field_F& x);

  m
    .def("mult_Field_Gn", (field_f_mult) &mult_Field_Gn)
    .def("mult_Field_Gd", (field_f_mult) &mult_Field_Gd)
    .def("multadd_Field_Gn", (field_f_multadd) &multadd_Field_Gn)
    .def("multadd_Field_Gd", (field_f_multadd) &multadd_Field_Gd)

    .def("mult_GM", (field_f_multgm) &mult_GM)
    .def("mult_iGM", (field_f_multgm) &mult_iGM)
    .def("mult_GMproj", (field_f_multgmproj) &mult_GMproj)
    .def("mult_GMproj2", (field_f_multgmproj2) &mult_GMproj2)
    ;

  py::class_<Field_G, Field>(m, "Field_G", R"pybdoc(
SU(N) gauge field
)pybdoc")
    .def(py::init<>(), REDIRECT)
    .def(py::init<int,int>(),
         py::arg("nvol"),
         py::arg("nex") =1)

    .def("reset", &Field_G::reset)
    .def("clone", &Field_G::clone)

    .def_property_readonly("nc", &Field_G::nc)

    .def("mat", [](Field_G& self, int site, int mn) {
        return self.mat(site, mn);
      },
      py::arg("site"), py::arg("mn")=0)
    .def("mat_dag", [](Field_G& self, int site, int mn) {
        return self.mat_dag(site, mn);
      },
      py::arg("site"), py::arg("mn")=0)
    .def("set_mat", &Field_G::set_mat)
    .def("add_mat", (void (Field_G::*)(const int site, const int mn, const Mat_SU_N& u)) &Field_G::add_mat)
    .def("add_mat", (void (Field_G::*)(const int site, const int mn, const Mat_SU_N& u, double prf)) &Field_G::add_mat)

    .def("xI", &Field_G::xI)
    .def("set_unit", &Field_G::set_unit)
    .def("set_random", (void (Field_G::*)(RandomNumbers *rand)) &Field_G::set_random)

    .def("reunit", &Field_G::reunit)
    ;

  m
    .def("mult_Field_Gnn", &mult_Field_Gnn)
    .def("mult_Field_Gdn", &mult_Field_Gdn)
    .def("mult_Field_Gnd", &mult_Field_Gnd)
    .def("mult_Field_Gdd", &mult_Field_Gdd)

    .def("multadd_Field_Gnn", &multadd_Field_Gnn)
    .def("multadd_Field_Gdn", &multadd_Field_Gdn)
    .def("multadd_Field_Gnd", &multadd_Field_Gnd)
    .def("multadd_Field_Gdd", &multadd_Field_Gdd)

    .def("at_Field_G", &at_Field_G)
    .def("ah_Field_G", &ah_Field_G)

    .def("mult_exp_Field_G", &mult_exp_Field_G)
    ;

#ifdef TRUNK
  py::class_<Field_F_1spinor, Field>(m, "Field_F_1spinor", R"pybdoc(
Staggered-type fermion field
)pybdoc")
    .def(py::init<>(), REDIRECT)
    .def(py::init<int,int>(),
         py::arg("nvol"),
         py::arg("nex") =1)

    .def("vec", &Field_F_1spinor::vec)
    .def("set_vec", &Field_F_1spinor::set_vec)
    .def("add_vec", &Field_F_1spinor::add_vec)
    .def("clear_vec", &Field_F_1spinor::clear_vec)
    ;

  typedef void (*field_f_1spinor_mult)(Field_F_1spinor& y, const int ex,
                                       const Field_G& u, int ex1,
                                       const Field_F_1spinor& x, int ex2);
  m
    .def("mult_Field_1spinor_Gn", (field_f_1spinor_mult) &mult_Field_Gn)
    .def("mult_Field_1spinor_Gd", (field_f_1spinor_mult) &mult_Field_Gd)
    ;  
#endif

}
