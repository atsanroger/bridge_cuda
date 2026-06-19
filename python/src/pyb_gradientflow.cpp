/*!
        @file    pyb_gradientflow.cpp

        @brief   python binding of gradient flow classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#include <tuple>
#include "pyb_defs.h"
#include "pyb_params.h"
#include "pyb_gradientflow.h"
#include "bridge.h"
#include "Field/field.h"
#include "Measurements/Gauge/gradientFlow.h"
#include "Parameters/parameters.h"
#include "pyb_wrap.h"

#ifdef TRUNK

class
__attribute__((visibility("hidden")))
PyFermionFlow : public FermionFlow
{
 public:
  PyFermionFlow(std::unique_ptr<Action>& action, py::object obj)
    : FermionFlow(action), obj_(obj)
  { ENTER;
    obj_.inc_ref();
  }
  ~PyFermionFlow()
  { ENTER;
    obj_.dec_ref();
  }

 private:
  py::object obj_;
};

class
__attribute__((visibility("hidden")))
PyFermionFlow_1pt_Function : public FermionFlow_1pt_Function
{
 public:
  PyFermionFlow_1pt_Function(
    const std::vector<Fprop*>& fprop,
    std::unique_ptr<Action>& action,
    std::vector<Parameters>& params,
    py::object action_)
    : FermionFlow_1pt_Function(fprop, action, &params_)
    , obj_(action_)
    , params_(params)
  { ENTER;
    obj_.inc_ref();
  }
  ~PyFermionFlow_1pt_Function()
  { ENTER;
    obj_.dec_ref();
  }

 private:
  py::object obj_;
  std::vector<Parameters> params_;
};

class
__attribute__((visibility("hidden")))
PyFermionFlow_2pt_Function : public FermionFlow_2pt_Function
{
 public:
  PyFermionFlow_2pt_Function(
    const std::vector<Fprop*>& fprop,
    std::unique_ptr<Action>& action,
    std::vector<Parameters>& params,
    py::object action_)
    : FermionFlow_2pt_Function(fprop, action, &params_)
    , obj_(action_)
    , params_(params)
  { ENTER;
    obj_.inc_ref();
  }
  ~PyFermionFlow_2pt_Function()
  { ENTER;
    obj_.dec_ref();
  }

 private:
  py::object obj_;
  std::vector<Parameters> params_;
};

#endif


void pybind_gradientflow(py::module m)
{
  //! GradientFlow
  //
  py::class_<GradientFlow>(m, "GradientFlow", R"pybdoc(
Gradient flow class

Parameters
----------
order_of_RungeKutta : int
step_size : float
order_of_approx_for_exp_iP : int
adaptive : {'true', 'false'}
tolerance : float, mandatory if adaptive=true
safety_factor : float, mandatory if adaptive=true
)pybdoc")
    .def(py::init([](py::object action, py::kwargs kwargs) -> GradientFlow*
                  {
                    auto obj = new PyWrap1<GradientFlow, Action>(action);

                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }), REDIRECT)

    .def("set_parameters", [](GradientFlow& self, py::kwargs kwargs)
         {
           self.set_parameters(make_parameters(kwargs));
         }, REDIRECT)

    .def("evolve", [](GradientFlow& self, double t, Field_G& U) -> std::tuple<double, double>
         {
           double tt = t;
           double r = self.evolve(tt, U);
           return std::tuple<double,double>(tt, r);
         }, REDIRECT)
    ;

#ifdef TRUNK

  //! FermionFlow
  //
  py::class_<FermionFlow>(m, "FermionFlow", R"pybdoc(
Fermion flow class

Parameters
----------
order_of_RungeKutta : int
step_size : float
order_of_approx_for_exp_iP : int
boundary_condition : List [int]
)pybdoc")
    .def(py::init([](py::object arg, py::kwargs kwargs) -> FermionFlow*
                  {
                    std::unique_ptr<Action> action_(arg.cast<Action*>());
                    auto obj = new PyFermionFlow(action_, arg);
                    action_.release();

                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }),
         py::arg("action"),
         REDIRECT)

    .def("laplacian",
         (void (FermionFlow::*)(Field_F&, const Field_F&, Field_G&)) &FermionFlow::laplacian,
         REDIRECT)
    .def("laplacian",
         (void (FermionFlow::*)(Field_F_1spinor&, const Field_F_1spinor&, Field_G&)) &FermionFlow::laplacian,
         REDIRECT)

    .def("del_symmetric",
         (void (FermionFlow::*)(int, Field_F&, const Field_F&, Field_G&)) &FermionFlow::del_symmetric,
         REDIRECT)
    .def("del_symmetric",
         (void (FermionFlow::*)(int, Field_F_1spinor&, const Field_F_1spinor&, Field_G&)) &FermionFlow::del_symmetric,
         REDIRECT)

    // .def("fermionFlow_1st", &FermionFlow::fermionFlow_1st, REDIRECT)
    // .def("fermionFlow_2nd", &FermionFlow::fermionFlow_2nd, REDIRECT)
    .def("fermionFlow_3rd",
         (void (FermionFlow::*)(Field_F&, Field_G&)) &FermionFlow::fermionFlow_3rd,
         REDIRECT)
    .def("fermionFlow_3rd",
         (void (FermionFlow::*)(Field_F_1spinor&, Field_G&)) &FermionFlow::fermionFlow_3rd,
         REDIRECT)
    // .def("fermionFlow_4th", &FermionFlow::fermionFlow_4th, REDIRECT)
    
    .def("fermionFlow_normal_order_3rd", &FermionFlow::fermionFlow_normal_order_3rd, REDIRECT)

    .def("evolve",
         (void (FermionFlow::*)(Field_F&, Field_G&)) &FermionFlow::evolve,
         REDIRECT)
    .def("evolve",
         (void (FermionFlow::*)(Field_F_1spinor&, Field_G&)) &FermionFlow::evolve,
         REDIRECT)
    .def("evolve",
         (void (FermionFlow::*)(int, Field_F&, Field_G&)) &FermionFlow::evolve,
         REDIRECT)

    .def("evolve2", &FermionFlow::evolve2, REDIRECT)

    .def("evolve_normal_order", &FermionFlow::evolve_normal_order, REDIRECT)
    ;

  //! FermionFlow 1pt function
  //
  py::class_<FermionFlow_1pt_Function>(m, "FermionFlow_1pt_Function", R"pybdoc(
Measures quark one-point function with fermion flow

Parameters
----------
number_of_noises : int
initial_tau : int
number_of_measurement_times : int
measurement_interval : int
gauge_store_interval : int
max_momentum : int
step_size : float
noise_type : string
filename_output : string, optional
)pybdoc")
    .def(py::init([](const std::vector<Fprop*>& fprop,
                     py::object action_,
                     std::vector<py::dict> params_,
                     py::dict dict_measurement,
                     py::dict dict_gflow,
                     py::dict dict_fflow,
                     py::dict dict_source_random
                  ) -> FermionFlow_1pt_Function*
                  {
                    std::vector<Parameters> params;
                    for (auto x : params_) {
                      params.push_back(make_parameters(x));
                    }

                    std::unique_ptr<Action> action(action_.cast<Action*>());
                    auto obj = new PyFermionFlow_1pt_Function(fprop, action, params, action_);
                    action.release();
                    
                    Parameters params_measurement = make_parameters(dict_measurement);
                    Parameters params_gflow = make_parameters(dict_gflow);
                    Parameters params_fflow = make_parameters(dict_fflow);
                    Parameters params_source_random = make_parameters(dict_source_random);

                    if (obj) obj->set_parameters(
                      params_measurement,
                      params_gflow,
                      params_fflow,
                      params_source_random);

                    return obj;
                  }),
         py::arg("fprops"),
         py::arg("action"),
         py::arg("params"),
         py::arg("params_measurement"),
         py::arg("params_gflow"),
         py::arg("params_fflow"),
         py::arg("params_source_random"),
         REDIRECT)

    .def("set_parameters", [](FermionFlow_1pt_Function& self,
                              py::dict dict_measurement,
                              py::dict dict_gflow,
                              py::dict dict_fflow,
                              py::dict dict_source_random)
         {
           Parameters params_measurement = make_parameters(dict_measurement);
           Parameters params_gflow = make_parameters(dict_gflow);
           Parameters params_fflow = make_parameters(dict_fflow);
           Parameters params_source_random = make_parameters(dict_source_random);
           return self.set_parameters(
             params_measurement,
             params_gflow,
             params_fflow,
             params_source_random);
         },
         py::arg("params_measurement"),
         py::arg("params_gradient_flow"),
         py::arg("params_fermion_flow"),
         py::arg("params_source_random"),
         REDIRECT)

#define DEFINE_METHOD(x) def(#x, &FermionFlow_1pt_Function:: x, REDIRECT)

    .DEFINE_METHOD(measure_disconnected)

    .DEFINE_METHOD(print_vev)
    .DEFINE_METHOD(print_correlator_t)
    .DEFINE_METHOD(print_correlator_t_FT)
    .DEFINE_METHOD(print_correlator_x)
    .DEFINE_METHOD(print_correlator_x_FT)

    // yet unimplemented
    // .DEFINE_METHOD(print_correlator_y)
    // .DEFINE_METHOD(print_correlator_y_FT)
    // .DEFINE_METHOD(print_correlator_z)
    // .DEFINE_METHOD(print_correlator_z_FT)

#undef DEFINE_METHOD
    ;

  //! FermionFlow 2pt function
  //
  using std::vector;
  
  py::class_<FermionFlow_2pt_Function>(m, "FermionFlow_2pt_Function", R"pybdoc(
Measures quark two-point function with fermion flow

Parameters
----------
initial_tau : int
number_of_measurement_times : int
measurement_interval : int
gauge_store_interval : int
max_momentum : int
step_size : float
filename_output : string, optional
)pybdoc")
    .def(py::init([](const std::vector<Fprop*>& fprop,
                     py::object action_,
                     std::vector<py::dict> params_,
                     py::dict dict_measurement,
                     py::dict dict_gflow,
                     py::dict dict_fflow
                  ) -> FermionFlow_2pt_Function*
                  {
                    std::vector<Parameters> params;
                    for (auto x : params_) {
                      params.push_back(make_parameters(x));
                    }

                    std::unique_ptr<Action> action(action_.cast<Action*>());
                    auto obj = new PyFermionFlow_2pt_Function(fprop, action, params, action_);
                    action.release();

                    Parameters params_measurement = make_parameters(dict_measurement);
                    Parameters params_gflow = make_parameters(dict_gflow);
                    Parameters params_fflow = make_parameters(dict_fflow);

                    if (obj) obj->set_parameters(
                      params_measurement,
                      params_gflow,
                      params_fflow);
                    
                    return obj;
                  }),
         py::arg("fprops"),
         py::arg("action"),
         py::arg("params"),
         py::arg("params_measurement"),
         py::arg("params_gflow"),
         py::arg("params_fflow"),
         REDIRECT)

    .def("set_parameters", [](FermionFlow_2pt_Function& self,
                              py::dict dict_measurement,
                              py::dict dict_gflow,
                              py::dict dict_fflow)
         {
           Parameters params_measurement = make_parameters(dict_measurement);
           Parameters params_gflow = make_parameters(dict_gflow);
           Parameters params_fflow = make_parameters(dict_fflow);
           return self.set_parameters(
             params_measurement,
             params_gflow,
             params_fflow);
         },
         py::arg("params_measurement"),
         py::arg("params_gradient_flow"),
         py::arg("params_fermion_flow"),
         REDIRECT)

    .def("measure_meson_correlator", &FermionFlow_2pt_Function::measure_meson_correlator, REDIRECT)

    .def("measure_EMT_correlator", &FermionFlow_2pt_Function::measure_EMT_correlator, REDIRECT)

    // method: print_meson_correlator_t
    //
    .def("print_meson_correlator_t",
         (double (FermionFlow_2pt_Function::*)
          (
            double,
            const vector<vector<Field_F> >&
          )
         )
         &FermionFlow_2pt_Function::print_meson_correlator_t,
         REDIRECT)

    .def("print_meson_correlator_t",
         (double (FermionFlow_2pt_Function::*)
          (
            double,
            const vector<vector<Field_F> >&,
            const vector<Field_F>&
          )
         )
         &FermionFlow_2pt_Function::print_meson_correlator_t,
         REDIRECT)

    // method: print_meson_correlator_t_FT
    //
    .def("print_meson_correlator_t_FT",
         (double (FermionFlow_2pt_Function::*)
          (
            double,
            const vector<vector<Field_F> >&
          )
         )
         &FermionFlow_2pt_Function::print_meson_correlator_t_FT,
         REDIRECT)

    .def("print_meson_correlator_t_FT",
         (double (FermionFlow_2pt_Function::*)
          (
            double,
            const vector<vector<Field_F> >&,
            const vector<Field_F>&
          )
         )
         &FermionFlow_2pt_Function::print_meson_correlator_t_FT,
         REDIRECT)

    // method: print_meson_correlator_x
    //
    .def("print_meson_correlator_x",
         (double (FermionFlow_2pt_Function::*)
          (
            double,
            const vector<vector<Field_F> >&
          )
         )
         &FermionFlow_2pt_Function::print_meson_correlator_x,
         REDIRECT)

    .def("print_meson_correlator_x",
         (double (FermionFlow_2pt_Function::*)
          (
            double,
            const vector<vector<Field_F> >&,
            const vector<Field_F>&
          )
         )
         &FermionFlow_2pt_Function::print_meson_correlator_x,
         REDIRECT)

    // method: print_meson_correlator_x_FT
    //
    .def("print_meson_correlator_x_FT",
         (double (FermionFlow_2pt_Function::*)
          (
            double,
            const vector<vector<Field_F> >&
          )
         )
         &FermionFlow_2pt_Function::print_meson_correlator_x_FT,
         REDIRECT)

    .def("print_meson_correlator_x_FT",
         (double (FermionFlow_2pt_Function::*)
          (
            double,
            const vector<vector<Field_F> >&,
            const vector<Field_F>&
          )
         )
         &FermionFlow_2pt_Function::print_meson_correlator_x_FT,
         REDIRECT)

#define DEFINE_METHOD(x) def(#x, &FermionFlow_2pt_Function:: x, REDIRECT)

    // method: print_O3O3a_correlator_{t,x}
    .DEFINE_METHOD(print_O3O3a_correlator_t)
    .DEFINE_METHOD(print_O3O3a_correlator_x)

    // method: print_O3O3a_correlator_{t,x}_FT
    .DEFINE_METHOD(print_O3O3a_correlator_t_FT)
    .DEFINE_METHOD(print_O3O3a_correlator_x_FT)

    // method: print_O3O3b_correlator_{t,x}
    .DEFINE_METHOD(print_O3O3b_correlator_t)
    .DEFINE_METHOD(print_O3O3b_correlator_x)

    // method: print_O3O3b_correlator_{t,x}_FT
    .DEFINE_METHOD(print_O3O3b_correlator_t_FT)
    .DEFINE_METHOD(print_O3O3b_correlator_x_FT)

    // method: print_O3O5_correlator_{t,x}
    .DEFINE_METHOD(print_O3O5_correlator_t)
    .DEFINE_METHOD(print_O3O5_correlator_x)

    // method: print_O3O5_correlator_{t,x}_FT
    .DEFINE_METHOD(print_O3O5_correlator_t_FT)
    .DEFINE_METHOD(print_O3O5_correlator_x_FT)

    // method: print_O5O3_correlator_{t,x}
    .DEFINE_METHOD(print_O5O3_correlator_t)
    .DEFINE_METHOD(print_O5O3_correlator_x)

    // method: print_O5O3_correlator_{t,x}_FT
    .DEFINE_METHOD(print_O5O3_correlator_t_FT)
    .DEFINE_METHOD(print_O5O3_correlator_x_FT)

#undef DEFINE_METHOD
    ;


#endif
}
