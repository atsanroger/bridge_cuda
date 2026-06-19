/*!
        @file    pyb_hmc.cpp

        @brief   python binding of HMC classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#include "pyb_defs.h"
#include "pyb_params.h"
#include "pyb_hmc.h"
#include "bridge.h"
#include "Field/field.h"
#include "HMC/hmc_General.h"
#include "HMC/integrator.h"
#include "Parameters/parameters.h"
#include "pyb_wrap.h"

// convert list of lists of actions to ActionList.
ActionList make_actionlist(const std::vector<std::vector<Action*> >& vlist)
{
  int nlevel = vlist.size();

  ActionList alist(nlevel);

  for (int lv = 0; lv < nlevel; ++lv) {
    for (auto& item: vlist[lv]) {
      alist.append(lv, item);
    }
  }

  return alist;
}

void pybind_hmc(py::module m)
{
  //! HMC General
  //
  py::class_<HMC_General>(m, "HMC_General", R"pybdoc(
General HMC update class

Parameters
----------
trajectory_length : float
Metropolis_test : {'true', 'false'}
)pybdoc")
    .def(py::init([](ActionList& actions,
                     std::vector<Director*> &directors,
                     PyHandle<Integrator> integp,
                     RandomNumbers *rand,
                     py::kwargs kwargs) -> HMC_General*
                  {
                    auto obj = new HMC_General(actions, directors, integp.get(), rand);
                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }),
         REDIRECT)

    .def(py::init([](std::vector<std::vector<Action*> > &actions,
                     std::vector<Director*> &directors,
                     PyHandle<Integrator> integp,
                     RandomNumbers *rand,
                     py::kwargs kwargs) -> HMC_General*
                  {
                    auto obj = new HMC_General(make_actionlist(actions), directors, integp.get(), rand);
                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }),
         REDIRECT)

    .def("update", &HMC_General::update, REDIRECT)

    .def("set_parameters", [](HMC_General& self, py::kwargs kwargs)
         {
           self.set_parameters(make_parameters(kwargs));
         }, REDIRECT)
    ;

  //! HMC Leapfrog
  //
  py::class_<HMC_Leapfrog>(m, "HMC_Leapfrog", R"pybdoc(
HMC with single level leapfrog integrator

Parameters
----------
number_of_steps : int
order_of_exp_iP : int
Metropolis_test : {'true', 'false'}

either or both of the following must be specified...
step_size : float
trajectory_length : float
)pybdoc")
    .def(py::init([](std::vector<Action*> &actions,
                     std::vector<Director*> &directors,
                     RandomNumbers *rand,
                     py::kwargs kwargs) -> HMC_Leapfrog*
                  {
                    auto obj = new HMC_Leapfrog(actions, directors, rand);
                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }),
         py::arg("actions"),
         py::arg("directors"),
         py::arg("rand"),
         REDIRECT)

    .def("update", &HMC_Leapfrog::update, REDIRECT)

    .def("set_parameters", [](HMC_Leapfrog& self, py::kwargs kwargs)
         {
           self.set_parameters(make_parameters(kwargs));
         }, REDIRECT)
    ;

  //! Integrator
  //
  py::class_<Integrator>(m, "Integrator", R"pybdoc(
Integrator classes
)pybdoc")
    .def(py::init([](py::kwargs kwargs) -> Integrator*
                  {
                    warn("invalid instantiation of integrator base class");
                    return nullptr;
                  }), REDIRECT)

    .def("set_parameters", [](Integrator& self, py::kwargs kwargs)
         {
           self.set_parameters(make_parameters(kwargs));
         }, REDIRECT)
    ;

  py::class_<Integrator_Leapfrog, Integrator>(m, "Integrator_Leapfrog")
    .def(py::init<Integrator*, Integrator*>(), REDIRECT)
    ;

  py::class_<Integrator_Omelyan, Integrator>(m, "Integrator_Omelyan")
    .def(py::init<Integrator*, Integrator*>(), REDIRECT)
    ;

  py::class_<Integrator_UpdateP, Integrator>(m, "Integrator_UpdateP")
    .def(py::init<std::vector<Action*> >(), REDIRECT)
    ;

  py::class_<Integrator_UpdateU, Integrator>(m, "Integrator_UpdateU")
    .def(py::init<std::vector<Director *> >(), REDIRECT)
    ;

  // nominal class definiton for pointer to integrator
  py::class_<PyHandle<Integrator> >(m, "IntegratorP")
    ;

  //! Integrator Builder
  //
  py::class_<Builder_Integrator>(m, "Builder_Integrator", R"pybdoc(
Builder of MD integrator for HMC

Parameters
----------
integrator : string, optinal
number_of_steps : int
order_of_exp_iP : int
lambda_Omelyan : float, mandatory when integrator is 'Omelyan'
)pybdoc")
    .def(py::init([](const ActionList& action_list,
                     std::vector<Director *> director,
                     py::kwargs kwargs)
                  {
                    auto obj = new Builder_Integrator(action_list, director);
                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }), REDIRECT)

    .def(py::init([](std::vector<std::vector<Action*> >& actions,
                     std::vector<Director *> director,
                     py::kwargs kwargs)
                  {
                    auto obj = new Builder_Integrator(make_actionlist(actions),
                                                      director);
                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }), REDIRECT)

    .def("set_parameters", [](Builder_Integrator& self, py::kwargs kwargs)
         {
           self.set_parameters(make_parameters(kwargs));
         }, REDIRECT)

    .def("build", [](Builder_Integrator& self) -> PyHandle<Integrator>
         {
           return PyHandle<Integrator>(self.build());
         }, REDIRECT)
    ;

  //! Action List
  //
  py::class_<ActionList>(m, "ActionList", R"pybdoc(
Lists of actions at respective integrator levels
)pybdoc")
    .def(py::init<int>(), REDIRECT)

    .def("append", [](ActionList& self, int level, Action *action)
         {
           self.append(level, action);
         }, REDIRECT)
    ;
  
  //! Langevin momentum
  //
  py::class_<Langevin_Momentum>(m, "Langevin_Momentum", R"pybdoc(
Langevin part of HMC for conjugate momentum to link variable
)pybdoc")
    .def(py::init([](py::object rand, py::kwargs kwargs)
                  {
                    auto obj = new PyWrap1<Langevin_Momentum,RandomNumbers>(rand);
                    // if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }),
         REDIRECT)

    .def("set_iP", &Langevin_Momentum::set_iP, REDIRECT)
    ;

}
