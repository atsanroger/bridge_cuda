/*!
        @file    pyb.cpp

        @brief   pybridge module definition

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
using std::vector;
#include <map>
#include <string>
using std::string;

#include <memory>
using std::unique_ptr;

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/iostream.h>
namespace py = pybind11;
// using namespace pybind11::literals;

#include "bridge.h"
using Bridge::vout;

#include "Tools/randomNumberManager.h"

#include "pyb_defs.h"
#include "pyb_utils.h"
#include "pyb_common.h"
#include "pyb_params.h"
#include "pyb_tools_su_n.h"
#include "pyb_field.h"
#include "pyb_index.h"
#include "pyb_fopr.h"
#include "pyb_solver.h"
#include "pyb_smear.h"
#include "pyb_shift_field.h"
#include "pyb_gauge_fixing.h"
#include "pyb_wilsonloop.h"
#include "pyb_source.h"
#include "pyb_fprop.h"
#include "pyb_gammamatrix.h"
#include "pyb_correlator.h"
#include "pyb_staple.h"
#include "pyb_action.h"
#include "pyb_force.h"
#include "pyb_hmc.h"
#include "pyb_rng.h"
#include "pyb_fft.h"
#include "pyb_gradientflow.h"
#include "pyb_topo.h"
#include "pyb_qsuscept.h"
#include "pyb_eigen.h"
#include "pyb_io.h"
#include "pyb_com.h"
#include "pyb_tools.h"

//! initialize bridge
//
// 1. bridge_initialize() is called automatically on importing the module.
//
// 2. setup, which wraps bridge_setup(), must be called by the user
//    with appropriate parameters.
//

void pyb_bridge_setup(const std::vector<int>& lattice_size,
                      const std::vector<int>& grid_size,
                      const int number_of_threads,
                      const int number_of_colors,
                      const std::string& logfile,
                      const std::string& ildg_logfile,
                      const std::string& verbose_level,
                      const std::string& random_number_type,
                      const unsigned long random_number_seed)
{
  bridge_setup(lattice_size,
               grid_size,
               number_of_threads,
               number_of_colors,
               logfile,
               ildg_logfile,
               verbose_level);

  // initialize random number manager
  RandomNumberManager::initialize(random_number_type, random_number_seed);

  return;
}


PYBIND11_MODULE(pybridge, m) {
  m.doc() = R"pybdoc(
    python binding for Bridge++ Lattice QCD codeset
  )pybdoc";

  py::add_ostream_redirect(m, "ostream_redirect");
  
  // call initializer
  {
    py::scoped_ostream_redirect sout(vout.getStream(), py::module::import("sys").attr("stdout"));
    bridge_initialize(nullptr, nullptr);
  }

  // setup lattice size (call by user)
  m.def("setup",
        &pyb_bridge_setup,
        R"pybdoc(
Setup basic parameters

Parameters
----------
lattice_size : array
    Size of the 4-dimensional lattice, given as an array of int.
grid_size : array, optional
    Number of division over parallel nodes, along each direction.
    It is valid when MPI is enabled. The default is [1,1,1,1].
    Automatic division is applied if an empty array [] is specified.
number_of_threads : int, optional
    Number of OpenMP threads. It is valid when multithreading is
    enabled. The default is 1.
number_of_colors : int, optional
    Number of colors, namely gauge group. The default is 3.
logfile : string, optional
    File name for log output.
ildg_logfile : string, optional
    File name for log output in ILDG metadata format.
verbose_level : string, optional
    Verbosity of log output. 'Crucial', 'General', 'Detailed',
    or 'Paranoiac'.
random_number_type : string, optional
    Algorithm for random number generator. The default is 'MT19937'.
random_number_seed : int, optional
    Seed for random number generator.

        )pybdoc",
        py::arg("lattice_size"),
        py::arg("grid_size") = std::vector<int>({ 1,1,1,1 }),
        py::arg("number_of_threads") = 1,
        py::arg("number_of_colors") = 3,
        py::arg("logfile") = "stdout",
        py::arg("ildg_logfile") = "stdout",
        py::arg("verbose_level") = "General",
        py::arg("random_number_type") = "MT19937",
        py::arg("random_number_seed") = 1234567UL,
        REDIRECT
  );
  
  // set finalizer
  auto atexit = py::module::import("atexit");
  atexit.attr("register")(py::cpp_function([]() {
        py::scoped_ostream_redirect sout(vout.getStream(), py::module::import("sys").attr("stdout"));

        RandomNumberManager::finalize();
        bridge_finalize();
      }));

  // setup classes
  pybind_common(m);
  pybind_params(m);
  pybind_tools_su_n(m);
  pybind_field(m);
  pybind_index(m);
  pybind_fopr(m);
  pybind_solver(m);
  pybind_smear(m);
  pybind_shift_field(m);
  pybind_gauge_fixing(m);
  pybind_wilsonloop(m);
  pybind_source(m);
  pybind_fprop(m);
  pybind_gammamatrix(m);
  pybind_correlator(m);
  pybind_staple(m);
  pybind_action(m);
  pybind_force(m);
  pybind_hmc(m);
  pybind_rng(m);
  pybind_fft(m);
  pybind_gradientflow(m);
  pybind_topo(m);
  pybind_measure_qsuscept(m);
  pybind_eigen(m);
  pybind_io(m);
  pybind_com(m);

  pybind_module_tools(m);
  

}

#ifdef MAIN

int main(int argc, char *argv[])
{
  return 0;
}

#endif
