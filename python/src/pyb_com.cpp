/*!
        @file    pyb_com.cpp

        @brief   interface to Communicator

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#include "pyb_defs.h"
#include "pyb_params.h"
#include "pyb_com.h"
#include "bridge.h"
#include "Communicator/communicator.h"
#ifdef USE_MPI
#include "Communicator/MPI/communicator_mpi.h"
#endif
#include "Parameters/parameters.h"

#ifdef USE_MPI
#ifdef USE_MPI4PY

#include <mpi.h>
#include <mpi4py/mpi4py.h>

// interoperability with mpi4py python MPI module:
// passing MPI communicator in bridge++ to python mpi4py
//
// cf. https://stackoverflow.com/questions/49259704/pybind11-possible-to-use-mpi4py
//

// wrapper for MPI_Comm
//
struct mpi4py_comm
{
  mpi4py_comm() =default;
  mpi4py_comm(MPI_Comm value_) : value(value_) {}

  operator MPI_Comm () { return value; }

  MPI_Comm value;
};

// type_caster for mpi4py_comm
//
namespace pybind11
{
  namespace detail
  {
    template<>
    struct type_caster<mpi4py_comm>
    {
     public:
      PYBIND11_TYPE_CASTER(mpi4py_comm, _("mpi4py_comm"));

      // Python -> C++
      bool load(handle src, bool)
      {
        PyObject *py_src = src.ptr();

        // check that an mpi4py communicator has been passed.
        if (PyObject_TypeCheck(py_src, &PyMPIComm_Type)) {

          // convert to regular MPI communicator
          value.value = *PyMPIComm_Get(py_src);
        } else {
          return false;
        }
        return !PyErr_Occurred();
      }

      // C++ -> Python
      static handle cast(mpi4py_comm src,
                         return_value_policy,
                         handle)
      {
        // Create an mpi4py handle
        return PyMPIComm_New(src.value);
      }
    };
  } // namespace detail
} // namespace type_caster

#endif
#endif


// A helper class that tranfers requests to communicator
//
class PyCommunicator
{
 public:
  static bool mpi_enabled()
  {
#ifdef USE_MPI
    return true;
#else
    return false;
#endif
  }

  static int self() { return Communicator::self(); }
  static int size() { return Communicator::size(); }

  static int ipe(const int idir) { return Communicator::ipe(idir); }
  static int npe(const int idir) { return Communicator::npe(idir); }
  
#ifdef USE_MPI
#ifdef USE_MPI4PY
  static mpi4py_comm get_comm()
  {
    return Communicator_impl::world();
  }
#endif
#endif
};


void pybind_com(py::module m)
{
#ifdef USE_MPI
#ifdef USE_MPI4PY
  // import mpi4py API
  //
  if (import_mpi4py() < 0) {
    throw std::runtime_error("Could not load mpi4py API.");
  }
#endif
#endif
  
  //! Communicator
  //
  py::class_<PyCommunicator>(m, "Communicator")
    .def_static("mpi_enabled", &PyCommunicator::mpi_enabled, REDIRECT)

    .def_static("self", &PyCommunicator::self, REDIRECT)
    .def_static("size", &PyCommunicator::size, REDIRECT)

    .def_static("ipe", &PyCommunicator::ipe, REDIRECT)
    .def_static("npe", &PyCommunicator::npe, REDIRECT)
    
#ifdef USE_MPI
#ifdef USE_MPI4PY
    .def_static("get_comm", &PyCommunicator::get_comm, REDIRECT)
#endif
#endif
    ;
}
