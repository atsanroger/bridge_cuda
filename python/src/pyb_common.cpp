/*!
        @file    pyb_common.cpp

        @brief   interface to CommonParameters

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#include <vector>
#include "pyb_common.h"
#include "bridge.h"

void pybind_common(py::module& m)
{
#define PROP(x) def_property_readonly_static(#x, [](py::object){ return CommonParameters:: x (); })

  py::enum_<Bridge::VerboseLevel>(m, "VerboseLevel")
    .value("Crucial", Bridge::CRUCIAL)
    .value("General", Bridge::GENERAL)
    .value("Detailed", Bridge::DETAILED)
    .value("Paranoiac", Bridge::PARANOIAC)
    .export_values()
    ;

  py::class_<CommonParameters>(m, "Common")
    .def_static("init", [](const std::vector<int>& vol
                ) {
                  CommonParameters::init(vol, {1,1,1,1});
                })
    .def_static("init", [](const std::vector<int>& vol,
                           const std::vector<int>& grid
                ) {
                  CommonParameters::init(vol, grid);
                })
    .def_static("init", [](const std::vector<int>& vol,
                           const std::vector<int>& grid,
                           const int nc
                ) {
                  CommonParameters::init(vol, grid, nc);
                })

    . PROP(Lx)
    . PROP(Ly)
    . PROP(Lz)
    . PROP(Lt)
    . PROP(Lvol)

    . PROP(NPEx)
    . PROP(NPEy)
    . PROP(NPEz)
    . PROP(NPEt)
    . PROP(NPE)

    . PROP(Nx)
    . PROP(Ny)
    . PROP(Nz)
    . PROP(Nt)
    . PROP(Nvol)

    . PROP(Nc)
    . PROP(Nd)
    . PROP(Ndim)

    . PROP(epsilon_criterion)
    . PROP(epsilon_criterion2)

    . PROP(Grid_map)

    . PROP(Vlevel)

    // number of threads from ThreadManager
    .def_property_readonly_static("Nthread", [](py::object) { return ThreadManager_OpenMP::get_num_threads_available(); })
    ;

#undef PROP
}

