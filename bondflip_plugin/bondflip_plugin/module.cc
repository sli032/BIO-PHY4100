// Copyright (c) 2009-2022 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

// Include the defined classes that are to be exported to python
#include "BondFlipUpdater.h"

#include <pybind11/pybind11.h>

using namespace hoomd::detail;

// specify the python module. Note that the name must explicitly match the PROJECT() name provided
// in CMakeLists (with an underscore in front)
PYBIND11_MODULE(_bondflip_plugin, m)
    {
    export_BondFlipUpdater(m);

#ifdef ENABLE_HIP
    export_BondFlipUpdaterGPU(m);
#endif
    }
