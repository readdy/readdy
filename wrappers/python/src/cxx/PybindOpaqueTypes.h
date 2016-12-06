/**
 * << detailed description >>
 *
 * @file PybindOpaqueTypes.h
 * @brief << brief description >>
 * @author clonker
 * @date 05.10.16
 */

#ifndef READDY_MAIN_PYBINDOPAQUETYPES_H
#define READDY_MAIN_PYBINDOPAQUETYPES_H

#include <pybind11/cast.h>

#include <vector>
#include <readdy/kernel/singlecpu/model/SCPUNeighborList.h>

PYBIND11_MAKE_OPAQUE(std::vector<double>);
PYBIND11_MAKE_OPAQUE(std::vector<unsigned long>);
PYBIND11_MAKE_OPAQUE(std::vector<readdy::kernel::singlecpu::model::Box>);
PYBIND11_MAKE_OPAQUE(std::vector<readdy::model::Vec3>);

#endif //READDY_MAIN_PYBINDOPAQUETYPES_H
