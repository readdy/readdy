/********************************************************************
 * Copyright © 2016 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * This file is part of ReaDDy.                                     *
 *                                                                  *
 * ReaDDy is free software: you can redistribute it and/or modify   *
 * it under the terms of the GNU Lesser General Public License as   *
 * published by the Free Software Foundation, either version 3 of   *
 * the License, or (at your option) any later version.              *
 *                                                                  *
 * This program is distributed in the hope that it will be useful,  *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of   *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
 * GNU Lesser General Public License for more details.              *
 *                                                                  *
 * You should have received a copy of the GNU Lesser General        *
 * Public License along with this program. If not, see              *
 * <http://www.gnu.org/licenses/>.                                  *
 ********************************************************************/


/**
 * << detailed description >>
 *
 * @file KernelWrap.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 03.08.16
 */

#include <readdy/model/Kernel.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace readdy {
namespace rpy {
struct KernelWrap : public readdy::model::Kernel {
    KernelWrap(const std::string &name) : Kernel(name) {}

    virtual readdy::model::actions::ActionFactory &getActionFactory() const override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_PURE(readdy::model::actions::ActionFactory &, readdy::model::Kernel, getActionFactory,);
    }

    virtual readdy::model::KernelStateModel &getKernelStateModel() const override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_PURE(readdy::model::KernelStateModel &, readdy::model::Kernel, getKernelStateModel,);
    }

    virtual readdy::model::KernelContext &getKernelContext() const override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_PURE(readdy::model::KernelContext &, readdy::model::Kernel, getKernelContext,);
    }

    virtual readdy::model::potentials::PotentialFactory &getPotentialFactory() const override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_PURE(
                readdy::model::potentials::PotentialFactory &, readdy::model::Kernel, getPotentialFactory,
        );
    }

    virtual readdy::model::reactions::ReactionFactory &getReactionFactory() const override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_PURE(readdy::model::reactions::ReactionFactory &, readdy::model::Kernel, getReactionFactory,)
    }

};

struct SCPUKernelWrap : public readdy::kernel::scpu::SCPUKernel {
    virtual kernel::scpu::SCPUStateModel &getKernelStateModel() const override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_NAME(kernel::scpu::SCPUStateModel &,
                               readdy::kernel::scpu::SCPUKernel, "get_kernel_state_model",
                               getKernelStateModel,);
    }

    virtual model::KernelContext &getKernelContext() const override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_NAME(model::KernelContext &, readdy::kernel::scpu::SCPUKernel, "get_kernel_context",
                               getKernelContext,);
    }

    virtual model::actions::ActionFactory &getActionFactory() const override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_NAME(model::actions::ActionFactory &, readdy::kernel::scpu::SCPUKernel,
                               "get_program_factory", getActionFactory,);
    }

    virtual std::vector<std::string> getAvailablePotentials() const override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_NAME(std::vector<std::string>, readdy::kernel::scpu::SCPUKernel,
                               "get_available_potentials", getAvailablePotentials,);
    }

    virtual model::potentials::PotentialFactory &getPotentialFactory() const override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_NAME(model::potentials::PotentialFactory &, readdy::kernel::scpu::SCPUKernel,
                               "get_potential_factory", getPotentialFactory,);
    }

    virtual model::reactions::ReactionFactory &getReactionFactory() const override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_NAME(model::reactions::ReactionFactory &, readdy::kernel::scpu::SCPUKernel,
                               "get_reaction_factory", getReactionFactory,);
    }

    virtual model::observables::ObservableFactory &getObservableFactory() const override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_NAME(model::observables::ObservableFactory &, readdy::kernel::scpu::SCPUKernel,
                               "get_observable_factory", getObservableFactory,);
    }
};
}
}
