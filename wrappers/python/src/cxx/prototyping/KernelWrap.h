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

    virtual readdy::model::actions::ActionFactory &getActionFactoryInternal() const override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_PURE(readdy::model::actions::ActionFactory &, readdy::model::Kernel, getActionFactoryInternal,);
    }

    virtual readdy::model::KernelStateModel &getKernelStateModelInternal() const override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_PURE(readdy::model::KernelStateModel &, readdy::model::Kernel, getKernelStateModelInternal,);
    }

    virtual readdy::model::KernelContext &getKernelContextInternal() const override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_PURE(readdy::model::KernelContext &, readdy::model::Kernel, getKernelContextInternal,);
    }

    virtual readdy::model::potentials::PotentialFactory &getPotentialFactoryInternal() const override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_PURE(
                readdy::model::potentials::PotentialFactory &, readdy::model::Kernel, getPotentialFactoryInternal,
        );
    }

    virtual readdy::model::reactions::ReactionFactory &getReactionFactoryInternal() const override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_PURE(readdy::model::reactions::ReactionFactory &, readdy::model::Kernel, getReactionFactoryInternal,)
    }

};

struct SCPUKernelWrap : public readdy::kernel::scpu::SCPUKernel {
    virtual kernel::scpu::SCPUStateModel &getKernelStateModelInternal() const override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_NAME(kernel::scpu::SCPUStateModel &,
                               readdy::kernel::scpu::SCPUKernel, "get_kernel_state_model",
                               getKernelStateModelInternal,);
    }

    virtual model::KernelContext &getKernelContextInternal() const override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_NAME(model::KernelContext &, readdy::kernel::scpu::SCPUKernel, "get_kernel_context",
                               getKernelContextInternal,);
    }

    virtual model::actions::ActionFactory &getActionFactoryInternal() const override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_NAME(model::actions::ActionFactory &, readdy::kernel::scpu::SCPUKernel,
                               "get_program_factory", getActionFactoryInternal,);
    }

    virtual std::vector<std::string> getAvailablePotentials() const override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_NAME(std::vector<std::string>, readdy::kernel::scpu::SCPUKernel,
                               "get_available_potentials", getAvailablePotentials,);
    }

    virtual model::potentials::PotentialFactory &getPotentialFactoryInternal() const override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_NAME(model::potentials::PotentialFactory &, readdy::kernel::scpu::SCPUKernel,
                               "get_potential_factory", getPotentialFactoryInternal,);
    }

    virtual model::reactions::ReactionFactory &getReactionFactoryInternal() const override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_NAME(model::reactions::ReactionFactory &, readdy::kernel::scpu::SCPUKernel,
                               "get_reaction_factory", getReactionFactoryInternal,);
    }

    virtual model::observables::ObservableFactory &getObservableFactoryInternal() const override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_NAME(model::observables::ObservableFactory &, readdy::kernel::scpu::SCPUKernel,
                               "get_observable_factory", getObservableFactoryInternal,);
    }
};
}
}
