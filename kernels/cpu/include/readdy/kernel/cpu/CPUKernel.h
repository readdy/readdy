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
 * @file CPUKernel.h
 * @brief << brief description >>
 * @author clonker
 * @date 23.06.16
 */


#pragma once

#include <readdy/model/Kernel.h>
#include <readdy/common/dll.h>
#include <readdy/common/thread/ctpl.h>
#include <readdy/common/thread/executor.h>
#include "CPUStateModel.h"

namespace readdy {
namespace kernel {
namespace cpu {


class CPUKernel : public readdy::model::Kernel {
public:
    static const std::string name;

    CPUKernel();

    ~CPUKernel() override;

    CPUKernel(const CPUKernel&) = delete;
    CPUKernel& operator=(const CPUKernel&) = delete;
    CPUKernel(CPUKernel&&) = delete;
    CPUKernel& operator=(CPUKernel&&) = delete;

    // factory method
    static readdy::model::Kernel* create();

    unsigned long getNThreads() const;

    void setNThreads(readdy::util::thread::Config::n_threads_type n);

    const CPUStateModel &getCPUKernelStateModel() const;
    CPUStateModel &getCPUKernelStateModel();

    const readdy::util::thread::Config& threadConfig() const;
    readdy::util::thread::Config& threadConfig();

    const readdy::util::thread::executor_base &executor() const;

    void initialize() override;

    void finalize() override;

protected:

    readdy::model::observables::ObservableFactory &getObservableFactoryInternal() const override;

    CPUStateModel &getKernelStateModelInternal() const override;

    readdy::model::KernelContext &getKernelContextInternal() const override;

    readdy::model::actions::ActionFactory &getActionFactoryInternal() const override;

    readdy::model::potentials::PotentialFactory &getPotentialFactoryInternal() const override;

    readdy::model::reactions::ReactionFactory &getReactionFactoryInternal() const override;

    readdy::model::compartments::CompartmentFactory &getCompartmentFactoryInternal() const override;

    readdy::model::top::TopologyActionFactory *getTopologyActionFactoryInternal() const override;

    struct Impl;
    std::unique_ptr<Impl> pimpl;
};

}
}
}

extern "C" const char* name();

extern "C" readdy::model::Kernel* createKernel();
