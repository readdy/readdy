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


//
// Created by clonker on 07.03.16.
//

#ifndef READDY_MAIN_SINGLECPUKERNEL_H
#define READDY_MAIN_SINGLECPUKERNEL_H

#include <readdy/model/RandomProvider.h>
#include <readdy/model/Kernel.h>
#include <readdy/kernel/singlecpu/SCPUStateModel.h>

namespace readdy {
namespace kernel {
namespace scpu {

class SCPUKernel : public readdy::model::Kernel {
public:

    static const std::string name;

    SCPUKernel();

    ~SCPUKernel();

    // move
    SCPUKernel(SCPUKernel &&rhs);

    SCPUKernel &operator=(SCPUKernel &&rhs);

    // factory method
    static std::unique_ptr<SCPUKernel> create();

    virtual SCPUStateModel &getKernelStateModel() const override;

    virtual readdy::model::KernelContext &getKernelContext() const override;

    virtual readdy::model::actions::ActionFactory &getActionFactory() const override;

    virtual std::vector<std::string> getAvailablePotentials() const override;

    virtual std::unique_ptr<readdy::model::potentials::Potential> createPotential(std::string &name) const override;

    virtual readdy::model::potentials::PotentialFactory &getPotentialFactory() const override;

    virtual readdy::model::reactions::ReactionFactory &getReactionFactory() const override;

    virtual readdy::model::observables::ObservableFactory &getObservableFactory() const override;

private:
    struct Impl;
    std::unique_ptr<readdy::kernel::scpu::SCPUKernel::Impl> pimpl;
};

}
}
}

#endif //READDY_MAIN_SINGLECPUKERNEL_H
