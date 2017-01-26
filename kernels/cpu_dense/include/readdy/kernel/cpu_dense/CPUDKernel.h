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
 * @file Kernel.h
 * @brief << brief description >>
 * @author clonker
 * @date 22.11.16
 */

#ifndef READDY_DENSE_KERNEL_H
#define READDY_DENSE_KERNEL_H

#include <readdy/model/Kernel.h>
#include <readdy/common/thread/Config.h>

#include "CPUDStateModel.h"

namespace readdy {
namespace kernel {
namespace cpu_dense {

class CPUDKernel : public readdy::model::Kernel {
public:
    static const std::string name;

    CPUDKernel();

    ~CPUDKernel();

    // factory method
    static readdy::model::Kernel* create();

    virtual readdy::model::actions::ActionFactory &getActionFactory() const override;

    virtual CPUDStateModel &getKernelStateModel() const override;

    virtual readdy::model::KernelContext &getKernelContext() const override;

    virtual readdy::model::potentials::PotentialFactory &getPotentialFactory() const override;

    virtual readdy::model::reactions::ReactionFactory &getReactionFactory() const override;

    virtual readdy::model::observables::ObservableFactory &getObservableFactory() const override;

    virtual readdy::model::compartments::CompartmentFactory &getCompartmentFactory() const override;

    virtual readdy::model::top::TopologyActionFactory *getTopologyActionFactory() const override;

    virtual readdy::model::compartments::CompartmentFactory &getCompartmentFactory() const override;

    unsigned long getNThreads() const;

    void setNThreads(readdy::util::thread::Config::n_threads_t n);

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl;
};

}
}
}

extern "C" const char* name();

extern "C" readdy::model::Kernel* createKernel();


#endif //READDY_DENSE_KERNEL_H
