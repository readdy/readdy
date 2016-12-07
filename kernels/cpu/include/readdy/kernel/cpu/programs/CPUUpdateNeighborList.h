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
 * @file UpdateNeighborList.h
 * @brief << brief description >>
 * @author clonker
 * @date 13.07.16
 */


#ifndef READDY_CPUKERNEL_UPDATENEIGHBORLIST_H
#define READDY_CPUKERNEL_UPDATENEIGHBORLIST_H

#include <readdy/model/programs/Programs.h>
#include <readdy/kernel/cpu/CPUKernel.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace programs {
class CPUUpdateNeighborList : public readdy::model::programs::UpdateNeighborList {
public:

    CPUUpdateNeighborList(CPUKernel *kernel) : kernel(kernel) {
    }

    virtual void execute() override {
        switch (action) {
            case create:
                kernel->getKernelStateModel().updateNeighborList();
                break;
            case clear:
                kernel->getKernelStateModel().clearNeighborList();
                break;
        }

    }

    virtual void setSkinSize(double skinSize) override {
        kernel->getKernelStateModel().getNeighborList()->setSkinSize(skinSize);
    }

    bool supportsSkin() const override {
        return true;
    }

private:
    CPUKernel *kernel;
};
}
}
}
}
#endif //READDY_CPUKERNEL_UPDATENEIGHBORLIST_H
