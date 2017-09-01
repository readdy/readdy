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


#pragma once
#include <readdy/model/actions/Actions.h>
#include <readdy/kernel/cpu/CPUKernel.h>
#include <readdy/kernel/cpu/nl/AdaptiveNeighborList.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace actions {
class CPUUpdateNeighborList : public readdy::model::actions::UpdateNeighborList {
    using super = readdy::model::actions::UpdateNeighborList;
    bool firstRun = true;
public:

    CPUUpdateNeighborList(CPUKernel *kernel, super::Operation op, readdy::scalar skin) : super(op, skin), kernel(kernel) {}

    void perform() override {
        if(firstRun) {
            if(skinSize >= 0) kernel->getCPUKernelStateModel().getNeighborList()->skin() = skinSize;
            firstRun = false;
        }
        switch (operation) {
            case create:
                kernel->getKernelStateModel().updateNeighborList();
                break;
            case clear:
                kernel->getKernelStateModel().clearNeighborList();
                break;
        }

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
