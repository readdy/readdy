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
 * @file Compartments.h
 * @brief Header file of CPU program Compartments
 * @author chrisfroe
 * @date 18.10.16
 */

#ifndef READDY_CPUKERNEL_COMPARTMENTS_H
#define READDY_CPUKERNEL_COMPARTMENTS_H

#include <readdy/kernel/singlecpu/actions/SCPUCompartments.h>
#include <readdy/kernel/cpu/CPUKernel.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace actions {

class CPUCompartments : public readdy::model::actions::Compartments {
public:
    using compartmentIdx_t = size_t;
    using particleType_t = unsigned int;

    CPUCompartments(CPUKernel const *const kernel);

    virtual void perform() override;

    virtual void registerCompartment(const std::function<bool(const readdy::model::Vec3)> fun) override;

    virtual void registerConversion(compartmentIdx_t compartmentIdx, std::string from, std::string to) override;

    virtual void registerConversion(compartmentIdx_t compartmentIdx, particleType_t from, particleType_t to);

protected:
    CPUKernel const *const kernel;
    std::vector<std::function<bool(readdy::model::Vec3)>> compartments;
    std::unordered_map<compartmentIdx_t, std::unordered_map<particleType_t, particleType_t>> conversions;
};

}
}
}
}

#endif //READDY_CPUKERNEL_COMPARTMENTS_H
