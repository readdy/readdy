/********************************************************************
 * Copyright © 2017 Computational Molecular Biology Group,          * 
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
 * @file CellLinkedList.h
 * @brief << brief description >>
 * @author clonker
 * @date 12.09.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <cstddef>
#include <readdy/common/Index.h>
#include <readdy/common/thread/Config.h>
#include <readdy/model/KernelContext.h>
#include <readdy/kernel/cpu/model/CPUParticleData.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace nl {

class CellLinkedList {
    CellLinkedList(model::CPUParticleData &data, const readdy::model::KernelContext &context,
                   const readdy::util::thread::Config &config);
private:

    scalar _skin {0};
    scalar _max_cutoff {0};
    scalar _max_cutoff_skin_squared {0};

    std::size_t _maxParticlesPerCell;
    util::Index3D _cellIndex;
    util::Index2D _cellNeighbors;
    util::Index2D _bins;

    std::reference_wrapper<model::CPUParticleData> _data;
    std::reference_wrapper<const readdy::model::KernelContext> _context;
    std::reference_wrapper<const readdy::util::thread::Config> _config;
};

}
}
}
}
