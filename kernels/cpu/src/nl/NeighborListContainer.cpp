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
 * @file NeighborListContainer.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 15.09.17
 * @copyright GNU Lesser General Public License v3.0
 */


#include <readdy/kernel/cpu/nl/NeighborListContainer.h>
#include <readdy/kernel/cpu/nl/NeighborList.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace nl {


NeighborListContainer::NeighborListContainer(const data_container_type &data, const thread_config_type &threadConfig)
        : _data(data), _config(threadConfig){}

NeighborListContainer::const_iterator NeighborListContainer::begin() const {
    return {_elements.begin(), _elements.end()};
}

NeighborListContainer::const_iterator NeighborListContainer::end() const {
    return const_iterator(_elements.end());
}

}
}
}
}