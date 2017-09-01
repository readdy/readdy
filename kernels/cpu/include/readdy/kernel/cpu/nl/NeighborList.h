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
 * @file NeighborList.h
 * @brief << brief description >>
 * @author clonker
 * @date 01.09.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <readdy/common/common.h>
#include <readdy/kernel/cpu/model/CPUParticleData.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace nl {

class NeighborList {
public:
    using data_t = readdy::kernel::cpu::model::CPUParticleData;
    using neighbors_t = data_t::neighbors_t;
    using iterator = data_t::neighbors_list_iterator;
    using const_iterator = data_t::neighbors_list_const_iterator;

    NeighborList(model::CPUParticleData &data, const readdy::model::KernelContext &context,
                 const readdy::util::thread::Config &config) : _data(data), _context(context), _config(config) {};

    virtual ~NeighborList() = default;

    virtual void set_up() = 0;

    virtual void update() = 0;

    virtual void clear() = 0;

    virtual void updateData(data_t::update_t &&update) = 0;

    const neighbors_t &neighbors_of(const data_t::index_t entry) const {
        const static neighbors_t no_neighbors{};
        if (_max_cutoff > 0) {
            return _data.get().neighbors_at(entry);
        }
        return no_neighbors;
    };

    scalar &skin() {
        return _skin;
    };

    const scalar &skin() const {
        return _skin;
    };

    iterator begin() {
        return _data.get().neighbors.begin();
    };

    iterator end() {
        return _data.get().neighbors.end();
    };

    const_iterator cbegin() const {
        return _data.get().neighbors.cbegin();
    };

    const_iterator cend() const {
        return _data.get().neighbors.cend();
    };

protected:
    scalar _skin {0};
    scalar _max_cutoff {0};
    scalar _max_cutoff_skin_squared {0};

    std::reference_wrapper<model::CPUParticleData> _data;
    std::reference_wrapper<const readdy::model::KernelContext> _context;
    std::reference_wrapper<const readdy::util::thread::Config> _config;
};

}
}
}
}