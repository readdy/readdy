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
 * @file NoNeighborList.h
 * @brief << brief description >>
 * @author clonker
 * @date 15.09.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include "NeighborList.h"
#include "../data/DefaultDataContainer.h"

namespace readdy {
namespace kernel {
namespace cpu {
namespace nl {

class NoNeighborList : public NeighborList {
public:
    NoNeighborList(const readdy::model::KernelContext &context, const readdy::util::thread::Config &config)
            : NeighborList(context, config), _data(context, config) {}

    const_iterator cbegin() const override {
        return NeighborListIterator{_neighbors.begin(), _neighbors.end(), true};
    }

    const_iterator cend() const override {
        return NeighborListIterator{_neighbors.end(), _neighbors.end(), true};
    }

    void set_up(const util::PerformanceNode &node) override {
        _neighbors.resize(_data.size());
    }

    void update(const util::PerformanceNode &node) override {
        _neighbors.resize(_data.size());
    }

    void clear(const util::PerformanceNode &node) override {}

    void updateData(DataUpdate &&update) override {
        _data.update(std::forward<DataUpdate>(update));
    }

    bool is_adaptive() const override {
        return false;
    }

    const data::EntryDataContainer *data() const override {
        return &_data;
    }

    std::size_t size() const override {
        return 0;
    }

    data::EntryDataContainer *data() override {
        return &_data;
    }

    const neighbors_type &neighbors_of(std::size_t entry) const override {
        return {};
    }

private:
    data::DefaultDataContainer _data;
    std::vector<std::vector<std::size_t>> _neighbors {};
};

}
}
}
}