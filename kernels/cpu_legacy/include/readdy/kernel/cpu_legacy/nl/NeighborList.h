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
#include <readdy/common/Timer.h>
#include <readdy/kernel/cpu_legacy/data/NLDataContainer.h>
#include "NeighborListIterator.h"

namespace readdy {
namespace kernel {
namespace cpu_legacy {
namespace nl {

class NeighborList {
public:
    using DataUpdate = std::tuple<std::vector<data::Entry>, std::vector<std::size_t>>;
    using neighbors_type = std::vector<std::size_t>;
    using const_iterator = NeighborListIterator;

    NeighborList(const readdy::model::Context &context,
                 const readdy::util::thread::Config &config) : _context(context), _config(config) {};

    virtual ~NeighborList() = default;

    const_iterator begin() const {
        return cbegin();
    };

    virtual const_iterator cbegin() const = 0;

    const_iterator end() const {
        return cend();
    }

    virtual const_iterator cend() const = 0;

    virtual std::size_t size() const = 0;

    virtual void set_up(const util::PerformanceNode &node) = 0;

    virtual void update(const util::PerformanceNode &node) = 0;

    virtual void clear(const util::PerformanceNode &node) = 0;

    virtual void updateData(DataUpdate &&update) = 0;

    virtual bool is_adaptive() const = 0;

    virtual const data::EntryDataContainer * data() const = 0;

    virtual data::EntryDataContainer *data() = 0;

    scalar &skin() {
        return _skin;
    };

    const scalar &skin() const {
        return _skin;
    };

protected:
    scalar _skin {0};
    scalar _max_cutoff {0};
    scalar _max_cutoff_skin_squared {0};

    std::reference_wrapper<const readdy::model::Context> _context;
    std::reference_wrapper<const readdy::util::thread::Config> _config;
};

}
}
}
}