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
 * @file NeighborListContainer.h
 * @brief << brief description >>
 * @author clonker
 * @date 15.09.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <vector>
#include <readdy/common/thread/Config.h>
#include <readdy/kernel/cpu_legacy/data/DefaultDataContainer.h>
#include "CellLinkedList.h"

namespace readdy {
namespace kernel {
namespace cpu {
namespace nl {

struct NLContainerConfig {
    using data_container_type = readdy::kernel::cpu::data::DefaultDataContainer;
    using thread_config_type = readdy::util::thread::Config;

    NLContainerConfig(const model::Context &context, const thread_config_type &threads,
                      const data_container_type &data);

    std::reference_wrapper<const model::Context> context;
    std::reference_wrapper<const data_container_type> data;
    std::reference_wrapper<const thread_config_type> threads;
};

class NeighborListContainer {

public:
    using value_type = std::size_t;
    using local_index_vector = std::vector<value_type>;
    using index_vector = std::vector<local_index_vector>;
    using indexer = std::vector<std::size_t>;

    using const_iterator = index_vector::const_iterator;

    explicit NeighborListContainer(NLContainerConfig config);

    index_vector &elements() {
        return _elements;
    }

    const index_vector &elements() const {
        return _elements;
    }

    const_iterator begin() const;

    const_iterator end() const;

    void clear();

    virtual void update(scalar cutoffSquared, const util::PerformanceNode &perf) = 0;

    virtual std::size_t size() const = 0;

protected:
    index_vector _elements;
    // todo: offsets array filled if "indexed==true"
    indexer _indexer;

    NLContainerConfig _config;
};

class ContiguousCLLNeighborListContainer : public NeighborListContainer {
public:
    ContiguousCLLNeighborListContainer(NLContainerConfig config, const ContiguousCellLinkedList &cll);

    void update(scalar cutoffSquared, const util::PerformanceNode &perf) override;

    std::size_t size() const override;

private:
    const ContiguousCellLinkedList &_cll;

    std::size_t _size {0};
};

class DynamicCLLNeighborListContainer : public NeighborListContainer {
public:
    DynamicCLLNeighborListContainer(NLContainerConfig config, const DynamicCellLinkedList &cll);

    void update(scalar cutoffSquared, const util::PerformanceNode &perf) override;

    std::size_t size() const override;

private:
    const DynamicCellLinkedList &_cll;

    std::size_t _size {0};
};

class CompactCLLNeighborListContainer : public NeighborListContainer {
public:
    CompactCLLNeighborListContainer(NLContainerConfig config, const CompactCellLinkedList &cll);

    void update(scalar cutoffSquared, const util::PerformanceNode &perf) override;

    void updateSerial(scalar cutoffSquared);

    void updateParallel(scalar cutoffSquared, const util::PerformanceNode &perf);

    std::size_t size() const override;

private:
    std::size_t _size {0};
    bool serialUpdate {false};
    const CompactCellLinkedList &_cll;
};

}
}
}
}