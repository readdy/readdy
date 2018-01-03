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
 * @file ContiguousCellLinkedList.h
 * @brief << brief description >>
 * @author clonker
 * @date 12/24/17
 */


#pragma once

#include "CellLinkedList.h"

namespace readdy {
namespace kernel {
namespace cpu {
namespace nl {

class ContiguousCellLinkedList : public CellLinkedList {
public:

    using count_type = std::size_t;
    using block_n_particles_type = std::vector<util::thread::copyable_atomic<std::size_t>>;

    ContiguousCellLinkedList(data_type &data, const readdy::model::Context &context, thread_pool &pool);

public:

    void update(const util::PerformanceNode &node) override {
        if (_max_cutoff > 0) {
            fillBins(node);
        }
    };

    void clear() override {
        _bins.clear();
    };

    const std::vector<std::size_t> &bins() const {
        return _bins;
    };

    const util::Index2D &binsIndex() const {
        return _binsIndex;
    };

    std::size_t *particlesBegin(std::size_t cellIndex) {
        return !_bins.empty() ? &_bins.at(_binsIndex(cellIndex, 0_z)) : nullptr;
    }

    const std::size_t *particlesBegin(std::size_t cellIndex) const {
        return !_bins.empty() ? &_bins.at(_binsIndex(cellIndex, 0_z)) : nullptr;
    };

    std::size_t *particlesEnd(std::size_t cellIndex) {
        return !_bins.empty() ? particlesBegin(cellIndex) + nParticles(cellIndex) : nullptr;
    };

    const std::size_t *particlesEnd(std::size_t cellIndex) const {
        return !_bins.empty() ? particlesBegin(cellIndex) + nParticles(cellIndex) : nullptr;
    };

    size_t nParticles(std::size_t cellIndex) const {
        return (*_blockNParticles.at(cellIndex)).load();
    };

    template<typename Function>
    void forEachNeighbor(std::size_t particle, const Function &function) const {
        forEachNeighbor(particle, cellOfParticle(particle), function);
    }

    template<typename Function>
    void forEachNeighbor(std::size_t particle, std::size_t cell, const Function& function) const {
        std::for_each(particlesBegin(cell), particlesEnd(cell), [&function, particle](auto x) {
            if(x != particle) function(x);
        });
        for(auto itNeighCell = neighborsBegin(cell); itNeighCell != neighborsEnd(cell); ++itNeighCell) {
            std::for_each(particlesBegin(*itNeighCell), particlesEnd(*itNeighCell), function);
        }
    }

protected:

    void setUpBins(const util::PerformanceNode &node) override { };

    void fillBins(const util::PerformanceNode &node);

private:

    count_type getMaxCounts(const util::PerformanceNode &node);

    util::Index2D _binsIndex;
    std::vector<std::size_t> _bins;
    block_n_particles_type _blockNParticles;
};

}
}
}
}
