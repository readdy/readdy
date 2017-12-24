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

    ContiguousCellLinkedList(data_type &data, const readdy::model::Context &context,
                             const util::thread::Config &config);

    void update(const util::PerformanceNode &node) override;

    void clear() override;

    const std::vector<std::size_t> &bins() const;

    const util::Index2D &binsIndex() const;

    std::size_t *particlesBegin(std::size_t cellIndex) override;

    const std::size_t *particlesBegin(std::size_t cellIndex) const override;

    std::size_t *particlesEnd(std::size_t cellIndex) override;

    const std::size_t *particlesEnd(std::size_t cellIndex) const override;

    size_t nParticles(std::size_t cellIndex) const override;

protected:
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
