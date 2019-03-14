/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
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

    void update() override {
        if (_cutoff > 0) {
            fillBins();
        }
    };

    void clear() override {
        _bins.clear();
        _isSetUp = false;
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

    void setUpBins() override { };

    void fillBins();

private:

    count_type getMaxCounts();

    util::Index2D _binsIndex;
    std::vector<std::size_t> _bins;
    block_n_particles_type _blockNParticles;
};

}
}
}
}
