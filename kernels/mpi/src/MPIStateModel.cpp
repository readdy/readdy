/********************************************************************
 * Copyright © 2019 Computational Molecular Biology Group,          *
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
 * « detailed description »
 *
 * @file MPIStateModel.cpp
 * @brief « brief description »
 * @author chrisfroe
 * @date 28.05.19
 */

#include <readdy/kernel/mpi/MPIStateModel.h>

namespace readdy::kernel::mpi {

MPIStateModel::MPIStateModel(Data &data, const readdy::model::Context &context) : _data(data), _context(context) {}

const std::vector<readdy::Vec3> MPIStateModel::getParticlePositions() const {
    throw std::runtime_error("impl");
}

const std::vector<MPIStateModel::Particle>
MPIStateModel::getParticles() const {
    throw std::runtime_error("impl");
}

void MPIStateModel::resetReactionCounts() {
    if (!reactionCounts().empty()) {
        for (auto &e : reactionCounts()) {
            e.second = 0;
        }
    } else {
        const auto &reactions = _context.get().reactions();
        for (const auto &entry : reactions.order1()) {
            for (auto reaction : entry.second) {
                reactionCounts()[reaction->id()] = 0;
            }
        }
        for (const auto &entry : reactions.order2()) {
            for (auto reaction : entry.second) {
                reactionCounts()[reaction->id()] = 0;
            }
        }
    }
}

void MPIStateModel::toDenseParticleIndices(
        std::vector<std::size_t>::iterator begin,
        std::vector<std::size_t>::iterator end) const {
    const auto &blanks = _data.get().blanks();
    std::transform(begin, end, begin, [&blanks](const std::size_t &ix) {
        auto result = ix;
        for (auto blankIx : blanks) {
            if (blankIx < ix) {
                --result;
            }
        }
        return result;
    });
}

void MPIStateModel::clear() {
    getParticleData()->clear();
    //topologies().clear();
    reactionRecords().clear();
    resetReactionCounts();
    virial() = {};
    energy() = 0;
}

}
