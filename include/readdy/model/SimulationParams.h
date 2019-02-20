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
 * @file SimulationParams.h
 * @brief Holds parameters that determine how the simulation is run.
 * @author chrisfroe
 * @date 20.02.2019
 */

#pragma once

#include <readdy/common/common.h>

namespace readdy::model {

class SimulationParams {
public:

    SimulationParams() = default;

    // no copies
    SimulationParams(const SimulationParams &) = delete;

    SimulationParams &operator=(const SimulationParams &) = delete;

    // move allowed
    SimulationParams(SimulationParams &&) = default;

    SimulationParams &operator=(SimulationParams &&) = default;

    virtual ~SimulationParams() = default;

    /**
     * Returns whether reactions with positions shall be recorded in the state model, then obtainable by
     * the readdy::model::observables::Reactions observable.
     * @return whether reactions shall be recorded in the state model, by default false
     */
    const bool &recordReactionsWithPositions() const {
        return _recordReactionsWithPositions;
    }

    /**
     * Returns whether reactions with positions shall be recorded in the state model, then obtainable by
     * the readdy::model::observables::Reactions observable.
     * @return whether reactions shall be recorded in the state model, by default false
     */
    bool &recordReactionsWithPositions() {
        return _recordReactionsWithPositions;
    }

    /**
     * Returns whether reaction counts shall be recorded in the state model (if it is supported). It is then obtainable
     * by the readdy::model::observables::ReactionCounts observable.
     * @return whether reaction counts shall be recorded
     */
    const bool &recordReactionCounts() const {
        return _recordReactionCounts;
    }

    /**
     * Returns whether reaction counts shall be recorded in the state model (if it is supported). It is then obtainable
     * by the readdy::model::observables::ReactionCounts observable.
     * @return whether reaction counts shall be recorded
     */
    bool &recordReactionCounts() {
        return _recordReactionCounts;
    }

    /**
     * Returns whether virial shall be recorded in the state model (if it is supported). It is then obtainable
     * by the readdy::model::observables::Virial observable.
     * @return whether virial shall be recorded
     */
    const bool &recordVirial() const {
        return _recordVirial;
    }

    /**
     * Returns whether virial shall be recorded in the state model (if it is supported). It is then obtainable
     * by the readdy::model::observables::Virial observable.
     * @return whether virial shall be recorded
     */
    bool &recordVirial() {
        return _recordVirial;
    }

protected:
    bool _recordReactionsWithPositions{false};
    bool _recordReactionCounts{false};
    bool _recordVirial{false};

    scalar neighborListCellWidth{-1.};
};

}
