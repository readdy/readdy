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
 * This header contains definitions for the structural topology reaction mode, i.e., its flags, as well as the
 * structural topology reaction itself.
 *
 * @file TopologyReaction.h
 * @brief Definition of everything belonging to the structural topology reaction.
 * @author clonker
 * @date 03.04.17
 * @copyright GPL-3
 */

#pragma once

#include <vector>
#include <memory>
#include <functional>
#include <bitset>

#include <readdy/common/macros.h>

#include "TopologyReactionAction.h"
#include "Operations.h"
#include "Recipe.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
class Kernel;
NAMESPACE_BEGIN(top)
class GraphTopology;
NAMESPACE_BEGIN(reactions)

/**
 * Struct holding the mode of the reaction:
 * - whether it should raise or roll back after a failed reaction
 * - whether it is expected to be still connected or if it should be "fissionated"
 */
struct Mode {
    /**
     * flag for raise or rollback
     */
    static constexpr std::size_t raise_or_rollback_flag = 0;
    /**
     * flag for expect connected or create children
     */
    static constexpr std::size_t expect_connected_or_create_children_flag = 1;
    /**
     * bitset over the flags
     */
    std::bitset<2> flags;

    /**
     * activate the 'raise' mode, automatically disables 'rollback'
     */
    void raise() {
        flags[raise_or_rollback_flag] = true;
    }

    /**
     * activate the 'rollback' mode, automatically disables 'raise'
     */
    void rollback() {
        flags[raise_or_rollback_flag] = false;
    }

    /**
     * activate the 'expect_connected' mode, automatically disables 'create_children'
     */
    void expect_connected() {
        flags[expect_connected_or_create_children_flag] = true;
    }

    /**
     * activate the 'create_children' mode, automatically disables 'expect_connected'
     */
    void create_children() {
        flags[expect_connected_or_create_children_flag] = false;
    }
};

class StructuralTopologyReaction {
public:
    /**
     * type of class holding the possible modes
     */
    using mode = Mode;
    /**
     * reaction recipe type
     */
    using reaction_recipe = Recipe;
    /**
     * reaction function type, yielding a recipe
     */
    using reaction_function = std::function<reaction_recipe(GraphTopology&)>;
    /**
     * rate function type, yielding a rate
     */
    using rate_function = std::function<scalar(const GraphTopology&)>;

    /**
     * creates a new instance by supplying a reaction function and a rate function
     * @param reaction_function the reaction function
     * @param rate_function the rate function
     */
    StructuralTopologyReaction(const reaction_function &reaction_function, const rate_function &rate_function);

    /**
     * creates a new instance by supplying a reaction function and a constant rate
     * @param reaction_function the reaction function
     * @param rate the rate
     */
    StructuralTopologyReaction(const reaction_function &reaction_function, const scalar &rate);

    /**
     * default copy
     */
    StructuralTopologyReaction(const StructuralTopologyReaction&) = default;

    /**
     * default copy assign
     */
    StructuralTopologyReaction& operator=(const StructuralTopologyReaction&) = default;

    /**
     * default move
     */
    StructuralTopologyReaction(StructuralTopologyReaction&&) = default;

    /**
     * default move assign
     */
    StructuralTopologyReaction& operator=(StructuralTopologyReaction&&) = default;

    /**
     * default destructor
     */
    ~StructuralTopologyReaction() = default;

    /**
     * Evaluates the rate of this reaction for a given topology.
     * @param topology the topology
     * @return the rate
     */
    scalar rate(const GraphTopology &topology) const {
        return _rate_function(topology);
    }

    /**
     * Yields a reaction recipe for a given topology.
     * @param topology the topology
     * @return the recipe
     */
    reaction_recipe operations(GraphTopology &topology) const {
        return _reaction_function(topology);
    }

    /**
     * checks if the reaction handler will raise if the post-reaction topology is invalid
     * @return true if the reaction handler raises upon invalid topology
     */
    const bool raises_if_invalid() const {
        return mode_.flags.test(mode::raise_or_rollback_flag);
    }

    /**
     * instruct the reaction handler to raise if this reaction yields invalid topologies, counter part to
     * `roll_back_if_invalid`
     */
    void raise_if_invalid() {
        mode_.raise();
    }

    /**
     * checks if the reaction handler triggers a roll-back if the post-reaction topology is invalid
     * @return true in case a roll back should be performed
     */
    const bool rolls_back_if_invalid() const {
        return !raises_if_invalid();
    }

    /**
     * instruct the reaction handler to roll back if this reaction yields invalid topologies, counter part to
     * `raise_if_invalid`
     */
    void roll_back_if_invalid() {
        mode_.rollback();
    }

    /**
     * checks if this reaction is expected to yield connected topology graphs only
     * @return true if only connected topology graphs should be yielded
     */
    const bool expects_connected_after_reaction() const {
        return mode_.flags.test(mode::expect_connected_or_create_children_flag);
    }

    /**
     * instruct the reaction handler to check whether to topology is still connected after the reaction, counter part
     * to create_child_topologies_after_reaction
     */
    void expect_connected_after_reaction() {
        mode_.expect_connected();
    }

    /**
     * checks if child topologies should be created if the topology reaction yielded more than one connected component
     * @return true if child topologies should be created
     */
    const bool creates_child_topologies_after_reaction() const {
        return !expects_connected_after_reaction();
    }

    /**
     * instruct the reaction handler to create child topologies if the reaction yields more than one connected
     * component in the graph - counter part to expect_connected_after_reaction
     */
    void create_child_topologies_after_reaction() {
        mode_.create_children();
    }

    /**
     * Executes the topology reaction on a topology and a kernel, possibly returns child topologies.
     * @param topology the topology
     * @param kernel the kernel
     * @return a vector of child topologies if they were created in the process
     */
    std::vector<GraphTopology> execute(GraphTopology &topology, const Kernel * kernel) const;

private:
    /**
     * the reaction function responsible of generating the reaction recipe out of a topology
     */
    reaction_function _reaction_function;
    /**
     * the rate function responsible of calculating a rate for a given topology
     */
    rate_function _rate_function;
    /**
     * the execution mode
     */
    mode mode_;
};

NAMESPACE_END(reactions)
NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
