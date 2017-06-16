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
 * @file TopologyReaction.h
 * @brief << brief description >>
 * @author clonker
 * @date 03.04.17
 * @copyright GNU Lesser General Public License v3.0
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
    static constexpr std::size_t raise_or_rollback_flag = 0;
    static constexpr std::size_t expect_connected_or_create_children_flag = 1;
    std::bitset<2> flags;

    void raise();

    void rollback();

    void expect_connected();

    void create_children();
};

class ReactionFunctionGenerator {
public:
    using reaction_recipe = Recipe;
    using reaction_function = std::function<reaction_recipe(void)>;

    virtual reaction_function generate(GraphTopology& topology) = 0;
};

class STDFunctionReactionFunctionGenerator : public ReactionFunctionGenerator {
public:
    using std_reaction_function = std::function<reaction_recipe(GraphTopology&)>;

    STDFunctionReactionFunctionGenerator(const std_reaction_function &fun);

    virtual reaction_function generate(GraphTopology &topology) override;

private:
    std_reaction_function fun;
};

class RateFunctionGenerator {
public:
    using rate_function = std::function<readdy::scalar(void)>;

    virtual rate_function generate(const GraphTopology& topology) = 0;
};

class STDFunctionRateFunctionGenerator : public RateFunctionGenerator {
public:
    using std_rate_function = std::function<scalar(const GraphTopology&)>;

    STDFunctionRateFunctionGenerator(const std_rate_function &fun);

    virtual rate_function generate(const GraphTopology &topology) override;

private:
    std_rate_function fun;
};

class TopologyReaction {
public:
    using mode = Mode;
    using reaction_recipe = ReactionFunctionGenerator::reaction_recipe;
    using reaction_function = STDFunctionReactionFunctionGenerator::std_reaction_function;
    using rate_function = STDFunctionRateFunctionGenerator::std_rate_function;

    TopologyReaction(const reaction_function &reaction_function, const rate_function &rate_function);

    TopologyReaction(const reaction_function &reaction_function, const double &rate);

    TopologyReaction(std::shared_ptr<ReactionFunctionGenerator> reaction_function_generator,
                     std::shared_ptr<RateFunctionGenerator> rate_function_generator);

    TopologyReaction(const TopologyReaction&) = default;

    TopologyReaction& operator=(const TopologyReaction&) = default;

    TopologyReaction(TopologyReaction&&) = default;

    TopologyReaction& operator=(TopologyReaction&&) = default;

    ~TopologyReaction() = default;

    double rate(const GraphTopology &topology) const;

    reaction_recipe operations(GraphTopology &topology) const;

    const bool raises_if_invalid() const;

    void raise_if_invalid();

    const bool rolls_back_if_invalid() const;

    void roll_back_if_invalid();

    const bool expects_connected_after_reaction() const;

    void expect_connected_after_reaction();

    const bool creates_child_topologies_after_reaction() const;

    void create_child_topologies_after_reaction();

    std::vector<GraphTopology> execute(GraphTopology &topology, const Kernel *const kernel);

private:
    std::shared_ptr<RateFunctionGenerator> rate_function_generator_;
    std::shared_ptr<ReactionFunctionGenerator> reaction_function_generator_;
    mode mode_;
};

NAMESPACE_END(reactions)
NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
