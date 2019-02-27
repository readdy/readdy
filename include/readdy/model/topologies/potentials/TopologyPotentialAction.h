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
 * @file TopologyAction.h
 * @brief << brief description >>
 * @author clonker
 * @date 30.01.17
 * @copyright BSD-3
 */

#pragma once
#include <readdy/common/macros.h>
#include <readdy/model/Context.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)
class Topology;
NAMESPACE_BEGIN(pot)

class TopologyPotentialAction {
public:
    explicit TopologyPotentialAction(const Context *const context) : context(context) {}
    TopologyPotentialAction(const TopologyPotentialAction&) = default;
    TopologyPotentialAction& operator=(const TopologyPotentialAction&) = delete;
    TopologyPotentialAction(TopologyPotentialAction&&) = default;
    TopologyPotentialAction& operator=(TopologyPotentialAction&&) = delete;
    virtual ~TopologyPotentialAction() = default;

protected:
    const Context* const context;
};

class EvaluatePotentialAction : public TopologyPotentialAction {
public:
    explicit EvaluatePotentialAction(const Context *const context) : TopologyPotentialAction(context) {}
    EvaluatePotentialAction(const EvaluatePotentialAction&) = default;
    EvaluatePotentialAction& operator=(const EvaluatePotentialAction&) = delete;
    EvaluatePotentialAction(EvaluatePotentialAction&&) = default;
    EvaluatePotentialAction& operator=(EvaluatePotentialAction&&) = delete;

    ~EvaluatePotentialAction() override = default;

    virtual scalar perform(const Topology*topology) = 0;
};

NAMESPACE_END(pot)
NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
