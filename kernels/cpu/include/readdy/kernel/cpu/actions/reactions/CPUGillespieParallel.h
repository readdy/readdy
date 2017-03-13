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
 * @file GillespieParallel.h
 * @brief << brief description >>
 * @author clonker
 * @date 20.10.16
 */

#pragma once

#include <readdy/model/actions/Actions.h>
#include <readdy/kernel/cpu/CPUKernel.h>
#include <readdy/kernel/singlecpu/actions/SCPUReactionImpls.h>
#include "ReactionUtils.h"

namespace readdy {
namespace kernel {
namespace cpu {
namespace actions {
namespace reactions {

class CPUGillespieParallel : public readdy::model::actions::reactions::GillespieParallel {
    using super = readdy::model::actions::reactions::GillespieParallel;
public:
    CPUGillespieParallel(kernel_t *const kernel, double timeStep);

    ~CPUGillespieParallel();

    virtual void perform() override;

    void clear();

    double getMaxReactionRadius() const;

    double getBoxWidth() const;

    unsigned int getLongestAxis() const;

    unsigned int getOtherAxis1() const;

    unsigned int getOtherAxis2() const;

    void setApproximateRate(bool approximateRate);

protected:
    kernel_t *const kernel;
    double maxReactionRadius = 0.0;
    double boxWidth = 0.0;
    unsigned int longestAxis;
    unsigned int otherAxis1, otherAxis2;

    struct SlicedBox {
        using particle_indices_t = std::vector<data_t::index_t>;
        particle_indices_t particleIndices{};
        unsigned int id = 0;
        vec_t lowerLeftVertex, upperRightVertex;
        double leftBoundary = 0;
        double rightBoundary = 0;
        particle_indices_t::size_type n_shells;
        unsigned int longestAxis;
        double boxWidth;
        double shellWidth = 0.0;

        long getShellIndex(const vec_t &pos) const;

        SlicedBox(unsigned int id, vec_t lowerLeftVertex, vec_t upperRightVertex, double maxReactionRadius,
                  unsigned int longestAxis);

        friend bool operator==(const SlicedBox &lhs, const SlicedBox &rhs) {
            return lhs.id == rhs.id;
        }

        friend bool operator!=(const SlicedBox &lhs, const SlicedBox &rhs) {
            return !(lhs == rhs);
        }

        bool isInBox(const vec_t &particle) const;
    };

    std::vector<SlicedBox> boxes;
    bool approximateRate = true;

    /**
     * look for the longest axis and divide it into n_threads parts, yielding halo boxes.
     */
    void setupBoxes();

    /**
     * Sorts particles into boxes and adds them to the problematic particles should they fall into
     * a halo region
     */
    void fillBoxes();

    /**
     * Executes the gillespie algorithm for each box and gives an update on problematic particles
     */
    virtual void handleBoxReactions();

    void findProblematicParticles(data_t::index_t entry, const SlicedBox &box, ctx_t ctx,
                                  const data_t& data, nl_t* nl, std::set<data_t::index_t> &problematic,
                                  const readdy::model::KernelContext::dist_squared_fun&) const;

};
}

}
}
}
}
