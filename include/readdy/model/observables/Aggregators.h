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
 * An Aggregator is the same as a CombinerObservable but the results that it produces can accumulate and depend on time.
 * This file contains the definitions for aggregators, these are:
 *   - MeanSquaredDisplacement
 *
 * @file Aggregators.h
 * @brief Definition of several aggregators.
 * @author chrisfroe
 * @date 07.11.16
 */

#pragma once

#include <readdy/model/observables/Observable.h>
#include <readdy/model/observables/Observables.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(observables)

class MeanSquaredDisplacement
        : public Combiner<std::pair<std::vector<time_step_type>, std::vector<scalar>>, Particles> {
public:
    MeanSquaredDisplacement(Kernel *kernel, stride_type stride, std::vector<std::string> typesToCount,
                            Particles *particlesObservable);

    MeanSquaredDisplacement(Kernel *kernel, stride_type stride, const std::vector<particle_type_type> &typesToCount,
                            Particles *particlesObservable);

    void evaluate() override = 0;

    std::string type() const override;

protected:
    std::vector<particle_type_type> typesToCount;
};

template<typename ParentObs>
class Trivial
        : public Combiner<std::pair<std::vector<time_step_type>, std::vector<typename ParentObs::result_type>>, ParentObs> {

    static_assert(
            std::is_base_of<readdy::model::observables::Observable<typename ParentObs::result_type>, ParentObs>::value,
            "ParentObs must extend readdy::model::observables::Observable");
public:

    using stride_type = typename ParentObs::stride_type;

    Trivial(Kernel *const kernel, stride_type stride, ParentObs *parentObs)
            : Combiner<std::pair<std::vector<time_step_type>, std::vector<typename ParentObs::result_type>>, ParentObs>(
            kernel, stride,
            parentObs) {}

    virtual void evaluate() {
        const auto &currentInput = std::get<0>(this->parentObservables)->getResult();
        auto &resultTimes = std::get<0>(this->result);
        auto &resultValues = std::get<1>(this->result);
        resultTimes.push_back(this->currentTimeStep());
        resultValues.push_back(currentInput);
    };

    std::string type() const override {
        return "Trivial";
    }
};

NAMESPACE_END(observables)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
