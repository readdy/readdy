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
 * This header file contains the definition of the ObservableFactory. Its purpose is to create observable of different
 * types in the form of unique_ptrs. The actual implementation of an observable can be changed by specializing the
 * dispatcher for its type and invoking a virtual (and then: overridden) method within the factory.
 *
 * @file ObservableFactory.h
 * @brief This header file contains the definition of the ObservableFactory.
 * @author clonker
 * @date 29.04.16
 */
#pragma once

#include <string>
#include <unordered_map>
#include <readdy/model/observables/Observable.h>
#include <readdy/common/Utils.h>
#include <readdy/model/observables/Observables.h>
#include "Aggregators.h"
#include "Positions.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
class Kernel;

NAMESPACE_BEGIN(observables)

class ObservableFactory {
public:
    ObservableFactory(Kernel *const kernel) : kernel(kernel) {};

    template<typename T, typename Obs1, typename Obs2>
    inline std::unique_ptr<T> create(Obs1 *obs1, Obs2 *obs2, unsigned int stride = 1) const {
        return std::make_unique<T>(kernel, obs1, obs2, stride);
    };

    template<typename R, typename... Args>
    inline std::unique_ptr<R> create(unsigned int stride, Args &&... args) const {
        return std::unique_ptr<R>(
                ObservableFactory::get_dispatcher<R, Args...>::impl(this, stride, std::forward<Args>(args)...));
    }

    virtual HistogramAlongAxis *
    createHistogramAlongAxis(unsigned int stride, std::vector<double> binBorders,
                             std::vector<std::string> typesToCount, unsigned int axis) const {
        throw std::runtime_error("should be overridden if a kernel supports this observable");
    }

    virtual NParticles *
    createNParticles(unsigned int stride, std::vector<std::string> typesToCount = {}) const {
        throw std::runtime_error("should be overridden if a kernel supports this observable");
    }

    virtual Forces *
    createForces(unsigned int stride, std::vector<std::string> typesToCount = {}) const {
        throw std::runtime_error("should be overridden if a kernel supports this observable");
    }

    virtual Positions *
    createPositions(unsigned int stride, std::vector<std::string> typesToCount = {}) const {
        throw std::runtime_error("should be overridden if a kernel supports this observable");
    }

    virtual RadialDistribution *
    createRadialDistribution(unsigned int stride, std::vector<double> binBorders, std::vector<std::string> typeCountFrom,
                             std::vector<std::string> typeCountTo,
                             double particleDensity) const {
        throw std::runtime_error("should be overridden if a kernel supports this observable");
    }

    virtual Particles *
    createParticles(unsigned int stride) const {
        throw std::runtime_error("should be overridden if a kernel supports this observable");
    }

    virtual MeanSquaredDisplacement *
    createMeanSquaredDisplacement(unsigned int stride, std::vector<std::string> typesToCount, Particles *particlesObservable) const {
        throw std::runtime_error("should be overridden if a kernel supports this observable");
    }

    virtual Reactions * createReactions(unsigned int stride) const {
        throw std::runtime_error("should be overridden if a kernel supports this observable");
    }

    virtual ReactionCounts* createReactionCounts(unsigned int stride) const {
        throw std::runtime_error("should be overridden if a kernel supports this observable");
    }

protected:
    Kernel *const kernel;

    template<typename T, typename... Args>
    struct get_dispatcher;

    template<typename T, typename... Args>
    struct get_dispatcher {
        static T *impl(const ObservableFactory *self, unsigned int stride, Args &&... args) {
            // this only invokes the normal constructor
            return new T(self->kernel, stride, std::forward<Args>(args)...);
        };
    };

};

READDY_CREATE_OBSERVABLE_FACTORY_DISPATCHER(HistogramAlongAxis)

READDY_CREATE_OBSERVABLE_FACTORY_DISPATCHER(NParticles)

READDY_CREATE_OBSERVABLE_FACTORY_DISPATCHER(Forces)

READDY_CREATE_OBSERVABLE_FACTORY_DISPATCHER(Positions)

READDY_CREATE_OBSERVABLE_FACTORY_DISPATCHER(RadialDistribution)

READDY_CREATE_OBSERVABLE_FACTORY_DISPATCHER(Particles)

READDY_CREATE_OBSERVABLE_FACTORY_DISPATCHER(MeanSquaredDisplacement)

READDY_CREATE_OBSERVABLE_FACTORY_DISPATCHER(Reactions)

READDY_CREATE_OBSERVABLE_FACTORY_DISPATCHER(ReactionCounts)

NAMESPACE_END(observables)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
