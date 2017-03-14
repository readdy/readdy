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
 * The PotentialFactory contains create-methods, that can be overridden by other kernels, in case
 * another kernel-specific implementation of the potential is needed. The dispatcher is responsible
 * for calling the correct create-method.
 *
 * @file PotentialFactory.h
 * @brief Header file containing the definition of the PotentialFactory.
 * @author clonker
 * @date 31.05.16
 */

#pragma once

#include <unordered_map>
#include <functional>

#include <readdy/model/potentials/PotentialsOrder1.h>
#include <readdy/model/potentials/PotentialsOrder2.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(potentials)

class PotentialFactory {
public:

    template<typename R, typename... Args>
    std::unique_ptr<R> createPotential(Args &&... args) const {
        return std::unique_ptr<R>(get_dispatcher<R, Args...>::impl(this, std::forward<Args>(args)...));
    }

    std::vector<std::string> getAvailablePotentials() const {
        return {getPotentialName<CubePotential>(), getPotentialName<SpherePotential>(),
                getPotentialName<HarmonicRepulsion>(), getPotentialName<WeakInteractionPiecewiseHarmonic>()};
    }

    virtual CubePotential *createCubePotential(const std::string &particleType, double forceConstant,
                                               const Vec3 &origin, const Vec3 &extent,
                                               bool considerParticleRadius) const {
        return new CubePotential(particleType, forceConstant, origin, extent, considerParticleRadius);
    };

    CubePotential *createCubePotential(const std::string &particleType, double forceConstant,
                                       const Vec3 &origin, const Vec3 &extent) const {
        return createCubePotential(particleType, forceConstant, origin, extent, true);
    }

    virtual SpherePotential *createSpherePotential(const std::string &particleType, double forceConstant,
                                                   const Vec3 &origin, double radius) const {
        return new SpherePotential(particleType, forceConstant, origin, radius);
    };

    virtual HarmonicRepulsion *
    createHarmonicRepulsion(const std::string &type1, const std::string &type2, double forceConstant) const {
        return new HarmonicRepulsion(type1, type2, forceConstant);
    };

    virtual WeakInteractionPiecewiseHarmonic *
    createWeakInteractionPiecewiseHarmonic(const std::string &type1, const std::string &type2,
                                           const double forceConstant,
                                           const WeakInteractionPiecewiseHarmonic::Configuration &config) const {
        return new WeakInteractionPiecewiseHarmonic(type1, type2, forceConstant, config);
    };

    WeakInteractionPiecewiseHarmonic *
    createWeakInteractionPiecewiseHarmonic(const std::string &type1, const std::string &type2,
                                           const double forceConstant, const double desiredDist,
                                           const double depth,
                                           const double cutoff) const {
        using config = WeakInteractionPiecewiseHarmonic::Configuration;
        return createWeakInteractionPiecewiseHarmonic(type1, type2, forceConstant, config{desiredDist, depth, cutoff});
    };

    LennardJones* createLennardJones(const std::string& type1, const std::string& type2, unsigned int m, unsigned int n, double cutoff, bool shift, double epsilon, double sigma) const {
        return new LennardJones(type1, type2, m, n, cutoff, shift, epsilon, sigma);
    }

protected:
    template<typename T, typename... Args>
    struct get_dispatcher;

    template<typename T, typename... Args>
    struct get_dispatcher {
        static T *impl(const PotentialFactory *self, Args &&... args) {
            return new T(std::forward<Args>(args)...);
        };
    };
};

READDY_CREATE_FACTORY_DISPATCHER(PotentialFactory, CubePotential)

READDY_CREATE_FACTORY_DISPATCHER(PotentialFactory, SpherePotential)

READDY_CREATE_FACTORY_DISPATCHER(PotentialFactory, HarmonicRepulsion)

READDY_CREATE_FACTORY_DISPATCHER(PotentialFactory, WeakInteractionPiecewiseHarmonic)

READDY_CREATE_FACTORY_DISPATCHER(PotentialFactory, LennardJones)

NAMESPACE_END(potentials)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
