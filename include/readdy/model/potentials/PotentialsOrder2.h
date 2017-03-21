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
 * This header contains the declarations of order 2 potentials. Currently:
 *   - Harmonic repulsion
 *   - Weak interaction piecewise harmonic
 *
 * @file PotentialsOrder2.h
 * @brief Contains the declaration of order 2 potentials.
 * @author clonker
 * @date 09.06.16
 */

#pragma once

#include <readdy/model/KernelContext.h>
#include <ostream>
#include "PotentialOrder2.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(potentials)

class HarmonicRepulsion : public PotentialOrder2 {
    using super = PotentialOrder2;
public:
    HarmonicRepulsion(const std::string &type1, const std::string &type2, double forceConstant);

    double getSumOfParticleRadii() const;

    double getSumOfParticleRadiiSquared() const;

    void describe(std::ostream &os) const override;

    double getForceConstant() const;

    virtual double getMaximalForce(double kbt) const noexcept override;

    double calculateEnergy(const Vec3 &x_ij) const override;

    void calculateForce(Vec3 &force, const Vec3 &x_ij) const override;

    void calculateForceAndEnergy(Vec3 &force, double &energy, const Vec3 &x_ij) const override;

    double getCutoffRadius() const override;

    double getCutoffRadiusSquared() const override;

protected:
    friend class readdy::model::KernelContext;

    void configureForTypes(const KernelContext *const ctx, particle_type_type type1, particle_type_type type2) override;

    double sumOfParticleRadii;
    double sumOfParticleRadiiSquared;
    const double forceConstant;
};

class WeakInteractionPiecewiseHarmonic : public PotentialOrder2 {
    using super = PotentialOrder2;
public:
    void describe(std::ostream &os) const override;

    class Configuration {
    public:
        Configuration(const double desiredParticleDistance, const double depthAtDesiredDistance,
                      const double noInteractionDistance);

        friend std::ostream &operator<<(std::ostream &os, const Configuration &configuration);

    private:
        friend class WeakInteractionPiecewiseHarmonic;

        const double desiredParticleDistance, depthAtDesiredDistance, noInteractionDistance, noInteractionDistanceSquared;
    };

    WeakInteractionPiecewiseHarmonic(const std::string &particleType1, const std::string &particleType2,
                                     const double forceConstant, const Configuration &config);

    virtual double getMaximalForce(double kbt) const noexcept override;

    double calculateEnergy(const Vec3 &x_ij) const override;

    void calculateForce(Vec3 &force, const Vec3 &x_ij) const override;

    void calculateForceAndEnergy(Vec3 &force, double &energy, const Vec3 &x_ij) const override;

    double getCutoffRadius() const override;

    double getCutoffRadiusSquared() const override;

protected:
    friend class readdy::model::KernelContext;

    void configureForTypes(const KernelContext *const ctx, particle_type_type type1, particle_type_type type2) override;

    const Configuration conf;
    const double forceConstant;
};

/**
 * Lennard-Jones potential class
 */
class LennardJones : public PotentialOrder2 {
    using super = PotentialOrder2;
public:
    /**
     * Constructs a Lennard-Jones-type potential between two particle types A and B (where possibly A = B) of the
     * form
     *
     * V_LJ(r) = k(epsilon, n, m) [ (sigma/r)^m - (sigma/r)^n ],
     *
     * where n,m are exponent 1 and 2, respectively, with m > n.
     * If shift == true, it will be defined as
     *
     * V_LJ_shifted(r) = V_LJ(r) - V_LJ(cutoffDistance)
     *
     * for r <= cutoffDistance, which makes a difference in energy, but not in force.
     *
     * @param particleType1 particle type A
     * @param particleType2 particle type B
     * @param m first exponent
     * @param n second exponent
     * @param cutoffDistance the cutoff distance
     * @param shift if it should be shifted or not
     * @param epsilon the well depth
     * @param sigma the distance at which the inter-particle potential is zero
     */
    LennardJones(const std::string &particleType1, const std::string &particleType2,
                 unsigned int m, unsigned int n, double cutoffDistance,
                 bool shift, double epsilon, double sigma);

    void describe(std::ostream &os) const override;

    virtual ~LennardJones();

    virtual double calculateEnergy(const Vec3 &x_ij) const override;

    virtual void calculateForce(Vec3 &force, const Vec3 &x_ij) const override;

    virtual void calculateForceAndEnergy(Vec3 &force, double &energy, const Vec3 &x_ij) const override;

    virtual double getCutoffRadius() const override;

    virtual double getCutoffRadiusSquared() const override;

    double getMaximalForce(double kbt) const noexcept override;

protected:
    friend class readdy::model::KernelContext;

    double energy(double r) const;

    virtual void configureForTypes(const KernelContext *const context, particle_type_type type1, particle_type_type type2) override;

    double m, n;
    double cutoffDistance, cutoffDistanceSquared;
    bool shift; // V_LJ_trunc = V_LJ(r) - V_LJ(r_cutoff)
    double epsilon; // depth
    double sigma;
    double k;
};

class ScreenedElectrostatics : public PotentialOrder2 {
    using super = PotentialOrder2;
public:
    ScreenedElectrostatics(const std::string &particleType1, const std::string &particleType2, double electrostaticStrength,
                           double inverseScreeningDepth, double repulsionStrength, double repulsionDistance, unsigned int exponent, double cutoff);

    virtual ~ScreenedElectrostatics();

    virtual double calculateEnergy(const Vec3 &x_ij) const override;

    virtual void calculateForce(Vec3 &force, const Vec3 &x_ij) const override;

    virtual void calculateForceAndEnergy(Vec3 &force, double &energy, const Vec3 &x_ij) const override;

    virtual double getCutoffRadius() const override;

    void describe(std::ostream &os) const override;

    virtual double getCutoffRadiusSquared() const override;

    double getMaximalForce(double kbt) const noexcept override;

protected:
    friend class readdy::model::KernelContext;

    virtual void configureForTypes(const KernelContext *const context, particle_type_type type1, particle_type_type type2) override;

    double electrostaticStrength;
    double inverseScreeningDepth;
    double repulsionStrength;
    double repulsionDistance;
    double exponent;
    double cutoff;
    double cutoffSquared;
};

template<typename T>
const std::string getPotentialName(typename std::enable_if<std::is_base_of<HarmonicRepulsion, T>::value>::type * = 0) {
    return "HarmonicRepulsion";
}

template<typename T>
const std::string
getPotentialName(typename std::enable_if<std::is_base_of<WeakInteractionPiecewiseHarmonic, T>::value>::type * = 0) {
    return "WeakInteractionPiecewiseHarmonic";
}

template<typename T>
const std::string
getPotentialName(typename std::enable_if<std::is_base_of<LennardJones, T>::value>::type * = 0) {
    return "LennardJones";
}

template<typename T>
const std::string
getPotentialName(typename std::enable_if<std::is_base_of<ScreenedElectrostatics, T>::value>::type * = 0) {
    return "ScreenedElectrostatics";
}
NAMESPACE_END(potentials)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
