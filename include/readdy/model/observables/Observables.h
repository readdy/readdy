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
 * Header file containing definitions for various observables. Currently:
 *  - Positions,
 *  - Particles,
 *  - RadialDistribution,
 *  - CenterOfMass,
 *  - HistogramAlongAxis,
 *  - NParticles,
 *  - TestCombiner,
 *  - Forces
 *
 * @file Observables.h
 * @brief Header file containing definitions for various observables.
 * @author clonker
 * @date 26.04.16
 * @todo for demonstration purposes, add a more meaningful combiner observable, such as velocity
 */

#ifndef READDY_MAIN_OBSERVABLES_H
#define READDY_MAIN_OBSERVABLES_H

#include <readdy/model/observables/Observable.h>
#include <readdy/model/Vec3.h>
#include <vector>
#include <iostream>
#include <set>
#include <readdy/model/Particle.h>
#include <readdy/model/reactions/ReactionRecord.h>

namespace readdy {
namespace model {
namespace observables {

class KernelContext;

class Positions : public Observable<std::vector<Vec3>> {
public:

    Positions(Kernel *const kernel, unsigned int stride = 1);

    void flush() override;

    Positions(Kernel *const kernel, unsigned int stride, std::vector<std::string> typesToCount);

    Positions(Kernel *const kernel, unsigned int stride, std::vector<unsigned int> typesToCount);

    virtual ~Positions();

protected:
    void initializeDataSet(io::File &file, const std::string &dataSetName, unsigned int flushStride) override;

    void append() override;

    std::vector<unsigned int> typesToCount;

    struct Impl;
    std::unique_ptr<Impl> pimpl;
};

class Particles
        : public Observable<std::tuple<std::vector<readdy::model::Particle::type_type>, std::vector<readdy::model::Particle::id_type>, std::vector<Vec3>>> {
public:
    Particles(Kernel *const kernel, unsigned int stride = 1);

    virtual ~Particles();

    void flush() override;

protected:
    void initializeDataSet(io::File &file, const std::string &dataSetName, unsigned int flushStride) override;

    void append() override;

    struct Impl;
    std::unique_ptr<Impl> pimpl;
};

class RadialDistribution : public Observable<std::pair<std::vector<double>, std::vector<double>>> {
public:
    RadialDistribution(Kernel *const kernel, unsigned int stride, std::vector<double> binBorders,
                       std::vector<unsigned int> typeCountFrom, std::vector<unsigned int> typeCountTo,
                       double particleToDensity);

    RadialDistribution(Kernel *const kernel, unsigned int stride, std::vector<double> binBorders,
                       const std::vector<std::string> &typeCountFrom, const std::vector<std::string> &typeCountTo,
                       double particleToDensity);

    virtual ~RadialDistribution();

    const std::vector<double> &getBinBorders() const;

    void evaluate() override;

    void flush() override;

protected:

    void setBinBorders(const std::vector<double> &binBorders);

    void initializeDataSet(io::File &file, const std::string &dataSetName, unsigned int flushStride) override;

    void append() override;

    struct Impl;
    std::unique_ptr<Impl> pimpl;
    std::vector<double> binBorders;
    std::vector<double> counts;
    std::vector<unsigned int> typeCountFrom, typeCountTo;
    double particleToDensity;
};

class CenterOfMass : public Observable<readdy::model::Vec3> {

public:
    CenterOfMass(Kernel *const kernel, unsigned int stride, unsigned int particleType);

    CenterOfMass(Kernel *const kernel, unsigned int stride, const std::vector<unsigned int> &particleTypes);

    CenterOfMass(Kernel *const kernel, unsigned int stride, const std::string &particleType);

    CenterOfMass(Kernel *const kernel, unsigned int stride, const std::vector<std::string> &particleType);

    virtual ~CenterOfMass();

    void flush() override;

    void evaluate() override;

protected:
    void initializeDataSet(io::File &file, const std::string &dataSetName, unsigned int flushStride) override;

    void append() override;

    struct Impl;
    std::unique_ptr<Impl> pimpl;
    std::set<unsigned int> particleTypes;
};

class HistogramAlongAxis : public Observable<std::vector<double>> {

public:
    HistogramAlongAxis(readdy::model::Kernel *const kernel, unsigned int stride,
                       std::vector<double> binBorders, std::set<unsigned int> typesToCount,
                       unsigned int axis);

    HistogramAlongAxis(Kernel *const kernel, unsigned int stride, std::vector<double> binBorders,
                       std::vector<std::string> typesToCount, unsigned int axis);

    void flush() override;

    virtual ~HistogramAlongAxis();

protected:
    struct Impl;
    std::unique_ptr<Impl> pimpl;

    void initializeDataSet(io::File &file, const std::string &dataSetName, unsigned int flushStride) override;

    void append() override;

    std::vector<double> binBorders;
    std::set<unsigned int> typesToCount;

    unsigned int axis;
};

class NParticles : public Observable<std::vector<unsigned long>> {

public:
    NParticles(Kernel *const kernel, unsigned int stride);

    NParticles(Kernel *const kernel, unsigned int stride, std::vector<std::string> typesToCount);

    NParticles(Kernel *const kernel, unsigned int stride, std::vector<unsigned int> typesToCount);

    void flush() override;

    virtual ~NParticles();

protected:
    struct Impl;
    std::unique_ptr<Impl> pimpl;

    void initializeDataSet(io::File &file, const std::string &dataSetName, unsigned int flushStride) override;

    void append() override;

    std::vector<unsigned int> typesToCount;
};

class Forces : public Observable<std::vector<readdy::model::Vec3>> {

public:
    Forces(Kernel *const kernel, unsigned int stride);

    Forces(Kernel *const kernel, unsigned int stride, std::vector<std::string> typesToCount);

    Forces(Kernel *const kernel, unsigned int stride, std::vector<unsigned int> typesToCount);

    virtual ~Forces();

    void flush() override;

protected:
    void initializeDataSet(io::File &file, const std::string &dataSetName, unsigned int flushStride) override;

    void append() override;

    struct Impl;
    std::unique_ptr<Impl> pimpl;

    std::vector<unsigned int> typesToCount;
};

class Reactions : public Observable<std::vector<reactions::ReactionRecord>> {
    using super = Observable<std::vector<reactions::ReactionRecord>>;

public:
    Reactions(Kernel *const kernel, unsigned int stride);
    virtual ~Reactions();

    virtual void flush() override;

protected:
    virtual void initialize(Kernel *const kernel) override;

    virtual void initializeDataSet(io::File &file, const std::string &dataSetName, unsigned int flushStride) override;

    virtual void append() override;

    struct Impl;
    std::unique_ptr<Impl> pimpl;
};

class ReactionCounts : public Observable<std::tuple<std::vector<std::size_t>, std::vector<std::size_t>>> {
public:
    ReactionCounts(Kernel *const kernel, unsigned int stride);

    virtual ~ReactionCounts();

    virtual void flush() override;

protected:
    virtual void initialize(Kernel *const kernel) override;

    virtual void initializeDataSet(io::File &file, const std::string &dataSetName, unsigned int flushStride) override;

    virtual void append() override;

    struct Impl;
    std::unique_ptr<Impl> pimpl;
};

}
}
}

#endif //READDY_MAIN_OBSERVABLES_H
