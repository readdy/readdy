/**
 * Header file containing definitions for various observables. Currently:
 *  - ParticlePositionObservable,
 *  - RadialDistributionObservable,
 *  - CenterOfMassObservable,
 *  - HistogramAlongAxisObservable,
 *  - NParticlesObservable,
 *  - TestCombinerObservable
 *
 * @file Observables.h
 * @brief Header file containing definitions for various observables.
 * @author clonker
 * @date 26.04.16
 * @todo for demonstration purposes, add a more meaningful combiner observable, such as velocity
 */

#ifndef READDY_MAIN_OBSERVABLES_H
#define READDY_MAIN_OBSERVABLES_H

#include <readdy/model/Observable.h>
#include <readdy/model/Vec3.h>
#include <vector>
#include <iostream>
#include <set>

namespace readdy {
namespace model {

class KernelContext;

class ParticlePositionObservable : public Observable<std::vector<Vec3>> {
public:
    ParticlePositionObservable(Kernel *const kernel, unsigned int stride = 1) : Observable(kernel, stride) {}

    ParticlePositionObservable(Kernel *const kernel, unsigned int stride, std::vector<std::string> typesToCount);

    ParticlePositionObservable(Kernel *const kernel, unsigned int stride, std::vector<unsigned int> typesToCount);

    virtual void evaluate() = 0;

protected:
    std::vector<unsigned int> typesToCount;
};

class RadialDistributionObservable : public Observable<std::pair<std::vector<double>, std::vector<double>>> {
public:
    RadialDistributionObservable(Kernel *const kernel, unsigned int stride, std::vector<double> binBorders,
                                 unsigned int typeCountFrom, unsigned int typeCountTo, double particleDensity);

    RadialDistributionObservable(Kernel *const kernel, unsigned int stride, std::vector<double> binBorders,
                                 const std::string &typeCountFrom, const std::string &typeCountTo,
                                 double particleDensity);

    virtual void evaluate() = 0;

    const std::vector<double> &getBinBorders() const;

    void setBinBorders(const std::vector<double> &binBorders);

protected:
    std::vector<double> binBorders;
    std::vector<double> counts;
    unsigned int typeCountFrom, typeCountTo;
    double particleDensity;
};

class CenterOfMassObservable : public Observable<readdy::model::Vec3> {

public:
    CenterOfMassObservable(Kernel *const kernel, unsigned int stride, unsigned int particleType);

    CenterOfMassObservable(Kernel *const kernel, unsigned int stride, const std::vector<unsigned int> &particleTypes);

    CenterOfMassObservable(Kernel *const kernel, unsigned int stride, const std::string &particleType);

    CenterOfMassObservable(Kernel *const kernel, unsigned int stride, const std::vector<std::string> &particleType);

    virtual void evaluate() override;

protected:
    std::set<unsigned int> particleTypes;
};

class HistogramAlongAxisObservable : public Observable<std::vector<double>> {

public:
    HistogramAlongAxisObservable(readdy::model::Kernel *const kernel, unsigned int stride,
                                 std::vector<double> binBorders, std::set<unsigned int> typesToCount,
                                 unsigned int axis);

    HistogramAlongAxisObservable(Kernel *const kernel, unsigned int stride, std::vector<double> binBorders,
                                 std::vector<std::string> typesToCount, unsigned int axis);

    virtual void evaluate() = 0;

protected:
    std::vector<double> binBorders;
    std::set<unsigned int> typesToCount;

    unsigned int axis;
};

class NParticlesObservable : public Observable<std::vector<unsigned long>> {

public:
    NParticlesObservable(Kernel *const kernel, unsigned int stride) : Observable(kernel, stride) {}

    NParticlesObservable(Kernel *const kernel, unsigned int stride, std::vector<std::string> typesToCount);

    NParticlesObservable(Kernel *const kernel, unsigned int stride, std::vector<unsigned int> typesToCount);

    virtual void evaluate() = 0;

protected:
    std::vector<unsigned int> typesToCount;
};

class TestCombinerObservable
        : public CombinerObservable<std::vector<double>, ParticlePositionObservable, ParticlePositionObservable> {
public:

    TestCombinerObservable(Kernel *const kernel, ParticlePositionObservable *obs1, ParticlePositionObservable *obs2,
                           unsigned int stride)
            : CombinerObservable(kernel, obs1, obs2, stride) {
    }

    virtual void evaluate() override;
};

class ForcesObservable : public Observable<std::vector<readdy::model::Vec3>> {

public:
    ForcesObservable(Kernel *const kernel, unsigned int stride) : Observable(kernel, stride), typesToCount({}) {}

    ForcesObservable(Kernel *const kernel, unsigned int stride, std::vector<std::string> typesToCount);

    ForcesObservable(Kernel *const kernel, unsigned int stride, std::vector<unsigned int> typesToCount);

    virtual void evaluate() = 0;

protected:
    std::vector<unsigned int> typesToCount;
};

}
}

#endif //READDY_MAIN_OBSERVABLES_H
