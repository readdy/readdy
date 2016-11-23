/**
 * << detailed description >>
 *
 * @file Observables.h
 * @brief << brief description >>
 * @author clonker
 * @date 22.11.16
 */

#ifndef READDY_KERNEL_CPU_DENSE_OBSERVABLES_H
#define READDY_KERNEL_CPU_DENSE_OBSERVABLES_H

#include <readdy/model/Observables.h>

namespace readdy {
namespace kernel {
namespace cpu_dense {
class Kernel;

namespace observables {

class ParticlePosition : public readdy::model::ParticlePositionObservable {
public:
    ParticlePosition(Kernel *const kernel, unsigned int stride, const std::vector<std::string> &typesToCount = {});

    virtual void evaluate() override;

protected:
    Kernel *const kernel;
};

class HistogramAlongAxis : public readdy::model::HistogramAlongAxisObservable {

public:
    HistogramAlongAxis(Kernel *const kernel, unsigned int stride,
                       const std::vector<double> &binBorders,
                       const std::vector<std::string> &typesToCount,
                       unsigned int axis);

    virtual void evaluate() override;

protected:
    Kernel *const kernel;
    size_t size;
};

class NParticles : public readdy::model::NParticlesObservable {
public:

    NParticles(Kernel *const kernel, unsigned int stride, std::vector<std::string> typesToCount = {});


    virtual void evaluate() override;

protected:
    Kernel *const kernel;
};

class Forces : public readdy::model::ForcesObservable {
public:
    Forces(Kernel *const kernel, unsigned int stride, std::vector<std::string> typesToCount = {});

    virtual ~Forces() {}

    virtual void evaluate() override;


protected:
    Kernel *const kernel;
};


}
}
}
}
#endif //READDY_KERNEL_CPU_DENSE_OBSERVABLES_H
