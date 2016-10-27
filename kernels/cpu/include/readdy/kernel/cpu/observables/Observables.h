/**
 * << detailed description >>
 *
 * @file Observables.h
 * @brief << brief description >>
 * @author clonker
 * @date 27.10.16
 */

#ifndef READDY_KERNEL_CPU_OBSERVABLES_H
#define READDY_KERNEL_CPU_OBSERVABLES_H

#include <readdy/model/Observables.h>

namespace readdy {
namespace kernel {
namespace cpu {
class CPUKernel;

namespace observables {

class ParticlePosition : public readdy::model::ParticlePositionObservable {
public:
    ParticlePosition(CPUKernel *const kernel, unsigned int stride, const std::vector<std::string> &typesToCount = {});

    virtual void evaluate() override;

protected:
    CPUKernel *const kernel;
};

class HistogramAlongAxis : public readdy::model::HistogramAlongAxisObservable {

public:
    HistogramAlongAxis(CPUKernel *const kernel, unsigned int stride,
                       const std::vector<double> &binBorders,
                       const std::vector<std::string> &typesToCount,
                       unsigned int axis);

    virtual void evaluate() override;

protected:
    CPUKernel *const kernel;
    size_t size;
};

class NParticles : public readdy::model::NParticlesObservable {
public:

    NParticles(CPUKernel *const kernel, unsigned int stride, std::vector<std::string> typesToCount = {});


    virtual void evaluate() override;

protected:
    CPUKernel *const kernel;
};

class Forces : public readdy::model::ForcesObservable {
public:
    Forces(CPUKernel *const kernel, unsigned int stride, std::vector<std::string> typesToCount = {});

    virtual ~Forces() {}

    virtual void evaluate() override;


protected:
    CPUKernel *const kernel;
};


}
}
}
}
#endif //READDY_KERNEL_CPU_OBSERVABLES_H
