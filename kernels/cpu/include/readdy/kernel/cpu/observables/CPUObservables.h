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

#include <readdy/model/observables/Observables.h>

namespace readdy {
namespace kernel {
namespace cpu {
class CPUKernel;

namespace observables {

class CPUParticlePosition : public readdy::model::observables::ParticlePosition {
public:
    CPUParticlePosition(CPUKernel *const kernel, unsigned int stride, const std::vector<std::string> &typesToCount = {});

    virtual void evaluate() override;

protected:
    CPUKernel *const kernel;
};

class CPUHistogramAlongAxis : public readdy::model::observables::HistogramAlongAxis {

public:
    CPUHistogramAlongAxis(CPUKernel *const kernel, unsigned int stride,
                       const std::vector<double> &binBorders,
                       const std::vector<std::string> &typesToCount,
                       unsigned int axis);

    virtual void evaluate() override;

protected:
    CPUKernel *const kernel;
    size_t size;
};

class CPUNParticles : public readdy::model::observables::NParticles {
public:

    CPUNParticles(CPUKernel *const kernel, unsigned int stride, std::vector<std::string> typesToCount = {});


    virtual void evaluate() override;

protected:
    CPUKernel *const kernel;
};

class CPUForces : public readdy::model::observables::Forces {
public:
    CPUForces(CPUKernel *const kernel, unsigned int stride, std::vector<std::string> typesToCount = {});

    virtual ~CPUForces() {}

    virtual void evaluate() override;


protected:
    CPUKernel *const kernel;
};


}
}
}
}
#endif //READDY_KERNEL_CPU_OBSERVABLES_H
