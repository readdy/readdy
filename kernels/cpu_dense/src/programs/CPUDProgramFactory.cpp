/**
 * << detailed description >>
 *
 * @file ProgramFactory.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.11.16
 */


#include <readdy/kernel/cpu_dense/programs/CPUDProgramFactory.h>
#include <readdy/kernel/cpu_dense/programs/CPUDEulerBDIntegrator.h>
#include <readdy/kernel/cpu_dense/programs/CPUDUpdateNeighborList.h>
#include <readdy/kernel/cpu_dense/programs/CPUDCalculateForces.h>
#include <readdy/kernel/cpu_dense/programs/CPUDCompartments.h>
#include <readdy/kernel/cpu_dense/programs/reactions/CPUDGillespie.h>
#include <readdy/kernel/cpu_dense/programs/reactions/CPUDUncontrolledApproximation.h>
#include <readdy/kernel/cpu_dense/programs/reactions/CPUDGillespieParallel.h>

namespace core_p = readdy::model::programs;

namespace readdy {
namespace kernel {
namespace cpu_dense {
namespace programs {
CPUDProgramFactory::CPUDProgramFactory(CPUDKernel *kernel) {
    factory[core_p::getProgramName<reactions::CPUDUncontrolledApproximation>()] = [kernel] {
        return new reactions::CPUDUncontrolledApproximation(kernel);
    };
    factory[core_p::getProgramName<CPUDEulerBDIntegrator>()] = [kernel] {
        return new CPUDEulerBDIntegrator(kernel);
    };
    factory[core_p::getProgramName<CPUDUpdateNeighborList>()] = [kernel] {
        return new CPUDUpdateNeighborList(kernel);
    };
    factory[core_p::getProgramName<CPUDCalculateForces>()] = [kernel] {
        return new CPUDCalculateForces(kernel);
    };
    factory[core_p::getProgramName<reactions::CPUDGillespie>()] = [kernel] {
        return new reactions::CPUDGillespie(kernel);
    };
    factory[core_p::getProgramName<reactions::CPUDGillespieParallel>()] = [kernel] {
        return new reactions::CPUDGillespieParallel(kernel);
    };
    factory[core_p::getProgramName<CPUDCompartments>()] = [kernel] {
        return new CPUDCompartments(kernel);
    };
}
}
}
}
}