/**
 * << detailed description >>
 *
 * @file CPUProgramFactory.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.06.16
 */

#include <readdy/kernel/cpu/programs/CPUProgramFactory.h>
#include <readdy/kernel/cpu/programs/CPUEulerBDIntegrator.h>
#include <readdy/kernel/cpu/programs/CPUUpdateNeighborList.h>
#include <readdy/kernel/cpu/programs/CPUCalculateForces.h>
#include <readdy/kernel/cpu/programs/CPUCompartments.h>
#include <readdy/kernel/cpu/programs/reactions/CPUGillespie.h>
#include <readdy/kernel/cpu/programs/reactions/CPUUncontrolledApproximation.h>
#include <readdy/kernel/cpu/programs/reactions/CPUGillespieParallel.h>
#include <readdy/kernel/cpu/programs/reactions/NextSubvolumesReactionScheduler.h>

namespace core_p = readdy::model::programs;

namespace readdy {
namespace kernel {
namespace cpu {
namespace programs {
CPUProgramFactory::CPUProgramFactory(CPUKernel *kernel) {
    factory[core_p::getProgramName<core_p::reactions::UncontrolledApproximation>()] = [kernel] {
        return new reactions::CPUUncontrolledApproximation(kernel);
    };
    factory[core_p::getProgramName<core_p::EulerBDIntegrator>()] = [kernel] {
        return new CPUEulerBDIntegrator(kernel);
    };
    factory[core_p::getProgramName<core_p::UpdateNeighborList>()] = [kernel] {
        return new CPUUpdateNeighborList(kernel);
    };
    factory[core_p::getProgramName<core_p::CalculateForces>()] = [kernel] {
        return new CPUCalculateForces(kernel);
    };
    factory[core_p::getProgramName<core_p::reactions::Gillespie>()] = [kernel] {
        return new reactions::CPUGillespie(kernel);
    };
    factory[core_p::getProgramName<core_p::reactions::GillespieParallel>()] = [kernel] {
        return new reactions::CPUGillespieParallel(kernel);
    };
    factory[core_p::getProgramName<core_p::reactions::NextSubvolumes>()] = [kernel] {
        return new reactions::NextSubvolumes(kernel);
    };
    factory[core_p::getProgramName<core_p::Compartments>()] = [kernel] {
        return new CPUCompartments(kernel);
    };
}
}
}
}
}