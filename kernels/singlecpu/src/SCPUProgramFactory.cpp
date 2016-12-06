//
// Created by clonker on 08.04.16.
//

#include <readdy/common/make_unique.h>
#include <readdy/kernel/singlecpu/programs/SCPUProgramFactory.h>
#include <readdy/kernel/singlecpu/programs/SCPUTestProgram.h>
#include <readdy/kernel/singlecpu/programs/SCPUAddParticle.h>
#include <readdy/kernel/singlecpu/programs/SCPUEulerBDIntegrator.h>
#include <readdy/kernel/singlecpu/programs/SCPUCalculateForces.h>
#include <readdy/kernel/singlecpu/programs/SCPUReactionImpls.h>
#include <readdy/kernel/singlecpu/programs/SCPUUpdateNeighborList.h>
#include <readdy/kernel/singlecpu/programs/SCPUCompartments.h>

namespace readdy {
namespace kernel {
namespace scpu {
namespace programs {
SCPUProgramFactory::SCPUProgramFactory(SCPUKernel *kernel) : kernel(kernel) {
    namespace core_p = readdy::model::programs;
    factory[core_p::getProgramName<SCPUTestProgram>()] = [] { return new SCPUTestProgram(); };
    factory[core_p::getProgramName<SCPUAddParticle>()] = [kernel] {
        return new SCPUAddParticle(kernel);
    };
    factory[core_p::getProgramName<SCPUEulerBDIntegrator>()] = [kernel] {
        return new SCPUEulerBDIntegrator(kernel);
    };
    factory[core_p::getProgramName<SCPUUpdateNeighborList>()] = [kernel] {
        return new SCPUUpdateNeighborList(kernel);
    };
    factory[core_p::getProgramName<SCPUCalculateForces>()] = [kernel] {
        return new SCPUCalculateForces(kernel);
    };
    factory[core_p::getProgramName<reactions::SCPUUncontrolledApproximation>()] = [kernel] {
        return new reactions::SCPUUncontrolledApproximation(kernel);
    };
    factory[core_p::getProgramName<core_p::reactions::Gillespie>()] = [kernel] {
        return new reactions::SCPUGillespie(kernel);
    };
    factory[core_p::getProgramName<SCPUCompartments>()] = [kernel] {
        return new SCPUCompartments(kernel);
    };
}
}
}
}
}


