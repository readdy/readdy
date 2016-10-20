//
// Created by clonker on 08.04.16.
//

#include <readdy/common/make_unique.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUProgramFactory.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUTestProgram.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUAddParticleProgram.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUEulerBDIntegrator.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUCalculateForces.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUReactionImpls.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUUpdateNeighborList.h>
#include <readdy/kernel/singlecpu/programs/Compartments.h>

namespace readdy {
namespace kernel {
namespace singlecpu {
namespace programs {
SingleCPUProgramFactory::SingleCPUProgramFactory(SingleCPUKernel *kernel) : kernel(kernel) {
    namespace core_p = readdy::model::programs;
    factory[core_p::getProgramName<SingleCPUTestProgram>()] = [] { return new SingleCPUTestProgram(); };
    factory[core_p::getProgramName<SingleCPUAddParticleProgram>()] = [kernel] {
        return new SingleCPUAddParticleProgram(kernel);
    };
    factory[core_p::getProgramName<SingleCPUEulerBDIntegrator>()] = [kernel] {
        return new SingleCPUEulerBDIntegrator(kernel);
    };
    factory[core_p::getProgramName<SingleCPUUpdateNeighborList>()] = [kernel] {
        return new SingleCPUUpdateNeighborList(kernel);
    };
    factory[core_p::getProgramName<SingleCPUCalculateForces>()] = [kernel] {
        return new SingleCPUCalculateForces(kernel);
    };
    factory[core_p::getProgramName<reactions::UncontrolledApproximation>()] = [kernel] {
        return new reactions::UncontrolledApproximation(kernel);
    };
    factory[core_p::getProgramName<core_p::reactions::Gillespie>()] = [kernel] {
        return new reactions::Gillespie(kernel);
    };
    factory[core_p::getProgramName<Compartments>()] = [kernel] {
        return new Compartments(kernel);
    };
}
}
}
}
}


