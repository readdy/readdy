/**
 * << detailed description >>
 *
 * @file ProgramFactory.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.11.16
 */


#include <readdy/kernel/cpu_dense/programs/ProgramFactory.h>
#include <readdy/kernel/cpu_dense/programs/EulerBDIntegrator.h>
#include <readdy/kernel/cpu_dense/programs/UpdateNeighborList.h>
#include <readdy/kernel/cpu_dense/programs/CalculateForces.h>
#include <readdy/kernel/cpu_dense/programs/Compartments.h>
#include <readdy/kernel/cpu_dense/programs/reactions/Gillespie.h>
#include <readdy/kernel/cpu_dense/programs/reactions/UncontrolledApproximation.h>
#include <readdy/kernel/cpu_dense/programs/reactions/GillespieParallel.h>

namespace core_p = readdy::model::programs;

namespace readdy {
namespace kernel {
namespace cpu_dense {
namespace programs {
ProgramFactory::ProgramFactory(Kernel *kernel) {
    factory[core_p::getProgramName<reactions::UncontrolledApproximation>()] = [kernel] {
        return new reactions::UncontrolledApproximation(kernel);
    };
    factory[core_p::getProgramName<EulerBDIntegrator>()] = [kernel] {
        return new EulerBDIntegrator(kernel);
    };
    factory[core_p::getProgramName<UpdateNeighborList>()] = [kernel] {
        return new UpdateNeighborList(kernel);
    };
    factory[core_p::getProgramName<CalculateForces>()] = [kernel] {
        return new CalculateForces(kernel);
    };
    factory[core_p::getProgramName<reactions::Gillespie>()] = [kernel] {
        return new reactions::Gillespie(kernel);
    };
    factory[core_p::getProgramName<reactions::GillespieParallel>()] = [kernel] {
        return new reactions::GillespieParallel(kernel);
    };
    factory[core_p::getProgramName<Compartments>()] = [kernel] {
        return new Compartments(kernel);
    };
}
}
}
}
}