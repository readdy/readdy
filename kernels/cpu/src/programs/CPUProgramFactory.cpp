/**
 * << detailed description >>
 *
 * @file CPUProgramFactory.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.06.16
 */

#include <readdy/kernel/cpu/programs/CPUProgramFactory.h>
#include <readdy/kernel/cpu/programs/Reactions.h>
#include <readdy/kernel/cpu/programs/CPUEulerBDIntegrator.h>
#include <readdy/kernel/cpu/programs/UpdateNeighborList.h>
#include <readdy/kernel/cpu/programs/CalculateForces.h>

namespace core_p = readdy::model::programs;

namespace readdy {
    namespace kernel {
        namespace cpu {
            namespace programs {
                CPUProgramFactory::CPUProgramFactory(CPUKernel *kernel) {
                    factory[core_p::getProgramName<core_p::reactions::UncontrolledApproximation>()] = [kernel] {
                        return new reactions::UncontrolledApproximation(kernel);
                    };
                    factory[core_p::getProgramName<core_p::EulerBDIntegrator>()] = [kernel] {
                        return new CPUEulerBDIntegrator(kernel);
                    };
                    factory[core_p::getProgramName<core_p::UpdateNeighborList>()] = [kernel] {
                        return new UpdateNeighborList(kernel);
                    };
                    factory[core_p::getProgramName<core_p::CalculateForces>()] = [kernel] {
                        return new CalculateForces(kernel);
                    };
                    factory[core_p::getProgramName<core_p::reactions::Gillespie>()] = [kernel] {
                        return new reactions::Gillespie(kernel);
                    };
                    factory[core_p::getProgramName<core_p::reactions::GillespieParallel>()] = [kernel] {
                        return new reactions::GillespieParallel(kernel);
                    };
                }
            }
        }
    }
}