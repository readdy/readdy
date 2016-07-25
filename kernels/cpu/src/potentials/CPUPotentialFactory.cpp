/**
 * << detailed description >>
 *
 * @file CPUPotentialFactory.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 13.07.16
 */

#include <readdy/kernel/cpu/potentials/CPUPotentialFactory.h>
#include <readdy/kernel/cpu/CPUKernel.h>
#include <readdy/model/potentials/PotentialsOrder1.h>
#include <readdy/kernel/cpu/potentials/CPUPotentials.h>

namespace readdy {
    namespace kernel {
        namespace cpu {
            namespace potentials {
                CPUPotentialFactory::CPUPotentialFactory(CPUKernel *const kernel) {
                    namespace p = readdy::model::potentials;
                    factory[p::getPotentialName<p::CubePotential>()] = [kernel] {
                        return new CPUCubePotential(kernel);
                    };
                    factory[p::getPotentialName<p::HarmonicRepulsion>()] = [kernel] {
                        return new CPUHarmonicRepulsion(kernel);
                    };
                }
            }
        }
    }
}