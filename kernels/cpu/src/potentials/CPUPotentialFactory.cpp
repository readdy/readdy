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
#include <readdy/kernel/singlecpu/potentials/PotentialsOrder1.h>
#include <readdy/kernel/singlecpu/potentials/PotentialsOrder2.h>

namespace readdy {
    namespace kernel {
        namespace cpu {
            namespace potentials {
                CPUPotentialFactory::CPUPotentialFactory(CPUKernel *const kernel) {
                    namespace p = readdy::model::potentials;
                    namespace singlecpu_pot = readdy::kernel::singlecpu::potentials;
                    factory[p::getPotentialName<p::CubePotential>()] = [kernel] {
                        return new singlecpu_pot::CubePotential(kernel);
                    };
                    factory[p::getPotentialName<p::HarmonicRepulsion>()] = [kernel] {
                        return new singlecpu_pot::HarmonicRepulsion(kernel);
                    };
                    factory[p::getPotentialName<p::WeakInteractionPiecewiseHarmonic>()] = [kernel] {
                        return new singlecpu_pot::WeakInteractionPiecewiseHarmonic(kernel);
                    };
                }
            }
        }
    }
}