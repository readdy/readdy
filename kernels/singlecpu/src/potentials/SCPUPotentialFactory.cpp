/**
 * @file SingleCPUPotentialFactory.cpp
 * @brief Fill the SingleCPUPotentialFactory with kernel-specific constructors for potentials.
 * @author clonker
 * @date 09.06.16
 */

#include <readdy/kernel/singlecpu/SCPUKernel.h>
#include <readdy/kernel/singlecpu/potentials/SCPUPotentialFactory.h>
#include <readdy/kernel/singlecpu/potentials/SCPUPotentialsOrder1.h>
#include <readdy/kernel/singlecpu/potentials/SCPUPotentialsOrder2.h>

namespace readdy {
namespace kernel {
namespace scpu {
namespace potentials {
SCPUPotentialFactory::SCPUPotentialFactory(SCPUKernel *const kernel)
        : readdy::model::potentials::PotentialFactory(), kernel(kernel) {
    namespace p = readdy::model::potentials;
    // order 1
    factory[p::getPotentialName<p::CubePotential>()] = [kernel] { return new SCPUCubePotential(kernel); };
    factory[p::getPotentialName<p::SpherePotential>()] = [kernel] { return new SCPUSpherePotential(kernel); };
    // order 2
    factory[p::getPotentialName<p::HarmonicRepulsion>()] = [kernel] { return new SCPUHarmonicRepulsion(kernel); };
    factory[p::getPotentialName<p::WeakInteractionPiecewiseHarmonic>()] = [kernel] {
        return new SCPUWeakInteractionPiecewiseHarmonic(kernel);
    };
}

}
}
}
}