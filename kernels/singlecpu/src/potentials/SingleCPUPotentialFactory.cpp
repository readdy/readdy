/**
 * << detailed description >>
 *
 * @file SingleCPUPotentialFactory.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 09.06.16
 */

#include <readdy/kernel/singlecpu/SingleCPUKernel.h>
#include <readdy/kernel/singlecpu/potentials/SingleCPUPotentialFactory.h>
#include <readdy/kernel/singlecpu/potentials/PotentialsOrder1.h>
#include <readdy/kernel/singlecpu/potentials/PotentialsOrder2.h>

namespace readdy {
namespace kernel {
namespace singlecpu {
namespace potentials {
SingleCPUPotentialFactory::SingleCPUPotentialFactory(SingleCPUKernel *const kernel)
        : readdy::model::potentials::PotentialFactory(), kernel(kernel) {
    namespace p = readdy::model::potentials;
    // order 1
    factory[p::getPotentialName<p::CubePotential>()] = [kernel] { return new CubePotential(kernel); };
    // order 2
    factory[p::getPotentialName<p::HarmonicRepulsion>()] = [kernel] { return new HarmonicRepulsion(kernel); };
    factory[p::getPotentialName<p::WeakInteractionPiecewiseHarmonic>()] = [kernel] {
        return new WeakInteractionPiecewiseHarmonic(kernel);
    };
}

}
}
}
}