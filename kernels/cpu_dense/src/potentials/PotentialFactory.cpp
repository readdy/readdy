/**
 * << detailed description >>
 *
 * @file PotentialFactory.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.11.16
 */

#include <readdy/kernel/cpu_dense/potentials/PotentialFactory.h>
#include <readdy/kernel/cpu_dense/Kernel.h>
#include <readdy/kernel/singlecpu/potentials/PotentialsOrder1.h>
#include <readdy/kernel/singlecpu/potentials/PotentialsOrder2.h>

namespace readdy {
namespace kernel {
namespace cpu_dense {
namespace potentials {
PotentialFactory::PotentialFactory(Kernel *const kernel) {
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