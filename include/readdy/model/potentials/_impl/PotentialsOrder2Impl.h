/**
 * << detailed description >>
 *
 * @file PotentialsOrder2Impl.h.h
 * @brief << brief description >>
 * @author clonker
 * @date 15.06.16
 */

namespace pot = readdy::model::potentials;

template<typename KernelType>
const std::string pot::_internal::PotentialName<pot::HarmonicRepulsion<KernelType>>::value = "HarmonicRepulsion";
