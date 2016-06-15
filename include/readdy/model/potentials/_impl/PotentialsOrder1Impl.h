/**
 * << detailed description >>
 *
 * @file PotentialsOrder1Impl.h
 * @brief << brief description >>
 * @author clonker
 * @date 15.06.16
 */

namespace pot = readdy::model::potentials;

template<typename KernelType>
const std::string pot::_internal::PotentialName<pot::CubePotential<KernelType>>::value = "Cube";
