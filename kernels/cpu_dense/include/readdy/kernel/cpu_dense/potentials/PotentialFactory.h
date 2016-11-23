/**
 * << detailed description >>
 *
 * @file CPUPotentialFactory.h
 * @brief << brief description >>
 * @author clonker
 * @date 22.11.16
 */


#ifndef READDY_DENSE_CPUPOTENTIALFACTORY_H
#define READDY_DENSE_CPUPOTENTIALFACTORY_H

#include <readdy/model/potentials/PotentialFactory.h>

namespace readdy {
namespace kernel {
namespace cpu_dense {
class Kernel;
namespace potentials {
class PotentialFactory : public readdy::model::potentials::PotentialFactory {
public:
    PotentialFactory(Kernel *const kernel);
};
}
}
}
}
#endif //READDY_DENSE_CPUPOTENTIALFACTORY_H
