/**
 * << detailed description >>
 *
 * @file P1Cube.h
 * @brief << brief description >>
 * @author clonker
 * @date 03.06.16
 */

#ifndef READDY_MAIN_P1CUBE_H
#define READDY_MAIN_P1CUBE_H

#include <string>
#include <readdy/model/potentials/PotentialOrder1.h>

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            namespace potentials {

                template<typename PotentialType>
                const std::string& getPotentialName() {
                    return readdy::model::potentials::_internal::PotentialName<PotentialType>::value;
                }

                class P1Cube : public readdy::model::potentials::PotentialOrder1 {
                public:
                    P1Cube();

                    virtual void evaluate() override;
                };

            }
        }
    }

    namespace model {
        namespace potentials {
            namespace _internal {
                namespace pot = readdy::kernel::singlecpu::potentials;
                template<> struct PotentialName<pot::P1Cube> { static const std::string value; };
            }
        }
    }
}

#endif //READDY_MAIN_P1CUBE_H