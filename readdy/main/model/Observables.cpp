/**
 * << detailed description >>
 *
 * @file Observables.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 26.04.16
 */

#include <readdy/model/Observables.h>
#include <readdy/model/Kernel.h>

namespace readdy {
    namespace model {

        ObservableBase::~ObservableBase() {
            kernel->deregisterObservable(this);
        }

        const std::string ObservableName<ParticlePositionObservable>::value = "ParticlePositionObservable";
        void ParticlePositionObservable::evaluate() {
            std::vector<Vec3> result = kernel->getKernelStateModel().getParticlePositions();
            *ParticlePositionObservable::result = result;
        }

        void TestCombinerObservable::evaluate() {
            std::vector<double> result;
            const auto &&r1 = obs1->getResult();
            const auto &&r2 = obs2->getResult();

            auto b1 = r1->begin();
            auto b2 = r2->begin();

            for (; b1 != (*r1).end();) {
                result.push_back((*b1) * (*b2));
                ++b1;
                ++b2;
            }

            *TestCombinerObservable::result = result;
        }
    }
}

