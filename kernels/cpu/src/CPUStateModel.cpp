/**
 * << detailed description >>
 *
 * @file CPUStateModel.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 13.07.16
 */

#include <readdy/kernel/cpu/CPUStateModel.h>

namespace readdy {
    namespace kernel {
        namespace cpu {
            void CPUStateModel::calculateForces() {

            }

            const std::vector<readdy::model::Vec3> CPUStateModel::getParticlePositions() const {
                return std::vector<readdy::model::Vec3>();
            }

            const std::vector<readdy::model::Particle> CPUStateModel::getParticles() const {
                return std::vector<readdy::model::Particle>();
            }

            void CPUStateModel::updateNeighborList() {

            }

            void CPUStateModel::addParticle(const readdy::model::Particle &p) {

            }

            void CPUStateModel::addParticles(const std::vector<readdy::model::Particle> &p) {

            }

            void CPUStateModel::removeParticle(const readdy::model::Particle &p) {

            }

            double CPUStateModel::getEnergy() const {
                return 0;
            }


        }
    }
}