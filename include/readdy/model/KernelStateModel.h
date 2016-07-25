/**
 * The KernelStateModel keeps information about the current state of the system, like particle positions and forces.
 * A listener can be attached, that fires when the time step changes.
 *
 * @file KernelStateModel.h
 * @brief Defines the KernelStateModel, which gives information about the system's current state.
 * @author clonker
 * @date 18/04/16
 */

#ifndef READDY_MAIN_KERNELSTATEMODEL_H
#define READDY_MAIN_KERNELSTATEMODEL_H

#include <vector>
#include <readdy/common/Types.h>
#include <boost/signals2/signal.hpp>
#include "Particle.h"
#include "Vec3.h"

namespace readdy {
    namespace model {
        class KernelStateModel {
        public:
            KernelStateModel();
            virtual ~KernelStateModel();

            // const accessor methods
            virtual const std::vector<Vec3> getParticlePositions() const = 0;
            virtual const std::vector<Particle> getParticles() const = 0;

            virtual void updateNeighborList() = 0;
            virtual void calculateForces() = 0;

            virtual void addParticle(const Particle &p) = 0;

            virtual void addParticles(const std::vector<Particle> &p) = 0;

            virtual void removeParticle(const Particle &p) = 0;

            virtual double getEnergy() const = 0;

            KernelStateModel(KernelStateModel&& rhs);
            KernelStateModel& operator=(KernelStateModel&& rhs);
        };
    }
}
#endif //READDY_MAIN_KERNELSTATEMODEL_H
