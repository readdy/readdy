/**
 * Header for the container of time dependent data in a kernel.
 *
 * @file KernelStateModel.h
 * @brief Defines KernelStateModel.
 * @author clonker
 * @date 18/04/16
 */

#ifndef READDY2_MAIN_KERNELSTATEMODEL_H
#define READDY2_MAIN_KERNELSTATEMODEL_H

#include <vector>
#include "Particle.h"
#include "Vec3.h"

namespace readdy {
    namespace model {
        class KernelStateModel {
        public:
            typedef unsigned long int time_step_type;

            virtual ~KernelStateModel();

            virtual void updateModel(time_step_type t, bool forces, bool distances) = 0;

            virtual void addParticle(const Particle &p) = 0;

            virtual void addParticles(const std::vector<Particle> &p) = 0;

            virtual std::vector<Vec3> getParticlePositions() = 0;
        };
    }
}
#endif //READDY2_MAIN_KERNELSTATEMODEL_H
