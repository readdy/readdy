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

namespace readdy {
    namespace model {
        class KernelStateModel {
        public:
            virtual ~KernelStateModel();

            virtual void updateModel(bool forces, bool distances) = 0;
            virtual void addParticle(const Particle &p) = 0;
            virtual void addParticles(const std::vector<Particle> &p) = 0;
        };
    }
}
#endif //READDY2_MAIN_KERNELSTATEMODEL_H
