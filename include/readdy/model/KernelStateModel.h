/**
 * Header for the container of time dependent data in a kernel.
 *
 * @file KernelStateModel.h
 * @brief Defines KernelStateModel.
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
        typedef boost::signals2::signal<void()> signal_t;
        public:
            KernelStateModel();
            virtual ~KernelStateModel();

            // const accessor methods
            virtual const time_step_type getCurrentTimeStep() const = 0;
            virtual const std::vector<Vec3> getParticlePositions() const = 0;
            virtual const std::vector<Particle> getParticles() const = 0;

            virtual void updateModel(time_step_type t, bool forces) = 0;

            virtual void addParticle(const Particle &p) = 0;

            virtual void addParticles(const std::vector<Particle> &p) = 0;

            virtual void removeParticle(const Particle &p) = 0;

            virtual double getEnergy() const = 0;

            boost::signals2::connection addListener(const signal_t::slot_type& l);

            KernelStateModel(KernelStateModel&& rhs);
            KernelStateModel& operator=(KernelStateModel&& rhs);

        protected:
            std::unique_ptr<signal_t> signal;

            void fireTimeStepChanged();
        };
    }
}
#endif //READDY_MAIN_KERNELSTATEMODEL_H
