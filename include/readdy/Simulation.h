//
// Created by Moritz Hoffmann on 18/02/16.
//

#ifndef READDY_SIMULATION_H
#define READDY_SIMULATION_H

#include <memory>

#include <boost/predef.h>
#include <vector>
#include <readdy/model/Vec3.h>
#include <readdy/model/KernelStateModel.h>
#include <functional>
#include <readdy/model/potentials/Potential.h>
#include <readdy/model/potentials/PotentialOrder2.h>
#include <readdy/model/Observable.h>
#include <readdy/model/potentials/PotentialOrder1.h>

#if BOOST_OS_MACOS
#include <array>
#endif

namespace readdy {
    /**
     * Simulation is the focus of the high-level C++ API of ReaDDy.
     * This is where the system is set up and run for a certain number of
     * timesteps.
     * Things like temperature, boxsize, reactions and potentials belong to the
     * context and are given to the kernel when run() is called.
     */
    class Simulation {
    public:
        Simulation();

        ~Simulation();

        // move
        Simulation(Simulation &&rhs);

        Simulation &operator=(Simulation &&rhs);

        double getKBT() const;

        void setKBT(double kBT);

        readdy::model::Vec3 getBoxSize() const;

        void setBoxSize(double dx, double dy, double dz);
        void setBoxSize(const readdy::model::Vec3& boxSize);

        std::array<bool, 3> getPeriodicBoundary() const;

        void setPeriodicBoundary(std::array<bool, 3> periodic);

        /**
         * Registers an observable with the kernel. A list of available observables can be obtained by
         *
         */
        template<typename T>
        boost::uuids::uuid registerObservable(unsigned int stride, std::function<void(typename T::result_t)>&& callbackFun);


        boost::uuids::uuid registerObservable(const std::string &name, unsigned int stride);
        boost::uuids::uuid registerObservable(readdy::model::ObservableBase& observable);
        std::vector<std::string> getAvailableObservables();
        void deregisterObservable(const boost::uuids::uuid uuid);


        void registerParticleType(const std::string &name, const double diffusionCoefficient, const double radius);

        void deregisterPotential(const boost::uuids::uuid& uuid);

        //const boost::uuids::uuid& registerPotentialOrder1(std::string name, const std::string &type);

        //----------------------
        // Order 1 potentials
        //----------------------
        boost::uuids::uuid registerBoxPotential(std::string particleType, double forceConstant, readdy::model::Vec3 origin, readdy::model::Vec3 extent, bool considerParticleRadius);
        void registerPotentialOrder1(readdy::model::potentials::PotentialOrder1 const* const ptr, const std::string &type);

        //----------------------
        // Order 2 potentials
        //----------------------
        boost::uuids::uuid registerHarmonicRepulsionPotential(std::string particleTypeA, std::string particleTypeB, double forceConstant);
        void registerPotentialOrder2(readdy::model::potentials::PotentialOrder2 const* const ptr, const std::string &type1, const std::string &type2);

        //void registerReaction(const Reaction& reaction);
        //void registerReactionByDescriptor(const std::string descriptor);

        void addParticle(double x, double y, double z, const std::string& type);

        const std::vector<readdy::model::Vec3> getParticlePositions() const;

        void setKernel(const std::string& kernel);

        bool isKernelSelected() const;

        const std::string& getSelectedKernelType() const;

        virtual void run(const readdy::model::time_step_type steps, const double timeStep);

    private:
        struct Impl;
        std::unique_ptr<readdy::Simulation::Impl> pimpl;

        void ensureKernelSelected() const;
    };

    class NoKernelSelectedException : public std::runtime_error {
    public:
        NoKernelSelectedException(const std::string &__arg);
    };

}

#endif //READDY_SIMULATION_H
