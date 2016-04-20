//
// Created by Moritz Hoffmann on 18/02/16.
//

#ifndef READDY2_SIMULATION_H
#define READDY2_SIMULATION_H

#include <memory>

#include <boost/predef.h>
#if BOOST_OS_MACOS
#include <array>
#endif

namespace readdy {
    /**
     * Simulation is the focus of the high-level C++ API of ReaDDy2.
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
        Simulation& operator=(Simulation &&rhs);
        // copy
        Simulation(const Simulation &rhs);
        Simulation& operator=(const Simulation &rhs);

        double getKBT() const;
        void setKBT(double kBT);
        std::array<double, 3> getBoxSize() const;
        void setBoxSize(double dx, double dy, double dz);
        std::array<bool, 3> getPeriodicBoundary() const;
        void setPeriodicBoundary(bool pb_x, bool pb_y, bool pb_z);
        std::string getKernel() const;
        void setKernel(const std::string name);
        void run(const unsigned long steps, const double timestep);
        //void registerParticleType(const std::string name, const double diffusionCoefficient);
        //void registerPotential(const Potential& potential);
        //void registerReaction(const Reaction& reaction);
        //void registerReactionByDescriptor(const std::string descriptor);

        void addParticle(double x, double y, double z, std::string type);
        void setKernel(const std::string kernel);
        bool isKernelSet() const;
        std::string getSelectedKernelType() const;
        virtual void run(const size_t steps, const double timeStep);

    private:
        struct Impl;
        std::unique_ptr<readdy::Simulation::Impl> pimpl;
    };

    class NoKernelSelectedException : public std::runtime_error {
    public:
        NoKernelSelectedException(const std::string &__arg);
    };

}

#endif //READDY2_SIMULATION_H
