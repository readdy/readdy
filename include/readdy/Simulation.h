//
// Created by Moritz Hoffmann on 18/02/16.
//

#ifndef READDY2_SIMULATION_H
#define READDY2_SIMULATION_H

#include <memory>

namespace readdy {
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

    private:
        struct Impl;
        std::unique_ptr<readdy::Simulation::Impl> pimpl;
    };

}

#endif //READDY2_SIMULATION_H
