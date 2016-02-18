//
// Created by Moritz Hoffmann on 18/02/16.
//

#ifndef READDY2_SIMULATION_H
#define READDY2_SIMULATION_H

#include <array>

namespace ReaDDy {
    class Simulation {
    protected:
        double kBT;
        std::array<double, 3> box_size;
        std::array<bool, 3> periodic_boundary;
    public:
        double getKBT() const;

        void setKBT(double kBT);
        void setBoxSize(double dx, double dy, double dz);
        void setPeriodicBoundary(bool pb_x, bool pb_y, bool pb_z);
    };
}

#endif //READDY2_SIMULATION_H
