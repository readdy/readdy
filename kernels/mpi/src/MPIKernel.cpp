/********************************************************************
 * Copyright © 2019 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 ********************************************************************/

/**
 * « detailed description »
 *
 * @file MPIKernel.cpp
 * @brief « brief description »
 * @author chrisfroe
 * @date 28.05.19
 */

#include <readdy/kernel/mpi/MPIKernel.h>
#include <mpi.h>

namespace readdy::kernel::mpi {

const std::string MPIKernel::name = "MPI";

readdy::model::Kernel *MPIKernel::create() {
    return new MPIKernel();
}

MPIKernel::MPIKernel() : Kernel(name), _stateModel(_data, _context), _actions(this), _observables(this) {
    // Since this kernel should be a drop-in replacement, we need to MPI_Init here
    // given the option that it is already initialized?
    int myWorldSize;
    int myRank;
    int nameLen;
    char myProcessorName[MPI_MAX_PROCESSOR_NAME];

    MPI_Comm_size(MPI_COMM_WORLD, &myWorldSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Get_processor_name(myProcessorName, &nameLen);

    rank = myRank;
    worldSize = myWorldSize;
    processorName = myProcessorName;

    readdy::log::console()->info("pid {} Rank {} / {} is on {}", static_cast<long>(getpid()), rank, worldSize,
                                 processorName);

}

bool MPIKernel::isValidDecomposition(const std::array<std::size_t, 3> nBoxesArr) {
    const auto cutoff = _context.calculateMaxCutoff();
    const auto periodic = _context.periodicBoundaryConditions();
    const auto box = _context.boxSize();
    const auto dx = box[0] / static_cast<scalar>(nBoxesArr[0]);
    const auto dy = box[1] / static_cast<scalar>(nBoxesArr[1]);
    const auto dz = box[2] / static_cast<scalar>(nBoxesArr[2]);
    if (dx < cutoff) {
        if (nBoxesArr[0] == 1 and not periodic[0]) {
            /* smaller than cutoff is ok, when there are no neighbors to be considered */
        } else {
            return false;
            throw std::logic_error(fmt::format("Resulting dx {} (nBoxesX {}, periodicX {}) of MPI box cannot be smaller than cutoff {}", dx, nBoxesArr[0], periodic[0], cutoff));
        }
    }
    if (dy < cutoff) {
        if (nBoxesArr[1] == 1 and not periodic[1]) {
            /* smaller than cutoff is ok, when there are no neighbors to be considered */
        } else {
            return false;
            throw std::logic_error(fmt::format("Resulting dy {} (nBoxesY {}, periodicY {}) of MPI box cannot be smaller than cutoff {}", dy, nBoxesArr[1], periodic[1], cutoff));
        }
    }
    if (dz < cutoff) {
        if (nBoxesArr[2] == 1 and not periodic[2]) {
            /* smaller than cutoff is ok, when there are no neighbors to be considered */
        } else {
            return false;
            throw std::logic_error(fmt::format("Resulting dz {} (nBoxesZ {}, periodicZ {}) of MPI box cannot be smaller than cutoff {}", dz, nBoxesArr[2], periodic[2], cutoff));
        }
    }
    return true;
}

void MPIKernel::initialize() {
    readdy::model::Kernel::initialize();

    // todo do any user configuration here

    {
        const auto conf = _context.kernelConfiguration();
        const auto minDx = conf.mpi.dx;
        const auto minDy = conf.mpi.dy;
        const auto minDz = conf.mpi.dz;
        const auto cutoff = _context.calculateMaxCutoff();

        const auto lx = _context.boxSize()[0];
        const auto ly = _context.boxSize()[1];
        const auto lz = _context.boxSize()[2];
        const auto periodicX = _context.periodicBoundaryConditions()[0];
        const auto periodicY = _context.periodicBoundaryConditions()[1];
        const auto periodicZ = _context.periodicBoundaryConditions()[2];

        auto dx = lx;
        std::size_t nBoxesX{0};
        while (dx > minDx) {
            nBoxesX++;
            dx = lx / static_cast<scalar>(nBoxesX);
        }
        if (nBoxesX > 1) {
            nBoxesX--;
            dx = lx / static_cast<scalar>(nBoxesX);
        } else {
            dx = lx;
        }

        auto dy = ly;
        std::size_t nBoxesY{0};
        while (dy > minDy) {
            nBoxesY++;
            dy = ly / static_cast<scalar>(nBoxesY);
        }
        if (nBoxesY > 1) {
            nBoxesY--;
            dy = ly / static_cast<scalar>(nBoxesY);
        } else {
            dy = ly;
        }

        auto dz = lz;
        std::size_t nBoxesZ{0};
        while (dz > minDz) {
            nBoxesZ++;
            dz = lz / static_cast<scalar>(nBoxesZ);
        }
        if (nBoxesZ > 1) {
            nBoxesZ--;
            dz = lz / static_cast<scalar>(nBoxesZ);
        } else {
            dz = lz;
        }
        
        // todo try +1 or -1 in some directions
        
        auto valid = isValidDecomposition({{nBoxesX, nBoxesY, nBoxesZ}});

        if (not valid) {
            for (intdelX)
        }

        
        auto nBoxes = nBoxesX * nBoxesY * nBoxesZ;
    }
    // todo domain decomposition
    // domains shall be as cubic as possible to optimize communication
    // find out if number of available ranks fits
    // first get user values for dx,dy,dz, find closest and if ranks do not fit, try + or - 1 box (around that) in x, y or z
    // a box must be at least as large as largest cutoff
    


    _stateModel.reactionRecords().clear();
    _stateModel.resetReactionCounts();
    _stateModel.virial() = Matrix33{{{0, 0, 0, 0, 0, 0, 0, 0, 0}}};
}

}

const char *name() {
    return readdy::kernel::mpi::MPIKernel::name.c_str();
}

readdy::model::Kernel *createKernel() {
    return readdy::kernel::mpi::MPIKernel::create();
}
