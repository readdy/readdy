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



void MPIKernel::initialize() {
    readdy::model::Kernel::initialize();

    // Spatial decomposition
    {
        const auto conf = _context.kernelConfiguration();
        std::array<scalar, 3> minDomainWidths {conf.mpi.dx, conf.mpi.dy, conf.mpi.dz};
        domain = std::make_unique<model::MPIDomain>(rank, worldSize, minDomainWidths, _context);

        // set up neighborList ? make halo regions clear for each kernel,
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
