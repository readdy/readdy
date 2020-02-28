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
}

MPIKernel::MPIKernel(readdy::model::Context ctx) : MPIKernel() {
    _context = std::move(ctx);
    initialize();
}


void MPIKernel::initialize() {
    readdy::model::Kernel::initialize();
    readdy::log::trace("MPIKernel::initialize");
    // Spatial decomposition
    {
        const auto conf = _context.kernelConfiguration();
        std::array<scalar, 3> minDomainWidths{conf.mpi.dx, conf.mpi.dy, conf.mpi.dz};
        _domain = std::make_unique<model::MPIDomain>(rank, worldSize, minDomainWidths, _context);
        _stateModel.setDomain(domain());
    }

    // Description of decomposition
    if (rank == 0) {
        std::string description;
        description += fmt::format("MPI Kernel uses domain decomposition:\n");
        description += fmt::format("--------------------------------\n");
        description += fmt::format(" - Number of domains on axes (x,y,z) ({},{},{})\n",
                                   _domain->nDomainsPerAxis()[0], _domain->nDomainsPerAxis()[1], _domain->nDomainsPerAxis()[2]);
        description += fmt::format(" - Domain widths on axes (x,y,z) ({},{},{})\n",
                                   _context.boxSize()[0] / _domain->nDomainsPerAxis()[0],
                                   _context.boxSize()[1] / _domain->nDomainsPerAxis()[1],
                                   _context.boxSize()[2] / _domain->nDomainsPerAxis()[2]);
        description += fmt::format(" - Used {} ranks of available worldSize {}\n", _domain->nUsedRanks(),
                                   _domain->worldSize());
        readdy::log::info(description);
        if (_domain->nUsedRanks() != _domain->worldSize()) {
            readdy::log::warn("! Number of used workers {} is not equal to what was allocated {} !",
                              _domain->nUsedRanks(), _domain->worldSize());
            readdy::log::warn("You should adapt your number of workers or tune the minimum domain widths");
        }
    }

    // Make a new communicator for only used processes
    if (_domain->nUsedRanks() != _domain->worldSize()) {
        // Create group of all in world
        MPI_Group worldGroup;
        MPI_Comm_group(MPI_COMM_WORLD, &worldGroup);

        // Remove all unnecessary ranks
        MPI_Group usedGroup;
        int removeRanges[1][3];
        removeRanges[0][0] = _domain->nUsedRanks();
        removeRanges[0][1] = _domain->worldSize() - 1;
        removeRanges[0][2] = 1;
        MPI_Group_range_excl(worldGroup, 1, removeRanges, &usedGroup);

        // the new communicator of used ranks -> use this in barriers and stuff
        MPI_Comm_create(MPI_COMM_WORLD, usedGroup, &_commUsedRanks);
    } else {
        _commUsedRanks = MPI_COMM_WORLD;
    }

    // propagate to other classes that need communicator and don't know the kernel
    _stateModel.commUsedRanks() = _commUsedRanks;

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
