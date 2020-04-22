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

MPIKernel::MPIKernel() : MPIKernel(readdy::model::Context{}) {}

// pay attention to order of initialization, which is defined by class hierarchy, then by order of declaration
MPIKernel::MPIKernel(const readdy::model::Context &ctx)
        : Kernel(name, ctx), _domain(_context), _data(&_domain), _actions(this), _observables(this),
          _stateModel(_data, _context, &_domain) {
    // Description of decomposition
    if (_domain.isMasterRank()) {
        readdy::log::info(_domain.describe());
        if (_domain.nUsedRanks() != _domain.worldSize()) {
            readdy::log::warn("! Number of used workers {} is not equal to what was allocated {} !",
                              _domain.nUsedRanks(), _domain.worldSize());
            readdy::log::warn("You should adapt your number of workers or tune the minimum domain widths");
        }
    }

    // Make a new communicator for only used processes
    if (_domain.nUsedRanks() != _domain.worldSize()) {
        // Create group of all in world
        MPI_Group worldGroup;
        MPI_Comm_group(MPI_COMM_WORLD, &worldGroup);

        // Remove all unnecessary ranks
        MPI_Group usedGroup;
        int removeRanges[1][3];
        removeRanges[0][0] = _domain.nUsedRanks();
        removeRanges[0][1] = _domain.worldSize() - 1;
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

const std::string MPIKernel::name = "MPI";

readdy::model::Kernel *MPIKernel::create(const readdy::model::Context &ctx) {
    return new MPIKernel(ctx);
}

readdy::model::Kernel *MPIKernel::create() {
    return new MPIKernel();
}

}

const char *name() {
    return readdy::kernel::mpi::MPIKernel::name.c_str();
}

readdy::model::Kernel *createKernel() {
    return readdy::kernel::mpi::MPIKernel::create();
}
