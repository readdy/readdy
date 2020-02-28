/********************************************************************
 * Copyright © 2020 Computational Molecular Biology Group,          *
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
 * << detailed description >>
 *
 * @file MPISession.h
 * @brief RAII wrapper for MPI_Init and MPI_Finalize and some utility
 * @author chrisfroe
 * @date 28.02.20
 */

#pragma once

#include <mpi.h>
#include <readdy/kernel/mpi/Timer.h>

/**
 * Tiny RAII wrapper for MPI Init and Finalize
 *  + a debugging lock
 *  + [@todo] and writing performance data
 */
class MPISession {
    int _worldSize;
    int _rank;
    int nameLen;
    char _processorName[MPI_MAX_PROCESSOR_NAME];

public:
    /** On most MPI implementations MPI_init */
    MPISession(int &argc, char **argv) {
        MPI_Init(&argc, &argv);

        MPI_Comm_size(MPI_COMM_WORLD, &_worldSize);
        MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
        MPI_Get_processor_name(_processorName, &nameLen);

        readdy::log::console()->set_level(spdlog::level::info);
        readdy::log::info("pid {} Rank {} / {} is on {}", static_cast<long>(getpid()), _rank, _worldSize, _processorName);
        waitForDebugger();
        readdy::log::console()->set_level(spdlog::level::warn);
    }

    std::string processorName() {
        return std::string(_processorName);
    }

    int rank() {
        return _rank;
    }

    int worldSize() {
        return _worldSize;
    }

    ~MPISession() {
        MPI_Finalize();
    }

    /**
     * Forces a certain processor with given rank to halt in a while loop and all others to wait at a barrier.
     * This allows gdb to attach to this pid and change `i`, which will continue the program.
     *
     * E.g. do the following on the commandline $ gdb -ex "attach $pid" -ex "set variable i=1" -ex "finish"
     *
     * To enable debugging set the environment variable READDY_MPI_DEBUG,
     * which can be exported to processes via `mpirun`.
     *
     * @param rank of the process calling this function
     * @param processorName name of the process calling this function
     */
    static void waitForDebugger() {
        if (getenv("READDY_MPI_DEBUG") != nullptr) {
            int rank;
            char processorName[MPI_MAX_PROCESSOR_NAME];
            int nameLen;
            MPI_Get_processor_name(processorName, &nameLen);
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);

            static int rankToDebug = 0;
            if (rank == rankToDebug) {
                volatile int i = 0;
                readdy::log::console()->warn("pid {} w/ rank {} on processor {} waiting for debugger",
                                             static_cast<unsigned long>(getpid()), rank, processorName);
                while (i == 0) { /* change ’i’ in the debugger */ }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
};
