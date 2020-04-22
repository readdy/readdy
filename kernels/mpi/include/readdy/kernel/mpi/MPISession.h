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
 * @file MPISession.h
 * @brief RAII wrapper for MPI_Init and MPI_Finalize and some utility
 * @author chrisfroe
 * @date 28.02.20
 */

#pragma once

#include <readdy/common/logging.h>

namespace readdy::kernel::mpi {

/** Tiny RAII wrapper for MPI Init and Finalize, and a debugging lock */
class MPISession {
    int _worldSize;
    int _rank;
    int nameLen;
    std::string _processorName;

public:
    MPISession(int &argc, char **argv);

    const std::string &processorName() {
        return _processorName;
    }

    int rank() {
        return _rank;
    }

    int worldSize() {
        return _worldSize;
    }

    ~MPISession();

    MPISession(const MPISession &) = delete;

    MPISession(MPISession &&) = delete;

    MPISession &operator=(const MPISession &) = delete;

    MPISession &operator=(MPISession &&) = delete;

    /* Create a barrier for all processors in the world communicator */
    static void barrier();

    /**
     * Forces a certain processor with given rank to halt in a while loop and all others to wait at a barrier.
     * This allows gdb to attach to this pid and change `i`, which will continue the program.
     *
     * E.g. do the following on the commandline $ gdb -ex "attach $pid" -ex "set variable i=1" -ex "finish"
     *
     * To enable debugging set the environment variable READDY_MPI_DEBUG,
     * which can be exported to processes via `mpirun`.
     */
    static void waitForDebugger();
};

}
