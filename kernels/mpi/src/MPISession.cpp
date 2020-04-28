/**
 * @file MPISession.cpp
 * @brief Implementation for MPISession
 * @author chrisfroe
 * @date 23.03.20
 */

#include <readdy/kernel/mpi/MPISession.h>
#include <mpi.h>

namespace readdy::kernel::mpi {

MPISession::MPISession(int &argc, char **argv) {
    char processorName[MPI_MAX_PROCESSOR_NAME];
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &_worldSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
    MPI_Get_processor_name(processorName, &nameLen);
    _processorName = std::string(processorName);

    readdy::log::info("pid {} Rank {} / {} is on {}", static_cast<long>(getpid()), _rank, _worldSize, _processorName);
    waitForDebugger();
}

MPISession::~MPISession() {
    MPI_Finalize();
}

void MPISession::barrier() {
    MPI_Barrier(MPI_COMM_WORLD);
}

void MPISession::waitForDebugger() {
    if (getenv("READDY_MPI_DEBUG") != nullptr) {
        int rank;
        char processorName[MPI_MAX_PROCESSOR_NAME];
        int nameLen;
        MPI_Get_processor_name(processorName, &nameLen);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        static int rankToDebug = 0;
        if (rank == rankToDebug) {
            volatile int i = 0;
            readdy::log::warn("pid {} w/ rank {} on processor {} waiting for debugger",
                              static_cast<unsigned long>(getpid()), rank, processorName);
            while (i == 0) { /* change ’i’ in the debugger, `set variable i=1` */ }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

}
