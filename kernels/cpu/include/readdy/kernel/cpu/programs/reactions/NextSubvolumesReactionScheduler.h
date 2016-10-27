/**
 * << detailed description >>
 *
 * @file SingleCPUNextSubvolumesReactionScheduler.h
 * @brief << brief description >>
 * @author clonker
 * @date 09.09.16
 */

#ifndef READDY_MAIN_SINGLECPUNEXTSUBVOLUMESREACTIONSCHEDULER_H
#define READDY_MAIN_SINGLECPUNEXTSUBVOLUMESREACTIONSCHEDULER_H

#include <readdy/model/programs/Programs.h>
#include "readdy/kernel/cpu/CPUKernel.h"

namespace readdy {
namespace kernel {
namespace cpu {
namespace programs {
namespace reactions {

class NextSubvolumes : public readdy::model::programs::reactions::NextSubvolumes {
using cell_index_t = unsigned int;
public:
    NextSubvolumes(const CPUKernel *const kernel);
    ~NextSubvolumes();

    virtual void execute() override;

    double getMaxReactionRadius() const;
private:
    struct ReactionEvent;
    struct GridCell;

    CPUKernel const* const kernel;

    // sets up a grid cell (rate, timestamp, next event)
    void setUpCell(GridCell& cell);
    // sets up the computational grid with spacing >= max(reaction radii)
    void setUpGrid();
    // assigns particles to the computational grid
    void assignParticles();
    // schedules the reactions
    void setUpEventQueue();
    // evaluates the collected events
    void evaluateReactions();
    // sets up the neighbor-linked-list structure
    void setUpNeighbors(GridCell& cell);
    GridCell * getCell(const readdy::model::Vec3& particlePosition);
    // fetches a cell at (i,j,k)
    GridCell * getCell(cell_index_t i, cell_index_t j, cell_index_t k);

    // array holding the number of boxes in each spatial direction
    std::array<unsigned int, 3> nCells;
    // size of each box
    readdy::model::Vec3 cellSize;

    std::vector<GridCell> cells;
    std::vector<GridCell*> eventQueue;
};

}
}
}
}
}

#endif //READDY_MAIN_SINGLECPUNEXTSUBVOLUMESREACTIONSCHEDULER_H
