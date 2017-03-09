/********************************************************************
 * Copyright © 2016 Computational Molecular Biology Group,          * 
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * This file is part of ReaDDy.                                     *
 *                                                                  *
 * ReaDDy is free software: you can redistribute it and/or modify   *
 * it under the terms of the GNU Lesser General Public License as   *
 * published by the Free Software Foundation, either version 3 of   *
 * the License, or (at your option) any later version.              *
 *                                                                  *
 * This program is distributed in the hope that it will be useful,  *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of   *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
 * GNU Lesser General Public License for more details.              *
 *                                                                  *
 * You should have received a copy of the GNU Lesser General        *
 * Public License along with this program. If not, see              *
 * <http://www.gnu.org/licenses/>.                                  *
 ********************************************************************/


/**
 * The observable Trajectory is a time series of all particles' positions. A single particle is represented by the
 * TrajectoryEntry. Unlike other observables, Trajectory must be constructed with a File object, because a trajectory
 * is not of much use except writing it to disk.
 *
 * @file Trajectory.h
 * @brief A trajectory keeps track of all particles' positions and saves them to a file.
 * @author chrisfroe
 * @author clonker
 * @date 12.12.16
 * @copyright GNU Lesser General Public License v3.0
 */
#ifndef READDY_MAIN_TRAJECTORY_H
#define READDY_MAIN_TRAJECTORY_H

#include <array>
#include <memory>
#include "readdy/model/Kernel.h"

namespace readdy {

namespace io { class File; }

namespace model {
namespace observables {

struct TrajectoryEntry {
    using pos_t = readdy::model::Vec3::value_t;

    TrajectoryEntry(const readdy::model::Particle &p, time_step_type t) : typeId(
            p.getType()), id(p.getId()), px(p.getPos()[0]), py(p.getPos()[1]), pz(p.getPos()[2]), t(t) {}

    TrajectoryEntry(const model::Particle::type_type typeId, const model::Particle::id_type id,
                    const model::Vec3 &position, time_step_type t)
            : typeId(typeId), id(id), px(position[0]), py(position[1]), pz(position[2]), t(t) {}

    readdy::model::Particle::type_type typeId;
    time_step_type t;
    readdy::model::Particle::id_type id;
    pos_t px, py, pz;

    friend std::ostream &operator<<(std::ostream &, const TrajectoryEntry&);
};

class Trajectory : public model::observables::Observable<std::vector<std::vector<TrajectoryEntry>>> {

    using observable_entry_t = std::vector<std::vector<TrajectoryEntry>>;
    using base_t = model::observables::Observable<observable_entry_t>;

public:

    const static std::string TRAJECTORY_GROUP_PATH;
    const static std::string TRAJECTORY_DATA_SET_NAME;

    Trajectory(model::Kernel *const kernel, unsigned int stride, unsigned int fStride, readdy::io::File &file);

    ~Trajectory();

    virtual void evaluate() override;

    virtual void append(observable_entry_t &);

    virtual void flush() override;

protected:
    void initializeDataSet(io::File &file, const std::string &dataSetName, unsigned int flushStride) override;

    void append() override;

protected:
    unsigned int count = 0;
    unsigned int flushStride = 0;
    struct Impl;
    std::unique_ptr<Impl> pimpl;
};

}
}
}

#endif //READDY_MAIN_TRAJECTORY_H
