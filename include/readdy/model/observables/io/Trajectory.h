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
 * << detailed description >>
 *
 * @file Trajectory.h
 * @brief << brief description >>
 * @author chrisfroe
 * @date 12.12.16
 * @copyright GNU Lesser General Public License v3.0
 */
#ifndef READDY_MAIN_TRAJECTORY_H
#define READDY_MAIN_TRAJECTORY_H

#include <array>
#include "readdy/io/File.h"
#include "readdy/model/Kernel.h"

namespace readdy {
namespace model {
namespace observables {

struct TrajectoryEntry {
    using pos_t = readdy::model::Vec3::entry_t;

    TrajectoryEntry(const readdy::model::Particle &p, readdy::model::observables::time_step_type t) : typeId(
            p.getType()), id(p.getId()), px(p.getPos()[0]), py(p.getPos()[1]), pz(p.getPos()[2]), t(t) {}

    TrajectoryEntry(const model::Particle::type_type typeId, const model::Particle::id_type id,
                    const model::Vec3 &position, readdy::model::observables::time_step_type t)
            : typeId(typeId), id(id), px(position[0]), py(position[1]), pz(position[2]), t(t) {}

    readdy::model::Particle::type_type typeId;
    readdy::model::observables::time_step_type t;
    readdy::model::Particle::id_type id;
    pos_t px, py, pz;
};

class Trajectory : public model::observables::Observable<std::vector<std::vector<TrajectoryEntry>>> {

    using observable_entry_t = std::vector<std::vector<TrajectoryEntry>>;
    using base_t = model::observables::Observable<observable_entry_t>;
    using blub = std::string;

public:

    const static std::string TRAJECTORY_GROUP_PATH;
    const static std::string TRAJECTORY_DATA_SET_NAME;

    Trajectory(model::Kernel *const kernel, unsigned int stride, unsigned int fs, readdy::io::File &file);

    ~Trajectory();

    virtual void evaluate();

    virtual void append(observable_entry_t &);

    virtual void flush();

    static readdy::io::h5::data_set_type_t getEntryTypeMemory();

    static readdy::io::h5::data_set_type_t getEntryTypeFile();

protected:
    unsigned int count = 0;
    unsigned int flushStride = 0;
    readdy::io::File &file;
    //Group trajectoryGroup;
    readdy::io::h5::handle_t dataSetHandle;
    readdy::io::h5::handle_t memorySpace;
    readdy::io::h5::handle_t entriesTypeMemory, entriesTypeFile;

};

}
}
}

#endif //READDY_MAIN_TRAJECTORY_H
