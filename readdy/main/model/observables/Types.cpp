/********************************************************************
 * Copyright © 2017 Computational Molecular Biology Group,          * 
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
 * @file Types.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 06.09.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/model/observables/io/Types.h>
#include <readdy/model/reactions/ReactionRecord.h>
#include <readdy/model/observables/io/TrajectoryEntry.h>

namespace readdy {
namespace model {
namespace observables {
namespace util {

using namespace h5rd;

CompoundH5Types getVec3Types(h5rd::Object::ParentFileRef ref) {
    NativeCompoundType nct = NativeCompoundTypeBuilder(sizeof(Vec3), std::move(ref))
            .insert<decltype(std::declval<Vec3>().x)>("x", offsetof(Vec3, x))
            .insert<decltype(std::declval<Vec3>().y)>("y", offsetof(Vec3, y))
            .insert<decltype(std::declval<Vec3>().z)>("z", offsetof(Vec3, z))
            .build();
    return std::make_tuple(nct, STDCompoundType(nct));
}

CompoundH5Types getReactionRecordTypes(h5rd::Object::ParentFileRef ref) {
    using entry = readdy::model::reactions::ReactionRecord;
    NativeCompoundType nct = NativeCompoundTypeBuilder(sizeof(entry), std::move(ref))
            .insert<decltype(std::declval<entry>().id)>("reaction_id", offsetof(entry, id))
            .insert<decltype(std::declval<entry>().type)>("reaction_type", offsetof(entry, type))
            .insertStdArray<decltype(std::declval<entry>().educts)>("educts", offsetof(entry, educts))
            .insertStdArray<decltype(std::declval<entry>().products)>("products", offsetof(entry, products))
            .insertArray<scalar, 3>("position", offsetof(entry, where))
            .insertStdArray<decltype(std::declval<entry>().types_from)>("types_from", offsetof(entry, types_from))
            .build();
    return std::make_tuple(nct, STDCompoundType(nct));
}

CompoundH5Types getTrajectoryEntryTypes(h5rd::Object::ParentFileRef ref) {
    using entry = readdy::model::observables::TrajectoryEntry;
    NativeCompoundType nct = NativeCompoundTypeBuilder(sizeof(entry), std::move(ref))
            .insert<decltype(std::declval<entry>().id)>("id", offsetof(entry, id))
            .insert<decltype(std::declval<entry>().typeId)>("typeId", offsetof(entry, typeId))
            .insert<decltype(std::declval<entry>().flavor)>("flavor", offsetof(entry, flavor))
            .insertArray<scalar, 3>("pos", offsetof(entry, pos))
            .build();
    return std::make_tuple(nct, STDCompoundType(nct));
}

}
}
}
}