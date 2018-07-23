/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
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
 * @file Types.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 06.09.17
 * @copyright GPL-3
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