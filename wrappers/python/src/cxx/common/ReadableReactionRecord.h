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
 * @file ReadableReactionRecord.h
 * @brief << brief description >>
 * @author clonker
 * @date 04.10.17
 * @copyright GPL-3
 */

#pragma once

#include <sstream>
#include <readdy/model/Particle.h>

namespace rpy {

struct ReadableReactionRecord {
    ReadableReactionRecord();

    ReadableReactionRecord(const ReadableReactionRecord &) = default;

    ReadableReactionRecord &operator=(const ReadableReactionRecord &) = default;

    ReadableReactionRecord(ReadableReactionRecord &&) = default;

    ReadableReactionRecord &operator=(ReadableReactionRecord &&) = default;

    ~ReadableReactionRecord() = default;

    std::string type;
    std::string reaction_label;
    std::vector<readdy::model::Particle::id_type> educts;
    std::vector<readdy::model::Particle::id_type> products;
    readdy::Vec3::data_arr where;

    friend std::ostream &operator<<(std::ostream &os, const ReadableReactionRecord &rrr);
};

ReadableReactionRecord convert(const readdy::model::reactions::ReactionRecord &reaction, const std::string &name);

std::string repr(const ReadableReactionRecord &rrr);

}
