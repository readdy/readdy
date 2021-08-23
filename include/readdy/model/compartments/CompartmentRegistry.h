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
 * @file CompartmentRegistry.h
 * @brief << brief description >>
 * @author clonker
 * @date 20.09.17
 * @copyright BSD-3
 */

#pragma once

#include <readdy/model/ParticleTypeRegistry.h>
#include <readdy/model/_internal/Util.h>
#include "Compartment.h"

namespace readdy::model { class Context; }

namespace readdy::model::compartments {

class CompartmentRegistry {
public:

    using CompartmentVector = std::vector<std::shared_ptr<readdy::model::compartments::Compartment>>;

    Compartment::id_type add(const std::shared_ptr<Compartment> &compartment);

    Compartment::id_type addCapsule(const Compartment::conversion_map &conversions, const std::string &uniqueName,
                                    Vec3 center, Vec3 direction, scalar length, scalar radius, bool inside);

    Compartment::id_type addCapsule(const Compartment::label_conversion_map &conversions, const std::string &uniqueName,
                                    Vec3 center, Vec3 direction, scalar length, scalar radius, bool inside) {
        return addCapsule(_internal::util::transformTypesMap(conversions, *_types), uniqueName, center,
                          direction, length, radius, inside);
    }

    Compartment::id_type addSphere(const Compartment::conversion_map &conversions, const std::string &uniqueName,
                                   const Vec3 &origin, scalar radius, bool largerOrLess);

    Compartment::id_type addSphere(const Compartment::label_conversion_map &conversions, const std::string &uniqueName,
                                   const Vec3 &origin, scalar radius, bool largerOrLess) {
        return addSphere(_internal::util::transformTypesMap(conversions, *_types), uniqueName, origin, radius,
                         !largerOrLess);
    }

    Compartment::id_type addPlane(const Compartment::conversion_map &conversions, const std::string &uniqueName,
                                  const Vec3 &normalCoefficients, scalar distance, bool largerOrLess);

    Compartment::id_type addPlane(const Compartment::label_conversion_map &conversions, const std::string &uniqueName,
                                  const Vec3 &normalCoefficients, scalar distance, bool largerOrLess) {
        return addPlane(_internal::util::transformTypesMap(conversions, *_types), uniqueName, normalCoefficients,
                        distance, largerOrLess);
    }


    [[nodiscard]] const CompartmentVector &get() const {
        return _compartments;
    }

    CompartmentVector &get() {
        return _compartments;
    }

private:
    CompartmentVector _compartments;

    const ParticleTypeRegistry *_types;

    friend readdy::model::Context;
};

}
