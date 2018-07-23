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
 * @file CompartmentRegistry.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 20.09.17
 * @copyright GPL-3
 */

#include <readdy/model/compartments/CompartmentRegistry.h>
#include <readdy/model/compartments/Compartments.h>
#include <readdy/model/_internal/Util.h>

namespace readdy {
namespace model {
namespace compartments {


Compartment::id_type
CompartmentRegistry::addSphere(const Compartment::conversion_map &conversions, const std::string &uniqueName,
                               const Vec3 &origin, scalar radius, bool largerOrLess) {
    _compartments.emplace_back(std::make_shared<Sphere>(conversions, uniqueName, origin, radius, largerOrLess));
    return _compartments.back()->getId();
}

Compartment::id_type
CompartmentRegistry::addPlane(const Compartment::conversion_map &conversions, const std::string &uniqueName,
                              const Vec3 &normalCoefficients, scalar distance, bool largerOrLess) {
    _compartments.emplace_back(std::make_shared<Plane>(conversions, uniqueName, normalCoefficients, distance,
                                                       largerOrLess));
    return _compartments.back()->getId();
}

CompartmentRegistry::CompartmentRegistry(const ParticleTypeRegistry &types) : _types(types) {}

}
}
}
