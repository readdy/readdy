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
 * @file Filter_detail.h
 * @brief << brief description >>
 * @author clonker
 * @date 06.09.17
 * @copyright BSD-3
 */

#pragma once

#include <H5Zpublic.h>

#include "../File.h"
#include "../Filter.h"
#include "../PropertyList.h"

namespace h5rd {

inline SZIPFilter::SZIPFilter(SZIPFilter::CodingMethod method, unsigned int pixelsPerBlock)
        : codingMethod(method), pixelsPerBlock(pixelsPerBlock) {}

inline bool SZIPFilter::available() const {
    return H5Zfilter_avail(H5Z_FILTER_SZIP) > 0;
}

inline void SZIPFilter::registerFilter() {}

inline void SZIPFilter::activate(PropertyList &plist) {
    if (available()) {
        switch (codingMethod) {
            case Entropy: {
                H5Pset_szip(plist.id(), H5_SZIP_EC_OPTION_MASK, pixelsPerBlock);
                break;
            }
            case NearestNeighbor: {
                H5Pset_szip(plist.id(), H5_SZIP_NN_OPTION_MASK, pixelsPerBlock);
                break;
            }
        }
    } else {
        throw std::runtime_error("Tried activating szip filter even though it was not available");
    }
}

inline bool NBITFilter::available() const {
    return H5Zfilter_avail(H5Z_FILTER_NBIT) > 0;
}

inline void NBITFilter::activate(PropertyList &plist) {
    if (!available()) {
        throw std::runtime_error("Tried activating nbit filter even though it was not available");
    }
    H5Pset_nbit(plist.id());
}

inline void NBITFilter::registerFilter() {}

inline ScaleOffsetFilter::ScaleOffsetFilter(ScaleType scaleType, int scaleFactor)
        : scaleType(scaleType), scaleFactor(scaleFactor) {}

inline bool ScaleOffsetFilter::available() const {
    return H5Zfilter_avail(H5Z_FILTER_SCALEOFFSET) > 0;
}

inline void ScaleOffsetFilter::registerFilter() {}

inline void ScaleOffsetFilter::activate(PropertyList &plist) {
    switch (scaleType) {
        case FloatingPointVariableMinBits: {
            H5Pset_scaleoffset(plist.id(), H5Z_SO_FLOAT_DSCALE, scaleFactor);
            break;
        }
        case FloatingPointFixedMinBits: {
            H5Pset_scaleoffset(plist.id(), H5Z_SO_FLOAT_ESCALE, scaleFactor);
            break;
        }
        case IntegerType: {
            H5Pset_scaleoffset(plist.id(), H5Z_SO_INT, scaleFactor);
            break;
        }
    }
}

inline bool ShuffleFilter::available() const {
    return H5Zfilter_avail(H5Z_FILTER_SHUFFLE) > 0;
}

inline void ShuffleFilter::activate(PropertyList &plist) {
    H5Pset_shuffle(plist.id());
}

inline void ShuffleFilter::registerFilter() {}

inline bool Fletcher32Filter::available() const {
    return H5Zfilter_avail(H5Z_FILTER_FLETCHER32) > 0;
}

inline void Fletcher32Filter::activate(PropertyList &plist) {
    H5Pset_fletcher32(plist.id());
}

inline void Fletcher32Filter::registerFilter() {}


}
