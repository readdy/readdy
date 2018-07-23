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
 * @file BloscFilter.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 06.09.17
 * @copyright GPL-3
 */


#include "readdy/io/BloscFilter.h"
#include "blosc_filter.h"

namespace readdy {
namespace io {


bool BloscFilter::available() const {
    return H5Zfilter_avail(FILTER_BLOSC) > 0;
}

void BloscFilter::activate(h5rd::PropertyList &plist) {
    registerFilter();
    unsigned int cd_values[7];
    // compression level 0-9 (0 no compression, 9 highest compression)
    cd_values[4] = compressionLevel;
    // 0: shuffle not active, 1: shuffle active
    cd_values[5] = shuffle ? 1 : 0;
    // the compressor to use
    switch (compressor) {
        case BloscLZ: {
            cd_values[6] = BLOSC_BLOSCLZ;
            break;
        }
        case LZ4: {
            cd_values[6] = BLOSC_LZ4;
            break;
        }
        case LZ4HC: {
            cd_values[6] = BLOSC_LZ4HC;
            break;
        }
        case SNAPPY: {
            cd_values[6] = BLOSC_SNAPPY;
            break;
        }
        case ZLIB: {
            cd_values[6] = BLOSC_ZLIB;
            break;
        }
        case ZSTD: {
            cd_values[6] = BLOSC_ZSTD;
            break;
        }
    }
    if (H5Pset_filter(plist.id(), FILTER_BLOSC, H5Z_FLAG_OPTIONAL, 7, cd_values) < 0) {
        H5Eprint(H5Eget_current_stack(), stderr);
        throw h5rd::Exception("Could not set blosc filter!");
    }
}

void BloscFilter::registerFilter() {
    static std::atomic_bool initialized{false};
    if (!initialized.load()) {
        char *version, *date;
        register_blosc(&version, &date);
        log::debug("registered blosc with version {} ({})", version, date);
        initialized = true;
    }
}

BloscFilter::BloscFilter(BloscFilter::Compressor compressor, unsigned int compressionLevel, bool shuffle) : compressor(
        compressor), compressionLevel(compressionLevel), shuffle(shuffle) {
    if (compressionLevel > 9) {
        throw std::invalid_argument("Blosc only allows compression levels ranging from 0 to 9.");
    }
}
}
}