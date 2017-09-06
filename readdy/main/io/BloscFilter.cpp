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
 * @file BloscFilter.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 06.09.17
 * @copyright GNU Lesser General Public License v3.0
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