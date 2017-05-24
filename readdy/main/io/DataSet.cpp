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
 * @file DataSet.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 24.05.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/io/DataSet.h>
#include "blosc_filter.h"

namespace readdy {
namespace io {
namespace blosc_compression {

void initialize() {
    static std::atomic_bool initialized{false};
    if (!initialized.load()) {
        initialized = true;
        char *version, *date;
        register_blosc(&version, &date);
        log::debug("registered blosc with version {} ({})", version, date);
    }

}

void activate(hid_t plist) {
    unsigned int cd_values[7];
    cd_values[4] = 4;       /* compression level */
    cd_values[5] = 1;       /* 0: shuffle not active, 1: shuffle active */
    cd_values[6] = BLOSC_ZSTD; /* the actual compressor to use */
    if (H5Pset_filter(plist, FILTER_BLOSC, H5Z_FLAG_OPTIONAL, 7, cd_values) < 0) {
        log::warn("could not set blosc filter!");
        H5Eprint(H5Eget_current_stack(), stderr);
    }
}

}
}
}
