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


//
// Created by Moritz Hoffmann on 19/02/16.
//

#include <iostream>
#include <readdy/io/IOUtils.h>
#include <hdf5.h>
#include <readdy/common/logging.h>

readdy::io::IOUtils::IOUtils() {
    log::console()->debug("ioutils instantiated");
    hsize_t dim2[] = {10};  /* Dimension size of the second dataset
                                       (in memory */
    auto x = H5S_UNLIMITED;
    std::cout << x << std::endl;
}
