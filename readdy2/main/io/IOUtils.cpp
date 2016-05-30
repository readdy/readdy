//
// Created by Moritz Hoffmann on 19/02/16.
//

#include <iostream>
#include <readdy/io/IOUtils.h>
#include <boost/log/trivial.hpp>
#include <hdf5.h>

readdy::io::IOUtils::IOUtils() {
    BOOST_LOG_TRIVIAL(debug) << "ioutils instantiated";
    hsize_t dim2[] = {10};  /* Dimension size of the second dataset
                                       (in memory */
    auto x = H5S_UNLIMITED;
    std::cout << x << std::endl;
}
