//
// Created by Moritz Hoffmann on 19/02/16.
//

#include <iostream>
#include <readdy/io/IOUtils.h>
#include <boost/log/trivial.hpp>

readdy::io::IOUtils::IOUtils() {
    BOOST_LOG_TRIVIAL(debug) << "ioutils instantiated";
}
