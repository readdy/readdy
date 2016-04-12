//
// Created by clonker on 08.03.16.
//

#include <boost/predef.h>
#include <readdy/common/Utils.h>

std::string readdy::utils::getOS() {
#if BOOST_OS_WINDOWS
#ifdef _WIN64
    return "win64";
#endif
    return "win32";
#elif BOOST_OS_MACOS
    return "osx";
#else
    return "unix";
#endif
}

bool readdy::utils::isWindows() {
#if BOOST_OS_WINDOWS
    return true;
#endif
    return false;
}
