//
// Created by clonker on 08.03.16.
//

#include <Utils.h>

std::string readdy::utils::getOS() {
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
#ifdef _WIN64
    return "win64";
#endif
    return "win32";
#elif defined(__APPLE__) && defined(TARGET_OS_MAC)
    return "osx";
#else
    return "unix";
#endif
}

bool readdy::utils::isWindows() {
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
    return true;
#endif
    return false;
}
