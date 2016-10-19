//
// Created by clonker on 08.03.16.
//

#include <vector>
#include <sstream>
#include <readdy/common/Utils.h>

namespace readdy {
namespace util {
bool isWindows() {
#if READDY_WINDOWS
    return true;
#endif
    return false;
}

std::string getOS() {
#if READDY_WINDOWS
#ifdef _WIN64
    return "win64";
#endif
    return "win32";
#elif READDY_WINDOWS
    return "osx";
#else
    return "unix";
#endif
}

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

}
}
