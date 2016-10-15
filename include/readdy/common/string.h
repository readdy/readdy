//
// Created by mho on 15/10/2016.
//

#ifndef READDY_MAIN_STRING_H
#define READDY_MAIN_STRING_H

#include <string>

namespace readdy {
namespace util {
namespace str {

inline bool has_suffix(const std::string &str, const std::string &suffix) {
    return str.size() >= suffix.size() && str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

}
}
}
#endif //READDY_MAIN_STRING_H
