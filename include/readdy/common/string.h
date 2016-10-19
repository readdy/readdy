/**
 * String utils (mainly to replace boost's string utils). So far only
 *     - has_suffix, which does exactly what the name says
 *
 * @file string.h
 * @brief Utility methods that deal with std::string's.
 * @author clonker
 * @date 15.10.16
 */

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
