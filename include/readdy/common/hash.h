/**
 * Defines utility functions for dealing with hashes.
 *
 * @file hash.h
 * @brief Utility functions for hashes.
 * @author clonker
 * @date 14.10.16
 */

#ifndef READDY_MAIN_HASH_H
#define READDY_MAIN_HASH_H

#include <cstddef>

namespace readdy {
namespace util {
namespace hash {
/**
 * Simplified version of boost hash combine.
 * @param seed the seed
 * @param v the value
 */
template<typename T>
void combine(std::size_t& seed, const T& v) {
    seed ^= v + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}
}
}
}
#endif //READDY_MAIN_HASH_H
