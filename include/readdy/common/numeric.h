/**
 * << detailed description >>
 *
 * @file numeric.h
 * @brief << brief description >>
 * @author clonker
 * @date 09.09.16
 */

#ifndef READDY_MAIN_NUMERIC_H
#define READDY_MAIN_NUMERIC_H

#include <type_traits>

namespace readdy {
namespace util {
namespace numeric {

template<typename T, typename D, typename std::enable_if<std::is_arithmetic<T>::value && std::is_arithmetic<D>::value, int>::type = 0>
inline typename std::make_unsigned<T>::type positive_modulo(T i, D n) {
    return static_cast<typename std::make_unsigned<T>::type>((i % n + n) % n);
}

}
}
}
#endif //READDY_MAIN_NUMERIC_H
