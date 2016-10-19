/**
 * Nodelete callable struct that can be supplied to, e.g., instances of unique_ptr such that the helt object is not
 * deleted along with the instance.
 *
 * @file nodelete.h
 * @brief Header file containing nodelete.
 * @author clonker
 * @date 18.10.16
 */

#ifndef READDY_MAIN_NODELETE_H
#define READDY_MAIN_NODELETE_H
namespace readdy {
namespace util {
struct nodelete {
    template<typename T>
    void operator()(T *) {}
};
}
}
#endif //READDY_MAIN_NODELETE_H
