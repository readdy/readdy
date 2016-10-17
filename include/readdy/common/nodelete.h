/**
 * << detailed description >>
 *
 * @file nodelete.h
 * @brief << brief description >>
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
