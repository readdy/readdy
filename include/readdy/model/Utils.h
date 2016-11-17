/**
 * << detailed description >>
 *
 * @file Utils.h
 * @brief << brief description >>
 * @author clonker
 * @date 17.11.16
 */

#ifndef READDY_MAIN_MODEL_UTILS_H
#define READDY_MAIN_MODEL_UTILS_H

#include <readdy/model/KernelContext.h>

namespace readdy {
namespace model {
namespace util {

double getRecommendedTimeStep(unsigned int N, KernelContext&);
double getMaximumDisplacement(KernelContext&);

}
}
}

#endif //READDY_MAIN_MODEL_UTILS_H
