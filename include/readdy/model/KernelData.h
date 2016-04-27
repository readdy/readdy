/**
 * << detailed description >>
 *
 * @file KernelData.h
 * @brief << brief description >>
 * @author clonker
 * @date 27.04.16
 */

#ifndef READDY2_MAIN_KERNELDATA_H
#define READDY2_MAIN_KERNELDATA_H

#include "KernelContext.h"
#include "KernelStateModel.h"

namespace readdy {
    namespace model {
        struct KernelData {
            virtual const KernelContext &getContext() const = 0;
            virtual const KernelStateModel &getModel() const = 0;
        };
    }
}

#endif //READDY2_MAIN_KERNELDATA_H
