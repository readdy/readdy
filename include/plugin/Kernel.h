//
// Created by clonker on 07.03.16.
//

#ifndef READDY2_MAIN_KERNEL_H
#define READDY2_MAIN_KERNEL_H

#include <string>
#include "Plugin.h"

namespace readdy {
    namespace plugin {
        class Kernel {
            virtual std::string getName() = 0;
        };

        class KernelFactory : PluginFactory<Kernel> {
            Kernel* create(std::string name);
        };

    }
}

#endif //READDY2_MAIN_KERNEL_H
