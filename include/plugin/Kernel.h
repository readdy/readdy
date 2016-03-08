//
// Created by clonker on 07.03.16.
//

#ifndef READDY2_MAIN_KERNEL_H
#define READDY2_MAIN_KERNEL_H

#include "Plugin.h"

namespace readdy {
    namespace plugin {
        class Kernel {
            virtual std::string getName() = 0;
        };

        class KernelProvider : PluginProvider<Kernel> {
        protected:
            KernelProvider();
        public:
            static KernelProvider & getInstance();
            std::unique_ptr<Kernel> get(std::string name);

            // prevent that copies can be created
            KernelProvider(KernelProvider const&) = delete;
            void operator=(KernelProvider const&) = delete;
        };

    }
}

#endif //READDY2_MAIN_KERNEL_H
