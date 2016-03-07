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

        class KernelFactory : PluginFactory<Kernel> {
        protected:
            KernelFactory();
        public:
            static KernelFactory& getInstance();
            std::unique_ptr<Kernel> create(std::string name);

            // prevent that copies can be created
            KernelFactory(KernelFactory const&) = delete;
            void operator=(KernelFactory const&) = delete;
        };

    }
}

#endif //READDY2_MAIN_KERNEL_H
