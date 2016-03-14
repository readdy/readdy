//
// Created by clonker on 14.03.16.
//

#ifndef READDY2_MAIN_KERNELPLUGINDECORATOR_H
#define READDY2_MAIN_KERNELPLUGINDECORATOR_H


#include <readdy/plugin/Kernel.h>

namespace readdy {
    namespace plugin {
        class KernelPluginDecorator : public readdy::plugin::Kernel {
        protected:
            const readdy::plugin::Kernel reference;

        public:
            KernelPluginDecorator(const readdy::plugin::Kernel &reference);
        };
    }
}


#endif //READDY2_MAIN_KERNELPLUGINDECORATOR_H
