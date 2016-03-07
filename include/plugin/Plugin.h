//
// Created by mho on 04.03.16.
//

#ifndef READDY2_MAIN_PLUGIN_H
#define READDY2_MAIN_PLUGIN_H

namespace readdy {
    namespace plugin {
        template<typename T>
        class PluginFactory {
            virtual T *create() = 0;
        };

        class Kernel {

        };

        typedef PluginFactory<Kernel> KernelFactory;
    }
}

#endif //READDY2_MAIN_PLUGIN_H
