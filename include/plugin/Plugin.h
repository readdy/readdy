//
// Created by mho on 04.03.16.
//

#ifndef READDY2_MAIN_PLUGIN_H
#define READDY2_MAIN_PLUGIN_H

#include <string>
#include <memory>

namespace readdy {
    namespace plugin {
        template<typename T>
        class PluginProvider {
            virtual std::unique_ptr<T> get(std::string name) = 0;
        };
    }
}

#endif //READDY2_MAIN_PLUGIN_H
