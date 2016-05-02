//
// Created by mho on 04.03.16.
//

#ifndef READDY2_MAIN_PLUGIN_H
#define READDY2_MAIN_PLUGIN_H

#include <string>
#include <memory>
#include <type_traits>
#include <sstream>
#include <unordered_map>

namespace readdy {
    namespace plugin {
        class NoSuchPluginException : public std::runtime_error {
        public:
            NoSuchPluginException(const std::string &__arg) : runtime_error(__arg) { }
        };

        class Plugin {
        public:
            virtual const std::string& getName() const = 0;
            virtual ~Plugin() {};
        };


    }
}
#endif //READDY2_MAIN_PLUGIN_H
