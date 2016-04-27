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

        template<typename T>
        class PluginProvider {
            static_assert(std::is_base_of<Plugin, T>::value, "T must extend readdy::plugin::Plugin");

        protected:
            std::unordered_map<std::string, const std::shared_ptr<T>> plugins;
        public:
            virtual const std::shared_ptr<T> get(std::string name) {
                auto it = plugins.find(name);
                if(it != plugins.end()) {
                    return it->second;
                }
                std::ostringstream os;
                os << "Could not load plugin with name \"" << name << "\"";
                throw NoSuchPluginException(os.str());
            }

            virtual void add(const std::string name, const std::shared_ptr<T>&& ptr) {
                plugins.emplace(std::make_pair(name, std::move(ptr)));
            }
        };
    }
}
#endif //READDY2_MAIN_PLUGIN_H
