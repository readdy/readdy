//
// Created by mho on 04.03.16.
//

#ifndef READDY2_MAIN_PLUGIN_H
#define READDY2_MAIN_PLUGIN_H

#include <string>
#include <memory>
#include <type_traits>
#include <map>
#include <sstream>

namespace readdy {
    namespace plugin {
        class NoSuchPluginException : public std::runtime_error {
        public:
            NoSuchPluginException(const std::string &__arg) : runtime_error(__arg) { }
        };

        class Plugin {
            virtual const std::string getName() const = 0;
        };

        template<typename T>
        class PluginProvider {
            static_assert(std::is_base_of<Plugin, T>::value, "T must extend readdy::plugin::Plugin");

        protected:
            std::map<std::string, std::shared_ptr<T>> plugins;
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
            virtual void add(std::shared_ptr<T> ptr) {
                plugins.emplace(ptr.get()->getName(), ptr);
            }
            virtual void add(T &t) {
                plugins.emplace(t.getName(), std::make_shared<T>(t));
            }
        };
    }
}
#endif //READDY2_MAIN_PLUGIN_H
