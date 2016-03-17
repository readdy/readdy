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
        public:
            virtual const std::string& getName() const = 0;
            virtual ~Plugin() {};
        };

        template<typename T>
        class PluginProvider {
            static_assert(std::is_base_of<Plugin, T>::value, "T must extend readdy::plugin::Plugin");

        protected:
            std::map<std::string, const std::shared_ptr<T>> plugins;
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

            virtual void add(const std::string name, std::shared_ptr<T>&& ptr) {
                std::cout << "before ("<< ptr.get()->getName() << "): " << ptr.use_count() << " .... ";
                plugins.emplace(std::make_pair(name, std::move(ptr)));
                std::cout << "after: " << ptr.use_count() << std::endl;
            }

            virtual void add(const std::string name, T&& t) {
                PluginProvider::add(name, std::move(std::make_shared<T>(std::move(t))));
            }
        };
    }
}
#endif //READDY2_MAIN_PLUGIN_H
