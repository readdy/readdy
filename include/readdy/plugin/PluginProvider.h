/**
 * << detailed description >>
 *
 * @file PluginProvider.h
 * @brief << brief description >>
 * @author clonker
 * @date 02.05.16
 */

#ifndef READDY2_MAIN_PLUGINPROVIDER_H
#define READDY2_MAIN_PLUGINPROVIDER_H

#include "Plugin.h"

namespace readdy {
    namespace plugin {
        template<typename T>
        class PluginProvider {
            static_assert(std::is_base_of<Plugin, T>::value, "T must extend readdy::plugin::Plugin");

        protected:
            std::unordered_map<std::string, const std::shared_ptr<T>> plugins;
        public:
            virtual const std::shared_ptr<T> get(std::string name) const {
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
#endif //READDY2_MAIN_PLUGINPROVIDER_H
