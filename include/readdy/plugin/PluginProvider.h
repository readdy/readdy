/**
 * << detailed description >>
 *
 * @file PluginProvider.h
 * @brief << brief description >>
 * @author clonker
 * @date 02.05.16
 */

#ifndef READDY_MAIN_PLUGINPROVIDER_H
#define READDY_MAIN_PLUGINPROVIDER_H

#include <readdy/model/Plugin.h>

namespace readdy {
    namespace plugin {

        class NoSuchPluginException : public std::runtime_error {
        public:
            NoSuchPluginException(const std::string &__arg) : runtime_error(__arg) { }
        };

        template<typename T>
        class PluginProvider {
            static_assert(std::is_base_of<readdy::model::Plugin, T>::value, "T must extend readdy::plugin::Plugin");

        protected:
            std::unordered_map<std::string, std::function<T *()>> factory;
        public:
            virtual std::unique_ptr<T> create(const std::string& name) const {
                auto it = factory.find(name);
                if(it != factory.end()) {
                    return std::unique_ptr<T>(it->second());
                }
                std::ostringstream os;
                os << "Could not load plugin with name \"" << name << "\"";
                throw NoSuchPluginException(os.str());
            }

            virtual void add(const std::string name, const std::function<T *()> creator) {
                factory.emplace(std::make_pair(name, creator));
            }
        };
    }
}
#endif //READDY_MAIN_PLUGINPROVIDER_H
