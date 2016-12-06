/********************************************************************
 * Copyright © 2016 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * This file is part of ReaDDy.                                     *
 *                                                                  *
 * ReaDDy is free software: you can redistribute it and/or modify   *
 * it under the terms of the GNU Lesser General Public License as   *
 * published by the Free Software Foundation, either version 3 of   *
 * the License, or (at your option) any later version.              *
 *                                                                  *
 * This program is distributed in the hope that it will be useful,  *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of   *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
 * GNU Lesser General Public License for more details.              *
 *                                                                  *
 * You should have received a copy of the GNU Lesser General        *
 * Public License along with this program. If not, see              *
 * <http://www.gnu.org/licenses/>.                                  *
 ********************************************************************/


/**
 * The PluginProvider is a factory superclass for plugin-typed objects. If no plugin with the requested name was found,
 * a NoSuchPluginException is thrown.
 *
 * @file PluginProvider.h
 * @brief This header contains the PluginProvider class and the NoSuchPluginException.
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
    NoSuchPluginException(const std::string &__arg) : runtime_error(__arg) {}
};

template<typename T>
class PluginProvider {
    static_assert(std::is_base_of<readdy::model::Plugin, T>::value, "T must extend readdy::plugin::Plugin");

protected:
    std::unordered_map<std::string, std::function<T *()>> factory;
public:
    virtual std::unique_ptr<T> create(const std::string &name) const {
        auto it = factory.find(name);
        if (it != factory.end()) {
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
