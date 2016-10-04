/**
 * Base class for all plugins.
 *
 * @file Plugin.h
 * @brief Header file containing the definitions for readdy::model::Plugin.
 * @author clonker
 * @date 04.03.16
 */

#ifndef READDY_MAIN_PLUGIN_H
#define READDY_MAIN_PLUGIN_H

#include <string>
#include <memory>
#include <type_traits>
#include <sstream>
#include <unordered_map>

namespace readdy {
namespace model {

class Plugin {
public:
    virtual const std::string &getName() const = 0;

    virtual ~Plugin() {};
};

}
}
#endif //READDY_MAIN_PLUGIN_H
