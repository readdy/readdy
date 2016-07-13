/**
 * << detailed description >>
 *
 * @file Utils.h
 * @brief << brief description >>
 * @author clonker
 * @date 13.07.16
 */

#ifndef READDY_MAIN_UTILS_H
#define READDY_MAIN_UTILS_H

#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/log/trivial.hpp>

namespace readdy {
    namespace testing {
        std::string getPluginsDirectory() {
            // test for several environment variables
            const std::string envs[] {"CONDA_ENV_PATH", "CONDA_PREFIX", "PREFIX"};
            const char *env = nullptr;
            for(auto&& key : envs) {
                env = std::getenv(key.c_str());
                if(env) {
                    BOOST_LOG_TRIVIAL(trace) << "Using env-variable for plugin dir prefix " << key << "=" << env;
                    break;
                }
            }
            std::string pluginDir = "lib/readdy_plugins";
            if (env) {
                auto _env = std::string(env);
                if (!boost::algorithm::ends_with(env, "/")) {
                    _env = _env.append("/");
                }
                pluginDir = _env.append(pluginDir);
            } else {
                BOOST_LOG_TRIVIAL(trace) << "no environment variables found that indicate plugins dir.";
            }
            return pluginDir;
        }
    }
}
#endif //READDY_MAIN_UTILS_H
