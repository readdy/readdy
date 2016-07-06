//
// Created by clonker on 08.03.16.
//

#include <boost/predef.h>
#include <readdy/common/Utils.h>
#include <boost/algorithm/string.hpp>
#include <boost/log/trivial.hpp>

namespace readdy {
    namespace utils {
        bool isWindows() {
            #if BOOST_OS_WINDOWS
            return true;
            #endif
            return false;
        }

        std::string getOS() {
            #if BOOST_OS_WINDOWS
            #ifdef _WIN64
            return "win64";
            #endif
            return "win32";
            #elif BOOST_OS_MACOS
            return "osx";
            #else
            return "unix";
            #endif
        }

        namespace testing {
            std::string getPluginsDirectory() {
                    // if we're in conda
                    const char *env = std::getenv("CONDA_ENV_PATH");
                    if(!env) {
                        env = std::getenv("PREFIX");
                        if(env) BOOST_LOG_TRIVIAL(trace) << "PREFIX=" << env;
                    } else {
                        BOOST_LOG_TRIVIAL(trace) << "CONDA_ENV_PATH=" << env;
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
}
