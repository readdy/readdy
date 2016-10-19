/**
 * Header file containing a console() logger which gives the spd console logger.
 *
 * @file logging.h
 * @brief Header responsible for the logging system.
 * @author clonker
 * @date 14.10.16
 */

#ifndef READDY_MAIN_LOGGING_H
#define READDY_MAIN_LOGGING_H

#include <spdlog/spdlog.h>

namespace readdy {
namespace log {
inline std::shared_ptr<spdlog::logger> console() {
    return spdlog::get("console");
}
}
}
#endif //READDY_MAIN_LOGGING_H
