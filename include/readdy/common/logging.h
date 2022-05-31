/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 ********************************************************************/


/**
 * Header file containing a console() logger which gives the spd console logger.
 *
 * @file logging.h
 * @brief Header responsible for the logging system.
 * @author clonker
 * @date 14.10.16
 */

#pragma once
#include <iostream>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include "macros.h"
namespace readdy::log {

READDY_API std::shared_ptr<spdlog::logger> console();

template<typename... Args>
void trace(spdlog::format_string_t<Args...> fmt, Args &&... args) {
    console()->trace(fmt, std::forward<Args>(args)...);
}

template<typename... Args>
void debug(spdlog::format_string_t<Args...> fmt, Args &&... args) {
    console()->debug(fmt, std::forward<Args>(args)...);
}

template<typename... Args>
void critical(spdlog::format_string_t<Args...> fmt, Args &&... args) {
    console()->critical(fmt, std::forward<Args>(args)...);
}

template<typename... Args>
void warn(spdlog::format_string_t<Args...> fmt, Args &&... args) {
    console()->warn(fmt, std::forward<Args>(args)...);
}

template<typename... Args>
constexpr void error(spdlog::format_string_t<Args...> fmt, Args &&... args) {
    console()->error(fmt, std::forward<Args>(args)...);
}

template<typename... Args>
void info(spdlog::format_string_t<Args...> fmt, Args &&... args) {
    console()->info(fmt, std::forward<Args>(args)...);
}

inline void set_level(spdlog::level::level_enum level) {
    console()->set_level(level);
}

class Level {
public:
    explicit Level(spdlog::level::level_enum newLevel) : oldLevel(console()->level()) {
        console()->set_level(newLevel);
    }

    Level(const Level &) = default;

    Level &operator=(const Level &) = default;

    Level(Level &&) = default;

    Level &operator=(Level &&) = default;

    ~Level() {
        console()->set_level(oldLevel);
    }

private:
    spdlog::level::level_enum oldLevel;
};

}
