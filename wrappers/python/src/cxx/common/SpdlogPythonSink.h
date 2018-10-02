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
 * << detailed description >>
 *
 * @file SpdlogPythonSink.h
 * @brief << brief description >>
 * @author clonker
 * @date 11.08.17
 * @copyright GPL-3
 */

#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/eval.h>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/base_sink.h>
#include <iostream>

namespace readdy {
namespace rpy {

class pysink : public spdlog::sinks::base_sink<std::mutex> {

protected:

    void flush_() override { /* no op */ }

    void sink_it_(const spdlog::details::log_msg &msg) override {
        if(should_log(msg.level)) {

            fmt::memory_buffer formatted;
            formatter_->format(msg, formatted);
            auto message = fmt::to_string(formatted);
            {
                // remove newline
                message.pop_back();
            }

            pybind11::gil_scoped_acquire gil;

            auto logging_module = pybind11::module::import("logging");
            std::string py_level;
            switch (msg.level) {
                case spdlog::level::trace: {
                    /* fall through */
                }
                case spdlog::level::debug: {
                    logging_module.attr("debug")(message);
                    break;
                }
                case spdlog::level::info: {
                    logging_module.attr("info")(message);
                    break;
                }
                case spdlog::level::warn: {
                    logging_module.attr("warning")(message);
                    break;
                }
                case spdlog::level::err: {
                    logging_module.attr("error")(message);
                    break;
                }
                case spdlog::level::critical: {
                    logging_module.attr("critical")(message);
                    break;
                }
                case spdlog::level::off:
                    break;
            }
        }
    }
};

}
}