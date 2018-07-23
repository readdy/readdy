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


//
// Created by mho on 10/08/16.
//

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>

#include <readdy/common/ReaDDyVec3.h>
#include <readdy/common/Timer.h>
#include <readdy/io/BloscFilter.h>
#include <pybind11/numpy.h>
#include "SpdlogPythonSink.h"
#include "ReadableParticle.h"

namespace py = pybind11;
using rvp = py::return_value_policy;

/**
 * Notice: Exporting classes here that are to be shared between prototyping and api module require the base
 * class to use be exported (preferably by the READDY_READDY_API macro defined in common/macros.h).
 */

void exportIO(py::module &);

void exportUtils(py::module& m);

void exportCommon(py::module& common) {
    using namespace pybind11::literals;
    common.def("set_logging_level", [](const std::string &level, bool python_console_out) -> void {
        /*spdlog::drop("console");
        spdlog::set_sync_mode();
        std::vector<spdlog::sink_ptr> sinks;
        auto sink_stdout = std::make_shared<spdlog::sinks::ansicolor_sink>(spdlog::sinks::stdout_sink_mt::instance());
        sinks.push_back(sink_stdout);
        if(python_console_out) {
            auto sink_pysink = std::make_shared<readdy::rpy::pysink>();
            sinks.push_back(sink_pysink);
        }*/

        auto l = [&level] {
            if (level == "trace") {
                return spdlog::level::trace;
            }
            if (level == "debug") {
                return spdlog::level::debug;
            }
            if (level == "info") {
                return spdlog::level::info;
            }
            if (level == "warn") {
                return spdlog::level::warn;
            }
            if (level == "err" || level == "error") {
                return spdlog::level::err;
            }
            if (level == "critical") {
                return spdlog::level::critical;
            }
            if (level == "off") {
                return spdlog::level::off;
            }
            readdy::log::warn("Did not select a valid logging level, setting to debug!");
            return spdlog::level::debug;
        }();

        readdy::log::set_level(l);

        /*auto logger =  spdlog::create("console", std::begin(sinks), std::end(sinks));
        logger->set_pattern("[          ] [%Y-%m-%d %H:%M:%S] [%t] [%l] %v");
        logger->set_level(l);*/

        /*if(python_console_out) {
            // update python loggers level
            py::gil_scoped_acquire gil;
            auto logging_module = pybind11::module::import("logging");
            auto args = py::dict("format"_a="%(message)s");
            switch(logger->level()) {
                case spdlog::level::trace:{
                }
                case spdlog::level::debug: {
                    args["level"] = "DEBUG";
                    break;
                }
                case spdlog::level::info: {
                    args["level"] = "INFO";
                    break;
                }
                case spdlog::level::warn: {
                    args["level"] = "WARNING";
                    break;
                }
                case spdlog::level::err: {
                    args["level"] = "ERROR";
                    break;
                }
                case spdlog::level::critical: {
                }
                case spdlog::level::off: {
                    args["level"] = "CRITICAL";
                    break;
                }
            }
            logging_module.attr("basicConfig")(**args);
        }*/
    }, "Function that sets the logging level. Possible arguments: \"trace\", \"debug\", \"info\", \"warn\", "
                       "\"err\", \"error\", \"critical\", \"off\".", "level"_a, "python_console_out"_a = true);
    common.def("register_blosc_hdf5_plugin", []() -> void {
        readdy::io::BloscFilter filter;
        filter.registerFilter();
    });
    {
        py::module io = common.def_submodule("io", "ReaDDy IO module");
        exportIO(io);
    }
    {
        py::module util = common.def_submodule("util", "ReaDDy util module");
        exportUtils(util);
    }
    {
        py::module perf = common.def_submodule("perf", "ReaDDy performance module");
        py::class_<readdy::util::PerformanceNode>(perf, "PerformanceNode")
                .def("__getitem__", [](const readdy::util::PerformanceNode &self, const std::string &label) -> const readdy::util::PerformanceNode& {
                    return self.child(label);
                }, rvp::reference)
                .def("__len__", [](const readdy::util::PerformanceNode &self) -> std::size_t {
                    return self.n_children();
                })
                .def("__repr__", [](const readdy::util::PerformanceNode &self) -> std::string {
                    return self.describe();
                })
                .def("keys", &readdy::util::PerformanceNode::keys)
                .def("clear", &readdy::util::PerformanceNode::clear)
                .def("time", [](const readdy::util::PerformanceNode &self) -> readdy::util::PerformanceData::time {
                    return self.data().cumulativeTime();
                })
                .def("count", [](const readdy::util::PerformanceNode &self) -> std::size_t {
                    return self.data().count();
                });
    }

    py::class_<rpy::ReadableParticle>(common, "Particle")
            .def_property_readonly("pos", &rpy::ReadableParticle::pos)
            .def_property_readonly("type", &rpy::ReadableParticle::type)
            .def_property_readonly("id", &rpy::ReadableParticle::id);

    py::class_<readdy::Vec3>(common, "Vec")
            .def(py::init<readdy::scalar, readdy::scalar, readdy::scalar>())
            .def(py::self + py::self)
            .def(py::self - py::self)
            .def(readdy::scalar() * py::self)
            .def(py::self / readdy::scalar())
            .def(py::self += py::self)
            .def(py::self *= readdy::scalar())
            .def(py::self == py::self)
            .def(py::self != py::self)
            .def(py::self * py::self)
            .def_property("x", [](const readdy::Vec3 &self) { return self.x; },
                          [](readdy::Vec3 &self, readdy::scalar x) { self.x = x; })
            .def_property("y", [](const readdy::Vec3 &self) { return self.y; },
                          [](readdy::Vec3 &self, readdy::scalar y) { self.y = y; })
            .def_property("z", [](const readdy::Vec3 &self) { return self.x; },
                          [](readdy::Vec3 &self, readdy::scalar z) { self.z = z; })
            .def("toarray", [](const readdy::Vec3 &self) { return self.data; })
            .def("__repr__", [](const readdy::Vec3 &self) {
                std::ostringstream stream;
                stream << self;
                return stream.str();
            })
            .def("__getitem__", [](const readdy::Vec3 &self, unsigned int i) {
                return self[i];
            });

    py::class_<readdy::Matrix33>(common, "Matrix33", py::buffer_protocol())
            .def_buffer([](readdy::Matrix33 &m) -> py::buffer_info {
                return py::buffer_info(
                        m.data().data(),
                        sizeof(readdy::Matrix33::data_arr::value_type),
                        py::format_descriptor<readdy::Matrix33::data_arr::value_type>::format(),
                        2,
                        {readdy::Matrix33::n(), readdy::Matrix33::m()},
                        { sizeof(readdy::Matrix33::data_arr::value_type) * readdy::Matrix33::n(),
                          sizeof(readdy::Matrix33::data_arr::value_type)}
                );
            })
            ;
    /**
     * .def(py::self + py::self)
            //.def(py::self - py::self)
            .def(readdy::scalar() * py::self)
            //.def(py::self / readdy::scalar())
            .def(py::self += py::self)
            .def(py::self *= readdy::scalar())
            .def(py::self == py::self)
            .def(py::self != py::self)
            .def("toarray", [](const readdy::Matrix33 &self) {
                np_array arr {3, 3};
                for(readdy::Matrix33::size_type i = 0; i < readdy::Matrix33::n(); ++i) {
                    for(readdy::Matrix33::size_type j = 0; j < readdy::Matrix33::m(); ++j) {
                        arr.mutable_at(i, j) = self.at(i, j);
                    }
                }
                return arr;
            })
     */
}
