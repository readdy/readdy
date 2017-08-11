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
 * << detailed description >>
 *
 * @file SpdlogPythonSink.h
 * @brief << brief description >>
 * @author clonker
 * @date 11.08.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <pybind11/pybind11.h>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/sink.h>
#include <iostream>

namespace readdy {
namespace rpy {

class pysink : public spdlog::sinks::base_sink<std::mutex> {

public:
    void flush() override {
        using namespace pybind11::literals;
        pybind11::gil_scoped_acquire gil;
        auto tup = pybind11::make_tuple("flush", true);
        pybind11::print("", *tup, "end"_a = "");
    }

protected:
    void _sink_it(const spdlog::details::log_msg &msg) override {
        using namespace pybind11::literals;
        pybind11::gil_scoped_acquire gil;
        pybind11::print(msg.formatted.str(), "end"_a = "");
    }
};

}
}