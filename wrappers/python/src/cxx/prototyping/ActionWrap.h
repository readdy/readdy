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
 * @file ActionWrap.h
 * @brief << brief description >>
 * @author clonker
 * @date 08.08.16
 */

#ifndef READDY_MAIN_ACTIONWRAP_H
#define READDY_MAIN_ACTIONWRAP_H

#include <readdy/model/actions/Action.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace readdy {
namespace rpy {
class PyAction : public readdy::model::actions::Action {
    using super = readdy::model::actions::Action;
public:
    using super::Action;

    virtual void perform(util::PerformanceNode &node = util::PerformanceNode::root()) override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_PURE(void, readdy::model::actions::Action, perform,)
    }
};
}
}

#endif //READDY_MAIN_ACTIONWRAP_H
