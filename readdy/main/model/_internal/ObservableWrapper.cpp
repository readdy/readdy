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


#include <readdy/model/_internal/ObservableWrapper.h>
#include <readdy/model/Kernel.h>

/**
 * << detailed description >>
 *
 * @file ObservableWrapper.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 10.05.16
 */


void readdy::model::observables::ObservableWrapper::operator()(observables::time_step_type t) {
    callback(t);
}

readdy::model::observables::ObservableWrapper::ObservableWrapper(readdy::model::Kernel *const kernel,
                                                                 const observables::observable_type &observable,
                                                                 unsigned int stride)
        : ObservableBase(kernel, stride), observable(observable) {
}

void readdy::model::observables::ObservableWrapper::evaluate() {
    observable(t_current);
}

void readdy::model::observables::ObservableWrapper::flush() {
    throw std::runtime_error("not supported");
}

void
readdy::model::observables::ObservableWrapper::initializeDataSet(readdy::io::File &file, const std::string &dataSetName,
                                                                 unsigned int flushStride) {
    throw std::runtime_error("not supported");
}

void readdy::model::observables::ObservableWrapper::append() {
    throw std::runtime_error("not supported");
}









