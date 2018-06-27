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
 * The observable handle is mainly used as a front for accessing an observable's configurational properties.
 * Mainly it enables easier access to the IO operations offered by the respective observable.
 *
 * @file ObservableHandle.h
 * @brief Definitions for the observable handle.
 * @author clonker
 * @date 06.01.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <readdy/model/observables/Observable.h>

NAMESPACE_BEGIN(readdy)

namespace detail {
template<typename T>
using is_observable_type = std::enable_if_t<std::is_base_of<model::observables::ObservableBase, T>::value>;
}

class ObservableHandle {
public:
    /**
     * Creates an empty observable handle.
     */
    explicit ObservableHandle(readdy::model::observables::ObservableBase *const obs) : _observable(obs) {};

    /**
     * This method sets up all required groups / data sets within the given hdf5 file for the observable to dump its
     * data in there and also instructs it to do so.
     *
     * @param file the file
     * @param dataSetName a custom data set name that will be used as postfix in the observable's group path
     * @param flushStride sets the hdf5 chunk size
     */
    void enableWriteToFile(File &file, const std::string &dataSetName, unsigned int flushStride) {
        if (_observable) {
            _observable->enableWriteToFile(file, dataSetName, flushStride);
        } else {
            log::warn("You just tried to enable write to file on a user-provided observable instance, "
                      "this is not supported!");
        }
    }

    std::string type() const {
        if(_observable) {
            return _observable->type();
        }
        throw std::runtime_error("No observable attached to this handle, therefore no type");
    }

    /**
     * Triggers a flush, i.e., everything that can be written will be written
     */
    void flush() {
        if(_observable) {
            _observable->flush();
        }
    }

    readdy::model::observables::ObservableBase *const get() const {
        return _observable;
    }

private:
    /**
     * a reference to the observable this handle belongs to
     */
    readdy::model::observables::ObservableBase *const _observable;
};

NAMESPACE_END(readdy)
