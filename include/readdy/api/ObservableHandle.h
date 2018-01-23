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
     * Type of an observable_id. Can be used to remove an observable from a Simulation instance.
     */
    using observable_id = std::size_t;

    /**
     * Creates an empty observable handle.
     */
    ObservableHandle();

    /**
     * Creates a new observable handle belonging to some existing observable instance
     * @param id the id of that observable instance
     * @param observable a pointer to the instance itself
     */
    ObservableHandle(observable_id id, model::observables::ObservableBase *observable);

    /**
     * This method sets up all required groups / data sets within the given hdf5 file for the observable to dump its
     * data in there and also instructs it to do so.
     *
     * @param file the file
     * @param dataSetName a custom data set name that will be used as postfix in the observable's group path
     * @param flushStride sets the hdf5 chunk size
     */
    void enableWriteToFile(File &file, const std::string &dataSetName, unsigned int flushStride);

    /**
     * Retrieves the id of this handle's observable
     * @return the id
     */
    observable_id getId() const;

    std::string getType() const;

    /**
     * Triggers a flush, i.e., everything that can be written will be written
     */
    void flush();

private:
    /**
     * the id
     */
    observable_id id;
    /**
     * a reference to the observable this handle belongs to
     */
    readdy::model::observables::ObservableBase *const observable;
};

NAMESPACE_END(readdy)

#include "bits/ObservableHandle_misc.h"
