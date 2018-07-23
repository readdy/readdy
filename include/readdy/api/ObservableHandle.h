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
 * The observable handle is mainly used as a front for accessing an observable's configurational properties.
 * Mainly it enables easier access to the IO operations offered by the respective observable.
 *
 * @file ObservableHandle.h
 * @brief Definitions for the observable handle.
 * @author clonker
 * @date 06.01.17
 * @copyright BSD-3
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
