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
 * @file SnapshotObservable.h
 * @brief << brief description >>
 * @author clonker
 * @date 29.08.16
 */
#pragma once

#include <readdy/model/observables/Observable.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(io)
class SnapshotWriter : public readdy::model::observables::Observable<void> {
public:
    SnapshotWriter(std::string fName,
                   model::Kernel *const kernel, unsigned int stride) : Observable(kernel, stride),
                                                                       fileName(std::move(fName)) {
        // todo write config (ie kernel context)
    }

    virtual void evaluate() override {
        // todo write positions (virtual)
        // todo write forces (virtual)
    }

protected:
    std::string fileName;
};
NAMESPACE_END(io)
NAMESPACE_END(readdy)
