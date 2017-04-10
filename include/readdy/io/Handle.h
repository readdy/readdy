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
 * @file Handle.h
 * @brief << brief description >>
 * @author clonker
 * @date 10.04.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <utility>
#include <stdexcept>
#include <H5Ipublic.h>
#include <readdy/common/common.h>
#include <ostream>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(io)

class Handle {
public:
    Handle() : hid_(H5I_INVALID_HID) {};

    Handle(hid_t hid) : hid_(hid) {};

    Handle(Handle &&rhs) : hid_(std::move(rhs.hid_)) {
        rhs.hid_ = H5I_INVALID_HID;
    }

    Handle &operator=(Handle &&rhs) {
        hid_ = std::move(rhs.hid_);
        rhs.hid_ = H5I_INVALID_HID;
        return *this;
    }

    Handle(const Handle &rhs) : hid_(rhs.hid_) {
        if (rhs.valid() && H5Iinc_ref(hid_) < 0) {
            throw std::runtime_error("failure in increasing reference counting of handle");
        }
    }

    Handle &operator=(const Handle &rhs) {
        hid_ = rhs.hid_;
        if (rhs.valid() && H5Iinc_ref(hid_) < 0) {
            throw std::runtime_error("failure in increasing reference counting of handle");
        }
        return *this;
    }

    virtual ~Handle() {
        if (valid() && H5Idec_ref(hid_) < 0) {
            log::warn("failure in decreasing reference counting of handle");
        }
    }

    bool valid() const {
        return hid_ != H5I_INVALID_HID && H5Iis_valid(hid_) != 0;
    }

    const hid_t hid() const {
        return hid_;
    }

    friend std::ostream &operator<<(std::ostream &os, const Handle &handle) {
        os << "Handle(hid: " << handle.hid_ << ")";
        return os;
    }

private:
    hid_t hid_;
};

NAMESPACE_END(io)
NAMESPACE_END(readdy)
