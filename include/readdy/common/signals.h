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
 * Strongly simplified signals implementation similar to the one that can be found in, e.g., boost.signals2 and halmd.
 * The implementation details are orientated on what can be found in halmd
 * (https://github.com/halmd-org/halmd/blob/testing/halmd/utility/signal.hpp), licensed under LGPL-3+.
 *
 * @file signals.h
 * @brief A signals implementation
 * @author clonker
 * @date 14.10.16
 */

#pragma once
#include <functional>
#include <list>
#include <memory>
#include <algorithm>
#include "make_unique.h"
#include "logging.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(signals)

class connection {
public:
    using slots_container_type = std::shared_ptr<void>;
    using disconnector_type = std::function<void(std::shared_ptr<void>)>;

    virtual ~connection() = default;

    bool connected() const {
        return !slots_container.expired();
    }

    void disconnect() {
        auto container = slots_container.lock();
        if (container) {
            disconnector(container);
            slots_container.reset();
        }
    }

    connection(connection const &) = delete;

    connection &operator=(connection const &) = delete;

    connection(connection &&rhs) : disconnector(std::move(rhs.disconnector)),
                                   slots_container(std::move(rhs.slots_container)) {
    }

    connection &operator=(connection &&rhs) {
        disconnector = std::move(rhs.disconnector);
        slots_container = std::move(rhs.slots_container);
        return *this;
    }

private:

    template<typename T>
    friend
    class slots_container;

    friend class scoped_connection;

    connection(slots_container_type slots_container, disconnector_type disconnector)
            : slots_container(slots_container), disconnector(disconnector) {
    }


    disconnector_type disconnector;
    std::weak_ptr<void> slots_container;
};

class scoped_connection {
public:
    virtual ~scoped_connection() {
        disconnect();
    }

    explicit scoped_connection(connection&& connection) : conn(std::move(connection)) {}

    void disconnect() {
        conn.disconnect();
    }

    scoped_connection(const scoped_connection &) = delete;

    scoped_connection &operator=(const scoped_connection &) = delete;

    scoped_connection(scoped_connection &&rhs) noexcept : conn(std::move(rhs.conn)) {}

    scoped_connection &operator=(scoped_connection &&rhs) noexcept {
        conn = std::move(rhs.conn);
        return *this;
    }

private:
    readdy::signals::connection conn;
};

template<typename T>
class slots_container {
    std::shared_ptr<std::list<T>> container;

    struct eraser {
        explicit eraser(typename std::list<T>::iterator it) : it(it) {}

        void operator()(std::shared_ptr<void> ptr) {
            auto slots = std::static_pointer_cast<std::list<T>>(ptr);
            slots->erase(it);
        }

    private:
        typename std::list<T>::iterator it;
    };

public:

    slots_container() : container(std::make_shared<std::list<T>>()) {}

    connection connect(const T &slot) {
        auto it = container->insert(container->end(), slot);
        return connection(container, eraser(it));
    }

    scoped_connection connect_scoped(const T &slot) {
        return scoped_connection(connect(slot));
    }

    const std::size_t n_slots() const {
        return container->size();
    }

    typename std::list<T>::iterator begin() {
        return container->begin();
    }

    typename std::list<T>::iterator end() {
        return container->end();
    }

};

template<typename F>
class signal;

template<typename... Args>
class signal<void(Args...)> : public slots_container<std::function<void(Args...)>> {

public:
    using slot_type = std::function<void(Args...)>;

    void fire_signal(Args... args) {
        for(auto &slot : *this) {
            slot(std::forward<Args>(args)...);
        }
    }

    void operator()(Args... args) {
        fire_signal(std::forward<Args>(args)...);
    }
};

NAMESPACE_END(signals)
NAMESPACE_END(readdy)
