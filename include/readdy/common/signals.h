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
#include <utility>
#include "logging.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(signals)

class connection {
public:
    /**
     * pointer to a slots container
     */
    using slots_container_type = std::shared_ptr<void>;
    /**
     * a disconnect function
     */
    using disconnector_type = std::function<void(std::shared_ptr<void>)>;

    /**
     * we can destroy a connection instance as usual
     */
    virtual ~connection() = default;

    /**
     * checks if the container is still alive
     * @return true if the container is alive
     */
    bool connected() const {
        return !slots_container.expired();
    }

    /**
     * disconnects from the container
     */
    void disconnect() {
        auto container = slots_container.lock();
        if (container) {
            disconnector(container);
            slots_container.reset();
        }
    }

    /**
     * no copying of connections
     */
    connection(connection const &) = delete;

    /**
     * no copy assign of connections
     */
    connection &operator=(connection const &) = delete;

    /**
     * connections can be moved
     * @param rhs the other connection
     */
    connection(connection &&rhs) noexcept : disconnector(std::move(rhs.disconnector)), 
                                            slots_container(std::move(rhs.slots_container)) { }

    /**
     * connections can be move-assigned
     * @param rhs the other connection
     * @return reference to this
     */
    connection &operator=(connection &&rhs) noexcept {
        disconnector = std::move(rhs.disconnector);
        slots_container = std::move(rhs.slots_container);
        return *this;
    }

private:

    template<typename T>
    friend class slots_container;

    friend class scoped_connection;

    /**
     * private constructor of connection, should only be called from slots_container or scoped_connection.
     * @param slots_container the respective container
     * @param disconnector a disconnector
     */
    connection(const slots_container_type &slots_container, disconnector_type disconnector)
            : slots_container(slots_container), disconnector(std::move(disconnector)) {
    }


    disconnector_type disconnector;
    std::weak_ptr<void> slots_container;
};

class scoped_connection {
public:
    /**
     * upon destruction, disconnect
     */
    virtual ~scoped_connection() {
        disconnect();
    }

    /**
     * constructs a new scoped connection based on an already existing connection
     * @param connection the already existing connection
     */
    explicit scoped_connection(connection&& connection) : conn(std::move(connection)) {}

    /**
     * disconnect
     */
    void disconnect() {
        conn.disconnect();
    }

    /**
     * no copy
     */
    scoped_connection(const scoped_connection &) = delete;

    /**
     * no copy
     */
    scoped_connection &operator=(const scoped_connection &) = delete;

    /**
     * the usual move
     * @param rhs
     */
    scoped_connection(scoped_connection &&rhs) noexcept : conn(std::move(rhs.conn)) {}

    /**
     * the usual move
     */
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
        /**
         * constructs a callable that will erase something of a linked list once called
         * @param it that something
         */
        explicit eraser(typename std::list<T>::iterator it) : it(it) {}

        /**
         * erase it
         * @param ptr pointer to the slots
         */
        void operator()(std::shared_ptr<void> ptr) {
            auto slots = std::static_pointer_cast<std::list<T>>(ptr);
            slots->erase(it);
        }

    private:
        typename std::list<T>::iterator it;
    };

public:

    /**
     * constructs a new slots container
     */
    slots_container() : container(std::make_shared<std::list<T>>()) {}

    /**
     * create a connection to a slot
     * @param slot the slot
     * @return a connection
     */
    connection connect(const T &slot) {
        auto it = container->insert(container->end(), slot);
        return connection(container, eraser(it));
    }

    /**
     * create a scoped connection to a slot
     * @param slot the slot
     * @return a scoped connection
     */
    scoped_connection connect_scoped(const T &slot) {
        return scoped_connection(connect(slot));
    }

    /**
     * yields the number of slots
     * @return n slots
     */
    const std::size_t n_slots() const {
        return container->size();
    }

    /**
     * begin of this container
     * @return iterator to the begin
     */
    typename std::list<T>::iterator begin() {
        return container->begin();
    }

    /**
     * end of this container
     * @return iterator to the end
     */
    typename std::list<T>::iterator end() {
        return container->end();
    }

};

template<typename F>
class signal;

template<typename... Args>
class signal<void(Args...)> : public slots_container<std::function<void(Args...)>> {

public:
    /**
     * the slot type to this signal
     */
    using slot_type = std::function<void(Args...)>;

    /**
     * fires a signal to all slots
     * @param args arguments for the slots
     */
    void fire_signal(Args... args) {
        for(auto &slot : *this) {
            slot(std::forward<Args>(args)...);
        }
    }

    /**
     * fires a signal to all slots
     * @param args arguments for the slots
     */
    void operator()(Args... args) {
        fire_signal(std::forward<Args>(args)...);
    }
};

NAMESPACE_END(signals)
NAMESPACE_END(readdy)
