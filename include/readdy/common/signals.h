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
