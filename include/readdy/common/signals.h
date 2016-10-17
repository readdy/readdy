/**
 * Strongly simplified signals implementation similar to the one that can be found in, e.g., boost.signals2.
 *
 * @file signals.h
 * @brief A signals implementation
 * @author clonker
 * @date 14.10.16
 */

#ifndef READDY_MAIN_SIGNALS_H
#define READDY_MAIN_SIGNALS_H

#include <functional>
#include <list>
#include <memory>
#include <algorithm>
#include "make_unique.h"
#include "logging.h"

namespace readdy {
namespace signals {

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
        if(container) {
            disconnector(container);
            slots_container.reset();
        }
    }
private:

    template<typename T>
    friend class slots_container;

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
    scoped_connection(connection connection) : connection(std::move(connection)) {}
    void disconnect() {
        connection.disconnect();
    }

private:
    connection connection;
};

template<typename T>
class slots_container {
    std::shared_ptr<std::list<T>> container;
    struct eraser {
        eraser(typename std::list<T>::iterator it) : it(it) { }
        void operator()(std::shared_ptr<void> ptr) {
            auto slots = std::static_pointer_cast<std::list<T>>(ptr);
            slots->erase(it);
        }
    private:
        typename std::list<T>::iterator it;
    };
public:

    slots_container() : container(std::make_shared<std::list<T>>()) { }

    connection connect(const T& slot) {
        auto it = container->insert(container->end(), slot);
        return connection(container, eraser(it));
    }

    scoped_connection connect_scoped(const T& slot) {
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
        std::for_each(this->begin(), this->end(), [&](slot_type& slot) { slot(std::forward<Args>(args)...); });
    }

    void operator()(Args... args) {
        fire_signal(std::forward<Args>(args)...);
    }
};

}
}

#endif //READDY_MAIN_SIGNALS_H
