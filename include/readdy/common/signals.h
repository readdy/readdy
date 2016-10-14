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

namespace readdy {
namespace signals {

class connection {
public:
    using slots_container_type = std::shared_ptr<void>;
    using disconnector_type = std::function<void(slots_container_type)>;

    virtual ~connection() = default;

    bool connected() const {
        return !slots_container.expired();
    }

    void disconnect() {
        auto container = slots_container.lock();
        if(container) {

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

class scoped_connection : public connection {
public:
    virtual ~scoped_connection() {
        disconnect();
    }

private:

    template<typename T>
    friend class slots_container;

    scoped_connection(const connection::slots_container_type &slots_container,
                      const connection::disconnector_type &disconnector)
            : connection(slots_container, disconnector) {}

};

template<typename T>
class slots_container {
public:

    slots_container() : container(std::make_shared<std::list<T>>()) { }

    connection connect(const T& slot) {
        auto it = container->insert(container->end(), slot);
        return connection(container, [it](connection::disconnector_type::argument_type ptr) {
            std::static_pointer_cast<std::list<T>>(ptr)->erase(it);
        });
    }

    scoped_connection connect_scoped(const T& slot) {
        auto it = container->insert(container->end(), slot);
        return scoped_connection(container, [it](connection::disconnector_type::argument_type ptr) {
            std::static_pointer_cast<std::list<T>>(ptr)->erase(it);
        });
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

private:
    std::shared_ptr<std::list<T>> container;
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
