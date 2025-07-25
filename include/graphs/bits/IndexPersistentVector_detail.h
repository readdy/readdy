//
// Created by mho on 11/6/19.
//

#pragma once

#include <vector>
#include <stack>
#include <algorithm>
#include <fmt/format.h>

namespace graphs {

struct PersistentIndex {
    std::size_t value;

    bool operator==(const PersistentIndex& rvi) const { return value == rvi.value; }
    bool operator!=(const PersistentIndex& rvi) const { return !(*this == rvi);}
    bool operator<(const PersistentIndex& rvi) const { return value < rvi.value; }
    bool operator >=(const PersistentIndex& rvi) const { return !(*this < rvi); }
    bool operator >(const PersistentIndex& rvi) const { return value > rvi.value; }
    bool operator <= (const PersistentIndex& rvi) const { return !(*this > rvi); }
};

}

namespace fmt {
template<>
struct formatter<graphs::PersistentIndex> {
    template <typename ParseContext>
    constexpr auto parse(ParseContext &ctx) { return ctx.begin(); }

    template <typename FormatContext>
    auto format(const graphs::PersistentIndex &v, FormatContext &ctx) const {
        return format_to(ctx.out(), "PersistentIndex[{}]", v.value);
    }
};
}

namespace graphs {

namespace detail {

template<typename T, typename = void>
struct can_be_deactivated : std::false_type {
};
template<typename T>
struct can_be_deactivated<T, std::void_t<decltype(std::declval<T>().deactivate())>> : std::true_type {
};
template<typename T, typename = void>
struct has_to_persistent : std::false_type {
};
template<typename T>
struct has_to_persistent<T, std::void_t<decltype(std::declval<T>().to_persistent())>> : std::true_type {
};
template<typename T, typename = void>
struct can_query_active : std::false_type {
};
template<typename T>
struct can_query_active<T, std::void_t<decltype(std::declval<T>().deactivated())>> : std::true_type {
};

template<typename IteratorL, typename IteratorR, typename std::enable_if<has_to_persistent<IteratorL>{} && has_to_persistent<IteratorR>{}, bool>::type = true>
bool operator==(const IteratorL &lhs, const IteratorR &rhs) { return lhs.to_persistent() == rhs.to_persistent(); }

template<typename IteratorL, typename IteratorR, typename std::enable_if<has_to_persistent<IteratorL>{} && has_to_persistent<IteratorR>{}, bool>::type = true>
bool operator!=(const IteratorL &lhs, const IteratorR &rhs) { return lhs.to_persistent() != rhs.to_persistent(); }

template<typename IteratorL, typename IteratorR, typename std::enable_if<has_to_persistent<IteratorL>{} && has_to_persistent<IteratorR>{}, bool>::type = true>
bool operator<(const IteratorL &lhs, const IteratorR &rhs) { return lhs.to_persistent() < rhs.to_persistent(); }

template<typename IteratorL, typename IteratorR, typename std::enable_if<has_to_persistent<IteratorL>{} && has_to_persistent<IteratorR>{}, bool>::type = true>
bool operator>(const IteratorL &lhs, const IteratorR &rhs) { return lhs.to_persistent() > rhs.to_persistent(); }

template<typename IteratorL, typename IteratorR, typename std::enable_if<has_to_persistent<IteratorL>{} && has_to_persistent<IteratorR>{}, bool>::type = true>
bool operator<=(const IteratorL &lhs, const IteratorR &rhs) { return lhs.to_persistent() <= rhs.to_persistent(); }

template<typename IteratorL, typename IteratorR, typename std::enable_if<has_to_persistent<IteratorL>{} && has_to_persistent<IteratorR>{}, bool>::type = true>
bool operator>=(const IteratorL &lhs, const IteratorR &rhs) { return lhs.to_persistent() >= rhs.to_persistent(); }

template<template<typename...> class BackingVector, typename T, typename... Rest>
class IndexPersistentContainer {
    static_assert(detail::can_be_deactivated<T>::value && detail::can_query_active<T>::value,
                  "IndexPersistentVector can only work with (ptr) element types which have a deactivate() method");
public:
    /**
     * the size type of this, inherited from the backing vector
     */
    using size_type = typename BackingVector<T, Rest...>::size_type;

    using persistent_index_t = PersistentIndex;

    /**
     * stack of blanks (indices) type
     */
    using BlanksList = std::vector<persistent_index_t>;

    /**
     * the difference type of this, inherited from the backing vector
     */
    using difference_type = typename BackingVector<T, Rest...>::difference_type;
    /**
     * the allocator type of this, inherited from the backing vector
     */
    using allocator_type = typename BackingVector<T, Rest...>::allocator_type;
    /**
     * the value type of this, inherited from the backing vector
     */
    using value_type = typename BackingVector<T, Rest...>::value_type;

    /**
     * the iterator type, same as backing vector's iterator
     */
    using persistent_iterator = typename BackingVector<T, Rest...>::iterator;
    /**
     * the const iterator type, same as backing vector's const iterator
     */
    using const_persistent_iterator = typename BackingVector<T, Rest...>::const_iterator;

    class const_active_iterator {
    public:
        using difference_type = typename const_persistent_iterator::difference_type;
        using value_type = typename const_persistent_iterator::value_type;
        using reference = typename const_persistent_iterator::reference;
        using pointer = typename const_persistent_iterator::pointer;
        using iterator_category = typename const_persistent_iterator::iterator_category;

        const_active_iterator() : parent(), begin(), end(), blanksPtr() {}

        const_active_iterator(const_persistent_iterator parent, const_persistent_iterator begin, const_persistent_iterator end,
                              const BlanksList *blanksPtr)
                : parent(parent), begin(begin), end(end), blanksPtr(blanksPtr) {
            skipBlanks();
        }

        const_active_iterator(const const_active_iterator &) = default;

        const_active_iterator(const_active_iterator &&) noexcept = default;

        ~const_active_iterator() = default;

        const_active_iterator &operator=(const const_active_iterator &) = default;

        const_active_iterator &operator=(const_active_iterator &&) noexcept = default;

        bool operator==(const const_active_iterator &rhs) const { return rhs.parent == parent; }

        bool operator!=(const const_active_iterator &rhs) const { return rhs.parent != parent; }

        bool operator<(const const_active_iterator &rhs) const { return parent < rhs.parent; }

        bool operator>(const const_active_iterator &rhs) const { return parent > rhs.parent; }

        bool operator<=(const const_active_iterator &rhs) const { return parent <= rhs.parent; }

        bool operator>=(const const_active_iterator &rhs) const { return parent >= rhs.parent; }

        const_active_iterator &operator++() {
            if (parent != end) {
                ++parent;
                skipBlanks();
            }
            return *this;
        }

        const_active_iterator operator++(int) {
            const_active_iterator copy(*this);
            ++(*this);
            return copy;
        }

        const_active_iterator &operator--() {
            --parent;
            while (parent >= begin && parent->deactivated()) {
                --parent;
            }
            return *this;
        }

        const_active_iterator operator--(int) {
            const_active_iterator copy(*this);
            copy--;
            return copy;
        }

        const_active_iterator &operator+=(size_type n) {
            auto pos = std::distance(begin, parent);
            auto targetPos = pos + n;
            auto it = std::lower_bound(blanksPtr->begin(), blanksPtr->end(), persistent_index_t{static_cast<std::size_t>(pos)});
            while (it != blanksPtr->end() && it->value <= targetPos) {
                ++targetPos;
                ++it;
            }
            parent += targetPos - pos;
            return *this;
        }

        const_active_iterator operator+(size_type n) const {
            const_active_iterator copy(*this);
            copy += n;
            return copy;
        }

        friend const_active_iterator operator+(size_type n, const const_active_iterator &it) {
            return it + n;
        }

        const_active_iterator &operator-=(size_type n) {
            auto pos = std::distance(begin, parent);
            auto targetPos = pos - n;
            auto it = std::lower_bound(blanksPtr->rbegin(), blanksPtr->rend(), persistent_index_t{static_cast<std::size_t>(pos)}, std::greater<>());
            while (it != blanksPtr->rend() && it->value >= targetPos) {
                --targetPos;
                ++it;
            }
            parent -= pos - targetPos;
            return *this;
        }

        const_active_iterator operator-(size_type n) const {
            auto copy = const_active_iterator(*this);
            copy -= n;
            return copy;
        }

        difference_type operator-(const_active_iterator rhs) const {
            auto dist = parent - rhs.parent;
            // find number of blanks in that range
            auto pos = std::distance(begin, parent);
            auto rhsPos = std::distance(begin, rhs.parent);
            auto nBlanksThis = std::distance(blanksPtr->begin(),
                                             std::lower_bound(blanksPtr->begin(), blanksPtr->end(), persistent_index_t{
                                                     static_cast<std::size_t>(pos)}));
            auto nBlanksThat = std::distance(blanksPtr->begin(),
                                             std::lower_bound(blanksPtr->begin(), blanksPtr->end(), persistent_index_t{
                                                     static_cast<std::size_t>(rhsPos)}));
            return dist - (nBlanksThis - nBlanksThat);
        }

        reference operator*() const {
            return *parent;
        }

        pointer operator->() const {
            return parent.operator->();
        }

        reference operator[](size_type n) const {
            return *(operator+(n));
        }

        [[nodiscard]] persistent_index_t persistent_index() const {
            return persistent_index_t{static_cast<std::size_t>(std::distance(begin, parent))};
        }

        const_persistent_iterator to_persistent() const {
            return parent;
        }

    private:
        void skipBlanks() {
            auto pos = std::distance(begin, parent);
            auto it = std::lower_bound(blanksPtr->begin(), blanksPtr->end(), persistent_index_t{static_cast<std::size_t>(pos)});
            while (it != blanksPtr->end() && parent != end && pos == it->value) {
                ++parent;
                ++pos;
                ++it;
            }
        }

        const_persistent_iterator parent;
        const_persistent_iterator begin;
        const_persistent_iterator end;
        const BlanksList *blanksPtr;
    };

    class active_iterator {
    public:
        using difference_type = typename persistent_iterator::difference_type;
        using value_type = typename persistent_iterator::value_type;
        using reference = typename persistent_iterator::reference;
        using pointer = typename persistent_iterator::pointer;
        using iterator_category = typename persistent_iterator::iterator_category;

        active_iterator() : parent(), begin(), end(), blanksPtr() {}

        active_iterator(persistent_iterator parent, persistent_iterator begin, persistent_iterator end, BlanksList *blanksPtr)
                : parent(parent), begin(begin), end(end), blanksPtr(blanksPtr) {
            skipBlanks();
        }

        reference operator*() const { return *parent; }

        pointer operator->() const { return parent.operator->(); }

        reference operator[](size_type i) const {
            return *(operator+(i));
        }

        active_iterator &operator++() {
            if (parent != end) {
                ++parent;
                skipBlanks();
            }
            return *this;
        }

        active_iterator operator++(int) {
            active_iterator copy(*this);
            ++(*this);
            return copy;
        }

        active_iterator &operator--() {
            --parent;
            while (parent >= begin && parent->deactivated()) {
                --parent;
            }
            return *this;
        }

        active_iterator operator--(int) {
            active_iterator copy(*this);
            --(*this);
            return copy;
        }

        active_iterator &operator+=(size_type n) {
            auto pos = std::distance(begin, parent);
            auto targetPos = pos + n;
            auto it = std::lower_bound(blanksPtr->begin(), blanksPtr->end(), persistent_index_t{static_cast<std::size_t>(pos)});
            while (it != blanksPtr->end() && it->value <= targetPos) {
                ++targetPos;
                ++it;
            }
            parent += targetPos - pos;
            return *this;
        }

        active_iterator operator+(size_type n) const {
            active_iterator copy(*this);
            copy += n;
            return copy;
        }

        friend active_iterator operator+(size_type n, const active_iterator &it) {
            return it + n;
        }

        active_iterator &operator-=(size_type n) {
            auto pos = std::distance(begin, parent);
            auto targetPos = pos - n;
            auto it = std::lower_bound(blanksPtr->rbegin(), blanksPtr->rend(), persistent_index_t{static_cast<std::size_t>(pos)}, std::greater<>());
            while (it != blanksPtr->rend() && it->value >= targetPos) {
                --targetPos;
                ++it;
            }
            parent -= pos - targetPos;
            return *this;
        }

        active_iterator operator-(size_type n) const {
            active_iterator copy(*this);
            copy -= n;
            return copy;
        }

        difference_type operator-(active_iterator rhs) const {
            auto dist = parent - rhs.parent;
            // find number of blanks in that range
            auto pos = std::distance(begin, parent);
            auto rhsPos = std::distance(begin, rhs.parent);
            auto nBlanksThis = std::distance(blanksPtr->begin(),
                                             std::lower_bound(blanksPtr->begin(), blanksPtr->end(), persistent_index_t{
                                                     static_cast<std::size_t>(pos)}));
            auto nBlanksThat = std::distance(blanksPtr->begin(),
                                             std::lower_bound(blanksPtr->begin(), blanksPtr->end(), persistent_index_t{
                                                     static_cast<std::size_t>(rhsPos)}));
            return dist - (nBlanksThis - nBlanksThat);
        }

        persistent_iterator to_persistent() const {
            return parent;
        }

        [[nodiscard]] persistent_index_t persistent_index() const {
            return persistent_index_t{static_cast<size_type>(std::distance(begin, parent))};
        }

        operator const_active_iterator() const {
            return {parent, begin, end, blanksPtr};
        }

    private:
        void skipBlanks() {
            auto pos = std::distance(begin, parent);
            auto it = std::lower_bound(blanksPtr->begin(), blanksPtr->end(), persistent_index_t{static_cast<std::size_t>(pos)});
            while (it != blanksPtr->end() && parent != end && pos == it->value) {
                ++parent;
                ++pos;
                ++it;
            }
        }

        persistent_iterator parent;
        persistent_iterator begin;
        persistent_iterator end;
        const BlanksList *blanksPtr;
    };

    using iterator = active_iterator;
    using const_iterator = const_active_iterator;

    /**
     * gives access to the backing vector
     * @return a reference to the backing vector
     */
    value_type *data_persistent() {
        return _backingVector.data();
    };

    /**
     * gives const access to the backing vector
     * @return a const reference to the backing vector
     */
    const value_type *data_persistent() const {
        return _backingVector.data();
    }

    /**_backingVector.at(index)
     * the size of this container, including blanks
     * @return the size
     */
    [[nodiscard]] size_type size_persistent() const {
        return _backingVector.size();
    }

    /**
     * returns whether the container is empty
     * @return true if it is empty, false otherwise
     */
    [[nodiscard]] bool empty_persistent() const {
        return _backingVector.empty();
    }

    [[nodiscard]] bool empty() const {
        return _backingVector.size() == _blanks.size();
    }

    /**
     * clears this container
     */
    void clear() {
        _backingVector.clear();
        _blanks.clear();
    }

    /**
     * Performs a push_back. If the blanks stack is empty, the element is simply pushed back to the backing vector,
     * otherwise it is inserted at the index the stack's first element is pointing to, which then is erased.
     * @param val the value to insert
     * @return an iterator pointing to the inserted element
     */
    iterator push_back(T &&val) {
        if (_blanks.empty()) {
            _backingVector.push_back(std::forward<T>(val));
            return iterator(std::prev(_backingVector.end()), _backingVector.begin(), _backingVector.end(), &_blanks);
        } else {
            const auto idx = _blanks.back();
            _blanks.pop_back();
            _backingVector.at(idx.value) = std::move(val);
            return {_backingVector.begin() + idx.value, std::begin(_backingVector), std::end(_backingVector),
                    &_blanks};
        }
    }

    /**
     * Performs a push_back. Same as the implementation with an r-value.
     * @param val the value to insert
     * @return an iterator pointing to the inserted element
     */
    iterator push_back(const T &val) {
        if (_blanks.empty()) {
            _backingVector.push_back(val);
            return {std::prev(_backingVector.end()), _backingVector.begin(), _backingVector.end(), &_blanks};
        } else {
            const auto idx = _blanks.back();
            _blanks.pop_back();
            _backingVector.at(idx.value) = val;
            return {_backingVector.begin() + idx.value, _backingVector.begin(), _backingVector.end(), &_blanks};
        }
    }

    /**
     * Performs an emplace_back. Works in the same way as push_back.
     * @tparam Args argument types
     * @param args arguments
     * @return an iterator pointed to the emplaced element
     */
    template<typename... Args>
    PersistentIndex emplace_back(Args &&... args) {
        if (_blanks.empty()) {
            _backingVector.emplace_back(std::forward<Args>(args)...);
            return {_backingVector.size() - 1};
        } else {
            const auto idx = _blanks.back();
            _blanks.pop_back();
            auto alloc = _backingVector.get_allocator();
            std::allocator_traits<decltype(alloc)>::construct(
                    alloc, &*_backingVector.begin() + idx.value, std::forward<Args>(args)...);
            return {idx};
        }
    }

    void erase(iterator pos) {
        deactivate(pos);
        insertBlank(pos.persistent_index());
    }

    void erase(persistent_iterator pos) {
        pos->deactivate();
        insertBlank(persistent_index_t{static_cast<std::size_t>(std::distance(_backingVector.begin(), pos))});
    }

    /**
     * Erases a range of elements.
     * @param start begin of the range, inclusive
     * @param end end of the range, exclusive
     */
    void erase(persistent_iterator start, const_persistent_iterator end) {
        auto offset = std::distance(_backingVector.begin(), start);
        for (auto it = start; it != end; ++it, ++offset) {
            deactivate(it);
            insertBlank(offset);
        }
    }

    void erase(iterator start, const_iterator end) {
        for (auto it = start; it != end; ++it) {
            deactivate(it);
            insertBlank(it.persistent_index());
        }
    }

    /**
     * Yields the number of deactivated elements, i.e., size() - n_deactivated() is the effective size of this
     * container.
     * @return the number of deactivated elements
     */
    [[nodiscard]] typename BlanksList::size_type n_deactivated() const {
        return _blanks.size();
    }

    /**
     * The effective size of this container
     * @return the effective size
     */
    [[nodiscard]] size_type size() const {
        return _backingVector.size() - n_deactivated();
    }

    /**
     * Yields a reference to the element at the requested index.
     * @param index the index
     * @return a reference to the element
     */
    T &at(size_type index) {
        return at((begin() + index).persistent_index());
    }

    /**
     * Yields a const reference to the element at the requested index.
     * @param index the index
     * @return a const reference to the element
     */
    const T &at(size_type index) const {
        return at((begin() + index).persistent_index());
    }

    T &at(PersistentIndex index) {
        auto &x = *(begin_persistent() + index.value);
        if (x.deactivated()) {
            throw std::invalid_argument(fmt::format("Requested deactivated element {}", index));
        }
        return x;
    }

    const T &at(PersistentIndex index) const {
        const auto &x = *(begin_persistent() + index.value);
        if (x.deactivated()) {
            throw std::invalid_argument(fmt::format("Requested deactivated element {}", index));
        }
        return x;
    }

    /**
     * Yields an iterator pointing to the begin of this container.
     * @return the iterator
     */
    persistent_iterator begin_persistent() noexcept {
        return _backingVector.begin();
    }

    iterator begin() noexcept {
        return {_backingVector.begin(), _backingVector.begin(), _backingVector.end(), &_blanks};
    }

    const_iterator begin() const noexcept {
        return cbegin();
    }

    const_iterator cbegin() const noexcept {
        return {_backingVector.cbegin(), _backingVector.cbegin(), _backingVector.cend(), &_blanks};
    }

    /**
     * Yields an iterator pointing to the end of this container.
     * @return the iterator
     */
    persistent_iterator end_persistent() noexcept {
        return _backingVector.end();
    }

    iterator end() noexcept {
        return {_backingVector.end(), _backingVector.begin(), _backingVector.end(), &_blanks};
    }

    const_iterator end() const noexcept {
        return cend();
    }

    const_iterator cend() const noexcept {
        return {_backingVector.end(), _backingVector.begin(), _backingVector.end(), &_blanks};
    }

    /**
     * Yields a const iterator pointing to the begin of this container.
     * @return the iterator
     */
    const_persistent_iterator begin_persistent() const noexcept {
        return _backingVector.begin();
    }

    /**
     * Yields a const iterator pointing to the begin of this container.
     * @return the iterator
     */
    const_persistent_iterator cbegin_persistent() const noexcept {
        return _backingVector.cbegin();
    }

    /**
     * Yields a const iterator pointing to the end of this container.
     * @return the iterator
     */
    const_persistent_iterator end_persistent() const noexcept {
        return _backingVector.end();
    }

    /**
     * Yields a const iterator pointing to the end of this container.
     * @return the iterator
     */
    const_persistent_iterator cend_persistent() const noexcept {
        return _backingVector.cend();
    }

    active_iterator persistent_to_active_iterator(persistent_iterator it) {
        return active_iterator(it, std::begin(_backingVector), std::end(_backingVector), &_blanks);
    }

    const_active_iterator persistent_to_active_iterator(const_persistent_iterator it) const {
        return cpersistent_to_active_iterator(it);
    }

    const_active_iterator cpersistent_to_active_iterator(const_persistent_iterator it) const {
        return const_active_iterator(it, std::begin(_backingVector), std::end(_backingVector), &_blanks);
    }

    [[nodiscard]] persistent_index_t persistentIndex(active_iterator it) const {
        return it.persistent_index();
    }

    [[nodiscard]] persistent_index_t persistentIndex(const_persistent_iterator it) const {
        auto d = std::distance(std::begin(_backingVector), it);
        if (d >= 0) {
            return {static_cast<std::size_t>(d)};
        } else {
            throw std::logic_error("Distance between begin and it was negative: d = " + std::to_string(d));
        }
    }

private:

    void deactivate(iterator it) {
        it->deactivate();
    }

    void insertBlank(typename BlanksList::value_type val) {
        auto it = std::lower_bound(_blanks.begin(), _blanks.end(), val, std::less<>());
        _blanks.insert(it, val);
    }

    BlanksList _blanks {};
    BackingVector<T, Rest...> _backingVector {};
};

}
}
