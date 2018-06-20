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
 * This header file contains definitions for various common utils. Currently:
 *   - getOS: returns a string corresponding to the executing operating system
 *   - isWindows: returns true if we are in windows
 *   - collections::hasKey: convenience method to check if a map contains a certain key
 *   - testing::getPluginsDirectory: Method that checks for some environment variables and then returns a potential
 *     directory in which the kernels are likely located.
 *
 * @file ObservableFactory.h
 * @brief Header file containing some common utils.
 * @author clonker
 * @date 08.03.16
 */

#pragma once
#include <string>
#include <iostream>
#include <memory>
#include <tuple>
#include <vector>
#include <algorithm>
#include <functional>
#include <utility>
#include <readdy/common/macros.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(util)
std::string getOS();

bool isWindows();

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);

std::vector<std::string> split(const std::string &s, char delim);

NAMESPACE_BEGIN(collections)

template<typename MapType, typename KeyType = std::string>
inline bool hasKey(const MapType &map, const KeyType &key) {
    return map.find(key) != map.end();
}

template<template<class, class, class...> class C, typename K, typename V, typename... Args>
inline const V &getOrDefault(const C<K, V, Args...> &m, const K &key, const V &defaultValue) {
    typename C<K, V, Args...>::const_iterator it = m.find(key);
    if (it == m.end()) {
        return defaultValue;
    }
    return it->second;
}

template<typename Collection, typename Fun>
inline void for_each_value(const Collection& collection, Fun f)  {
    for(auto&& e : collection) {
        for(auto&& inner : e.second) {
            f(e.first, inner);
        }
    }
}

template<typename Collection, typename Fun>
inline void for_each_value_ref(const Collection& collection, Fun &f)  {
    for(auto&& e : collection) {
        for(auto&& inner : e.second) {
            f(inner);
        }
    }
}

template <typename map, typename Fun>
inline void for_each_value_in_map(map &m, Fun f) {
    for (auto&& entry : m) {
        f(entry.second);
    }
}

template<typename order_iterator, typename value_iterator>
void reorder_destructive(order_iterator order_begin, order_iterator order_end, value_iterator v) {
    typedef typename std::iterator_traits<value_iterator>::value_type value_t;
    typedef typename std::iterator_traits<order_iterator>::value_type index_t;
    typedef typename std::iterator_traits<order_iterator>::difference_type diff_t;

    diff_t remaining = order_end - 1 - order_begin;
    for (index_t s = index_t(); remaining > 0; ++s) {
        index_t d = order_begin[s];
        if (d == (diff_t) -1) continue;
        --remaining;
        value_t temp = std::move(v[s]);
        for (index_t d2; d != s; d = d2) {
            std::swap(temp, v[d]);
            std::swap(order_begin[d], d2 = (diff_t) -1);
            --remaining;
        }
        v[s] = std::move(temp);
    }
}

template<typename T, typename Predicate>
typename std::vector<T>::iterator insert_sorted(std::vector<T> &vec, T const &item, Predicate pred) {
    return vec.insert(std::upper_bound(vec.begin(), vec.end(), item, pred), item);
}

NAMESPACE_END(collections)

template<typename Derived, typename Base, typename Del>
std::unique_ptr<Derived, Del> static_unique_ptr_cast( std::unique_ptr<Base, Del>&& p ) {
    auto d = static_cast<Derived *>(p.release());
    return std::unique_ptr<Derived, Del>(d, std::move(p.get_deleter()));
}

template<typename Derived, typename Base>
std::unique_ptr<Derived> static_unique_ptr_cast_no_del( std::unique_ptr<Base>&& p ) {
    auto d = static_cast<Derived *>(p.release());
    return std::unique_ptr<Derived>(d);
}

template<typename Derived, typename Base, typename Del>
std::unique_ptr<Derived, Del> dynamic_unique_ptr_cast( std::unique_ptr<Base, Del>&& p ) {
    if(auto *result = dynamic_cast<Derived *>(p.get())) {
        p.release();
        return std::unique_ptr<Derived, Del>(result, std::move(p.get_deleter()));
    }
    return std::unique_ptr<Derived, Del>(nullptr, p.get_deleter());
}

NAMESPACE_END(util)
NAMESPACE_END(readdy)
