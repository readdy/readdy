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
#include <readdy/common/common.h>

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
