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


//
// Created by mho on 10/10/2016.
//



#ifndef READDY_MAIN_MACROS_H
#define READDY_MAIN_MACROS_H


/**
 * Defines for the current OS
 */
#ifdef _WIN32
#define READDY_WINDOWS true
#define READDY_OSX false
#define READDY_LINUX false
#elif __APPLE__
#define READDY_WINDOWS false
#define READDY_OSX true
#define READDY_LINUX false
#elif __linux__ || __unix__
#define READDY_WINDOWS false
#define READDY_OSX false
#define READDY_LINUX true
#else
#error "Unknown compiler"
#endif

/**
 * Export symbols / change their visibility
 */
#if READDY_WINDOWS
#  define READDY_API __declspec(dllexport)
#else
#  define READDY_API __attribute__ ((visibility("default")))
#endif

/**
 * Utilities
 */
#define READDY_CREATE_FACTORY_DISPATCHER(FACTORY, TYPE) template<typename... Args> \
struct FACTORY::get_dispatcher<TYPE, Args...> { \
    static TYPE *impl(const FACTORY *self, Args &&... args) { \
        return self->create##TYPE(std::forward<Args>(args)...); \
    } \
};
#define READDY_CREATE_FACTORY_DISPATCHER2(FACTORY, NS, TYPE) template<typename... Args> \
struct FACTORY::get_dispatcher<NS::TYPE, Args...> { \
    static NS::TYPE *impl(const FACTORY *self, Args &&... args) { \
        return self->create##TYPE(std::forward<Args>(args)...); \
    } \
};
#define READDY_CREATE_OBSERVABLE_FACTORY_DISPATCHER(TYPE) template<typename... Args> \
struct ObservableFactory::get_dispatcher<TYPE, Args...> { \
    static TYPE *impl(const ObservableFactory *self, unsigned int stride, Args &&... args) { \
        return self->create##TYPE(stride, std::forward<Args>(args)...); \
    } \
};
#define READDY_CREATE_COMPARTMENT_FACTORY_DISPATCHER(TYPE) template<typename... Args> \
struct CompartmentFactory::get_dispatcher<TYPE, Args...> { \
    static TYPE *impl(const CompartmentFactory *self, const std::unordered_map<particleType_t, particleType_t> &convMap, Args &&... args) { \
        return self->create##TYPE(convMap, std::forward<Args>(args)...); \
    } \
};

#if !defined(NAMESPACE_BEGIN)
#  define NAMESPACE_BEGIN(name) namespace name {
#endif
#if !defined(NAMESPACE_END)
#  define NAMESPACE_END(name) }
#endif

#endif //READDY_MAIN_MACROS_H
