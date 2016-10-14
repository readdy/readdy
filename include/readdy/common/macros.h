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
#  define READDY_EXPORT __declspec(dllexport)
#else
#  define READDY_EXPORT __attribute__ ((visibility("default")))
#endif

#endif //READDY_MAIN_MACROS_H
