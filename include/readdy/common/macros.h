//
// Created by mho on 10/10/2016.
//

#ifndef READDY_MAIN_MACROS_H
#define READDY_MAIN_MACROS_H
#ifdef _WIN32
#  define READDY_EXPORT __declspec(dllexport)
#else
#  define READDY_EXPORT __attribute__ ((visibility("default")))
#endif
#endif //READDY_MAIN_MACROS_H
