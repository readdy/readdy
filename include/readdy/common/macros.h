//
// Created by clonker on 9/19/2021.
//

#pragma once

#ifdef WIN32
#    ifdef LIBRARY_EXPORTS
#        define READDY_API __declspec(dllexport)
#    else
#        define READDY_API __declspec(dllimport)
#    endif
#elif
#    define READDY_API
#endif
