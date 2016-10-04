/**
 * This file contains the declaration of the base class of all programs. They
 * have a name and potentially a templatized ProgramName struct.
 *
 * @file Program.h
 * @brief Declaration of the program base class.
 * @author clonker
 * @date 08.04.16
 */

#ifndef READDY_MAIN_PROGRAM_H
#define READDY_MAIN_PROGRAM_H

#include <memory>

#include <boost/predef.h>

#if BOOST_OS_MACOS
#include <string>
#endif

namespace readdy {
namespace model {
namespace programs {

class Program {
public:
    Program() {};

    virtual void execute() = 0;

    virtual ~Program() = default;

};

namespace _internal {

template<typename T>
struct ProgramName {
    static const std::string value;
};
}

}
}
}


#endif //READDY_MAIN_PROGRAM_H
