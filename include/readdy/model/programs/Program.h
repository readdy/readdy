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

#if READDY_OSX
#include <string>
#endif

namespace readdy {
namespace model {
namespace programs {

class Program {
public:
    virtual ~Program() = default;

    virtual void execute() = 0;
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
