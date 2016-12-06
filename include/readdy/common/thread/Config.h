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
 * << detailed description >>
 *
 * @file Config.h
 * @brief << brief description >>
 * @author clonker
 * @date 05.09.16
 */

#ifndef READDY_MAIN_CONFIG_H
#define READDY_MAIN_CONFIG_H

namespace readdy {
namespace util {
namespace thread {

struct Config {
    using n_threads_t = decltype(std::thread::hardware_concurrency());

    Config();
    n_threads_t nThreads() const;
    void setNThreads(const n_threads_t);
private:
    n_threads_t m_nThreads;
};

}
}
}
#endif //READDY_MAIN_CONFIG_H
