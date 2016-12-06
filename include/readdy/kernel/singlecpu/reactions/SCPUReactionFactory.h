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
 * @file SingleCPUReactionFactory.h
 * @brief << brief description >>
 * @author clonker
 * @date 21.06.16
 */

#include <readdy/model/reactions/ReactionFactory.h>

#ifndef READDY_MAIN_SINGLECPUREACTIONFACTORY_H
#define READDY_MAIN_SINGLECPUREACTIONFACTORY_H
namespace readdy {
namespace kernel {
namespace scpu {

class SCPUKernel;
namespace reactions {

class SCPUReactionFactory : public readdy::model::reactions::ReactionFactory {

public:
    SCPUReactionFactory(SCPUKernel const *const kernel);

protected:
    virtual readdy::model::reactions::Conversion *
    createConversion(const std::string &name, unsigned int from, unsigned int to,
                     const double rate) const override;

    virtual readdy::model::reactions::Enzymatic *
    createEnzymatic(const std::string &name, unsigned int catalyst, unsigned int from, unsigned int to,
                    const double rate, const double eductDistance) const override;

    virtual readdy::model::reactions::Fission *
    createFission(const std::string &name, unsigned int from, unsigned int to1, unsigned int to2,
                  const double rate, const double productDistance, const double weight1 = 0.5,
                  const double weight2 = 0.5) const override;

    virtual readdy::model::reactions::Fusion *
    createFusion(const std::string &name, unsigned int from1, unsigned int from2, unsigned int to,
                 const double rate, const double eductDistance, const double weight1 = 0.5,
                 const double weight2 = 0.5) const override;

    SCPUKernel const *const kernel;
};

}
}
}
}
#endif //READDY_MAIN_SINGLECPUREACTIONFACTORY_H
