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
 * The KernelPluginDecorator class wraps a loaded kernel instance and the corresponding shared_library
 * instance by using the decorator pattern.
 *
 * @file KernelPluginDecorator.h
 * @brief This file contains the KernelPluginDecorator class definitions.
 * @author clonker
 * @date 14.03.16
 */

#ifndef READDY_MAIN_KERNELPLUGINDECORATOR_H
#define READDY_MAIN_KERNELPLUGINDECORATOR_H

#include <readdy/model/Kernel.h>
#include <readdy/common/dll.h>

namespace readdy {
namespace plugin {
namespace _internal {

class KernelPluginDecorator : public readdy::model::Kernel {
protected:
    std::unique_ptr<readdy::model::Kernel> reference;
    std::unique_ptr<readdy::util::dll::shared_library> lib;

public:

    KernelPluginDecorator(const std::string& sharedLib);

    virtual ~KernelPluginDecorator();

    virtual readdy::model::top::TopologyActionFactory *getTopologyActionFactory() const override;

    virtual readdy::model::KernelStateModel &getKernelStateModel() const override;

    virtual readdy::model::KernelContext &getKernelContext() const override;

    virtual const std::string &getName() const override;

    virtual std::vector<std::string> getAvailablePotentials() const override;

    virtual readdy::model::actions::ActionFactory &getActionFactory() const override;

    virtual readdy::model::reactions::ReactionFactory &getReactionFactory() const override;

    virtual readdy::model::compartments::CompartmentFactory &getCompartmentFactory() const override;

    virtual readdy::signals::scoped_connection
    connectObservable(model::observables::ObservableBase *const observable) override;

    virtual void evaluateObservables(readdy::model::observables::time_step_type t) override;

    virtual std::tuple<std::unique_ptr<readdy::model::observables::ObservableWrapper>, readdy::signals::scoped_connection>
    registerObservable(const model::observables::observable_type &observable, unsigned int stride) override;

    virtual std::vector<std::string> getAvailableActions() const override;

    virtual void addParticle(const std::string &type, const model::Vec3 &pos) override;

    virtual unsigned int getTypeId(const std::string &string) const override;


protected:
    virtual readdy::model::observables::ObservableFactory &getObservableFactory() const override;

    virtual readdy::model::potentials::PotentialFactory &getPotentialFactory() const override;


};

class InvalidPluginException : public std::runtime_error {
public:
    InvalidPluginException(const std::string &__arg);
};

const std::string loadKernelName(const std::string &sharedLib);
}
}
}


#endif //READDY_MAIN_KERNELPLUGINDECORATOR_H
