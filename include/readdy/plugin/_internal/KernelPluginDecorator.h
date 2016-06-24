//
// Created by clonker on 14.03.16.
//

#ifndef READDY_MAIN_KERNELPLUGINDECORATOR_H
#define READDY_MAIN_KERNELPLUGINDECORATOR_H


#include <readdy/model/Kernel.h>
#include <boost/dll/shared_library.hpp>

namespace readdy {
    namespace plugin {
        namespace _internal {
            class KernelPluginDecorator : public readdy::model::Kernel {
            protected:
                std::unique_ptr<readdy::model::Kernel> reference;
                boost::dll::shared_library lib;

            public:
                KernelPluginDecorator(const boost::filesystem::path sharedLib);
                virtual ~KernelPluginDecorator();

                virtual readdy::model::KernelStateModel& getKernelStateModel() const override;

                virtual readdy::model::KernelContext& getKernelContext() const override;

                virtual const std::string &getName() const override;

                virtual std::unique_ptr<readdy::model::potentials::Potential> createPotential(std::string &name) const override;

                virtual std::vector<std::string> getAvailablePotentials() const override;

                virtual readdy::model::programs::ProgramFactory &getProgramFactory() const override;

                virtual readdy::model::reactions::ReactionFactory &getReactionFactory() const override;

            protected:
                virtual readdy::model::potentials::PotentialFactory &getPotentialFactory() const override;


            };

            class InvalidPluginException : public std::runtime_error {
            public:
                InvalidPluginException(const std::string &__arg);
            };

            const std::string loadKernelName(const boost::filesystem::path& sharedLib);
        }
    }
}


#endif //READDY_MAIN_KERNELPLUGINDECORATOR_H
