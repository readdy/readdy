//
// Created by clonker on 14.03.16.
//

#ifndef READDY2_MAIN_KERNELPLUGINDECORATOR_H
#define READDY2_MAIN_KERNELPLUGINDECORATOR_H


#include <readdy/plugin/Kernel.h>
#include <boost/dll/shared_library.hpp>

namespace readdy {
    namespace plugin {
        namespace _internal {
            class KernelPluginDecorator : public readdy::plugin::Kernel {
            protected:
                std::shared_ptr<readdy::plugin::Kernel> reference;
                boost::dll::shared_library lib;

            public:
                KernelPluginDecorator(const boost::filesystem::path sharedLib);
                virtual ~KernelPluginDecorator();

                virtual std::unique_ptr<readdy::plugin::Program> createProgram(const std::string& name) const override;

                virtual std::vector<std::string> getAvailablePrograms() const override;

                virtual readdy::model::KernelStateModel& getKernelStateModel() const override;

                virtual readdy::model::KernelContext& getKernelContext() const override;

                virtual const std::string &getName() const override;
            };

            class InvalidPluginException : public std::runtime_error {
            public:
                InvalidPluginException(const std::string &__arg);
            };
        }
    }
}


#endif //READDY2_MAIN_KERNELPLUGINDECORATOR_H
