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

                virtual std::shared_ptr<readdy::plugin::Program> createProgram(std::string name) override;

                virtual std::vector<std::string> getAvailablePrograms() override;

                virtual std::shared_ptr<readdy::model::KernelStateModel> getKernelStateModel() override;

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
