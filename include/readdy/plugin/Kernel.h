//
// Created by clonker on 07.03.16.
//

#ifndef READDY2_MAIN_KERNEL_H
#define READDY2_MAIN_KERNEL_H

#include <map>
#include <iostream>
#include <boost/log/trivial.hpp>
#include <boost/filesystem.hpp>
#include <readdy/plugin/Plugin.h>


namespace readdy {
    namespace plugin {
        class Kernel : public Plugin {
        protected:
            const std::string name;
        public:
            Kernel(const std::string name);
            virtual ~Kernel();

            virtual const std::string& getName() const override;
        };

        class KernelProvider : public PluginProvider<Kernel> {
        protected:
            // cannot instantiate directly
            KernelProvider();

            ~KernelProvider() {
                std::cout << "destroying kernel provider" << std::endl;
            }

            bool isSharedLibrary(const boost::filesystem::path &path);

        public:
            static KernelProvider &getInstance();

            void loadKernelsFromDirectory(const std::string& directory);

            template<class D>
            void addAs(D &&kernel) {
                const std::string name = kernel.getName();
                std::shared_ptr<Kernel> shared = std::make_shared<Kernel>(std::move(kernel));
                PluginProvider::add(name, std::move(shared));
            }

            void add(const boost::filesystem::path &sharedLib);
            void add(Kernel &k);

            const std::string getDefaultKernelDirectory();

        private:
            // prevent that copies can be created
            KernelProvider(KernelProvider const &) = delete;
            void operator=(KernelProvider const &) = delete;
        };

    }
}

#endif //READDY2_MAIN_KERNEL_H
