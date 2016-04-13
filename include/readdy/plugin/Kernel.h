/**
 * This file contains the class definitions for Kernel and KernelProvider.
 * A Kernel is used to execute Programs, i.e., instances of readdy::plugin::Program.
 * The kernels can be built in or provided by shared libs in directories, which are loaded by the KernelProvider.
 * Each Kernel has a readdy::plugin::Kernel::name by which it can be accessed in the KernelProvider.
 *
 * @file Kernel.h
 * @brief File containing the definitions of Kernel and KernelProvider.
 * @author clonker
 * @date 07.03.16
 */

#ifndef READDY2_MAIN_KERNEL_H
#define READDY2_MAIN_KERNEL_H

#include <map>
#include <iostream>
#include <boost/log/trivial.hpp>
#include <boost/filesystem.hpp>
#include <readdy/plugin/Plugin.h>
#include <boost/log/sources/logger.hpp>
#include "Program.h"


namespace readdy {
    namespace plugin {
        /**
         * Base class of kernels.
         * A Kernel is used to execute Programs, i.e., instances of readdy::plugin::Program.
         * The kernels can be built in or provided by shared libs in directories, which are loaded by the KernelProvider.
         * Each Kernel has a #name by which it can be accessed in the KernelProvider.
         */
        class Kernel : public Plugin {
        protected:
            /**
             * The name of the kernel.
             */
            std::string name;
        public:
            /**
             * Constructs a kernel with a given name.
             */
            Kernel(const std::string name);
            /**
             * The kernel destructor.
             */
            virtual ~Kernel();
            /**
             * This method returns the name of the kernel.
             *
             * @return The name
             */
            virtual const std::string &getName() const override;
            /**
             * Create a program that can be executed on this kernel.
             * If the requested program is not available on the kernel, a nullptr is returned.
             *
             * @param name the name of the program
             * @see getAvailablePrograms()
             * @return The program if it was available, otherwise nullptr
             */
            virtual std::shared_ptr<readdy::plugin::Program> createProgram(std::string name);
            /**
             * Returns a vector containing all available program names for this specific kernel instance.
             *
             * @see createProgram(name)
             * @return The program names.
             */
            virtual std::vector<std::string> getAvailablePrograms();
        };
        /**
         * The KernelProvider is a singleton which can be accessed by getInstance()
         * and provides Kernels that can be added directly or loaded from directories.
         * If loaded from directories (with loadKernelsFromDirectory(string), the
         * specified directory will be scanned for shared libraries with the required
         * symbols, i.e., with an implementation of the Kernel class.
         */
        class KernelProvider : public PluginProvider<Kernel> {
        protected:
            /**
             * The constructor of KernelProvider. As it is a singleton, it is protected.
             */
            KernelProvider();
            /**
             * The destructor of KernelProvider.
             */
            ~KernelProvider() {
                BOOST_LOG_TRIVIAL(debug) << "destroying kernel provider";
            }
            /**
             * A protected method that determines if a boost::filesystem::path points to a shared library.
             *
             * @param path the path
             * @return True if the path points to a shared library, otherwise false.
             */
            bool isSharedLibrary(const boost::filesystem::path &path);

        public:
            /**
             * Method that returns the singleton KernelProvider.
             *
             * @return The KernelProvider.
             */
            static KernelProvider &getInstance();

            /**
             * Method to load kernels (non-recursively) in shared libraries from a directory.
             *
             * @param directory the directory in which the shared libraries are located.
             */
            void loadKernelsFromDirectory(const std::string &directory);

            /**
             * @todo write documentation
             */
            template<class D>
            void addAs(D &&kernel) {
                const std::string name = kernel.getName();
                std::shared_ptr<Kernel> shared = std::make_shared<Kernel>(std::move(kernel));
                PluginProvider::add(name, std::move(shared));
            }

            /**
             * @todo write documentation
             */
            void add(const boost::filesystem::path &sharedLib);

            /**
             * @todo write documentation
             */
            void add(Kernel &k);

            /**
             * @todo write documentation
             */
            const std::string getDefaultKernelDirectory();

        private:
            // prevent that copies can be created
            KernelProvider(KernelProvider const &) = delete;

            void operator=(KernelProvider const &) = delete;

        };

    }
}

#endif //READDY2_MAIN_KERNEL_H
