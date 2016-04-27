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
#include <readdy/plugin/Observables.h>
#include <readdy/plugin/Program.h>
#include <readdy/model/KernelStateModel.h>
#include <readdy/model/KernelContext.h>
#include <boost/signals2/signal.hpp>
#include "Observable.h"


namespace readdy {
    namespace plugin {
        /**
         * Base class of kernels.
         * A Kernel is used to execute Programs, i.e., instances of readdy::plugin::Program.
         * The kernels can be built in or provided by shared libs in directories, which are loaded by the KernelProvider.
         * Each Kernel has a #name by which it can be accessed in the KernelProvider.
         */
        class Kernel : public Plugin {
        typedef boost::signals2::signal<void(const std::shared_ptr<readdy::model::KernelContext>, const std::shared_ptr<readdy::model::KernelStateModel>)> signal_t;
        protected:
            /**
             * The name of the kernel.
             */
            std::string name;
            /**
             * todo
             */
            signal_t signal;
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

            typedef signal_t::slot_type ObservableType;

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
             * Get a vector of the registered predefined observable names, which can be created by createObservable(name).
             *
             * @return the vector of available observable names
             */
            virtual std::vector<std::string> getAvailableObservables();
            /**
             * Creates an observable that is already available as part of the kernel implementation. The list of observables can be obtained by getAvailableObservables().
             *
             * @return a shared pointer to the created observable
             */
            virtual std::shared_ptr<Observable> createObservable(std::string name);

            /**
             * Registers an observable to the kernel signal.
             */
            boost::signals2::connection registerObservable(const Observable &observable);
            /**
             * Registers an observable to the kernel signal.
             */
            boost::signals2::connection registerObservable(const ObservableType &observable, unsigned int stride);
            /**
             * @todo document this
             */
            template<class T>
            std::shared_ptr<T> createProgramAs(std::string name) {
                return std::dynamic_pointer_cast<T>(createProgram(name));
            }

            /**
             * Returns a vector containing all available program names for this specific kernel instance.
             *
             * @see createProgram(name)
             * @return The program names.
             */
            virtual std::vector<std::string> getAvailablePrograms();

            /**
             * @todo implement this properly
             */
            virtual std::shared_ptr<readdy::model::KernelStateModel> getKernelStateModel();

            /**
             * @todo implement & document this properly
             */
            virtual std::shared_ptr<readdy::model::KernelContext> getKernelContext();
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
             * Method that allows to move a kernel into the KernelProvider and thus make it available.
             *
             * @param kernel the kernel that should be moved
             */
            void add(const Kernel &&kernel);

            /**
             * Method that allows to add a kernel to the KernelProvider by providing a path to a shared lib (containing an implementation of a kernel).
             *
             * @param sharedLib the path to the shared lib
             */
            void add(const boost::filesystem::path &sharedLib);

            /**
             * Method that gives the default kernel directory, i.e., where the kernel implementations are usually to be found.
             * First it is checked, if the environment variable 'READDY_PLUGIN_DIR' is set. In that case, the default kernel directory is the contents
             * of that environment variable.
             * Otherwise, the default kernel directory on unix systems is
             * \code{.unparsed}
             * /usr/local/readdy/lib/readdy_plugins
             * \endcode
             * and the default kernel directory on windows systems is
             * \code{.unparsed}
             * C:\\Program Files\ReaDDy2\lib\readdy_plugins
             * \endcode
             *
             * @return the default kernel directory.
             */
            static const std::string getDefaultKernelDirectory();

        private:
            // prevent that copies can be created
            KernelProvider(KernelProvider const &) = delete;

            // prevent that copies can be created
            void operator=(KernelProvider const &) = delete;

        };

    }
}

#endif //READDY2_MAIN_KERNEL_H
