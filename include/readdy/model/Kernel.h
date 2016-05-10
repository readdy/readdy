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
#include <readdy/model/Plugin.h>
#include <boost/log/sources/logger.hpp>
#include <readdy/model/Program.h>
#include <readdy/model/KernelStateModel.h>
#include <readdy/model/KernelContext.h>
#include <readdy/model/_internal/ObservableFactory.h>
#include <boost/signals2/signal.hpp>
#include <boost/signals2/shared_connection_block.hpp>
#include <readdy/model/_internal/ObservableWrapper.h>

namespace readdy {
    namespace model {
        /**
         * Base class of kernels.
         * A Kernel is used to execute Programs, i.e., instances of readdy::plugin::Program.
         * The kernels can be built in or provided by shared libs in directories, which are loaded by the KernelProvider.
         * Each Kernel has a #name by which it can be accessed in the KernelProvider.
         */
        class Kernel : public Plugin {
        public:
            /**
             * Constructs a kernel with a given name.
             */
            Kernel(const std::string &name);

            /**
             * The kernel destructor.
             */
            virtual ~Kernel();

            Kernel(const Kernel &rhs) = delete;

            Kernel &operator=(const Kernel &rhs) = delete;

            Kernel(Kernel &&rhs);

            Kernel &operator=(Kernel &&rhs);

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
            virtual std::unique_ptr<Program> createProgram(const std::string &name) const;

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
            virtual std::unique_ptr<ObservableBase> createObservable(const std::string &name);

            template<typename T>
            std::unique_ptr<T> createObservable() {
                return getObservableFactory().create<T>();
            }

            template<typename T, typename Obs1, typename Obs2>
            std::unique_ptr<T> createObservable(Obs1 *obs1, Obs2 *obs2, unsigned int stride = 1) {
                return getObservableFactory().create<T>(obs1, obs2, stride);
            };

            // todo test this
            template<typename T>
            std::tuple<std::unique_ptr<T>, boost::signals2::connection> createAndRegisterObservable(unsigned int stride) {
                auto &&obs = createObservable<T>();
                obs->setStride(stride);
                auto &&connection = registerObservable(obs.get());
                return std::make_tuple(std::move(obs), connection);
            };

            // todo test this
            std::tuple<std::unique_ptr<ObservableBase>, boost::signals2::connection> createAndRegisterObservable(const std::string &name, unsigned int stride);

            /**
             * Registers an observable to the kernel signal.
             */
            boost::signals2::connection registerObservable(ObservableBase *const observable);

            //todo
            void evaluateObservablesAutomatically(bool evaluate);

            void evaluateObservables();

            void evaluateAllObservables();

            /**
             * Registers an observable to the kernel signal.
             */
            std::tuple<boost::signals2::connection, ObservableWrapper> registerObservable(const ObservableType &observable, unsigned int stride);

            /**
             * Creates a specified program and returns a pointer to it in the templated type.
             *
             * @return a pointer to the created program.
             */
            template<class T>
            std::unique_ptr<T> createProgramAs(std::string name) {
                return std::dynamic_pointer_cast<T>(createProgram(name));
            }

            /**
             * Returns a vector containing all available program names for this specific kernel instance.
             *
             * @see createProgram(name)
             * @return The program names.
             */
            virtual std::vector<std::string> getAvailablePrograms() const;

            /**
             * @todo implement this properly
             */
            virtual readdy::model::KernelStateModel &getKernelStateModel() const;

            /**
             * @todo implement & document this properly
             */
            virtual readdy::model::KernelContext &getKernelContext() const;

        protected:
            struct Impl;
            std::unique_ptr<Impl> pimpl;

            _internal::ObservableFactory &getObservableFactory() const;
        };


    }
}

#endif //READDY2_MAIN_KERNEL_H
