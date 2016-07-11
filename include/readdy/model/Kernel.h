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

#ifndef READDY_MAIN_KERNEL_H
#define READDY_MAIN_KERNEL_H

#include <map> 
#include <iostream>
#include <boost/filesystem.hpp>
#include <readdy/model/Plugin.h>
#include <boost/log/sources/logger.hpp>
#include <readdy/model/programs/Program.h>
#include <readdy/model/KernelStateModel.h>
#include <readdy/model/KernelContext.h>
#include <readdy/model/_internal/ObservableFactory.h>
#include <boost/signals2/signal.hpp>
#include <boost/signals2/shared_connection_block.hpp>
#include <readdy/model/_internal/ObservableWrapper.h>
#include <readdy/model/potentials/PotentialFactory.h>
#include <readdy/model/programs/ProgramFactory.h>
#include <readdy/model/reactions/ReactionFactory.h>

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
            virtual std::unique_ptr<programs::Program> createProgram(const std::string &name) const {
                return getProgramFactory().createProgram(name);
            };


            template<typename ProgramType>
            std::unique_ptr<ProgramType> createProgram() const {
                return getProgramFactory().createProgramAs<ProgramType>(programs::getProgramName<ProgramType>());
            }

            template<typename T, typename... Args>
            std::unique_ptr<T> createObservable(unsigned int stride, Args... args) {
                return getObservableFactory().create<T>(stride, std::forward<Args>(args)...);
            }

            template<typename T, typename Obs1, typename Obs2>
            std::unique_ptr<T> createObservable(Obs1 *obs1, Obs2 *obs2, unsigned int stride = 1) {
                return getObservableFactory().create<T>(obs1, obs2, stride);
            };

            /**
             * Creates an observable and connects it to the kernel.
             *
             * @return a tuple of the created observable and a scoped_connection object.
             */
            template<typename T>
            std::tuple<std::unique_ptr<T>, boost::signals2::scoped_connection> createAndConnectObservable(unsigned int stride) {
                auto &&obs = createObservable<T>(stride);
                auto &&connection = connectObservable(obs.get());
                return std::make_tuple(std::move(obs), connection);
            };

            /**
             * Connects an observable to the kernel signal.
             *
             * @return A scoped_connection object that, once deleted, releases the connection of the observable.
             */
            boost::signals2::scoped_connection connectObservable(ObservableBase *const observable);

            /**
             * If set to true, all (for the current timestep unblocked) observables
             * will be evaluated once the model gets advanced in time.
             */
            void evaluateObservablesAutomatically(bool evaluate);

            /**
             * Evaluates all unblocked observables.
             */
            void evaluateObservables(readdy::model::time_step_type t);

            /**
             * Evaluates all observables, regardless if they are blocked or not.
             */
            void evaluateAllObservables(readdy::model::time_step_type t);

            /**
             * Deconnects the observable of the signal, deletes the
             * corresponding connection block object.
             */
            void deconnectObservable(ObservableBase *const observable);

            /**
             * Registers an observable to the kernel signal.
             */
            std::tuple<std::unique_ptr<ObservableWrapper>, boost::signals2::scoped_connection> registerObservable(const ObservableType &observable, unsigned int stride);

            virtual readdy::model::programs::ProgramFactory &getProgramFactory() const;

            /**
             * Returns a vector containing all available program names for this specific kernel instance.
             *
             * @see createProgram(name)
             * @return The program names.
             */
            std::vector<std::string> getAvailablePrograms() const {
                return getProgramFactory().getAvailablePrograms();
            }

            /**
             * Adds a particle of the type "type" at position "pos".
             */
            void addParticle(const std::string& type, const Vec3 &pos);

            /**
             * @todo implement this properly
             */
            virtual readdy::model::KernelStateModel &getKernelStateModel() const;

            /**
             * @todo implement & document this properly
             */
            virtual readdy::model::KernelContext &getKernelContext() const;

            virtual std::vector<std::string> getAvailablePotentials() const;

            virtual std::unique_ptr<readdy::model::potentials::Potential> createPotential(std::string &name) const;

            template<typename T>
            std::unique_ptr<T> createPotentialAs() const {
                return createPotentialAs<T>(readdy::model::potentials::getPotentialName<T>());
            }

            template<typename T>
            std::unique_ptr<T> createPotentialAs(const std::string &name) const {
                return getPotentialFactory().createPotentialAs<T>(name);
            }

            virtual readdy::model::potentials::PotentialFactory &getPotentialFactory() const;

            virtual readdy::model::reactions::ReactionFactory &getReactionFactory() const;

            virtual readdy::model::_internal::ObservableFactory &getObservableFactory() const;
        protected:
            struct Impl;
            std::unique_ptr<Impl> pimpl;


        };


    }
}

#endif //READDY_MAIN_KERNEL_H
