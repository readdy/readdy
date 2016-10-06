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

namespace detail {
template<typename T, typename... Args>
struct get_reaction_dispatcher;
}

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
    std::unique_ptr<T> createObservable(unsigned int stride, Args &&... args) {
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
    template<typename T, typename... Args>
    std::tuple<std::unique_ptr<T>, boost::signals2::scoped_connection>
    createAndConnectObservable(unsigned int stride, Args &&... args) {
        auto &&obs = createObservable<T>(stride, std::forward<Args>(args)...);
        auto &&connection = connectObservable(obs.get());
        return std::make_tuple(std::move(obs), std::move(connection));
    };

    /**
     * Connects an observable to the kernel signal.
     *
     * @return A scoped_connection object that, once deleted, releases the connection of the observable.
     */
    virtual boost::signals2::scoped_connection connectObservable(ObservableBase *const observable);

    /**
     * Evaluates all unblocked observables.
     */
    virtual void evaluateObservables(readdy::model::time_step_type t);

    /**
     * Evaluates all observables, regardless if they are blocked or not.
     */
    virtual void evaluateAllObservables(readdy::model::time_step_type t);

    /**
     * Deconnects the observable of the signal, deletes the
     * corresponding connection block object.
     */
    virtual void disconnectObservable(ObservableBase *const observable);

    /**
     * Registers an observable to the kernel signal.
     */
    virtual std::tuple<std::unique_ptr<ObservableWrapper>, boost::signals2::scoped_connection>
    registerObservable(const ObservableType &observable, unsigned int stride);

    virtual readdy::model::programs::ProgramFactory &getProgramFactory() const = 0;

    /**
     * Returns a vector containing all available program names for this specific kernel instance.
     *
     * @see createProgram(name)
     * @return The program names.
     */
    virtual std::vector<std::string> getAvailablePrograms() const {
        return getProgramFactory().getAvailablePrograms();
    }

    /**
     * Adds a particle of the type "type" at position "pos".
     */
    virtual void addParticle(const std::string &type, const Vec3 &pos);

    /**
     * @todo implement this properly
     */
    virtual readdy::model::KernelStateModel &getKernelStateModel() const = 0;

    /**
     * @todo implement & document this properly
     */
    virtual readdy::model::KernelContext &getKernelContext() const = 0;

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

    template<typename T, typename... Args>
    const short registerReaction(Args &&... args) {
        return getKernelContext().registerReaction(
                detail::get_reaction_dispatcher<T, Args...>::impl(this, std::forward<Args>(args)...));
    };

    std::unique_ptr<reactions::Reaction<1>>
    createConversionReaction(const std::string &name, const std::string &from, const std::string &to,
                             const double rate) const;

    std::unique_ptr<reactions::Reaction<2>>
    createEnzymaticReaction(const std::string &name, const std::string &catalyst, const std::string &from,
                            const std::string &to, const double rate, const double eductDistance) const;

    std::unique_ptr<reactions::Reaction<1>>
    createFissionReaction(const std::string &name, const std::string &from, const std::string &to1,
                          const std::string &to2, const double rate, const double productDistance,
                          const double weight1 = 0.5, const double weight2 = 0.5) const;

    std::unique_ptr<reactions::Reaction<2>>
    createFusionReaction(const std::string &name, const std::string &from1, const std::string &from2,
                         const std::string &to, const double rate, const double eductDistance,
                         const double weight1 = 0.5, const double weight2 = 0.5) const;

    std::unique_ptr<reactions::Reaction<1>>
    createDecayReaction(const std::string &name, const std::string &type, const double rate) const;

    virtual readdy::model::potentials::PotentialFactory &getPotentialFactory() const = 0;

    virtual readdy::model::reactions::ReactionFactory &getReactionFactory() const = 0;

    virtual readdy::model::_internal::ObservableFactory &getObservableFactory() const;

    virtual unsigned int getTypeId(const std::string &) const;

protected:
    struct Impl;
    std::unique_ptr<Impl> pimpl;
};
namespace detail {
template<typename... Args>
struct get_reaction_dispatcher<readdy::model::reactions::Conversion, Args...> {
    static std::unique_ptr<readdy::model::reactions::Reaction<1>>
    impl(const readdy::model::Kernel *const self, Args &&... args) {
        return self->createConversionReaction(std::forward<Args>(args)...);
    }
};

template<typename... Args>
struct get_reaction_dispatcher<readdy::model::reactions::Enzymatic, Args...> {
    static std::unique_ptr<readdy::model::reactions::Reaction<2>>
    impl(const readdy::model::Kernel *const self, Args &&... args) {
        return self->createEnzymaticReaction(std::forward<Args>(args)...);
    }
};

template<typename... Args>
struct get_reaction_dispatcher<readdy::model::reactions::Fission, Args...> {
    static std::unique_ptr<readdy::model::reactions::Reaction<1>>
    impl(const readdy::model::Kernel *const self, Args &&... args) {
        return self->createFissionReaction(std::forward<Args>(args)...);
    }
};

template<typename... Args>
struct get_reaction_dispatcher<readdy::model::reactions::Fusion, Args...> {
    static std::unique_ptr<readdy::model::reactions::Reaction<2>>
    impl(const readdy::model::Kernel *const self, Args &&... args) {
        return self->createFusionReaction(std::forward<Args>(args)...);
    }
};

template<typename... Args>
struct get_reaction_dispatcher<readdy::model::reactions::Decay, Args...> {
    static std::unique_ptr<readdy::model::reactions::Reaction<1>>
    impl(const readdy::model::Kernel *const self, Args &&... args) {
        return self->createDecayReaction(std::forward<Args>(args)...);
    }
};
}

}
}

#endif //READDY_MAIN_KERNEL_H
