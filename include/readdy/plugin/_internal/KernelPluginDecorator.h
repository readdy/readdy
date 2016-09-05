/**
 * The KernelPluginDecorator class wraps a loaded kernel instance and the corresponding boost::dll::shared_library
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

                virtual boost::signals2::scoped_connection
                connectObservable(model::ObservableBase *const observable) override;

                virtual void disconnectObservable(model::ObservableBase *const observable) override;

                virtual std::unique_ptr<model::programs::Program> createProgram(const std::string &name) const override;

                virtual void evaluateObservables(readdy::model::time_step_type t) override;

                virtual void evaluateAllObservables(readdy::model::time_step_type t) override;

                virtual std::tuple<std::unique_ptr<model::ObservableWrapper>, boost::signals2::scoped_connection>
                registerObservable(const model::ObservableType &observable, unsigned int stride) override;

                virtual std::vector<std::string> getAvailablePrograms() const override;

                virtual void addParticle(const std::string &type, const model::Vec3 &pos) override;

                virtual unsigned int getTypeId(const std::string &string) const override;


            protected:
                virtual readdy::model::_internal::ObservableFactory &getObservableFactory() const override;
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
