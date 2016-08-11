/**
 * << detailed description >>
 *
 * @file KernelWrap.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 03.08.16
 */

#include <boost/python.hpp>
#include <readdy/model/Kernel.h>

namespace bpy = boost::python;

using scpu_kernel_t = readdy::kernel::singlecpu::SingleCPUKernel;
namespace readdy {
    namespace py {
        struct KernelWrap : public readdy::model::Kernel, bpy::wrapper<readdy::model::Kernel> {
            KernelWrap(const std::string &name) : Kernel(name) {}

            virtual readdy::model::programs::ProgramFactory &getProgramFactory() const override {
                return this->get_override("getProgramFactory")();
            }

            virtual readdy::model::KernelStateModel &getKernelStateModel() const override {
                return this->get_override("getKernelStateModel")();
            }

            virtual readdy::model::KernelContext &getKernelContext() const override {
                return this->get_override("getKernelContext")();
            }

            virtual readdy::model::potentials::PotentialFactory &getPotentialFactory() const override {
                return this->get_override("getPotentialFactory")();
            }

            virtual readdy::model::reactions::ReactionFactory &getReactionFactory() const override {
                return this->get_override("getReactionFactory")();
            }

        };

        struct SingleCPUKernelWrap : public scpu_kernel_t, bpy::wrapper<scpu_kernel_t>  {
            virtual kernel::singlecpu::SingleCPUKernelStateModel &getKernelStateModel() const override {
                if(auto f = this->get_override("get_kernel_state_model")) return f();
                return scpu_kernel_t::getKernelStateModel();
            }

            virtual model::KernelContext &getKernelContext() const override {
                if(auto f = this->get_override("get_kernel_context")) return f();
                return scpu_kernel_t::getKernelContext();
            }

            virtual model::programs::ProgramFactory &getProgramFactory() const override {
                if(auto f = this->get_override("get_program_factory")) return f();
                return scpu_kernel_t::getProgramFactory();
            }

            virtual std::vector<std::string> getAvailablePotentials() const override {
                if(auto f = this->get_override("get_available_potentials")) return f();
                return scpu_kernel_t::getAvailablePotentials();
            }

            virtual model::potentials::PotentialFactory &getPotentialFactory() const override {
                if(auto f = this->get_override("get_potential_factory")) return f();
                return scpu_kernel_t::getPotentialFactory();
            }

            virtual model::reactions::ReactionFactory &getReactionFactory() const override {
                if(auto f = this->get_override("get_reaction_factory")) return f();
                return scpu_kernel_t::getReactionFactory();
            }

            virtual model::_internal::ObservableFactory &getObservableFactory() const override {
                if(auto f = this->get_override("get_observable_factory")) return f();
                return scpu_kernel_t::getObservableFactory();
            }
        };
    }
}
