/**
 * << detailed description >>
 *
 * @file KernelWrap.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 03.08.16
 */

#include <readdy/model/Kernel.h>

namespace readdy {
namespace py {
struct KernelWrap : public readdy::model::Kernel {
    KernelWrap(const std::string &name) : Kernel(name) {}

    virtual readdy::model::programs::ProgramFactory &getProgramFactory() const override {
        PYBIND11_OVERLOAD_PURE(readdy::model::programs::ProgramFactory &, readdy::model::Kernel, getProgramFactory,);
    }

    virtual readdy::model::KernelStateModel &getKernelStateModel() const override {
        PYBIND11_OVERLOAD_PURE(readdy::model::KernelStateModel &, readdy::model::Kernel, getKernelStateModel,);
    }

    virtual readdy::model::KernelContext &getKernelContext() const override {
        PYBIND11_OVERLOAD_PURE(readdy::model::KernelContext &, readdy::model::Kernel, getKernelContext,);
    }

    virtual readdy::model::potentials::PotentialFactory &getPotentialFactory() const override {
        PYBIND11_OVERLOAD_PURE(
                readdy::model::potentials::PotentialFactory &, readdy::model::Kernel, getPotentialFactory,
        );
    }

    virtual readdy::model::reactions::ReactionFactory &getReactionFactory() const override {
        PYBIND11_OVERLOAD_PURE(readdy::model::reactions::ReactionFactory &, readdy::model::Kernel, getReactionFactory,)
    }

};

struct SingleCPUKernelWrap : public readdy::kernel::singlecpu::SingleCPUKernel {
    virtual kernel::singlecpu::SingleCPUKernelStateModel &getKernelStateModel() const override {
        PYBIND11_OVERLOAD_NAME(kernel::singlecpu::SingleCPUKernelStateModel &,
                               readdy::kernel::singlecpu::SingleCPUKernel, "get_kernel_state_model",
                               getKernelStateModel,);
    }

    virtual model::KernelContext &getKernelContext() const override {
        PYBIND11_OVERLOAD_NAME(model::KernelContext &, readdy::kernel::singlecpu::SingleCPUKernel, "get_kernel_context",
                               getKernelContext,);
    }

    virtual model::programs::ProgramFactory &getProgramFactory() const override {
        PYBIND11_OVERLOAD_NAME(model::programs::ProgramFactory &, readdy::kernel::singlecpu::SingleCPUKernel,
                               "get_program_factory", getProgramFactory,);
    }

    virtual std::vector<std::string> getAvailablePotentials() const override {
        PYBIND11_OVERLOAD_NAME(std::vector<std::string>, readdy::kernel::singlecpu::SingleCPUKernel,
                               "get_available_potentials", getAvailablePotentials,);
    }

    virtual model::potentials::PotentialFactory &getPotentialFactory() const override {
        PYBIND11_OVERLOAD_NAME(model::potentials::PotentialFactory &, readdy::kernel::singlecpu::SingleCPUKernel,
                               "get_potential_factory", getPotentialFactory,);
    }

    virtual model::reactions::ReactionFactory &getReactionFactory() const override {
        PYBIND11_OVERLOAD_NAME(model::reactions::ReactionFactory &, readdy::kernel::singlecpu::SingleCPUKernel,
                               "get_reaction_factory", getReactionFactory,);
    }

    virtual model::_internal::ObservableFactory &getObservableFactory() const override {
        PYBIND11_OVERLOAD_NAME(model::_internal::ObservableFactory &, readdy::kernel::singlecpu::SingleCPUKernel,
                               "get_observable_factory", getObservableFactory,);
    }
};
}
}
