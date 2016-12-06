/**
 * << detailed description >>
 *
 * @file Kernel.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 22.11.16
 */

#include "readdy/kernel/cpu_dense/CPUDKernel.h"
#include <readdy/kernel/cpu_dense/observables/CPUDObservableFactory.h>
#include <readdy/kernel/cpu_dense/potentials/CPUDPotentialFactory.h>
#include <readdy/kernel/cpu_dense/programs/CPUDProgramFactory.h>

namespace readdy {
namespace kernel {
namespace cpu_dense {

const std::string CPUDKernel::name = "CPU_Dense";

struct CPUDKernel::Impl {
    std::unique_ptr<programs::CPUDProgramFactory> programFactory;
    std::unique_ptr<potentials::CPUDPotentialFactory> potentialFactory;
    std::unique_ptr<readdy::model::reactions::ReactionFactory> reactionFactory;
    std::unique_ptr<observables::CPUDObservableFactory> observableFactory;
    std::unique_ptr<CPUDStateModel> stateModel;
    std::unique_ptr<readdy::model::KernelContext> context;
    std::unique_ptr<readdy::util::thread::Config> config;
};

readdy::model::Kernel* CPUDKernel::create() {
    return new CPUDKernel();
}

readdy::model::programs::ProgramFactory &CPUDKernel::getProgramFactory() const {
    return *pimpl->programFactory;
}

CPUDKernel::CPUDKernel() : readdy::model::Kernel(name), pimpl(std::make_unique<Impl>()) {
    pimpl->config = std::make_unique<readdy::util::thread::Config>();
    pimpl->reactionFactory = std::make_unique<readdy::model::reactions::ReactionFactory>();
    pimpl->context = std::make_unique<readdy::model::KernelContext>();
    pimpl->programFactory = std::make_unique<programs::CPUDProgramFactory>(this);
    pimpl->stateModel = std::make_unique<CPUDStateModel>(pimpl->context.get(), pimpl->config.get());
    pimpl->potentialFactory = std::make_unique<potentials::CPUDPotentialFactory>(this);
    pimpl->observableFactory = std::make_unique<observables::CPUDObservableFactory>(this);
}

CPUDStateModel &CPUDKernel::getKernelStateModel() const {
    return *pimpl->stateModel;
}

readdy::model::KernelContext &CPUDKernel::getKernelContext() const {
    return *pimpl->context;
}

readdy::model::potentials::PotentialFactory &CPUDKernel::getPotentialFactory() const {
    return *pimpl->potentialFactory;
}

readdy::model::reactions::ReactionFactory &CPUDKernel::getReactionFactory() const {
    return *pimpl->reactionFactory;
}

readdy::model::_internal::ObservableFactory &CPUDKernel::getObservableFactory() const {
    return *pimpl->observableFactory;
}

unsigned long CPUDKernel::getNThreads() const {
    return pimpl->config->nThreads();
}

void CPUDKernel::setNThreads(readdy::util::thread::Config::n_threads_t n) {
    if (n > 0 && n <= std::thread::hardware_concurrency()) {
        pimpl->config->setNThreads(n);
    } else {
        log::console()->error("Tried to set number of threads to {}, but there are only {} hardware thread contexts "
                                      "available.", n, std::thread::hardware_concurrency());
    }
}

CPUDKernel::~CPUDKernel() = default;

}
}
}


const char* name()  {
    return readdy::kernel::cpu_dense::CPUDKernel::name.c_str();
}

readdy::model::Kernel* createKernel() {
    return readdy::kernel::cpu_dense::CPUDKernel::create();
}
