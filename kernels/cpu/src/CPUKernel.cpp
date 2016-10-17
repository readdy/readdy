/**
 * << detailed description >>
 *
 * @file CPUKernel.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.06.16
 */

#include <readdy/kernel/cpu/CPUKernel.h>
#include <readdy/kernel/cpu/programs/CPUProgramFactory.h>
#include <readdy/kernel/cpu/potentials/CPUPotentialFactory.h>
#include <readdy/kernel/cpu/reactions/CPUReactionFactory.h>
#include <readdy/kernel/cpu/observables/ObservableFactory.h>
#include <thread>
#include <readdy/kernel/cpu/util/Config.h>

namespace readdy {
namespace kernel {
namespace cpu {
const std::string CPUKernel::name = "CPU";

struct CPUKernel::Impl {
    std::unique_ptr<programs::CPUProgramFactory> programFactory;
    std::unique_ptr<potentials::CPUPotentialFactory> potentialFactory;
    std::unique_ptr<reactions::CPUReactionFactory> reactionFactory;
    std::unique_ptr<observables::ObservableFactory> observableFactory;
    std::unique_ptr<CPUStateModel> stateModel;
    std::unique_ptr<readdy::model::KernelContext> context;
    std::unique_ptr<util::Config> config;
};

readdy::model::Kernel* CPUKernel::create() {
    return new CPUKernel();
}

readdy::model::programs::ProgramFactory &CPUKernel::getProgramFactory() const {
    return *pimpl->programFactory;
}

CPUKernel::CPUKernel() : readdy::model::Kernel(name), pimpl(std::make_unique<Impl>()) {
    pimpl->config = std::make_unique<util::Config>();
    pimpl->reactionFactory = std::make_unique<reactions::CPUReactionFactory>(this);
    pimpl->context = std::make_unique<readdy::model::KernelContext>();
    pimpl->programFactory = std::make_unique<programs::CPUProgramFactory>(this);
    pimpl->stateModel = std::make_unique<CPUStateModel>(pimpl->context.get(), pimpl->config.get());
    pimpl->potentialFactory = std::make_unique<potentials::CPUPotentialFactory>(this);
    pimpl->observableFactory = std::make_unique<observables::ObservableFactory>(this);
}

CPUStateModel &CPUKernel::getKernelStateModel() const {
    return *pimpl->stateModel;
}

readdy::model::KernelContext &CPUKernel::getKernelContext() const {
    return *pimpl->context;
}

readdy::model::potentials::PotentialFactory &CPUKernel::getPotentialFactory() const {
    return *pimpl->potentialFactory;
}

readdy::model::reactions::ReactionFactory &CPUKernel::getReactionFactory() const {
    return *pimpl->reactionFactory;
}

readdy::model::_internal::ObservableFactory &CPUKernel::getObservableFactory() const {
    return *pimpl->observableFactory;
}

unsigned long CPUKernel::getNThreads() const {
    return pimpl->config->nThreads;
}

void CPUKernel::setNThreads(unsigned long n) {
    if (n > 0 && n < std::thread::hardware_concurrency()) {
        pimpl->config->nThreads = n;
    }
}

CPUKernel::~CPUKernel() = default;


}
}
}


const char* name()  {
    return readdy::kernel::cpu::CPUKernel::name.c_str();
}

readdy::model::Kernel* createKernel() {
    return readdy::kernel::cpu::CPUKernel::create();
}
