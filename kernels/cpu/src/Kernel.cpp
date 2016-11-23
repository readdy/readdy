/**
 * << detailed description >>
 *
 * @file CPUKernel.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.06.16
 */

#include <readdy/kernel/cpu/Kernel.h>
#include <readdy/kernel/cpu/programs/ProgramFactory.h>
#include <readdy/kernel/cpu/potentials/PotentialFactory.h>
#include <readdy/kernel/cpu/reactions/ReactionFactory.h>
#include <readdy/kernel/cpu/observables/ObservableFactory.h>

namespace readdy {
namespace kernel {
namespace cpu {
const std::string Kernel::name = "CPU";

struct Kernel::Impl {
    std::unique_ptr<programs::ProgramFactory> programFactory;
    std::unique_ptr<potentials::PotentialFactory> potentialFactory;
    std::unique_ptr<reactions::ReactionFactory> reactionFactory;
    std::unique_ptr<observables::ObservableFactory> observableFactory;
    std::unique_ptr<StateModel> stateModel;
    std::unique_ptr<readdy::model::KernelContext> context;
    std::unique_ptr<readdy::util::thread::Config> config;
};

readdy::model::Kernel* Kernel::create() {
    return new Kernel();
}

readdy::model::programs::ProgramFactory &Kernel::getProgramFactory() const {
    return *pimpl->programFactory;
}

Kernel::Kernel() : readdy::model::Kernel(name), pimpl(std::make_unique<Impl>()) {
    pimpl->config = std::make_unique<readdy::util::thread::Config>();
    pimpl->reactionFactory = std::make_unique<reactions::ReactionFactory>(this);
    pimpl->context = std::make_unique<readdy::model::KernelContext>();
    pimpl->programFactory = std::make_unique<programs::ProgramFactory>(this);
    pimpl->stateModel = std::make_unique<StateModel>(pimpl->context.get(), pimpl->config.get());
    pimpl->potentialFactory = std::make_unique<potentials::PotentialFactory>(this);
    pimpl->observableFactory = std::make_unique<observables::ObservableFactory>(this);
}

StateModel &Kernel::getKernelStateModel() const {
    return *pimpl->stateModel;
}

readdy::model::KernelContext &Kernel::getKernelContext() const {
    return *pimpl->context;
}

readdy::model::potentials::PotentialFactory &Kernel::getPotentialFactory() const {
    return *pimpl->potentialFactory;
}

readdy::model::reactions::ReactionFactory &Kernel::getReactionFactory() const {
    return *pimpl->reactionFactory;
}

readdy::model::_internal::ObservableFactory &Kernel::getObservableFactory() const {
    return *pimpl->observableFactory;
}

unsigned long Kernel::getNThreads() const {
    return pimpl->config->nThreads();
}

void Kernel::setNThreads(readdy::util::thread::Config::n_threads_t n) {
    if (n > 0 && n <= std::thread::hardware_concurrency()) {
        pimpl->config->setNThreads(n);
    } else {
        log::console()->error("Tried to set number of threads to {}, but there are only {} hardware thread contexts "
                                      "available.", n, std::thread::hardware_concurrency());
    }
}

Kernel::~Kernel() = default;


}
}
}


const char* name()  {
    return readdy::kernel::cpu::Kernel::name.c_str();
}

readdy::model::Kernel* createKernel() {
    return readdy::kernel::cpu::Kernel::create();
}
