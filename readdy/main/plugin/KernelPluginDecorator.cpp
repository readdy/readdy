//
// Created by clonker on 14.03.16.
//

#include <readdy/plugin/_internal/KernelPluginDecorator.h>
#include <readdy/common/nodelete.h>

namespace plug = readdy::plugin::_internal;
namespace dll = readdy::util::dll;

const std::string &plug::KernelPluginDecorator::getName() const {
    return reference.get()->getName();
}

readdy::model::KernelStateModel &readdy::plugin::_internal::KernelPluginDecorator::getKernelStateModel() const {
    return (*reference).getKernelStateModel();
}

readdy::plugin::_internal::KernelPluginDecorator::KernelPluginDecorator(const std::string& sharedLib)
        : readdy::model::Kernel::Kernel(sharedLib) {
    // load the library
    lib = std::make_unique<dll::shared_library>(sharedLib, RTLD_LAZY | RTLD_GLOBAL);
    // check if library has required symbol
    if (!lib->has_symbol("createKernel")) {
        log::console()->error("skipping, since it had no createKernel symbol");
        throw plug::InvalidPluginException("library " + sharedLib + " had no createKernel symbol.");
    }
    // load the kernel
    {
        readdy::model::Kernel* kernel = lib->load<readdy::model::Kernel*()>("createKernel")();
        log::console()->debug("loaded kernel with name {}", kernel->getName());
        reference = std::unique_ptr<readdy::model::Kernel>(std::move(kernel));
    }
}

readdy::plugin::_internal::KernelPluginDecorator::~KernelPluginDecorator() {
    // release and delete reference kernel
    reference.reset(nullptr);
}

readdy::model::KernelContext &readdy::plugin::_internal::KernelPluginDecorator::getKernelContext() const {
    return reference->getKernelContext();
}

std::unique_ptr<readdy::model::potentials::Potential>
readdy::plugin::_internal::KernelPluginDecorator::createPotential(std::string &name) const {
    return reference->createPotential(name);
}

std::vector<std::string> readdy::plugin::_internal::KernelPluginDecorator::getAvailablePotentials() const {
    return reference->getAvailablePotentials();
}

readdy::model::potentials::PotentialFactory &
readdy::plugin::_internal::KernelPluginDecorator::getPotentialFactory() const {
    return reference->getPotentialFactory();
}

readdy::model::programs::ProgramFactory &readdy::plugin::_internal::KernelPluginDecorator::getProgramFactory() const {
    return reference->getProgramFactory();
}

readdy::model::reactions::ReactionFactory &
readdy::plugin::_internal::KernelPluginDecorator::getReactionFactory() const {
    return reference->getReactionFactory();
}

readdy::model::observables::ObservableFactory &
readdy::plugin::_internal::KernelPluginDecorator::getObservableFactory() const {
    return reference->getObservableFactory();
}

readdy::signals::scoped_connection
readdy::plugin::_internal::KernelPluginDecorator::connectObservable(readdy::model::observables::ObservableBase *const observable) {
    return reference->connectObservable(observable);
}

std::unique_ptr<readdy::model::programs::Program>
readdy::plugin::_internal::KernelPluginDecorator::createProgram(const std::string &name) const {
    return reference->createProgram(name);
}

void readdy::plugin::_internal::KernelPluginDecorator::evaluateObservables(model::observables::time_step_type t) {
    reference->evaluateObservables(t);
}

std::tuple<std::unique_ptr<readdy::model::observables::ObservableWrapper>, readdy::signals::scoped_connection>
readdy::plugin::_internal::KernelPluginDecorator::registerObservable(
        const model::observables::observable_type &observable, unsigned int stride
) {
    return reference->registerObservable(observable, stride);
}

std::vector<std::string> readdy::plugin::_internal::KernelPluginDecorator::getAvailablePrograms() const {
    return reference->getAvailablePrograms();
}

void
readdy::plugin::_internal::KernelPluginDecorator::addParticle(const std::string &type, const readdy::model::Vec3 &pos) {
    reference->addParticle(type, pos);
}

unsigned int readdy::plugin::_internal::KernelPluginDecorator::getTypeId(const std::string &string) const {
    return reference->getTypeId(string);
}


plug::InvalidPluginException::InvalidPluginException(const std::string &__arg) : runtime_error(__arg) {}

const std::string readdy::plugin::_internal::loadKernelName(const std::string &sharedLib) {
    auto lib = dll::shared_library(sharedLib, RTLD_LAZY | RTLD_GLOBAL);
    if (!lib.has_symbol("name")) {
        throw plug::InvalidPluginException("library " + sharedLib + " had no name() symbol");
    } else {
        auto fun = lib.load<const char*()>("name");
        return fun();
    }
}