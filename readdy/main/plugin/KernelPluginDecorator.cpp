//
// Created by clonker on 14.03.16.
//

#include <readdy/plugin/_internal/KernelPluginDecorator.h>
#include <boost/dll.hpp>

namespace plug = readdy::plugin::_internal;
namespace dll = boost::dll;

const std::string &plug::KernelPluginDecorator::getName() const {
    return reference.get()->getName();
}

readdy::plugin::_internal::KernelPluginDecorator::KernelPluginDecorator(const boost::filesystem::path sharedLib)
        : readdy::model::Kernel::Kernel(sharedLib.string()) {
    // load the library
    lib = dll::shared_library(sharedLib, dll::load_mode::rtld_lazy | dll::load_mode::rtld_global);
    // check if library has required symbol
    if (!lib.has("createKernel")) {
        BOOST_LOG_TRIVIAL(error) << "... skipping, since it had no createKernel symbol";
        if (lib.is_loaded()) lib.unload();
        std::string errMsg = std::string("library ").append(sharedLib.string()).append(" had no createKernel symbol.");
        throw plug::InvalidPluginException(errMsg);
    }
    // load the kernel
    {
        typedef std::unique_ptr<readdy::model::Kernel> (kernel_t)();
        boost::function<kernel_t> factory = dll::import_alias<kernel_t>(lib, "createKernel");
        reference = factory();
    }
}

readdy::model::KernelStateModel& readdy::plugin::_internal::KernelPluginDecorator::getKernelStateModel() const {
    return (*reference).getKernelStateModel();
}

readdy::plugin::_internal::KernelPluginDecorator::~KernelPluginDecorator() {
    reference.reset();
}

readdy::model::KernelContext& readdy::plugin::_internal::KernelPluginDecorator::getKernelContext() const {
    return (*reference).getKernelContext();
}

std::unique_ptr<readdy::model::potentials::Potential> readdy::plugin::_internal::KernelPluginDecorator::createPotential(std::string &name) const {
    return reference->createPotential(name);
}

std::vector<std::string> readdy::plugin::_internal::KernelPluginDecorator::getAvailablePotentials() const {
    return reference->getAvailablePotentials();
}

readdy::model::potentials::PotentialFactory &readdy::plugin::_internal::KernelPluginDecorator::getPotentialFactory() const {
    return reference->getPotentialFactory();
}

readdy::model::programs::ProgramFactory &readdy::plugin::_internal::KernelPluginDecorator::getProgramFactory() const {
    return reference->getProgramFactory();
}

readdy::model::reactions::ReactionFactory &readdy::plugin::_internal::KernelPluginDecorator::getReactionFactory() const {
    return reference->getReactionFactory();
}

readdy::model::_internal::ObservableFactory &readdy::plugin::_internal::KernelPluginDecorator::getObservableFactory() const {
    return reference->getObservableFactory();
}


plug::InvalidPluginException::InvalidPluginException(const std::string &__arg) : runtime_error(__arg) { }

const std::string readdy::plugin::_internal::loadKernelName(const boost::filesystem::path& sharedLib) {
    auto lib = dll::shared_library(sharedLib, dll::load_mode::rtld_lazy | dll::load_mode::rtld_global);
    if(!lib.has("name")) {
        if (lib.is_loaded()) lib.unload();
        std::string errMsg = std::string("library ").append(sharedLib.string()).append(" had no name() symbol.");
        throw plug::InvalidPluginException(errMsg);
    } else {
        boost::shared_ptr<std::string> name = dll::import_alias<std::string>(lib, "name");
        return *name;
    }
}