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
        : readdy::plugin::Kernel::Kernel(sharedLib.string()) {
    BOOST_LOG_TRIVIAL(debug) << "... loading kernel " << sharedLib.string();
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
        typedef std::shared_ptr<readdy::plugin::Kernel> (kernel_t)();
        boost::function<kernel_t> factory = dll::import_alias<kernel_t>(lib, "createKernel");
        reference = factory();
        BOOST_LOG_TRIVIAL(debug) << "loaded.";
    }
}

readdy::model::KernelStateModel& readdy::plugin::_internal::KernelPluginDecorator::getKernelStateModel() const {
    return (*reference).getKernelStateModel();
}

std::unique_ptr<readdy::plugin::Program> readdy::plugin::_internal::KernelPluginDecorator::createProgram(const std::string& name) const {
    return (*reference).createProgram(name);
}

std::vector<std::string> readdy::plugin::_internal::KernelPluginDecorator::getAvailablePrograms() const {
    return (*reference).getAvailablePrograms();
}

readdy::plugin::_internal::KernelPluginDecorator::~KernelPluginDecorator() {
    BOOST_LOG_TRIVIAL(debug) << "destroying decorator of " << getName();
    BOOST_LOG_TRIVIAL(debug) << "use count: " << reference.use_count();
    reference.reset();
}

readdy::model::KernelContext& readdy::plugin::_internal::KernelPluginDecorator::getKernelContext() const {
    return (*reference).getKernelContext();
}

plug::InvalidPluginException::InvalidPluginException(const std::string &__arg) : runtime_error(__arg) { }
