//
// Created by clonker on 14.03.16.
//

#include <readdy/plugin/_internal/KernelPluginDecorator.h>

namespace plug = readdy::plugin::_internal;

plug::KernelPluginDecorator::KernelPluginDecorator(const readdy::plugin::Kernel &reference, boost::dll::shared_library &&lib)
        : reference(reference), lib(lib), readdy::plugin::Kernel(reference.getName()) {
    BOOST_LOG_TRIVIAL(debug) << "Wrapping kernel " << reference.getName() << " with kernel decorator.";
}

const std::string plug::KernelPluginDecorator::getName() const {
    return reference.getName();
}
