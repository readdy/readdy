//
// Created by clonker on 14.03.16.
//

#include <readdy/plugin/_internal/KernelPluginDecorator.h>

namespace plug = readdy::plugin;

readdy::plugin::KernelPluginDecorator::KernelPluginDecorator(const readdy::plugin::Kernel &reference)
        : reference(reference), readdy::plugin::Kernel(reference.getName()) {
    BOOST_LOG_TRIVIAL(debug) << "Wrapping kernel " << reference.getName() << " with kernel decorator.";
}
