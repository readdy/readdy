//
// Created by clonker on 14.03.16.
//

#include <readdy/plugin/_internal/KernelPluginDecorator.h>

namespace plug = readdy::plugin::_internal;

const std::string plug::KernelPluginDecorator::getName() const {
    return reference.getName();
}
// TODO write this properly:
// declare copy, move constructors (see http://stackoverflow.com/a/18300974 | http://stackoverflow.com/questions/8114276/how-do-i-pass-a-unique-ptr-argument-to-a-constructor-or-a-function) | http://stackoverflow.com/questions/3279543/what-is-the-copy-and-swap-idiom | http://stackoverflow.com/questions/3106110/what-are-move-semantics
// declare constructor that loads the lib via dll (!)
// rename: Decorator -> Adapter (or something something, as we then bridge low level api w/ high level api)
readdy::plugin::_internal::KernelPluginDecorator::KernelPluginDecorator(readdy::plugin::Kernel &&reference, boost::dll::shared_library &&lib)
        : reference(std::move(reference)),
          lib(std::move(lib)),
          readdy::plugin::Kernel(reference.getName()) {
}
