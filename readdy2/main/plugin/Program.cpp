//
// Created by clonker on 08.04.16.
//


#include <readdy/plugin/Program.h>
#include <boost/make_unique.hpp>

struct readdy::plugin::Program::Impl {
    std::string name;
};

// constructor initializing impl
readdy::plugin::Program::Program(const std::string name) : impl_ptr{boost::make_unique<readdy::plugin::Program::Impl>()} {
    (*impl_ptr).name = name;
}

const std::string readdy::plugin::Program::getName() const {
    return (*impl_ptr).name;
}

