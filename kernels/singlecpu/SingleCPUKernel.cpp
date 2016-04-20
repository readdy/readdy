//
// Created by clonker on 07.03.16.
//

#include "SingleCPUKernel.h"
#include "SingleCPUProgramFactory.h"
#include "programs/SingleCPUTestProgram.h"
#include "programs/SingleCPUAddParticleProgram.h"
#include "programs/SingleCPUDiffuseProgram.h"
#include <boost/make_unique.hpp>
#include <unordered_map>

namespace kern = readdy::kernel::singlecpu;
struct readdy::kernel::singlecpu::SingleCPUKernel::Impl {
    std::unordered_map<std::string, std::shared_ptr<kern::SingleCPUProgramFactory>> programFactories {};
    std::shared_ptr<kern::SingleCPUKernelStateModel> model = std::make_shared<kern::SingleCPUKernelStateModel>();
    std::shared_ptr<readdy::model::KernelContext> context = std::make_shared<readdy::model::KernelContext>();
    std::shared_ptr<readdy::utils::RandomProvider> rand = std::make_shared<readdy::utils::RandomProvider>();
};
kern::SingleCPUKernel:: SingleCPUKernel() : readdy::plugin::Kernel("SingleCPU"), pimpl(boost::make_unique<kern::SingleCPUKernel::Impl>()){
    BOOST_LOG_TRIVIAL(debug) << "Single CPU Kernel instantiated, registering program factories...";
    using factory_ptr_type = std::shared_ptr<kern::SingleCPUProgramFactory>;

    factory_ptr_type ptr = std::make_shared<kern::SingleCPUProgramFactory>(pimpl->context, pimpl->model, pimpl->rand);
    (*pimpl).programFactories.emplace(kern::programs::SingleCPUTestProgram::getName(), ptr);
    (*pimpl).programFactories.emplace(kern::programs::SingleCPUAddParticleProgram::getName(), ptr);
    (*pimpl).programFactories.emplace(kern::programs::SingleCPUDiffuseProgram::getName(), ptr);
    BOOST_LOG_TRIVIAL(debug) << "...done";
}

/**
 * factory method
 */
std::shared_ptr<kern::SingleCPUKernel> kern::SingleCPUKernel::create() {
    return std::make_shared<kern::SingleCPUKernel>();
}
/**
 * Destructor: default
 */
readdy::kernel::singlecpu::SingleCPUKernel::~SingleCPUKernel() = default;
/**
 * Copy operations
 */
/*kern::SingleCPUKernel::SingleCPUKernel(const kern::SingleCPUKernel &rhs) : readdy::plugin::Kernel(rhs), pimpl(boost::make_unique<kern::SingleCPUKernel::Impl>()) {};
kern::SingleCPUKernel &kern::SingleCPUKernel::operator=(kern::SingleCPUKernel &rhs) {
    *pimpl = *rhs.pimpl;
    return *this;
}*/

std::shared_ptr<readdy::plugin::Program> readdy::kernel::singlecpu::SingleCPUKernel::createProgram(std::string name) {
    auto it = (*pimpl).programFactories.find(name);
    if(it != (*pimpl).programFactories.end()) {
        return (*it->second).createProgram(name);
    }
    return nullptr;
}

std::vector<std::string> readdy::kernel::singlecpu::SingleCPUKernel::getAvailablePrograms() {
    std::vector<std::string> keys;
    for(auto&& entry : (*pimpl).programFactories) {
        keys.push_back(entry.first);
    }
    return keys;
}

std::shared_ptr<readdy::model::KernelStateModel> readdy::kernel::singlecpu::SingleCPUKernel::getKernelStateModel() {
    return (*pimpl).model;
}

std::shared_ptr<readdy::utils::RandomProvider> readdy::kernel::singlecpu::SingleCPUKernel::getRandomProvider() const {
    return pimpl->rand;
}

std::shared_ptr<readdy::model::KernelContext> readdy::kernel::singlecpu::SingleCPUKernel::getKernelContext() {
    return pimpl->context;
}

/**
 * Move operations default
 */
kern::SingleCPUKernel &kern::SingleCPUKernel::operator=(kern::SingleCPUKernel&& rhs) = default;
kern::SingleCPUKernel::SingleCPUKernel(kern::SingleCPUKernel&& rhs) = default;



