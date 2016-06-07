//
// Created by clonker on 07.03.16.
//

#include <readdy/kernel/singlecpu/SingleCPUKernel.h>
#include <readdy/kernel/singlecpu/SingleCPUProgramFactory.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUTestProgram.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUAddParticleProgram.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUDiffuseProgram.h>

namespace kern = readdy::kernel::singlecpu;

const std::string kern::SingleCPUKernel::name = "SingleCPU";

struct readdy::kernel::singlecpu::SingleCPUKernel::Impl {
    std::unordered_map<std::string, std::shared_ptr<kern::SingleCPUProgramFactory>> programFactories {};
    std::unique_ptr<kern::SingleCPUKernelStateModel> model = std::make_unique<kern::SingleCPUKernelStateModel>();
    std::shared_ptr<readdy::model::KernelContext> context = std::make_shared<readdy::model::KernelContext>();
    std::shared_ptr<readdy::model::RandomProvider> rand = std::make_shared<readdy::model::RandomProvider>();
};
kern::SingleCPUKernel:: SingleCPUKernel() : readdy::model::Kernel(name), pimpl(std::make_unique<kern::SingleCPUKernel::Impl>()){
    BOOST_LOG_TRIVIAL(debug) << "Single CPU Kernel instantiated, registering program factories...";
    using factory_ptr_type = std::shared_ptr<kern::SingleCPUProgramFactory>;

    factory_ptr_type ptr = std::make_shared<kern::SingleCPUProgramFactory>(this);
    (*pimpl).programFactories.emplace(kern::programs::SingleCPUTestProgram::getName(), ptr);
    (*pimpl).programFactories.emplace(kern::programs::SingleCPUAddParticleProgram::getName(), ptr);
    (*pimpl).programFactories.emplace(kern::programs::SingleCPUDiffuseProgram::getName(), ptr);
    BOOST_LOG_TRIVIAL(debug) << "...done";
}

/**
 * factory method
 */
std::unique_ptr<kern::SingleCPUKernel> kern::SingleCPUKernel::create() {
    return std::make_unique<kern::SingleCPUKernel>();
}
/**
 * Destructor: default
 */
readdy::kernel::singlecpu::SingleCPUKernel::~SingleCPUKernel() = default;

std::unique_ptr<readdy::model::Program> readdy::kernel::singlecpu::SingleCPUKernel::createProgram(const std::string& name) const {
    const auto&& it = (*pimpl).programFactories.find(name);
    if(it != (*pimpl).programFactories.end()) {
        return (*it->second).createProgram(name);
    }
    return nullptr;
}

std::vector<std::string> readdy::kernel::singlecpu::SingleCPUKernel::getAvailablePrograms() const {
    std::vector<std::string> keys;
    for(auto&& entry : (*pimpl).programFactories) {
        keys.push_back(entry.first);
    }
    return keys;
}

readdy::model::KernelStateModel& readdy::kernel::singlecpu::SingleCPUKernel::getKernelStateModel() const {
    return *pimpl->model;
}

readdy::model::RandomProvider& readdy::kernel::singlecpu::SingleCPUKernel::getRandomProvider() const {
    return *pimpl->rand;
}

readdy::model::KernelContext& readdy::kernel::singlecpu::SingleCPUKernel::getKernelContext() const {
    return *pimpl->context;
}

readdy::kernel::singlecpu::SingleCPUKernelStateModel &readdy::kernel::singlecpu::SingleCPUKernel::getKernelStateModelSingleCPU() const {
    return *pimpl->model;
}

kern::SingleCPUKernel &kern::SingleCPUKernel::operator=(kern::SingleCPUKernel&& rhs) = default;
kern::SingleCPUKernel::SingleCPUKernel(kern::SingleCPUKernel&& rhs) = default;


