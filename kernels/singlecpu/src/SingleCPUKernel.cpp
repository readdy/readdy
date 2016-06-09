//
// Created by clonker on 07.03.16.
//

#include <readdy/kernel/singlecpu/SingleCPUKernel.h>
#include <readdy/kernel/singlecpu/SingleCPUProgramFactory.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUTestProgram.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUAddParticleProgram.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUDiffuseProgram.h>
#include <readdy/kernel/singlecpu/potentials/SingleCPUPotentialFactory.h>


namespace readdy {
    namespace kernel {
        namespace singlecpu {
            const std::string SingleCPUKernel::name = "SingleCPU";
            struct SingleCPUKernel::Impl {
                std::unordered_map<std::string, std::shared_ptr<SingleCPUProgramFactory>> programFactories{};
                std::unique_ptr<SingleCPUKernelStateModel> model = std::make_unique<SingleCPUKernelStateModel>();
                std::unique_ptr<readdy::model::KernelContext> context = std::make_unique<readdy::model::KernelContext>();
                std::unique_ptr<readdy::model::RandomProvider> rand = std::make_unique<readdy::model::RandomProvider>();
                std::unique_ptr<potentials::SingleCPUPotentialFactory> potentials;
            };

            SingleCPUKernel::SingleCPUKernel() : readdy::model::Kernel(name), pimpl(std::make_unique<SingleCPUKernel::Impl>()) {
                using factory_ptr_type = std::shared_ptr<SingleCPUProgramFactory>;

                factory_ptr_type ptr = std::make_shared<SingleCPUProgramFactory>(this);
                (*pimpl).programFactories.emplace(programs::SingleCPUTestProgram::getName(), ptr);
                (*pimpl).programFactories.emplace(programs::SingleCPUAddParticleProgram::getName(), ptr);
                (*pimpl).programFactories.emplace(programs::SingleCPUDiffuseProgram::getName(), ptr);


                (*pimpl).potentials = std::make_unique<potentials::SingleCPUPotentialFactory>(this);
            }

            /**
             * factory method
             */
            std::unique_ptr<SingleCPUKernel> SingleCPUKernel::create() {
                return std::make_unique<SingleCPUKernel>();
            }

            /**
             * Destructor: default
             */
            SingleCPUKernel::~SingleCPUKernel() = default;

            std::unique_ptr<readdy::model::Program> SingleCPUKernel::createProgram(const std::string &name) const {
                const auto &&it = (*pimpl).programFactories.find(name);
                if (it != (*pimpl).programFactories.end()) {
                    return (*it->second).createProgram(name);
                }
                return nullptr;
            }

            std::vector<std::string> SingleCPUKernel::getAvailablePrograms() const {
                std::vector<std::string> keys;
                for (auto &&entry : (*pimpl).programFactories) {
                    keys.push_back(entry.first);
                }
                return keys;
            }

            readdy::model::KernelStateModel &SingleCPUKernel::getKernelStateModel() const {
                return *pimpl->model;
            }

            readdy::model::RandomProvider &SingleCPUKernel::getRandomProvider() const {
                return *pimpl->rand;
            }

            readdy::model::KernelContext &SingleCPUKernel::getKernelContext() const {
                return *pimpl->context;
            }

            SingleCPUKernelStateModel &SingleCPUKernel::getKernelStateModelSingleCPU() const {
                return *pimpl->model;
            }

            std::vector<std::string> SingleCPUKernel::getAvailablePotentials() const {
                return pimpl->potentials->getAvailablePotentials();
            }

            std::unique_ptr<readdy::model::potentials::Potential> SingleCPUKernel::createPotential(std::string &name) const {
                return pimpl->potentials->createPotential(name);
            }

            readdy::model::potentials::PotentialFactory &SingleCPUKernel::getPotentialFactory() const {
                return *pimpl->potentials;
            }


            SingleCPUKernel &SingleCPUKernel::operator=(SingleCPUKernel &&rhs) = default;

            SingleCPUKernel::SingleCPUKernel(SingleCPUKernel &&rhs) = default;

        }
    }
}




