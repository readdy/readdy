//
// Created by mho on 7/2/18.
//

#include <pybind11/pybind11.h>

#include <readdy/kernel/cpu/CPUKernel.h>
#include <readdy/model/actions/UserDefinedAction.h>

namespace py = pybind11;
namespace rnd = readdy::model::rnd;

class EBDIntegrator : public readdy::model::actions::UserDefinedAction {
public:

    explicit EBDIntegrator(readdy::scalar timeStep) : UserDefinedAction(timeStep) {}

    void perform(const readdy::util::PerformanceNode &node) override {
        auto cpuKernel = dynamic_cast<readdy::kernel::cpu::CPUKernel*>(kernel());
        auto data = cpuKernel->getCPUKernelStateModel().getParticleData();

        const auto size = data->size();
        const auto &context = cpuKernel->context();
        const auto dt = timeStep();

        auto worker = [&context, data, dt](std::size_t, std::size_t beginIdx, auto entry_begin, auto entry_end)  {
            const auto kbt = context.kBT();
            std::size_t idx = beginIdx;
            const auto &box = context.boxSize().data();
            const auto &pbc = context.periodicBoundaryConditions().data();
            for (auto it = entry_begin; it != entry_end; ++it, ++idx) {
                if(!it->deactivated) {
                    const readdy::scalar D = context.particleTypes().diffusionConstantOf(it->type);
                    const auto randomDisplacement = std::sqrt(2. * D * dt) * rnd::normal3<readdy::scalar>(0, 1);
                    const auto deterministicDisplacement = it->force * dt * D / kbt;
                    it->pos += randomDisplacement + deterministicDisplacement;
                    readdy::bcs::fixPosition(it->pos, box, pbc);
                }
            }
        };

        std::vector<readdy::util::thread::joining_future<void>> waitingFutures;
        waitingFutures.reserve(cpuKernel->getNThreads());
        auto &pool  = cpuKernel->pool();
        {
            auto it = data->begin();
            std::vector<std::function<void(std::size_t)>> executables;
            executables.reserve(cpuKernel->getNThreads());

            auto granularity = cpuKernel->getNThreads();
            const std::size_t grainSize = size / granularity;

            std::size_t idx = 0;
            for (std::size_t i = 0; i < granularity-1; ++i) {
                auto itNext = it + grainSize;
                if(it != itNext) {
                    waitingFutures.emplace_back(pool.push(worker, idx, it, itNext));
                }
                it = itNext;
                idx += grainSize;
            }
            if(it != data->end()) {
                waitingFutures.emplace_back(pool.push(worker, idx, it, data->end()));
            }
        }
    }
};


PYBIND11_MODULE (custom_integrator_example, m) {
    py::module::import("readdy");
    py::class_<EBDIntegrator, readdy::model::actions::UserDefinedAction, std::shared_ptr<EBDIntegrator>>(m, "EBDIntegrator").def(py::init<readdy::scalar>());
}
