/**
 * << detailed description >>
 *
 * @file CPUEulerBDIntegrator.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 07.07.16
 */

#include <readdy/kernel/cpu/programs/CPUEulerBDIntegrator.h>
#include <thread>
#include <readdy/model/RandomProvider.h>

namespace readdy {
    namespace kernel {
        namespace cpu {
            namespace programs {
                void CPUEulerBDIntegrator::execute() {
                    const auto &&pd = kernel->getKernelStateModel().getParticleData();
                    const auto size = pd->size();
                    std::vector<std::thread> threads(kernel->getNCores());
                    const std::size_t grainSize = size / kernel->getNCores();

                    const auto &context = kernel->getKernelContext();
                    using _it_vec3_t = std::vector<readdy::model::Vec3>::iterator;
                    using _it_uint_t = std::vector<unsigned int>::iterator;

                    auto worker = [&context](_it_vec3_t pos0, _it_vec3_t posEnd, _it_vec3_t forces0,
                                             _it_uint_t types0) {
                        readdy::model::RandomProvider provider;
                        const auto &fixPos = context.getFixPositionFun();
                        const auto &kbt = context.getKBT();
                        const auto &&dt = context.getTimeStep();
                        for (auto it = pos0; it != posEnd; ++it) {
                            const double D = context.getDiffusionConstant(*types0);
                            const auto randomDisplacement = sqrt(2. * D * dt) * (provider.getNormal3());
                            *it += randomDisplacement;
                            const auto deterministicDisplacement = *forces0 * dt * D / kbt;
                            *it += deterministicDisplacement;
                            fixPos(*it);
                            ++forces0;
                            ++types0;
                        }

                    };
                    auto work_iter = pd->begin_positions();
                    std::size_t pos = 0;
                    for (auto it = std::begin(threads); it != std::end(threads) - 1; ++it) {
                        *it = std::thread(worker, work_iter, work_iter + grainSize, pd->begin_forces() + pos,
                                          pd->begin_types() + pos);
                        work_iter += grainSize;
                        pos += grainSize;
                    }
                    threads.back() = std::thread(worker, work_iter, pd->end_positions(), pd->begin_forces() + pos,
                                                 pd->begin_types() + pos);

                    for (auto &&i : threads) {
                        i.join();
                    }


                }

                CPUEulerBDIntegrator::CPUEulerBDIntegrator(CPUKernel *kernel) : kernel(kernel) { }

            }
        }
    }
}