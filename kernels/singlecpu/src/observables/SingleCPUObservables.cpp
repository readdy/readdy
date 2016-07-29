/**
 * << detailed description >>
 *
 * @file SingleCPUObservables.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 30.06.16
 */
#include <readdy/kernel/singlecpu/SingleCPUKernel.h>
#include <readdy/kernel/singlecpu/observables/SingleCPUObservables.h>

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            namespace observables {
                SingleCPUHistogramAlongAxisObservable::SingleCPUHistogramAlongAxisObservable(
                        readdy::model::Kernel *const kernel, unsigned int stride, const std::vector<double> &binBorders,
                        const std::vector<std::string> &typesToCount, unsigned int axis)
                        : readdy::model::HistogramAlongAxisObservable(kernel, stride, binBorders, typesToCount, axis),
                          singleCPUKernel(dynamic_cast<SingleCPUKernel *>(kernel)) {
                    size = result.size();
                }

                void SingleCPUHistogramAlongAxisObservable::evaluate() {
                    std::fill(result.begin(), result.end(), 0);

                    const auto &model = singleCPUKernel->getKernelStateModel();
                    const auto data = model.getParticleData();

                    auto it_pos = data->begin_positions();
                    auto it_types = data->begin_types();

                    while (it_pos != data->end_positions()) {
                        if (typesToCount.find(*it_types) != typesToCount.end()) {
                            const auto &vec = *it_pos;
                            auto upperBound = std::upper_bound(binBorders.begin(), binBorders.end(), vec[axis]);
                            if (upperBound != binBorders.end()) {
                                unsigned long binBordersIdx = upperBound - binBorders.begin();
                                if (binBordersIdx > 1) {
                                    ++result[binBordersIdx - 1];
                                }
                            }
                        }
                        ++it_pos;
                        ++it_types;
                    }
                }


            }


        }
    }
}


