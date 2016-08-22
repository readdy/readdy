/**
 * << detailed description >>
 *
 * @file SingleCPUReactionFactory.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 21.06.16
 */

#include <readdy/kernel/singlecpu/SingleCPUKernel.h>
#include <readdy/kernel/singlecpu/reactions/SingleCPUReactionFactory.h>
#include <readdy/kernel/singlecpu/reactions/SingleCPUReactions.h>

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            namespace reactions {
                SingleCPUReactionFactory::SingleCPUReactionFactory(SingleCPUKernel const *const kernel) : kernel(
                        kernel) {

                }

                readdy::model::reactions::Conversion *
                SingleCPUReactionFactory::createConversion(const std::string &name, unsigned int from, unsigned int to,
                                                           const double &rate) const {
                    return new Conversion(name, from, to, rate);
                }

                readdy::model::reactions::Enzymatic *
                SingleCPUReactionFactory::createEnzymatic(const std::string &name, unsigned int catalyst,
                                                          unsigned int from, unsigned int to, const double &rate,
                                                          const double &eductDistance) const {
                    return new Enzymatic(name, catalyst, from, to, rate, eductDistance);
                }

                readdy::model::reactions::Fission *
                SingleCPUReactionFactory::createFission(const std::string &name, unsigned int from, unsigned int to1,
                                                        unsigned int to2,
                                                        const double &rate, const double productDistance,
                                                        const double &weight1, const double &weight2) const {
                    return new Fission(name, from, to1, to2, rate, productDistance, weight1, weight2);
                }

                readdy::model::reactions::Fusion *
                SingleCPUReactionFactory::createFusion(const std::string &name, unsigned int from1, unsigned int from2,
                                                       unsigned int to, const double &rate,
                                                       const double &eductDistance, const double &weight1,
                                                       const double &weight2) const {
                    return new Fusion(name, from1, from2, to, rate, eductDistance, weight1, weight2);
                }

            }
        }
    }
}