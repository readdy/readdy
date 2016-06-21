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

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            namespace reactions {
                SingleCPUReactionFactory::SingleCPUReactionFactory(SingleCPUKernel const *const kernel) : kernel(kernel){

                }

                readdy::model::reactions::Conversion *SingleCPUReactionFactory::createConversion(const std::string &name, unsigned int from, unsigned int to, const double &rate) const {
                    return nullptr;
                }

                readdy::model::reactions::Enzymatic *SingleCPUReactionFactory::createEnzymatic(const std::string &name, unsigned int catalyst, unsigned int from, unsigned int to, const double &rate,
                                                                                       const double &eductDistance) const {
                    return nullptr;
                }

                readdy::model::reactions::Fission *SingleCPUReactionFactory::createFission(const std::string &name, unsigned int from, unsigned int to1, unsigned int to2, const double productDistance,
                                                                                   const double &rate) const {
                    return nullptr;
                }

                readdy::model::reactions::Fusion *SingleCPUReactionFactory::createFusion(const std::string &name, unsigned int from1, unsigned int from2, unsigned int to, const double &rate,
                                                                                 const double &eductDistance) const {
                    return nullptr;
                }

                readdy::model::reactions::Death *SingleCPUReactionFactory::createDeath(const std::string &name, unsigned int typeFrom, const double &rate) const {
                    return nullptr;
                }


            }
        }
    }
}