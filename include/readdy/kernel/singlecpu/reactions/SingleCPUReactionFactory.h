/**
 * << detailed description >>
 *
 * @file SingleCPUReactionFactory.h
 * @brief << brief description >>
 * @author clonker
 * @date 21.06.16
 */

#include <readdy/model/reactions/ReactionFactory.h>

#ifndef READDY_MAIN_SINGLECPUREACTIONFACTORY_H
#define READDY_MAIN_SINGLECPUREACTIONFACTORY_H
namespace readdy {
    namespace kernel {
        namespace singlecpu {
            class SingleCPUKernel;
            namespace reactions {
                class SingleCPUReactionFactory : public readdy::model::reactions::ReactionFactory {

                public:
                    SingleCPUReactionFactory(SingleCPUKernel const *const kernel);

                protected:
                    virtual readdy::model::reactions::Conversion *
                    createConversion(const std::string &name, unsigned int from, unsigned int to,
                                     const double rate) const override;

                    virtual readdy::model::reactions::Enzymatic *
                    createEnzymatic(const std::string &name, unsigned int catalyst, unsigned int from, unsigned int to,
                                    const double rate, const double eductDistance) const override;

                    virtual readdy::model::reactions::Fission *
                    createFission(const std::string &name, unsigned int from, unsigned int to1, unsigned int to2,
                                  const double rate, const double productDistance, const double weight1 = 0.5,
                                  const double weight2 = 0.5) const override;

                    virtual readdy::model::reactions::Fusion *
                    createFusion(const std::string &name, unsigned int from1, unsigned int from2, unsigned int to,
                                 const double rate, const double eductDistance, const double weight1 = 0.5,
                                 const double weight2 = 0.5) const override;

                    SingleCPUKernel const *const kernel;
                };
            }
        }
    }
}
#endif //READDY_MAIN_SINGLECPUREACTIONFACTORY_H
