/**
 * << detailed description >>
 *
 * @file SingleCPUReactions.h
 * @brief << brief description >>
 * @author clonker
 * @date 21.06.16
 */

#ifndef READDY_MAIN_SINGLECPUREACTIONS_H
#define READDY_MAIN_SINGLECPUREACTIONS_H

#include <readdy/model/reactions/Conversion.h>
#include <readdy/model/Particle.h>
#include <readdy/model/reactions/Enzymatic.h>
#include <readdy/model/reactions/Fission.h>
#include <readdy/model/reactions/Fusion.h>
#include <readdy/model/RandomProvider.h>
#include <readdy/common/make_unique.h>

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            namespace reactions {
                class SingleCPUConversion : public readdy::model::reactions::Conversion {

                public:
                    SingleCPUConversion(const std::string &name, unsigned int typeFrom, unsigned int typeTo, const double &rate) : Conversion(name, typeFrom, typeTo, rate) { }

                    virtual void perform(const readdy::model::Particle &p1_in, const readdy::model::Particle &p2_in, readdy::model::Particle &p1_out, readdy::model::Particle &p2_out) const override;
                };

                class SingleCPUEnzymatic : public readdy::model::reactions::Enzymatic {

                public:
                    SingleCPUEnzymatic(const std::string &name, unsigned int catalyst, unsigned int from, unsigned int to, const double &rate, const double &eductDistance)
                            : Enzymatic(name, catalyst, from, to, rate, eductDistance) { }

                    virtual void perform(const readdy::model::Particle &p1_in, const readdy::model::Particle &p2_in, readdy::model::Particle &p1_out, readdy::model::Particle &p2_out) const override;
                };

                class SingleCPUFission : public readdy::model::reactions::Fission {

                public:
                    SingleCPUFission(const std::string &name, unsigned int from, unsigned int to1, unsigned int to2, const double productDistance, const double &rate, const double& weight1, const double &weight2)
                            : Fission(name, from, to1, to2, productDistance, rate, weight1, weight2) { }

                    virtual void perform(const readdy::model::Particle &p1_in, const readdy::model::Particle &p2_in, readdy::model::Particle &p1_out, readdy::model::Particle &p2_out) const override;

                protected:
                    std::unique_ptr<readdy::model::RandomProvider> rand = std::make_unique<readdy::model::RandomProvider>();
                };

                class SingleCPUFusion : public readdy::model::reactions::Fusion {

                public:
                    SingleCPUFusion(const std::string &name, unsigned int from1, unsigned int from2, unsigned int to, const double &rate, const double &eductDistance, const double &weight1 = 0.5, const double &weight2 = 0.5)
                            : Fusion(name, from1, from2, to, rate, eductDistance, weight1, weight2) { }

                    virtual void perform(const readdy::model::Particle &p1_in, const readdy::model::Particle &p2_in, readdy::model::Particle &p1_out, readdy::model::Particle &p2_out) const override;

                };

            }
        }
    }
}


#endif //READDY_MAIN_SINGLECPUREACTIONS_H
