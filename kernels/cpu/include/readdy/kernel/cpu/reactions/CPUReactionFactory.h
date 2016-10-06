/**
 * << detailed description >>
 *
 * @file CPUReactionFactory.h
 * @brief << brief description >>
 * @author clonker
 * @date 13.07.16
 */


#ifndef READDY_MAIN_CPUREACTIONFACTORY_H
#define READDY_MAIN_CPUREACTIONFACTORY_H

#include <readdy/model/reactions/ReactionFactory.h>
#include <readdy/model/RandomProvider.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace reactions {

struct Conversion : public readdy::model::reactions::Conversion {

    Conversion(const std::string &name, unsigned int typeFrom, unsigned int typeTo, const double rate)
            : readdy::model::reactions::Conversion(name, typeFrom, typeTo, rate) {}

    virtual void perform(const readdy::model::Particle &p1_in, const readdy::model::Particle &p2_in,
                         readdy::model::Particle &p1_out, readdy::model::Particle &p2_out,
                         const rnd_ptr &rnd) const override {
        p1_out.setPos(p1_in.getPos());
        p1_out.setType(getTypeTo());
        p1_out.setId(p1_in.getId());
    }

};

struct Enzymatic : public readdy::model::reactions::Enzymatic {

    Enzymatic(const std::string &name, unsigned int catalyst, unsigned int from, unsigned int to,
              const double rate, const double eductDistance) : readdy::model::reactions::Enzymatic(
            name, catalyst, from, to, rate, eductDistance) {

    }

    virtual void perform(const readdy::model::Particle &p1_in, const readdy::model::Particle &p2_in,
                         readdy::model::Particle &p1_out, readdy::model::Particle &p2_out,
                         const rnd_ptr &rnd) const override {
        if (p1_in.getType() == getCatalyst()) {
            // p1 is the catalyst
            p1_out.setType(getCatalyst());
            p1_out.setPos(p1_in.getPos());
            p1_out.setId(p1_in.getId());
            p2_out.setType(getTo());
            p2_out.setPos(p2_in.getPos());
        } else {
            // p2 is the catalyst
            p1_out.setType(getCatalyst());
            p1_out.setPos(p2_in.getPos());
            p1_out.setId(p2_in.getId());
            p2_out.setType(getTo());
            p2_out.setPos(p1_in.getPos());
        }
    }

};

struct Fission : public readdy::model::reactions::Fission {

    Fission(const std::string &name, unsigned int from, unsigned int to1, unsigned int to2,
            const double rate, const double productDistance, const double weight1,
            const double weight2) : readdy::model::reactions::Fission(name, from, to1, to2,
                                                                      rate, productDistance, weight1,
                                                                      weight2) {}

    virtual void perform(const readdy::model::Particle &p1_in, const readdy::model::Particle &p2_in,
                         readdy::model::Particle &p1_out, readdy::model::Particle &p2_out,
                         const rnd_ptr &rnd) const override {
        // as long as the orientation is uniform, it does not matter of which type p1_in and p2_in are.
        auto n3 = rnd->getNormal3();
        n3 /= sqrt(n3 * n3);
        p1_out.setType(getTo1());
        p1_out.setPos(p1_in.getPos() + getWeight1() * getProductDistance() * n3);

        p2_out.setType(getTo2());
        p2_out.setPos(p1_in.getPos() - getWeight2() * getProductDistance() * n3);
    }

};

struct Fusion : public readdy::model::reactions::Fusion {

    Fusion(const std::string &name, unsigned int from1, unsigned int from2, unsigned int to,
           const double rate, const double eductDistance, const double weight1,
           const double weight2) : readdy::model::reactions::Fusion(name, from1, from2, to, rate,
                                                                    eductDistance, weight1,
                                                                    weight2) {}

    virtual void perform(const readdy::model::Particle &p1_in, const readdy::model::Particle &p2_in,
                         readdy::model::Particle &p1_out, readdy::model::Particle &p2_out,
                         const rnd_ptr &rnd) const override {
        p1_out.setType(getTo());
        if (getFrom1() == p1_in.getType()) {
            p1_out.setPos(p1_in.getPos() + getWeight1() * (p2_in.getPos() - p1_in.getPos()));
        } else {
            p1_out.setPos(p2_in.getPos() + getWeight1() * (p1_in.getPos() - p2_in.getPos()));
        }
    }

};

class CPUReactionFactory : public readdy::model::reactions::ReactionFactory {
public:
    CPUReactionFactory(CPUKernel *const kernel) : kernel(kernel) {

    }

protected:
    virtual readdy::model::reactions::Conversion *createConversion(const std::string &name,
                                                                   unsigned int from,
                                                                   unsigned int to,
                                                                   const double rate) const override {
        return new Conversion(name, from, to, rate);
    }

    virtual readdy::model::reactions::Enzymatic *createEnzymatic(const std::string &name,
                                                                 unsigned int catalyst,
                                                                 unsigned int from, unsigned int to,
                                                                 const double rate,
                                                                 const double eductDistance) const override {
        return new Enzymatic(name, catalyst, from, to, rate, eductDistance);
    }

    virtual readdy::model::reactions::Fission *createFission(const std::string &name, unsigned int from,
                                                             unsigned int to1, unsigned int to2,
                                                             const double rate,
                                                             const double productDistance,
                                                             const double weight1,
                                                             const double weight2) const override {
        return new Fission(name, from, to1, to2, rate, productDistance, weight1, weight2);
    }

    virtual readdy::model::reactions::Fusion *createFusion(const std::string &name, unsigned int from1,
                                                           unsigned int from2, unsigned int to,
                                                           const double rate,
                                                           const double eductDistance,
                                                           const double weight1,
                                                           const double weight2) const override {
        return new Fusion(name, from1, from2, to, rate, eductDistance, weight1,
                          weight2);
    }


protected:
    CPUKernel *const kernel;
};


}
}
}
}
#endif //READDY_MAIN_CPUREACTIONFACTORY_H
