/**
 * << detailed description >>
 *
 * @file CPUStateModel.h
 * @brief << brief description >>
 * @author clonker
 * @date 13.07.16
 */

#ifndef READDY_CPUKERNEL_CPUSTATEMODEL_H
#define READDY_CPUKERNEL_CPUSTATEMODEL_H


#include <readdy/model/KernelStateModel.h>
#include <readdy/model/KernelContext.h>
#include <readdy/kernel/cpu/model/ParticleIndexPair.h>
#include <readdy/kernel/singlecpu/model/ParticleData.h>
#include <readdy/kernel/cpu/model/NeighborList.h>
#include <readdy/kernel/cpu/util/Config.h>

namespace readdy {
namespace kernel {
namespace cpu {
class CPUStateModel : public readdy::model::KernelStateModel {

public:
    CPUStateModel(readdy::model::KernelContext *const context, util::Config const *const config);

    ~CPUStateModel();

    virtual const std::vector<readdy::model::Vec3> getParticlePositions() const override;

    virtual const std::vector<readdy::model::Particle> getParticles() const override;

    virtual void updateNeighborList() override;

    virtual void calculateForces() override;

    virtual void addParticle(const readdy::model::Particle &p) override;

    virtual void addParticles(const std::vector<readdy::model::Particle> &p) override;

    virtual void removeParticle(const readdy::model::Particle &p) override;

    virtual void removeAllParticles() override;

    virtual double getEnergy() const override;

    readdy::kernel::singlecpu::model::ParticleData *const getParticleData() const;

    model::NeighborList *const getNeighborList() const;

    virtual void clearNeighborList() override;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl;
    util::Config const *const config;
};
}
}
}

#endif //READDY_CPUKERNEL_CPUSTATEMODEL_H
