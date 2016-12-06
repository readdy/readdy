/**
 * @file Compartments.h
 * @brief Header file of CPU program Compartments
 * @author chrisfroe
 * @date 18.10.16
 */

#ifndef READDY_CPUKERNEL_COMPARTMENTS_H
#define READDY_CPUKERNEL_COMPARTMENTS_H

#include <readdy/kernel/singlecpu/programs/SCPUCompartments.h>
#include <readdy/kernel/cpu/Kernel.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace programs {

class Compartments : public readdy::model::programs::Compartments {
public:
    using compartmentIdx_t = size_t;
    using particleType_t = unsigned int;

    Compartments(Kernel const *const kernel);

    virtual void execute() override;

    virtual void registerCompartment(const std::function<bool(const readdy::model::Vec3)> fun) override;

    virtual void registerConversion(compartmentIdx_t compartmentIdx, std::string from, std::string to) override;

    virtual void registerConversion(compartmentIdx_t compartmentIdx, particleType_t from, particleType_t to);

protected:
    Kernel const *const kernel;
    std::vector<std::function<bool(readdy::model::Vec3)>> compartments;
    std::unordered_map<compartmentIdx_t, std::unordered_map<particleType_t, particleType_t>> conversions;
};

}
}
}
}

#endif //READDY_CPUKERNEL_COMPARTMENTS_H
