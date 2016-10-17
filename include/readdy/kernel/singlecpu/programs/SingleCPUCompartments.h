/**
 * The SingleCPUCompartments program defines compartments via characteristic functions that map from Vec3 to bool. For every compartment
 * one can then define conversions that should take place as soon as a particle enters the compartment. Note that the user is responsible for
 * keeping the compartments disjoint.
 *
 * @file SingleCPUCompartments.cpp
 * @brief SingleCPUCompartments defines compartments via characteristic functions and performs conversions of particle species on these compartments.
 * @author chrisfroe
 * @date 13.10.16
 */

#ifndef READDY_MAIN_SINGLECPUCOMPARTMENTS_H
#define READDY_MAIN_SINGLECPUCOMPARTMENTS_H

#include <readdy/model/programs/Programs.h>
#include <readdy/kernel/singlecpu/SingleCPUKernel.h>

namespace readdy {
namespace kernel {
namespace singlecpu {
namespace programs {

class Compartments : readdy::model::programs::Compartments {
public:
    using compartmentIdx_t = size_t;
    using particleType_t = unsigned int;

    Compartments(SingleCPUKernel const *const kernel);

    virtual void execute() override;

    virtual void registerCompartment(const std::function<bool(const readdy::model::Vec3)> fun) override;

    virtual void registerConversion(compartmentIdx_t compartmentIdx, std::string from, std::string to) override;

    virtual void registerConversion(compartmentIdx_t compartmentIdx, particleType_t from, particleType_t to);

protected:
    SingleCPUKernel const *const kernel;
    std::vector<std::function<bool(readdy::model::Vec3)>> compartments;
    std::unordered_map<compartmentIdx_t, std::unordered_map<particleType_t, particleType_t>> conversions;
};

}
}
}
}
#endif //READDY_MAIN_SINGLECPUCOMPARTMENTS_H
