/**
 * @file Compartments.cpp
 * @brief Implementation of SingleCPU program Compartments
 * @author chrisfroe
 * @date 13.10.16
 */

#include <readdy/kernel/singlecpu/programs/SCPUCompartments.h>

namespace readdy {
namespace kernel {
namespace scpu {
namespace programs {

SCPUCompartments::SCPUCompartments(SCPUKernel const *const kernel) : kernel(kernel) {}

void SCPUCompartments::execute() {
    const auto &ctx = kernel->getKernelContext();
    auto data = kernel->getKernelStateModel().getParticleData();
    auto posIt = data->begin_positions();
    auto typesIt = data->begin_types();
    while (posIt != data->end_positions()) {
        for (auto i=0; i<compartments.size(); ++i) {
            if (compartments[i](*posIt)) {
                if (conversions[i].find(*typesIt) != conversions[i].end()) {
                    const auto targetType = conversions[i][*typesIt];
                    *typesIt = targetType;
                }
            }
        }
        ++posIt;
        ++typesIt;
    }
}

void SCPUCompartments::registerCompartment(const std::function<bool(const readdy::model::Vec3)> fun) {
    compartments.push_back(std::move(fun));
}

void SCPUCompartments::registerConversion(compartmentIdx_t compartmentIdx, particleType_t from, particleType_t to) {
    if (compartmentIdx >= compartments.size()) {
        throw std::runtime_error("Given compartment does not exist. Register it first.");
    }
    if (conversions.find(compartmentIdx) == conversions.end()) {
        conversions.emplace(compartmentIdx, std::unordered_map<particleType_t, particleType_t>());
    }
    conversions[compartmentIdx].emplace(from, to);
}

void SCPUCompartments::registerConversion(compartmentIdx_t compartmentIdx, std::string from, std::string to) {
    // Since this program is not part of the default readdy functionality it shall not be able to
    // create particleTypes, i.e. if 'from' or 'to' do not exist the conversion cannot be registered
    const auto typeMapping = kernel->getKernelContext().getTypeMapping();
    auto findFrom = typeMapping.find(from);
    auto findTo = typeMapping.find(to);
    if (findFrom == typeMapping.end() || findTo == typeMapping.end()) {
        throw readdy::model::UnknownParticleType("Particle type is unknown to context (and shall not be registered)");
    }
    registerConversion(compartmentIdx, findFrom->second, findTo->second);
}

}
}
}
}

