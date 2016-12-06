/**
 * << detailed description >>
 *
 * @file Compartments.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.11.16
 */


#include <readdy/kernel/cpu_dense/programs/CPUDCompartments.h>

namespace readdy {
namespace kernel {
namespace cpu_dense {
namespace programs {

CPUDCompartments::CPUDCompartments(const CPUDKernel *const kernel) : kernel(kernel) {}

void CPUDCompartments::execute() {
    const auto &ctx = kernel->getKernelContext();
    long long idx = 0;
    for(auto& e : *kernel->getKernelStateModel().getParticleData()) {
        for (auto i = 0; i < compartments.size(); ++i) {
            if (compartments[i](e.position())) {
                if (conversions[i].find(e.type) != conversions[i].end()) {
                    e.type = conversions[i][e.type];
                }
            }
        }
        ++idx;
    }
}

void CPUDCompartments::registerCompartment(const std::function<bool(const readdy::model::Vec3)> fun) {
    compartments.push_back(std::move(fun));
}

void CPUDCompartments::registerConversion(compartmentIdx_t compartmentIdx, particleType_t from, particleType_t to) {
    if (compartmentIdx >= compartments.size()) {
        throw std::runtime_error("Given compartment does not exist. Register it first.");
    }
    if (conversions.find(compartmentIdx) == conversions.end()) {
        conversions.emplace(compartmentIdx, std::unordered_map<particleType_t, particleType_t>());
    }
    conversions[compartmentIdx].emplace(from, to);
}

void CPUDCompartments::registerConversion(compartmentIdx_t compartmentIdx, std::string from, std::string to) {
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
