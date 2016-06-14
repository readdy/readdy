/**
 * << detailed description >>
 *
 * @file KernelContext.cpp
 * @brief Implementation file of the KernelContext.
 * @author clonker
 * @date 18.04.16
 * @todo make proper reference to KernelContext.h, is kBT really indepdendent of t?
 */

#include <readdy/model/KernelContext.h>
#include <readdy/common/make_unique.h>
#include <unordered_map>
#include <boost/log/trivial.hpp>
#include <readdy/common/Utils.h>
#include <readdy/model/_internal/ParticleTypePair.h>

namespace readdy {
    namespace model {
        struct ParticleTypePairHasher {
            std::size_t operator()(const readdy::model::_internal::ParticleTypePair &k) const {
                return hash_value(k);
            }
        };

        struct readdy::model::KernelContext::Impl {
            uint typeCounter;
            std::unordered_map<std::string, uint> typeMapping;
            double kBT = 0;
            std::array<double, 3> box_size{};
            std::array<bool, 3> periodic_boundary{};
            std::unordered_map<uint, double> diffusionConstants{};
            std::unordered_map<uint, double> particleRadii{};
            std::unordered_map<_internal::ParticleTypePair, std::vector<potentials::Potential *>, readdy::model::ParticleTypePairHasher> potentialRegistry{};
            double timeStep;
        };


        double KernelContext::getKBT() const {
            return (*pimpl).kBT;
        }

        void KernelContext::setKBT(double kBT) {
            (*pimpl).kBT = kBT;
        }

        void KernelContext::setBoxSize(double dx, double dy, double dz) {
            (*pimpl).box_size = {dx, dy, dz};
        }

        void KernelContext::setPeriodicBoundary(bool pb_x, bool pb_y, bool pb_z) {
            (*pimpl).periodic_boundary = {pb_x, pb_y, pb_z};
        }

        KernelContext::KernelContext() : pimpl(std::make_unique<KernelContext::Impl>()) { }

        std::array<double, 3> &KernelContext::getBoxSize() const {
            return pimpl->box_size;
        }

        std::array<bool, 3> &KernelContext::getPeriodicBoundary() const {
            return pimpl->periodic_boundary;
        }

        KernelContext::KernelContext(const KernelContext &rhs) : pimpl(std::make_unique<KernelContext::Impl>(*rhs.pimpl)) { }

        KernelContext &KernelContext::operator=(const KernelContext &rhs) {
            *pimpl = *rhs.pimpl;
            return *this;
        }

        double KernelContext::getDiffusionConstant(const std::string &particleType) const {
            return pimpl->diffusionConstants[pimpl->typeMapping[particleType]];
        }

        void KernelContext::setDiffusionConstant(const std::string &particleType, double D) {
            pimpl->diffusionConstants[getOrCreateTypeId(particleType)] = D;
        }

        double KernelContext::getTimeStep() const {
            return pimpl->timeStep;
        }

        void KernelContext::setTimeStep(double dt) {
            pimpl->timeStep = dt;
        }

        unsigned int KernelContext::getParticleTypeID(const std::string &name) const {
            return pimpl->typeMapping[name];
        }

        double KernelContext::getDiffusionConstant(uint particleType) const {
            return pimpl->diffusionConstants[particleType];
        }

        double KernelContext::getParticleRadius(const std::string &type) const {
            return getParticleRadius(pimpl->typeMapping[type]);
        }

        double KernelContext::getParticleRadius(const unsigned int &type) const {
            return pimpl->particleRadii[type];
        }

        void KernelContext::setParticleRadius(const std::string &particleType, const double r) {
            pimpl->particleRadii[getOrCreateTypeId(particleType)] = r;
        }

        unsigned int KernelContext::getOrCreateTypeId(const std::string &particleType) {
            uint t_id;
            if (pimpl->typeMapping.find(particleType) != pimpl->typeMapping.end()) {
                t_id = pimpl->typeMapping[particleType];
            } else {
                t_id = ++(pimpl->typeCounter);
                pimpl->typeMapping.emplace(particleType, t_id);
            }
            return t_id;
        }

        void KernelContext::registerOrder2Potential(potentials::Potential &potential, const std::string &type1, const std::string &type2) {
            // wlog: type1 <= type2
            auto type1Id = pimpl->typeMapping[type1];
            auto type2Id = pimpl->typeMapping[type2];
            _internal::ParticleTypePair pp {type1Id, type2Id};
            if(pimpl->potentialRegistry.find(pp) == pimpl->potentialRegistry.end()) {
                pimpl->potentialRegistry.emplace(pp, std::vector<potentials::Potential*>());
            }
            pimpl->potentialRegistry[pp].push_back(&potential);
        }

        std::vector<potentials::Potential *> KernelContext::getOrder2Potentials(const std::string &type1, const std::string &type2) const {
            _internal::ParticleTypePair pp {pimpl->typeMapping[type1], pimpl->typeMapping[type2]};
            return pimpl->potentialRegistry[pp];
        }

        std::vector<potentials::Potential *> KernelContext::getOrder2Potentials(const unsigned int type1, const unsigned int type2) const {
            return pimpl->potentialRegistry[{type1, type2}];
        }

        std::unordered_set<std::tuple<unsigned int, unsigned int>> KernelContext::getAllOrder2RegisteredPotentialTypes() const {
            std::unordered_set<std::tuple<unsigned int, unsigned int>> result {};
            for(auto it = pimpl->potentialRegistry.begin(); it != pimpl->potentialRegistry.end(); ++it) {
                result.insert(std::make_tuple(it->first.t1, it->first.t2));
            }
            return result;
        }

        void KernelContext::configure() {
            // TODO implement this: Configure observables (computable params)
        }


        KernelContext &KernelContext::operator=(KernelContext &&rhs) = default;

        KernelContext::KernelContext(KernelContext &&rhs) = default;

        KernelContext::~KernelContext() = default;
    }
}






