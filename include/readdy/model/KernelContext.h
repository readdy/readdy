/**
 * todo
 *
 * @file KernelContext.h
 * @brief Container class for time independent information of the KernelContext.
 * @author clonker
 * @date 18.04.16
 * @todo write docs, is kbt really time indep (or should be treated as such)?
 */

#ifndef READDY_MAIN_KERNELCONTEXT_H
#define READDY_MAIN_KERNELCONTEXT_H

#include <array>
#include <memory>
#include <boost/predef.h>
#include <vector>
#include <unordered_set>
#include <readdy/model/potentials/PotentialOrder1.h>
#include <readdy/model/potentials/PotentialOrder2.h>

#if BOOST_OS_MAC
#include <string>
#endif

namespace readdy {
    namespace model {
        class ParticleTypePairHasher;
        class KernelContext {
        public:
            double getKBT() const;

            void setKBT(double kBT);

            std::array<double, 3>& getBoxSize() const;

            void setBoxSize(double dx, double dy, double dz);

            std::array<bool, 3>& getPeriodicBoundary() const;
            unsigned int getParticleTypeID(const std::string& name) const;
            void setPeriodicBoundary(bool pb_x, bool pb_y, bool pb_z);

            double getDiffusionConstant(const std::string& particleType) const;
            double getDiffusionConstant(const unsigned int particleType) const;
            void setDiffusionConstant(const std::string& particleType, double D);

            double getParticleRadius(const std::string& type) const;
            double getParticleRadius(const unsigned int& type) const;
            void setParticleRadius(const std::string& type, const double r);

            double getTimeStep() const;
            void setTimeStep(double dt);

            void deregisterPotential(const boost::uuids::uuid &potential);

            const boost::uuids::uuid& registerOrder1Potential(potentials::PotentialOrder1 const* const potential, const std::string &type);
            const std::vector<std::unique_ptr<potentials::PotentialOrder1>>& getOrder1Potentials(const std::string& type) const;
            const std::vector<std::unique_ptr<potentials::PotentialOrder1>>& getOrder1Potentials(const unsigned int type) const;
            std::unordered_set<unsigned int> getAllOrder1RegisteredPotentialTypes() const;

            const boost::uuids::uuid& registerOrder2Potential(potentials::PotentialOrder2 const* const potential, const std::string &type1, const std::string &type2);
            const std::vector<std::unique_ptr<potentials::PotentialOrder2>>& getOrder2Potentials(const std::string &type1, const std::string &type2) const;
            const std::vector<std::unique_ptr<potentials::PotentialOrder2>>& getOrder2Potentials(const unsigned int type1, const unsigned int type2) const;
            std::unordered_set<std::tuple<unsigned int, unsigned int>, readdy::model::ParticleTypePairHasher> getAllOrder2RegisteredPotentialTypes() const;

            // ctor and dtor
            KernelContext();

            ~KernelContext();

            // move
            KernelContext(KernelContext &&rhs);

            KernelContext &operator=(KernelContext &&rhs);

            // copy
            KernelContext(const KernelContext &rhs) = delete;

            KernelContext &operator=(const KernelContext &rhs) = delete;

        private:
            struct Impl;
            std::unique_ptr<readdy::model::KernelContext::Impl> pimpl;
            unsigned int getOrCreateTypeId(const std::string &name);
        };

    }
}
#endif //READDY_MAIN_KERNELCONTEXT_H
