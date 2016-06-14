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
#include <readdy/model/potentials/Potential.h>
#include <vector>
#include <unordered_set>

#if BOOST_OS_MAC
#include <string>
#endif

namespace readdy {
    namespace model {
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

            void registerOrder2Potential(potentials::Potential &potential, const std::string &type1, const std::string &type2);
            std::vector<potentials::Potential*> getOrder2Potentials(const std::string &type1, const std::string &type2) const;
            std::vector<potentials::Potential*> getOrder2Potentials(const unsigned int type1, const unsigned int type2) const;
            std::unordered_set<std::tuple<unsigned int, unsigned int>> getAllOrder2RegisteredPotentialTypes() const;

            void configure();

            // ctor and dtor
            KernelContext();

            ~KernelContext();

            // move
            KernelContext(KernelContext &&rhs);

            KernelContext &operator=(KernelContext &&rhs);

            // copy
            KernelContext(const KernelContext &rhs);

            KernelContext &operator=(const KernelContext &rhs);

        private:
            struct Impl;
            std::unique_ptr<readdy::model::KernelContext::Impl> pimpl;
            unsigned int getOrCreateTypeId(const std::string &name);
        };

    }
}
#endif //READDY_MAIN_KERNELCONTEXT_H
