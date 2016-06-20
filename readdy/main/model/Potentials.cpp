/**
 * << detailed description >>
 *
 * @file Potentials.cpp.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 20.06.16
 */

#include <readdy/model/Kernel.h>
#include <readdy/model/potentials/PotentialsOrder1.h>
#include <readdy/model/potentials/PotentialsOrder2.h>

namespace readdy {
    namespace model {
        namespace potentials {

            /////////////////////////////////////////////////////////////////////////////
            //
            // Potentials order 1
            //
            /////////////////////////////////////////////////////////////////////////////

            /**
             * Cube potential
             */

            CubePotential::CubePotential(const Kernel *const kernel) : PotentialOrder1(getPotentialName<CubePotential>()), kernel(kernel) { }

            void CubePotential::configureForType(const unsigned int &type) {
                particleRadius = kernel->getKernelContext().getParticleRadius(type);
                for (auto i = 0; i < 3; i++) {
                    if (origin[i] < origin[i] + extent[i]) {
                        min[i] = origin[i];
                        max[i] = origin[i] + extent[i];
                    } else {
                        min[i] = origin[i] + extent[i];
                        max[i] = origin[i];
                    }
                }
            }
            const Vec3 &CubePotential::getOrigin() const { return origin; }
            void CubePotential::setOrigin(const Vec3 &origin) { CubePotential::origin = origin; }
            const Vec3 &CubePotential::getExtent() const { return extent; }
            void CubePotential::setExtent(const Vec3 &extent) { CubePotential::extent = extent; }
            double CubePotential::getForceConstant() const { return forceConstant; }
            void CubePotential::setForceConstant(double forceConstant) { CubePotential::forceConstant = forceConstant; }
            bool CubePotential::isConsiderParticleRadius() const { return considerParticleRadius; }
            void CubePotential::setConsiderParticleRadius(bool considerParticleRadius) { CubePotential::considerParticleRadius = considerParticleRadius; }
            double CubePotential::getParticleRadius() const { return particleRadius; }




            /////////////////////////////////////////////////////////////////////////////
            //
            // Potentials order 2
            //
            /////////////////////////////////////////////////////////////////////////////

            /*
             * Harmonic repulsion
             */

            void HarmonicRepulsion::configureForTypes(unsigned int type1, unsigned int type2) {
                auto r1 = kernel->getKernelContext().getParticleRadius(type1);
                auto r2 = kernel->getKernelContext().getParticleRadius(type2);
                sumOfParticleRadii = r1 + r2;
                sumOfParticleRadiiSquared = sumOfParticleRadii * sumOfParticleRadii;
            }
            HarmonicRepulsion::HarmonicRepulsion(const Kernel *const kernel) : PotentialOrder2(_internal::PotentialName<HarmonicRepulsion>::value), kernel(kernel) { }
            double HarmonicRepulsion::getSumOfParticleRadii() const {
                return sumOfParticleRadii;
            }
            double HarmonicRepulsion::getSumOfParticleRadiiSquared() const {
                return sumOfParticleRadiiSquared;
            }
            double HarmonicRepulsion::getForceConstant() const {
                return forceConstant;
            }
            void HarmonicRepulsion::setForceConstant(double forceConstant) {
                HarmonicRepulsion::forceConstant = forceConstant;
            }


        }
    }
}