/**
 * << detailed description >>
 *
 * @file Observables.h
 * @brief << brief description >>
 * @author clonker
 * @date 26.04.16
 */

#ifndef READDY_MAIN_OBSERVABLES_H
#define READDY_MAIN_OBSERVABLES_H

#include <readdy/model/Observable.h>
#include <readdy/model/Vec3.h>
#include <vector>
#include <readdy/common/Types.h>
#include <iostream>
#include <set>

namespace readdy {
    namespace model {

        class KernelContext;

        class ParticlePositionObservable : public Observable<std::vector<Vec3>> {
        public:


            ParticlePositionObservable(Kernel *const kernel, unsigned int stride = 1) : Observable(kernel, stride) { }

            virtual ~ParticlePositionObservable() {
            }

            virtual void evaluate() override;
        };

        class RadialDistributionObservable : public Observable<std::pair<std::vector<double>, std::vector<double>>> {
        public:
            RadialDistributionObservable(Kernel *const kernel, unsigned int stride, std::vector<double> binBorders, unsigned int typeCountFrom, unsigned int typeCountTo, double particleDensity);
            RadialDistributionObservable(Kernel *const kernel, unsigned int stride, std::vector<double> binBorders, const std::string& typeCountFrom, const std::string& typeCountTo, double particleDensity);

            virtual void evaluate() override;

            const std::vector<double> &getBinBorders() const;
            void setBinBorders(const std::vector<double> &binBorders);
        protected:
            std::vector<double> binBorders;
            std::vector<double> counts;
            unsigned int typeCountFrom, typeCountTo;
            double particleDensity;
        };

        class CenterOfMassObservable : public Observable<readdy::model::Vec3> {

        public:
            CenterOfMassObservable(Kernel *const kernel, unsigned int stride, unsigned int particleType);
            CenterOfMassObservable(Kernel *const kernel, unsigned int stride, const std::vector<unsigned int>& particleTypes);
            CenterOfMassObservable(Kernel *const kernel, unsigned int stride, const std::string& particleType);
            CenterOfMassObservable(Kernel *const kernel, unsigned int stride, const std::vector<std::string>& particleType);

            virtual void evaluate() override;

        protected:
            std::set<unsigned int> particleTypes;
        };

        class HistogramAlongAxisObservable : public Observable<std::vector<double>> {

        public:
            HistogramAlongAxisObservable(readdy::model::Kernel *const kernel, unsigned int stride,
                                                     std::vector<double> binBorders, std::set<unsigned int> typesToCount, unsigned int axis);
            HistogramAlongAxisObservable(Kernel *const kernel, unsigned int stride, std::vector<double> binBorders, std::vector<std::string> typesToCount, unsigned int axis);

            virtual void evaluate() = 0;

        protected:
            std::vector<double> binBorders;
            std::set<unsigned int> typesToCount;
            std::set<unsigned int> transformTypes(std::vector<std::string> types, const readdy::model::KernelContext &ctx);
            unsigned int axis;
        };

        class NParticlesObservable : public Observable<long> {

        public:
            NParticlesObservable(Kernel *const kernel, unsigned int stride) : Observable(kernel, stride) { }
            virtual void evaluate() = 0;
        };

        class TestCombinerObservable : public CombinerObservable<std::vector<double>, ParticlePositionObservable, ParticlePositionObservable> {
        public:

            TestCombinerObservable(Kernel *const kernel, ParticlePositionObservable * obs1, ParticlePositionObservable * obs2, unsigned int stride)
                    : CombinerObservable(kernel, obs1, obs2, stride) {
            }

            virtual void evaluate() override;
        };
    }
}

#endif //READDY_MAIN_OBSERVABLES_H
