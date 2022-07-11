/********************************************************************
 * Copyright © 2019 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 ********************************************************************/

/**
 * @file TestIntegration.cpp
 * @brief « brief description »
 * @author chrisfroe
 * @date 03.06.19
 */

#include <catch2/catch_test_macros.hpp>
#include <readdy/model/Kernel.h>
#include <readdy/kernel/mpi/MPIKernel.h>
#include <readdy/kernel/singlecpu/SCPUKernel.h>
#include <readdy/model/RandomProvider.h>

namespace rnd = readdy::model::rnd;
namespace rmo = readdy::model::observables;
namespace rkmu = readdy::kernel::mpi::util;

using ParticlePODSet = std::unordered_set<rkmu::ParticlePOD, rkmu::HashPOD>;
struct CompareVec3 {
    bool operator() (const readdy::Vec3 &v1, const readdy::Vec3 &v2) const {return v1.x <  v2.x; }
};

/// Static in contrast to dynamic (and thus stochastic) results like reaction counts,
/// i.e. StaticResults only depends on the current particle configuration
struct StaticResults {
    rmo::Virial::result_type virial;
    rmo::Energy::result_type energy;
    rmo::Positions::result_type positions;
    rmo::Forces::result_type forces;
    rmo::Particles::result_type particles;
    rmo::NParticles::result_type nParticles;
    // rmo::RadialDistribution::result_type radialDistribution; // not yet implemented for MPI
    rmo::HistogramAlongAxis::result_type histogramAlongAxis;

    bool operator==(const StaticResults &other) const {
        // virial, approx the same since order of computation might differ -> rounding errors
        {
            const auto &v1 = virial.data();
            const auto &v2 = other.virial.data();
            for (std::size_t i = 0; i < v1.size(); ++i) {
                if (v1[i] != Approx(v2[i])) {
                    return false;
                }
            }
        }

        // energy, approx the same since order of computation might differ -> rounding errors
        if (energy != Approx(other.energy)) {
            return false;
        }

        // particles, ignore unique ids, ignore order
        {
            ParticlePODSet set1;
            const auto&[types1, ids1, pos1] = particles;

            for (std::size_t i = 0; i < types1.size(); ++i) {
                set1.emplace(pos1[i], types1[i]);
            }

            ParticlePODSet set2;
            const auto&[types2, ids2, pos2] = other.particles;
            for (std::size_t i = 0; i < types2.size(); ++i) {
                set2.emplace(pos2[i], types2[i]);
            }

            if (set1 != set2) {
                return false;
            }
        }

        // positions, ignore order, abuse the ParticlePODSet with a fixed type
        {
            ParticlePODSet set1;
            for (const auto &p : positions) {
                set1.emplace(p, 0);
            }

            ParticlePODSet set2;
            for (const auto &p : other.positions) {
                set2.emplace(p, 0);
            }

            if (set1 != set2) {
                return false;
            }
        }

        // forces, sort by x, then must be approximately same
        {
            std::vector<readdy::Vec3> f1;
            for (const auto &f : forces) {
                f1.push_back(f);
            }
            std::sort(f1.begin(), f1.end(), CompareVec3{});

            std::vector<readdy::Vec3> f2;
            for (const auto &f : other.forces) {
                f2.push_back(f);
            }
            std::sort(f2.begin(), f2.end(), CompareVec3{});

            for (std::size_t i=0; i < f1.size(); ++i) {
                const auto &v1 = f1[i];
                const auto &v2 = f2[i];
                for (std::size_t j = 0; j<3; ++j) {
                    if (not (v1[j] == Approx(v2[j]))) {
                        return false;
                    }
                }
            }
        }

        // nParticles
        {
            if (nParticles.size() != other.nParticles.size()) {
                return false;
            }
            if (not std::equal(nParticles.begin(), nParticles.end(), other.nParticles.begin())) {
                return false;
            }
        }

        // histogramAlongAxis, should be approx same
        {
            if (histogramAlongAxis.size() != other.histogramAlongAxis.size()) {
                return false;
            }
            if (not std::equal(histogramAlongAxis.begin(), histogramAlongAxis.end(), other.histogramAlongAxis.begin(),
                               [](const auto &x1, const auto &x2) -> bool {
                                   return x1 == Approx(x2);
                               })) {
                return false;
            }
        }

        // all good
        return true;
    }
};

readdy::model::Context getContext() {
    readdy::model::Context ctx;
    ctx.boxSize() = {10., 10., 10.};
    ctx.kernelConfiguration().mpi.dx = 4.9;
    ctx.kernelConfiguration().mpi.dy = 4.9;
    ctx.kernelConfiguration().mpi.dz = 4.9;
    ctx.particleTypes().add("A", 0.5);
    ctx.particleTypes().add("B", 1.0);
    ctx.reactions().add("fusion: B +(1) B -> A", 100.);
    ctx.potentials().addHarmonicRepulsion("B", "B", 1., 2.);
    return ctx;
}

StaticResults observeState(const std::vector<readdy::model::Particle> &initParticles, readdy::model::Kernel *kernel) {
    auto virial = kernel->observe().virial(1);
    auto energy = kernel->observe().energy(1);
    auto positions = kernel->observe().positions(1);
    auto forces = kernel->observe().forces(1);
    auto particles = kernel->observe().particles(1);
    auto nParticles = kernel->observe().nParticles(1);
    auto histogramAlongAxis = kernel->observe().histogramAlongAxis(1, {-1, 0.5, 0.5, 1}, {"B"}, 0);

    kernel->actions().addParticles(initParticles)->perform();
    kernel->actions().initializeKernel()->perform();
    kernel->actions().createNeighborList(kernel->context().calculateMaxCutoff())->perform();
    kernel->actions().updateNeighborList()->perform();
    kernel->actions().calculateForces()->perform();

    virial->call(0);
    energy->call(0);
    positions->call(0);
    forces->call(0);
    particles->call(0);
    nParticles->call(0);
    histogramAlongAxis->call(0);

    return {
            .virial = virial->getResult(),
            .energy = energy->getResult(),
            .positions = positions->getResult(),
            .forces = forces->getResult(),
            .particles = particles->getResult(),
            .nParticles = nParticles->getResult(),
            .histogramAlongAxis = histogramAlongAxis->getResult(),
    };
}

TEST_CASE("Integration test state observables compared to SCPU", "[mpi]") {
    auto ctx = getContext();
    double volumeOccupation{0.7};
    auto numberParticles = static_cast<std::size_t>(
            volumeOccupation * ctx.boxVolume() / (
                    4. / 3. * readdy::util::numeric::pi<double>() * std::pow(ctx.calculateMaxCutoff() / 2., 3)
            )
    );

    const auto &box = ctx.boxSize();
    auto idB = ctx.particleTypes().idOf("B");
    std::vector<readdy::model::Particle> initParticles;
    for (std::size_t i = 0; i < numberParticles; ++i) {
        readdy::Vec3 p{rnd::uniform_real() * box[0] - 0.5 * box[0],
                       rnd::uniform_real() * box[1] - 0.5 * box[1],
                       rnd::uniform_real() * box[2] - 0.5 * box[2]};
        initParticles.emplace_back(p, idB);
    }

    readdy::kernel::mpi::MPIKernel mpiKernel(ctx);

    if (not mpiKernel.domain().isIdleRank()) {
        StaticResults mpiResults = observeState(initParticles, &mpiKernel);
        if (mpiKernel.domain().isMasterRank()) {
            readdy::kernel::scpu::SCPUKernel scpuKernel;
            scpuKernel.context() = ctx;
            StaticResults scpuResults = observeState(initParticles, &scpuKernel);
            REQUIRE(mpiResults == scpuResults);
        }
    }
}
