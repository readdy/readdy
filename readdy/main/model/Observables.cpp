/********************************************************************
 * Copyright © 2016 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * This file is part of ReaDDy.                                     *
 *                                                                  *
 * ReaDDy is free software: you can redistribute it and/or modify   *
 * it under the terms of the GNU Lesser General Public License as   *
 * published by the Free Software Foundation, either version 3 of   *
 * the License, or (at your option) any later version.              *
 *                                                                  *
 * This program is distributed in the hope that it will be useful,  *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of   *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
 * GNU Lesser General Public License for more details.              *
 *                                                                  *
 * You should have received a copy of the GNU Lesser General        *
 * Public License along with this program. If not, see              *
 * <http://www.gnu.org/licenses/>.                                  *
 ********************************************************************/


/**
 * Core library implementation of observables. Since every kernel usually has its own implementation, there are mostly constructors here.
 *
 * @file Observables.cpp
 * @brief Implementation of observables
 * @author clonker
 * @author chrisfroe
 * @date 26.04.16
 */

#include <readdy/model/observables/Observables.h>
#include <readdy/model/Kernel.h>
#include <readdy/model/_internal/Util.h>
#include <readdy/common/numeric.h>
#include <readdy/io/DataSet.h>

const std::string OBSERVABLES_GROUP_PATH = "/readdy/observables";

namespace readdy {
namespace model {
namespace observables {

struct Vec3POD {
    using entry_t = readdy::model::Vec3::entry_t;

    Vec3POD(Vec3::entry_t x, Vec3::entry_t y, Vec3::entry_t z) : x(x), y(y), z(z) {}

    Vec3POD(const Vec3 &v) : x(v[0]), y(v[1]), z(v[2]) {}

    entry_t x, y, z;

    friend std::ostream &operator<<(std::ostream &os, const Vec3POD &v) {
        os << "Vec3POD(" << v.x << ", " << v.y << ", " << v.z << ")";
        return os;
    }
};

class Vec3MemoryType : public readdy::io::DataSetType {
public:
    Vec3MemoryType() {
        using entry_t = Vec3POD;
        tid = []() -> hid_t {
            hid_t entryTypeMemory = H5Tcreate(H5T_COMPOUND, sizeof(entry_t));
            // init vec pod
            io::NativeDataSetType<entry_t::entry_t> posType{};
            H5Tinsert(entryTypeMemory, "x", HOFFSET(entry_t, x), posType.tid);
            H5Tinsert(entryTypeMemory, "y", HOFFSET(entry_t, y), posType.tid);
            H5Tinsert(entryTypeMemory, "z", HOFFSET(entry_t, z), posType.tid);
            return entryTypeMemory;
        }();
    }
};

class Vec3FileType : public readdy::io::DataSetType {
public:
    Vec3FileType() {
        using entry_t = Vec3POD;
        tid = []() -> hid_t {
            std::size_t size = 3 * sizeof(entry_t::entry_t);
            hid_t entryTypeFile = H5Tcreate(H5T_COMPOUND, size);
            // data types
            io::STDDataSetType<entry_t::entry_t> posType{};
            // init trajectory pod
            std::size_t cumsize = 0;
            H5Tinsert(entryTypeFile, "x", cumsize, posType.tid);
            cumsize += sizeof(entry_t::entry_t);
            H5Tinsert(entryTypeFile, "y", cumsize, posType.tid);
            cumsize += sizeof(entry_t::entry_t);
            H5Tinsert(entryTypeFile, "z", cumsize, posType.tid);
            return entryTypeFile;
        }();
    }
};


struct Positions::Impl {
    using writer_t = io::DataSet<Vec3POD, true>;
    std::unique_ptr<writer_t> writer;
};

Positions::Positions(Kernel *const kernel, unsigned int stride,
                     std::vector<std::string> typesToCount) :
        Positions(kernel, stride,
                  _internal::util::transformTypes2(typesToCount, kernel->getKernelContext())) {}

Positions::Positions(Kernel *const kernel, unsigned int stride,
                     std::vector<unsigned int> typesToCount) :
        Observable(kernel, stride), typesToCount(typesToCount), pimpl(std::make_unique<Impl>()) {}

void Positions::append() {
    std::vector<Vec3POD> podVec(result.begin(), result.end());
    pimpl->writer->append({1}, &podVec);
}

Positions::Positions(Kernel *const kernel, unsigned int stride) : Observable(kernel, stride) {}

void Positions::initializeDataSet(io::File &file, const std::string &dataSetName, unsigned int flushStride) {
    if (!pimpl->writer) {
        std::vector<readdy::io::h5::dims_t> fs = {flushStride};
        std::vector<readdy::io::h5::dims_t> dims = {readdy::io::h5::UNLIMITED_DIMS};
        auto dataSet = std::make_unique<io::DataSet<Vec3POD, true>>(
                dataSetName, file.createGroup(OBSERVABLES_GROUP_PATH), fs, dims,
                Vec3MemoryType(), Vec3FileType()
        );
        pimpl->writer = std::move(dataSet);
    }
}

void Positions::flush() {
    pimpl->writer->flush();
}

Positions::~Positions() = default;

struct RadialDistribution::Impl {
    using writer_t = io::DataSet<double, false>;
    std::unique_ptr<writer_t> writerRadialDistribution;
};

RadialDistribution::RadialDistribution(Kernel *const kernel, unsigned int stride,
                                       std::vector<double> binBorders, std::vector<unsigned int> typeCountFrom,
                                       std::vector<unsigned int> typeCountTo, double particleToDensity)
        : Observable(kernel, stride), typeCountFrom(typeCountFrom), typeCountTo(typeCountTo),
          particleToDensity(particleToDensity), pimpl(std::make_unique<Impl>()) {
    setBinBorders(binBorders);
}

void RadialDistribution::evaluate() {
    if (binBorders.size() > 1) {
        std::fill(counts.begin(), counts.end(), 0);
        const auto particles = kernel->getKernelStateModel().getParticles();
        auto isInCollection = [](const readdy::model::Particle &p, const std::vector<unsigned int> &collection) {
            return std::find(collection.begin(), collection.end(), p.getType()) != collection.end();
        };
        const auto nFromParticles = std::count_if(particles.begin(), particles.end(),
                                                  [this, isInCollection](const readdy::model::Particle &p) {
                                                      return isInCollection(p, typeCountFrom);
                                                  });
        {
            const auto &distSquared = kernel->getKernelContext().getDistSquaredFun();
            for (auto &&pFrom : particles) {
                if (isInCollection(pFrom, typeCountFrom)) {
                    for (auto &&pTo : particles) {
                        if (isInCollection(pTo, typeCountTo) && pFrom.getId() != pTo.getId()) {
                            const auto dist = sqrt(distSquared(pFrom.getPos(), pTo.getPos()));
                            auto upperBound = std::upper_bound(binBorders.begin(), binBorders.end(), dist);
                            if (upperBound != binBorders.end()) {
                                const auto binBordersIdx = upperBound - binBorders.begin();
                                counts[binBordersIdx - 1]++;
                            }
                        }
                    }
                }
            }
        }

        auto &radialDistribution = std::get<1>(result);
        {
            const auto &binCenters = std::get<0>(result);
            auto &&it_centers = binCenters.begin();
            auto &&it_distribution = radialDistribution.begin();
            for (auto &&it_counts = counts.begin(); it_counts != counts.end(); ++it_counts) {
                const auto idx = it_centers - binCenters.begin();
                const auto lowerRadius = binBorders[idx];
                const auto upperRadius = binBorders[idx + 1];
                *it_distribution =
                        (*it_counts) /
                        (4 / 3 * util::numeric::pi() * (std::pow(upperRadius, 3) - std::pow(lowerRadius, 3)) *
                         nFromParticles * particleToDensity);
                ++it_distribution;
                ++it_centers;
            }
        }
    }
}

const std::vector<double> &RadialDistribution::getBinBorders() const {
    return binBorders;
}

void RadialDistribution::setBinBorders(const std::vector<double> &binBorders) {
    if (binBorders.size() > 1) {
        RadialDistribution::binBorders = binBorders;
        auto nCenters = binBorders.size() - 1;
        result = std::make_pair(std::vector<double>(nCenters), std::vector<double>(nCenters));
        counts = std::vector<double>(nCenters);
        auto &binCenters = std::get<0>(result);
        auto it_begin = binBorders.begin();
        auto it_begin_next = it_begin + 1;
        size_t idx = 0;
        while (it_begin_next != binBorders.end()) {
            binCenters[idx++] = (*it_begin + *it_begin_next) / 2;
            ++it_begin;
            ++it_begin_next;
        }
    } else {
        log::warn("Argument bin borders' size should be at least two to make sense.");
    }

}

RadialDistribution::RadialDistribution(Kernel *const kernel, unsigned int stride,
                                       std::vector<double> binBorders,
                                       const std::vector<std::string> &typeCountFrom,
                                       const std::vector<std::string> &typeCountTo, double particleToDensity)
        : RadialDistribution(kernel, stride, binBorders,
                             _internal::util::transformTypes2(typeCountFrom, kernel->getKernelContext()),
                             _internal::util::transformTypes2(typeCountTo, kernel->getKernelContext()),
                             particleToDensity
) {}

void RadialDistribution::initializeDataSet(io::File &file, const std::string &dataSetName, unsigned int flushStride) {
    if (!pimpl->writerRadialDistribution) {
        auto &centers = std::get<0>(result);
        std::vector<readdy::io::h5::dims_t> fs = {flushStride, centers.size()};
        std::vector<readdy::io::h5::dims_t> dims = {readdy::io::h5::UNLIMITED_DIMS, centers.size()};
        const auto path = OBSERVABLES_GROUP_PATH + "/" + dataSetName;
        auto group = file.createGroup(path);
        log::debug("created group with path {}", path);
        group.write("bin_centers", centers);
        auto dataSet = std::make_unique<Impl::writer_t>(
                "distribution", group, fs, dims
        );
        pimpl->writerRadialDistribution = std::move(dataSet);
    }
}

void RadialDistribution::append() {
    auto &dist = std::get<1>(result);
    pimpl->writerRadialDistribution->append({1, dist.size()}, dist.data());
}

void RadialDistribution::flush() {
    if (pimpl->writerRadialDistribution) pimpl->writerRadialDistribution->flush();
}

RadialDistribution::~RadialDistribution() = default;

struct CenterOfMass::Impl {
    using writer_t = io::DataSet<Vec3POD, false>;
    std::unique_ptr<writer_t> ds;
};

CenterOfMass::CenterOfMass(readdy::model::Kernel *const kernel, unsigned int stride,
                           unsigned int particleType)
        : Observable(kernel, stride), particleTypes({particleType}), pimpl(std::make_unique<Impl>()) {
}

CenterOfMass::CenterOfMass(Kernel *const kernel, unsigned int stride,
                           const std::string &particleType)
        : CenterOfMass(kernel, stride, kernel->getKernelContext().getParticleTypeID(particleType)) {}

void CenterOfMass::evaluate() {
    readdy::model::Vec3 com{0, 0, 0};
    unsigned long n_particles = 0;
    for (auto &&p : kernel->getKernelStateModel().getParticles()) {
        if (particleTypes.find(p.getType()) != particleTypes.end()) {
            ++n_particles;
            com += p.getPos();
        }
    }
    com /= n_particles;
    result = com;
}

CenterOfMass::CenterOfMass(Kernel *const kernel, unsigned int stride,
                           const std::vector<unsigned int> &particleTypes)
        : Observable(kernel, stride), particleTypes(particleTypes.begin(), particleTypes.end()),
          pimpl(std::make_unique<Impl>()) {
}

CenterOfMass::CenterOfMass(Kernel *const kernel, unsigned int stride, const std::vector<std::string> &particleType)
        : Observable(kernel, stride), pimpl(std::make_unique<Impl>()), particleTypes() {
    for (auto &&pt : particleType) {
        particleTypes.emplace(kernel->getKernelContext().getParticleTypeID(pt));
    }

}

void CenterOfMass::initializeDataSet(io::File &file, const std::string &dataSetName, unsigned int flushStride) {
    if (!pimpl->ds) {
        std::vector<readdy::io::h5::dims_t> fs = {flushStride};
        std::vector<readdy::io::h5::dims_t> dims = {readdy::io::h5::UNLIMITED_DIMS};
        auto group = file.createGroup(OBSERVABLES_GROUP_PATH);
        auto dataSet = std::make_unique<Impl::writer_t>(
                dataSetName, group, fs, dims, Vec3MemoryType(), Vec3FileType()
        );
        log::debug("created data set with path {}", OBSERVABLES_GROUP_PATH + "/" + dataSetName);
        pimpl->ds = std::move(dataSet);
    }
}

void CenterOfMass::append() {
    auto pod = Vec3POD(result);
    pimpl->ds->append({1}, &pod);
}

void CenterOfMass::flush() {
    if (pimpl->ds) pimpl->ds->flush();
}

CenterOfMass::~CenterOfMass() = default;

NParticles::NParticles(Kernel *const kernel, unsigned int stride,
                       std::vector<std::string> typesToCount)
        : NParticles(kernel, stride,
                     _internal::util::transformTypes2(typesToCount, kernel->getKernelContext())) {

}

NParticles::NParticles(Kernel *const kernel, unsigned int stride,
                       std::vector<unsigned int> typesToCount)
        : Observable(kernel, stride), typesToCount(typesToCount), pimpl(std::make_unique<Impl>()) {
}

NParticles::~NParticles() = default;

struct NParticles::Impl {
    using data_set_t = io::DataSet<unsigned long, false>;
    std::unique_ptr<data_set_t> ds;
};

NParticles::NParticles(Kernel *const kernel, unsigned int stride)
        : Observable(kernel, stride), pimpl(std::make_unique<Impl>()) {}

void NParticles::initializeDataSet(io::File &file, const std::string &dataSetName, unsigned int flushStride) {
    if (!pimpl->ds) {
        const auto size = typesToCount.empty() ? 1 : typesToCount.size();
        std::vector<readdy::io::h5::dims_t> fs = {flushStride, size};
        std::vector<readdy::io::h5::dims_t> dims = {readdy::io::h5::UNLIMITED_DIMS, size};
        auto group = file.createGroup(OBSERVABLES_GROUP_PATH);
        auto dataSet = std::make_unique<Impl::data_set_t>(
                dataSetName, group, fs, dims
        );
        log::debug("created data set with path {}", OBSERVABLES_GROUP_PATH + "/" + dataSetName);
        pimpl->ds = std::move(dataSet);
    }
}

void NParticles::append() {
    pimpl->ds->append({1, result.size()}, result.data());
}

void NParticles::flush() {
    if (pimpl->ds) pimpl->ds->flush();
}

struct Forces::Impl {
    using data_set_t = io::DataSet<Vec3POD, true>;
    std::unique_ptr<data_set_t> dataSet;
};

Forces::Forces(Kernel *const kernel, unsigned int stride, std::vector<std::string> typesToCount)
        : Forces(kernel, stride,
                 _internal::util::transformTypes2(typesToCount, kernel->getKernelContext())) {}

Forces::Forces(Kernel *const kernel, unsigned int stride, std::vector<unsigned int> typesToCount)
        : Observable(kernel, stride), typesToCount(typesToCount), pimpl(std::make_unique<Impl>()) {}

Forces::Forces(Kernel *const kernel, unsigned int stride) : Observable(kernel, stride), typesToCount({}),
                                                            pimpl(std::make_unique<Impl>()) {}

void Forces::flush() {
    if (pimpl->dataSet) pimpl->dataSet->flush();
}

void Forces::initializeDataSet(io::File &file, const std::string &dataSetName, unsigned int flushStride) {
    if (!pimpl->dataSet) {
        std::vector<readdy::io::h5::dims_t> fs = {flushStride};
        std::vector<readdy::io::h5::dims_t> dims = {readdy::io::h5::UNLIMITED_DIMS};
        auto dataSet = std::make_unique<io::DataSet<Vec3POD, true>>(
                dataSetName, file.createGroup(OBSERVABLES_GROUP_PATH), fs, dims,
                Vec3MemoryType(), Vec3FileType()
        );
        pimpl->dataSet = std::move(dataSet);
    }
}

void Forces::append() {
    std::vector<Vec3POD> podVec(result.begin(), result.end());
    pimpl->dataSet->append({1}, &podVec);
}

Forces::~Forces() = default;

struct HistogramAlongAxis::Impl {
    using data_set_t = io::DataSet<double, false>;
    std::unique_ptr<data_set_t> dataSet;
};

HistogramAlongAxis::HistogramAlongAxis(readdy::model::Kernel *const kernel, unsigned int stride,
                                       std::vector<double> binBorders, std::set<unsigned int> typesToCount,
                                       unsigned int axis)
        : Observable(kernel, stride), binBorders(binBorders), typesToCount(typesToCount), axis(axis),
          pimpl(std::make_unique<Impl>()) {
    auto nCenters = binBorders.size() - 1;
    result = std::vector<double>(nCenters);
}


HistogramAlongAxis::HistogramAlongAxis(Kernel *const kernel, unsigned int stride,
                                       std::vector<double> binBorders,
                                       std::vector<std::string> typesToCount,
                                       unsigned int axis)
        : HistogramAlongAxis(kernel, stride, binBorders,
                             _internal::util::transformTypes(typesToCount, kernel->getKernelContext()),
                             axis) {

}

void HistogramAlongAxis::initializeDataSet(io::File &file, const std::string &dataSetName, unsigned int flushStride) {
    if (!pimpl->dataSet) {
        const auto size = result.size();
        std::vector<readdy::io::h5::dims_t> fs = {flushStride, size};
        std::vector<readdy::io::h5::dims_t> dims = {readdy::io::h5::UNLIMITED_DIMS, size};
        const auto path = OBSERVABLES_GROUP_PATH + "/" + dataSetName;
        auto group = file.createGroup(OBSERVABLES_GROUP_PATH);
        log::debug("created data set with path {}", path);
        auto dataSet = std::make_unique<Impl::data_set_t>(dataSetName, group, fs, dims);
        pimpl->dataSet = std::move(dataSet);
    }
}

void HistogramAlongAxis::append() {
    pimpl->dataSet->append({1, result.size()}, result.data());
}

void HistogramAlongAxis::flush() {
    if (pimpl->dataSet) pimpl->dataSet->flush();
}

HistogramAlongAxis::~HistogramAlongAxis() = default;

struct Particles::Impl {
    using particle_t = readdy::model::Particle;
    using types_writer_t = io::DataSet<particle_t::type_type, true>;
    using ids_writer_t = io::DataSet<particle_t::id_type, true>;
    using pos_writer_t = io::DataSet<Vec3POD, true>;
    std::unique_ptr<types_writer_t> dataSetTypes;
    std::unique_ptr<ids_writer_t> dataSetIds;
    std::unique_ptr<pos_writer_t> dataSetPositions;
};

Particles::Particles(Kernel *const kernel, unsigned int stride) : Observable(kernel, stride),
                                                                  pimpl(std::make_unique<Impl>()) {}

void Particles::initializeDataSet(io::File &file, const std::string &dataSetName, unsigned int flushStride) {
    if (!pimpl->dataSetTypes) {
        std::vector<readdy::io::h5::dims_t> fs = {flushStride};
        std::vector<readdy::io::h5::dims_t> dims = {readdy::io::h5::UNLIMITED_DIMS};
        auto group = file.createGroup(OBSERVABLES_GROUP_PATH + "/" + dataSetName);
        {
            auto dataSetTypes = std::make_unique<Impl::types_writer_t>("types", group, fs, dims);
            pimpl->dataSetTypes = std::move(dataSetTypes);
        }
        {
            auto dataSetIds = std::make_unique<Impl::ids_writer_t>("ids", group, fs, dims);
            pimpl->dataSetIds = std::move(dataSetIds);
        }
        {
            auto dataSetPositions = std::make_unique<Impl::pos_writer_t>(
                    "positions", group, fs, dims, Vec3MemoryType(), Vec3FileType()
            );
            pimpl->dataSetPositions = std::move(dataSetPositions);
        }
    }
}

void Particles::append() {
    {
        auto &types = std::get<0>(result);
        pimpl->dataSetTypes->append({1}, &types);
    }
    {
        auto &ids = std::get<1>(result);
        pimpl->dataSetIds->append({1}, &ids);
    }
    {
        const auto &positions = std::get<2>(result);
        std::vector<Vec3POD> podVec(positions.begin(), positions.end());
        pimpl->dataSetPositions->append({1}, &podVec);
    }
}

void Particles::flush() {
    if (pimpl->dataSetTypes) pimpl->dataSetTypes->flush();
    if (pimpl->dataSetIds) pimpl->dataSetIds->flush();
    if (pimpl->dataSetPositions) pimpl->dataSetPositions->flush();
}

Particles::~Particles() = default;

struct ReactionRecordPOD {
    int reactionType;
    time_step_type when;
    Particle::id_type educts[2];
    Particle::id_type products[2];
    Vec3::entry_t position[3];

    ReactionRecordPOD() = default;

    ReactionRecordPOD(const readdy::model::reactions::ReactionRecord &record, readdy::time_step_type t)
            : reactionType(0), when(0), educts{0, 0}, products{0, 0}, position{0, 0, 0} {
        reactionType = static_cast<int>(record.type);
        when = t;
        educts[0] = record.educts[0];
        educts[1] = record.educts[1];
        products[0] = record.products[0];
        products[1] = record.products[1];
        const auto &data = record.where.data();
        position[0] = data[0];
        position[1] = data[1];
        position[2] = data[2];
    }
};

class ReactionRecordPODMemoryType : public readdy::io::DataSetType {
public:
    ReactionRecordPODMemoryType() {
        using entry_t = ReactionRecordPOD;
        log::debug("attempting to build memory type");
        tid = []() -> hid_t {
            auto nativeType = io::NativeCompoundTypeBuilder(sizeof(entry_t))
                    .insert<decltype(std::declval<ReactionRecordPOD>().reactionType)>("reaction_type", offsetof(entry_t, reactionType))
                    .insert<decltype(std::declval<ReactionRecordPOD>().when)>("when", offsetof(entry_t, when))
                    .insertArray<decltype(std::declval<ReactionRecordPOD>().educts), 2>("educts", offsetof(entry_t, educts))
                    .insertArray<decltype(std::declval<ReactionRecordPOD>().products), 2>("products", offsetof(entry_t, products))
                    .insertArray<decltype(std::declval<ReactionRecordPOD>().reactionType), 3>("position", offsetof(entry_t, position))
                    .build();
            return nativeType.tid;
        }();
    }
};

class ReactionRecordPODFileType : public readdy::io::DataSetType {
public:
    ReactionRecordPODFileType() {
        using entry_t = ReactionRecordPOD;
        tid = []() -> hid_t {
            ReactionRecordPODMemoryType memType{};
            auto file_type = H5Tcopy(memType.tid);
            H5Tpack(file_type);
            return file_type;
        }();
    }
};


struct Reactions::Impl {
    using reactions_record_t = ReactionRecordPOD;
    using reactions_writer_t = io::DataSet<reactions_record_t, true>;
    std::unique_ptr<reactions_writer_t> writer;
};

Reactions::Reactions(Kernel *const kernel, unsigned int stride, bool withPositions)
        : super(kernel, stride), pimpl(std::make_unique<Impl>()), withPositions_(withPositions) {
}

const bool &Reactions::withPositions() const {
    return withPositions_;
}

void Reactions::flush() {
    if (pimpl->writer) pimpl->writer->flush();
}

void Reactions::initializeDataSet(io::File &file, const std::string &dataSetName, unsigned int flushStride) {
    if (!pimpl->writer) {
        std::vector<readdy::io::h5::dims_t> fs = {flushStride};
        std::vector<readdy::io::h5::dims_t> dims = {readdy::io::h5::UNLIMITED_DIMS};
        auto group = file.createGroup(OBSERVABLES_GROUP_PATH);
        {
            auto dataSet = std::make_unique<Impl::reactions_writer_t>(
                    dataSetName, group, fs, dims, ReactionRecordPODMemoryType(), ReactionRecordPODFileType()
            );
            pimpl->writer = std::move(dataSet);
        }
    }
}

void Reactions::append() {
    std::vector<ReactionRecordPOD> pods{};
    pods.reserve(result.size());
    for (const auto &rr : result) {
        auto pod = ReactionRecordPOD(rr, t_current);
        pods.push_back(std::move(pod));
    }
    pimpl->writer->append({1}, &pods);
}

void Reactions::initialize(Kernel *const kernel) {
    kernel->getKernelContext().recordReactionsWithPositions() = true;
}

Reactions::~Reactions() = default;


}
}
}
