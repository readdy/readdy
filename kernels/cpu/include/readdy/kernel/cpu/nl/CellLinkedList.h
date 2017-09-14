/********************************************************************
 * Copyright © 2017 Computational Molecular Biology Group,          * 
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
 * << detailed description >>
 *
 * @file CellLinkedList.h
 * @brief << brief description >>
 * @author clonker
 * @date 12.09.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <cstddef>
#include <readdy/common/Index.h>
#include <readdy/common/thread/Config.h>
#include <readdy/model/KernelContext.h>
#include <readdy/kernel/cpu/model/CPUParticleData.h>
#include <readdy/common/Timer.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace nl {

class CellLinkedList {
public:
    CellLinkedList(model::CPUParticleData &data, const readdy::model::KernelContext &context,
                   const readdy::util::thread::Config &config);

    void setUp(scalar skin, std::uint8_t radius, const util::PerformanceNode &node);

    virtual void update(const util::PerformanceNode &node) = 0;
    
    virtual void clear() = 0;
    
    const util::Index3D &cellIndex() const;
    
    const util::Index2D &neighborIndex() const;
    
    const std::vector<std::size_t> &neighbors() const;
    
    virtual std::size_t *particlesBegin(std::size_t cellIndex) = 0;
    virtual const std::size_t *particlesBegin(std::size_t cellIndex) const = 0;
    
    virtual std::size_t *particlesEnd(std::size_t cellIndex) = 0;
    virtual const std::size_t *particlesEnd(std::size_t cellIndex) const = 0;

    virtual std::size_t nParticles(std::size_t cellIndex) const = 0;

    std::size_t *neighborsBegin(std::size_t cellIndex);
    const std::size_t *neighborsBegin(std::size_t cellIndex) const;

    std::size_t *neighborsEnd(std::size_t cellIndex);
    const std::size_t *neighborsEnd(std::size_t cellIndex) const;

    std::size_t nNeighbors(std::size_t cellIndex) const;

protected:
    virtual void setUpBins(const util::PerformanceNode &node) = 0;

    bool _is_set_up {false};

    scalar _skin {0};
    scalar _max_cutoff {0};
    scalar _max_cutoff_skin_squared {0};
    std::uint8_t _radius;

    Vec3 _cellSize {0, 0, 0};

    util::Index3D _cellIndex;
    // index of size (n_cells x (1 + nAdjacentCells)), where the first element tells how many adj cells are stored
    util::Index2D _cellNeighbors;
    // backing vector of _cellNeighbors index of size (n_cells x (1 + nAdjacentCells))
    std::vector<std::size_t> _cellNeighborsContent;

    std::reference_wrapper<model::CPUParticleData> _data;
    std::reference_wrapper<const readdy::model::KernelContext> _context;
    std::reference_wrapper<const readdy::util::thread::Config> _config;
};

class ContiguousCellLinkedList : public CellLinkedList {
public:
    ContiguousCellLinkedList(model::CPUParticleData &data, const readdy::model::KernelContext &context,
                             const util::thread::Config &config);
    void update(const util::PerformanceNode &node) override;

    void clear() override;
    
    const std::vector<std::size_t> &bins() const;
    
    const util::Index2D &binsIndex() const;

    std::size_t *particlesBegin(std::size_t cellIndex) override;

    const std::size_t *particlesBegin(std::size_t cellIndex) const override;

    std::size_t *particlesEnd(std::size_t cellIndex) override;

    const std::size_t *particlesEnd(std::size_t cellIndex) const override;

    size_t nParticles(std::size_t cellIndex) const override;

protected:
    void initializeBinsStructure();

    void fillBins(const util::PerformanceNode &node);

    void setUpBins(const util::PerformanceNode &node) override;

private:
    std::size_t _maxParticlesPerCell {0};
    util::Index2D _binsIndex;
    std::vector<std::size_t> _bins;
};

class CompactCellLinkedList : public CellLinkedList {
public:
    CompactCellLinkedList(model::CPUParticleData &data, const readdy::model::KernelContext &context,
                          const util::thread::Config &config);

    void update(const util::PerformanceNode &node) override;

    void clear() override;

    size_t *particlesBegin(std::size_t cellIndex) override;

    const size_t *particlesBegin(std::size_t cellIndex) const override;

    size_t *particlesEnd(std::size_t cellIndex) override;

    const size_t *particlesEnd(std::size_t cellIndex) const override;

    size_t nParticles(std::size_t cellIndex) const override;

    const std::vector<std::size_t> &head() const;

    const std::vector<std::size_t> &list() const;

protected:
    void setUpBins(const util::PerformanceNode &node) override;

    void fillBins(const util::PerformanceNode &node);

    std::vector<std::size_t> _head;
    // particles, 1-indexed
    std::vector<std::size_t> _list;
};

class DynamicCellLinkedList : public CellLinkedList {
public:
    using Bins = std::vector<std::vector<std::size_t>>;

    DynamicCellLinkedList(model::CPUParticleData &data, const readdy::model::KernelContext &context,
                          const util::thread::Config &config);

    void update(const util::PerformanceNode &node) override;

    void clear() override;

    std::size_t *particlesBegin(std::size_t cellIndex) override;

    const std::size_t *particlesBegin(std::size_t cellIndex) const override;

    std::size_t *particlesEnd(std::size_t cellIndex) override;

    const std::size_t *particlesEnd(std::size_t cellIndex) const override;

    std::size_t nParticles(std::size_t cellIndex) const override;

protected:
    void setUpBins(const util::PerformanceNode &node) override;

    void fillBins(const util::PerformanceNode &node);

private:
    Bins _bins;
};

}
}
}
}
