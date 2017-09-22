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
 * @file CLLNeighborList.h
 * @brief << brief description >>
 * @author clonker
 * @date 12.09.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include "NeighborList.h"
#include "CellLinkedList.h"
#include "NeighborListContainer.h"

namespace readdy {
namespace kernel {
namespace cpu {
namespace nl {

class ContiguousCLLNeighborList : public NeighborList {
public:
    ContiguousCLLNeighborList(std::uint8_t cll_radius, const readdy::model::Context &context,
                              const readdy::util::thread::Config &config);

    ContiguousCLLNeighborList(data::EntryDataContainer *data, std::uint8_t cll_radius,
                              const readdy::model::Context &context, const readdy::util::thread::Config &config);

    bool is_adaptive() const override;

    void set_up(const util::PerformanceNode &node) override;

    void fill_verlet_list(const util::PerformanceNode &node);

    void update(const util::PerformanceNode &node) override;

    void clear(const util::PerformanceNode &node) override;

    void updateData(DataUpdate &&update) override;

    virtual const_iterator cbegin() const override;

    virtual const_iterator cend() const override;

    const data::EntryDataContainer *data() const override;

    data::EntryDataContainer *data() override;

    std::size_t size() const override;

private:

    std::uint8_t cll_radius;
    data::DefaultDataContainer _data;

    ContiguousCellLinkedList ccll;
    ContiguousCLLNeighborListContainer ccllContainer;
    bool _is_set_up {false};
};

class DynamicCLLNeighborList : public NeighborList {
public:
    DynamicCLLNeighborList(const readdy::model::Context &context,
                           const readdy::util::thread::Config &config);

    DynamicCLLNeighborList(data::EntryDataContainer *data, const readdy::model::Context &context,
                           const readdy::util::thread::Config &config);

    bool is_adaptive() const override;

    void set_up(const util::PerformanceNode &node) override;

    void fill_verlet_list(const util::PerformanceNode &node);

    void update(const util::PerformanceNode &node) override;

    void clear(const util::PerformanceNode &node) override;

    void updateData(DataUpdate &&update) override;

    virtual const_iterator cbegin() const override;

    virtual const_iterator cend() const override;

    const data::EntryDataContainer *data() const override;

    data::EntryDataContainer *data() override;

    std::size_t size() const override;

private:
    data::DefaultDataContainer _data;
    DynamicCellLinkedList dcll;
    DynamicCLLNeighborListContainer dcllContainer;
    bool _is_set_up {false};
};

class CompactCLLNeighborList : public NeighborList {
public:
    CompactCLLNeighborList(std::uint8_t cll_radius, const readdy::model::Context &context, const util::thread::Config &config);

    CompactCLLNeighborList(data::EntryDataContainer *data, std::uint8_t cll_radius,
                           const model::Context &context, const util::thread::Config &config);

    bool is_adaptive() const override;

    void set_up(const util::PerformanceNode &node) override;

    void update(const util::PerformanceNode &node) override;

    void clear(const util::PerformanceNode &node) override;

    void updateData(DataUpdate &&update) override;

    void fill_verlet_list(const util::PerformanceNode &node);

    virtual const_iterator cbegin() const override;

    virtual const_iterator cend() const override;

    const data::EntryDataContainer *data() const override;

    data::EntryDataContainer *data() override;

    const CompactCLLNeighborListContainer &container() const;

    const CompactCellLinkedList &cellLinkedList() const;

    CompactCellLinkedList &cellLinkedList();

    std::size_t size() const override;

private:
    std::uint8_t cll_radius;
    data::DefaultDataContainer _data;
    CompactCellLinkedList ccll;
    CompactCLLNeighborListContainer ccllContainer;
    bool _is_set_up {false};
};

}
}
}
}
