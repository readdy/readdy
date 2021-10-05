//
// Created by mho on 10/28/19.
//

#pragma once

#include <list>
#include <ostream>
#include <vector>
#include <sstream>

#include <fmt/format.h>
#include "IndexPersistentVector.h"

namespace graphs {

template<typename... T>
class Vertex {
public:
    using data_type = std::tuple<T...>;
    using NeighborList = std::vector<PersistentIndex>;
    using size_type = typename std::vector<Vertex<T...>>::size_type;

    Vertex(T&&... data);

    Vertex(data_type data);

    Vertex(const Vertex &);

    Vertex &operator=(const Vertex &);

    Vertex(Vertex && other) noexcept;

    Vertex &operator=(Vertex && rhs) noexcept;

    virtual ~Vertex();

    template<typename... T2>
    friend std::ostream &operator<<(std::ostream &os, const Vertex<T2...> &vertex);

    const NeighborList &neighbors() const;

    NeighborList &neighbors();

    void addNeighbor(NeighborList::value_type neighbor);

    void removeNeighbor(NeighborList::value_type neighbor);

    const auto &data() const;

    void setData(data_type data);

    void deactivate();

    bool deactivated() const;

    auto operator->() {
        return &std::get<0>(_data);
    }

    auto operator->() const {
        return &std::get<0>(_data);
    }

private:
    NeighborList _neighbors{};
    data_type _data;
    bool _deactivated {false};
};

}

#include "bits/Vertex_detail.h"
