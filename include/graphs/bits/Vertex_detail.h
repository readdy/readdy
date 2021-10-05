//
// Created by mho on 11/3/19.
//

#pragma once

#include "../Vertex.h"

namespace graphs {

template<typename... T>
inline Vertex<T...>::Vertex(T&&... data) : _data(std::tuple<T...>(std::forward<T>(data)...)) {}

template<typename... T>
inline Vertex<T...>::Vertex(Vertex::data_type data) : _data(std::move(data)) {}

template<typename... T>
inline std::ostream &operator<<(std::ostream &os, const Vertex<T...> &vertex) {
    os << fmt::format("{}", vertex);
    return os;
}

template<typename... T>
inline const typename Vertex<T...>::NeighborList &Vertex<T...>::neighbors() const {
    return _neighbors;
}

template<typename... T>
inline typename Vertex<T...>::NeighborList &Vertex<T...>::neighbors() {
    return _neighbors;
}

template<typename... T>
inline const auto &Vertex<T...>::data() const {
    if constexpr (std::tuple_size_v<data_type> == 1) {
        return std::get<0>(_data);
    } else {
        return _data;
    }
}

template<typename... T>
inline void Vertex<T...>::setData(data_type data) {
    _data = std::move(data);
}

template<typename... T>
inline void Vertex<T...>::addNeighbor(NeighborList::value_type neighbor) {
    auto neighborIt = std::find(_neighbors.begin(), _neighbors.end(), neighbor);
    auto newNeighbor = neighborIt == _neighbors.end();
    if(newNeighbor) {
        _neighbors.push_back(neighbor);
    }
}

template<typename... T>
inline void Vertex<T...>::removeNeighbor(NeighborList::value_type neighbor) {
    auto it = std::find(_neighbors.begin(), _neighbors.end(), neighbor);
    if (it != _neighbors.end()) {
        _neighbors.erase(it);
    }
}

template<typename... T>
inline void Vertex<T...>::deactivate() {
    _deactivated = true;
}

template<typename... T>
inline bool Vertex<T...>::deactivated() const {
    return _deactivated;
}

template<typename... T>
inline Vertex<T...>::~Vertex() = default;

template<typename... T>
inline Vertex< T...> &Vertex<T...>::operator=(Vertex &&rhs) noexcept = default;

template<typename... T>
inline Vertex<T...>::Vertex(Vertex &&other) noexcept = default;

template<typename... T>
inline Vertex<T...>::Vertex(const graphs::Vertex<T...> &) = default;

template<typename... T>
inline Vertex<T...>& Vertex<T...>::operator=(const Vertex&) = default;

}

namespace fmt {
template<typename... T>
struct formatter<graphs::Vertex<T...>> {
    template <typename ParseContext>
    constexpr auto parse(ParseContext &ctx) { return ctx.begin(); }

    template <typename FormatContext>
    auto format(const graphs::Vertex<T...> &v, FormatContext &ctx) {
        std::stringstream ss;
        bool first {true};
        for (const auto neighbor : v.neighbors()) {
            if(!first) {
                ss << ",";
            }
            ss << neighbor.value;
            first = false;
        }
        return format_to(ctx.out(), "Vertex[{}, neighbors=[{}]]", v.data(), ss.str());
    }
};
}
