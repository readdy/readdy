//
// Created by mho on 11/11/19.
//

#pragma once

#include <cstddef>

#include <graphs/graphs.h>

#include "../../common/common.h"

namespace readdy::model::top {

struct VertexData {
    using ParticleIndex = std::size_t;

    ParticleIndex particleIndex {};
    std::string data {};
};

using Vertex = graphs::Vertex<VertexData>;
using Graph = graphs::Graph<Vertex>;

}
