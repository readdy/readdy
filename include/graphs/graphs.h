//
// Created by mho on 11/4/19.
//

#pragma once

#include "Vertex.h"
#include "Graph.h"

namespace graphs{
    using DefaultVertex = Vertex<std::size_t>;
    using DefaultGraph = graphs::Graph<graphs::IndexPersistentVector, DefaultVertex>;
}
