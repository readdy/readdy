/**
 * This file contains the index_persistent_vector, a vector structure together with a stack of 'blanks'. Removal of
 * elements will result in a push back onto the stack of their respective indices, rendering them 'blank'. This handling
 * potentially increases the memory requirements but avoids the shift of access indices.
 *
 * @file index_persistent_vector.h
 * @brief Definitions for the index_persistent_vector
 * @author clonker
 * @date 09.06.17
 * @copyright BSD-3
 */

#pragma once

#include <functional>

#include "bits/IndexPersistentVector_detail.h"

namespace graphs {

struct PersistentIndex;

template<typename T>
using IndexPersistentVector = detail::IndexPersistentContainer<std::vector, T>;

}
