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
 * << detailed description >>
 *
 * @file AccumulativeWriter_bits.h
 * @brief << brief description >>
 * @author clonker
 * @date 06.01.17
 * @copyright GNU Lesser General Public License v3.0
 */

#ifndef READDY_MAIN_ACCUMULATIVEWRITER_BITS_H
#define READDY_MAIN_ACCUMULATIVEWRITER_BITS_H

#include "../AccumulativeWriter.h"

namespace readdy {
namespace model {
namespace observables {

template<typename DataSetType, typename AppendType>
inline AccumulativeWriter<DataSetType, AppendType>::AccumulativeWriter(unsigned int flushStride,
                                                                       std::unique_ptr<DataSetType> &&dataSet)
        : flushStride(flushStride), count(0), dataSet(std::move(dataSet)), appends{} {
    appends.reserve(flushStride);
}

template<typename DataSetType, typename AppendType>
inline void AccumulativeWriter<DataSetType, AppendType>::append(AppendType &data) {
    if(flushStride == 0 || flushStride == 1) {
        dataSet->append({1}, &data);
    } else {
        appends.push_back(data);
        ++count;
        if (count % flushStride == 0) {
            count = 0;
            dataSet->append({appends.size()}, appends.data());
            appends.clear();
        }
    }
}

}
}
}
#endif //READDY_MAIN_ACCUMULATIVEWRITER_BITS_H
