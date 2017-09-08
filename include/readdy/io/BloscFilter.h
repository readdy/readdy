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
 * @file BloscFilter.h
 * @brief << brief description >>
 * @author clonker
 * @date 06.09.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <h5rd/h5rd.h>

#include <readdy/common/common.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(io)

class BloscFilter : public h5rd::Filter {
public:

    enum Compressor {
        BloscLZ, LZ4, LZ4HC, SNAPPY, ZLIB, ZSTD
    };

    explicit BloscFilter(Compressor compressor = Compressor::BloscLZ, unsigned int compressionLevel = 9,
                         bool shuffle = true);

    ~BloscFilter() override = default;

    bool available() const override;

    void activate(h5rd::PropertyList &plist) override;

    void registerFilter() override;

private:
    bool shuffle;
    Compressor compressor;
    // 0 - no compression; 9 - maximal compression
    unsigned int compressionLevel;
};

NAMESPACE_END(io)
NAMESPACE_END(readdy)