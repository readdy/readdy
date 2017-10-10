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
 * This header file contains the definitions of the hdf5 blosc filter plugin for the h5rd wrapper api.
 * Blosc is a meta compressor that speeds up the application of compression algorithms. More information can be found
 * under http://blosc.org. This filter implementation is based on hdf5-blosc (https://github.com/Blosc/hdf5-blosc).
 *
 * @file BloscFilter.h
 * @brief Definitions of the hdf5 blosc filter plugin for h5rd
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

    /**
     * The available compressors.
     */
    enum Compressor {
        BloscLZ, LZ4, LZ4HC, SNAPPY, ZLIB, ZSTD
    };

    /**
     * Creates a new BloscFilter instance.
     * @param compressor the backing compressor to use, by default the blosc internal LZ4 implementation
     * @param compressionLevel the compression level (0 is the lowest and 9 is the highest compression level)
     * @param shuffle whether to perform bitshuffle
     */
    explicit BloscFilter(Compressor compressor = Compressor::BloscLZ, unsigned int compressionLevel = 9,
                         bool shuffle = true);

    /**
     * default destructor
     */
    ~BloscFilter() override = default;

    /**
     * Checks whether this filter is available as hdf5 filter plugin
     * @return true if it is available
     */
    bool available() const override;

    /**
     * Activates the filter for a data set and its corresponding property list upon creation.
     * @param plist the property list
     */
    void activate(h5rd::PropertyList &plist) override;

    /**
     * Registers the filter plugin with hdf5
     */
    void registerFilter() override;

private:
    /**
     * whether to perform bit shuffle
     */
    bool shuffle;
    /**
     * the compressor to use
     */
    Compressor compressor;
    /**
     * 0 - no compression; 9 - maximal compression
     */
    unsigned int compressionLevel;
};

NAMESPACE_END(io)
NAMESPACE_END(readdy)