/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
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
 * @copyright BSD-3
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