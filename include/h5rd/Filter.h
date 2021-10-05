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
 * << detailed description >>
 *
 * @file Filter.h
 * @brief << brief description >>
 * @author clonker
 * @date 06.09.17
 * @copyright BSD-3
 */

#pragma once

#include "common.h"

namespace h5rd {

class PropertyList;

class Filter {
public:

    /**
     * default destructor for filter
     */
    virtual ~Filter() = default;

    /**
     * checks availability
     * @return true if available
     */
    virtual bool available() const = 0;

    /**
     * Activate filter for a data set
     * @param plist corresponding property list
     */
    virtual void activate(PropertyList &plist) = 0;

    /**
     * register filter, only relevant for external ones
     */
    virtual void registerFilter() = 0;
};

/**
 * SZIP compression cannot be applied to compound datatypes, array datatypes, variable-length datatypes, enumerations,
 * or any other user-defined datatypes. If an SZIP filter is set in a dataset creation property list used to
 * create a dataset containing a non-allowed datatype, the call to H5Dcreate will fail; the conflict can be detected
 * only when the property list is used.
 */
class SZIPFilter : public Filter {
public:

    enum CodingMethod {
        Entropy, NearestNeighbor
    };

    SZIPFilter(CodingMethod method, unsigned int pixelsPerBlock);

    ~SZIPFilter() override = default;

    bool available() const override;

    void activate(PropertyList &plist) override;

    void registerFilter() override;

private:
    CodingMethod codingMethod;
    unsigned int pixelsPerBlock;
};

class NBITFilter : public Filter {
public:

    ~NBITFilter() override = default;

    bool available() const override;

    void activate(PropertyList &plist) override;

    void registerFilter() override;
};

class ScaleOffsetFilter : public Filter {
public:

    enum ScaleType {
        FloatingPointVariableMinBits, FloatingPointFixedMinBits, IntegerType
    };

    ScaleOffsetFilter(ScaleType scaleType, int scaleFactor);

    ~ScaleOffsetFilter() override = default;

    bool available() const override;

    void activate(PropertyList &plist) override;

    void registerFilter() override;

private:
    ScaleType scaleType;
    int scaleFactor;
};

class ShuffleFilter : public Filter {
public:

    ~ShuffleFilter() override = default;

    bool available() const override;

    void activate(PropertyList &plist) override;

    void registerFilter() override;
};

class Fletcher32Filter : public Filter {
public:

    ~Fletcher32Filter() override = default;

    bool available() const override;

    void activate(PropertyList &plist) override;

    void registerFilter() override;
};

}

#include "detail/Filter_detail.h"
