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
 * @file Exception.h
 * @brief << brief description >>
 * @author clonker
 * @date 05.09.17
 * @copyright BSD-3
 */

#pragma once

#include <utility>
#include <sstream>
#include <exception>
#include <H5Epublic.h>
#include "common.h"

namespace h5rd {


class Exception : public std::runtime_error {
public:
    explicit Exception(const std::string &msg) : std::runtime_error(h5stack(msg)) {}

    static std::string h5stack(const std::string &msg) {
        handle_id stackId = H5Eget_current_stack();
        if (stackId >= 0) {

            H5E_walk2_t walker = [](unsigned int n, const H5E_error2_t *desc, void *clientData) -> herr_t {

                auto *ss = static_cast<std::stringstream *>(clientData);

#if H5_VERSION_GE(2,0,0)
                char major_buf[256], minor_buf[256];
                H5Eget_msg(desc->maj_num, nullptr, major_buf, sizeof(major_buf));
                H5Eget_msg(desc->min_num, nullptr, minor_buf, sizeof(minor_buf));
#else
                char *major_buf = H5Eget_major(desc->maj_num);
                char *minor_buf = H5Eget_minor(desc->min_num);
#endif

                std::string err_string("(");
                err_string += major_buf;
                err_string += ") ";
                err_string += minor_buf;

#if !H5_VERSION_GE(2,0,0)
                free(major_buf);
                free(minor_buf);
#endif

                *ss << err_string << std::endl;
                return 0;
            };

            std::stringstream ss;
            H5Ewalk2(stackId, H5E_WALK_UPWARD, walker, &ss);
            return msg + ": " + ss.str();
        }
        return msg + ": Unknown hdf5 error!";
    }
};
}