/**
 * Blosc for HDF5 - An HDF5 filter that uses the Blosc compressor.
 *
 * Copyright (C) 2009-2015 Francesc Alted <francesc@blosc.org>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

#ifndef FILTER_BLOSC_H
#define FILTER_BLOSC_H

#ifdef __cplusplus
extern "C" {
#endif

#include "blosc.h"

/* Filter revision number, starting at 1 */
/* #define FILTER_BLOSC_VERSION 1 */
#define FILTER_BLOSC_VERSION 2	/* multiple compressors since Blosc 1.3 */

/* Filter ID registered with the HDF Group */
#define FILTER_BLOSC 32001

/* Registers the filter with the HDF5 library. */
#if defined(_MSC_VER)
__declspec(dllexport)
#endif	/* defined(_MSC_VER) */
int register_blosc(char **version, char **date);

#ifdef __cplusplus
}
#endif

#endif
