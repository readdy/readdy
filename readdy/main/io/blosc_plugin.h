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

/*
 * Dynamically loaded filter plugin for HDF5 blosc filter.
 *
 * Author: Kiyoshi Masui <kiyo@physics.ubc.ca>
 * Created: 2014
 *
 *
 * Header file
 * -----------
 *
 * This provides dynamically loaded HDF5 filter functionality (introduced
 * in HDF5-1.8.11, May 2013) to the blosc HDF5 filter.
 *
 * Usage: compile as a shared library and install either to the default
 * search location for HDF5 filter plugins (on Linux
 * /usr/local/hdf5/lib/plugin) or to a location pointed to by the
 * HDF5_PLUGIN_PATH environment variable.
 *
 */


#ifndef PLUGIN_BLOSC_H
#define PLUGIN_BLOSC_H

#include "H5PLextern.h"


H5PL_type_t H5PLget_plugin_type(void);


const void* H5PLget_plugin_info(void);


#endif    // PLUGIN_BLOSC_H

