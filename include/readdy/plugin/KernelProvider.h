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
 * The KernelProvider gives means to access kernels that are in shared libs. The can be indexed by calling
 * loadKernelsFromDirectory(directory).
 *
 * @file KernelProvider.h
 * @brief This file contains the KernelProvider class. It is used to make kernels that reside in shared libs accessible.
 * @author clonker
 * @date 02.05.16
 */

#pragma once

#include <memory>

#include <readdy/model/Kernel.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(util)
NAMESPACE_BEGIN(dll)
class shared_library;
NAMESPACE_END(dll)
NAMESPACE_END(util)

NAMESPACE_BEGIN(plugin)
class KernelDeleter {
    std::shared_ptr<readdy::util::dll::shared_library> ptr;
public:
    /**
     * Creates a KernelDeleter for an internal kernel (which is already included in the core library)
     */
    KernelDeleter();

    /**
     * Creates a KernelDeleter for an external kernel (i.e., from a shared library)
     * @param libPtr pointer to the library
     */
    explicit KernelDeleter(const std::shared_ptr<readdy::util::dll::shared_library> &libPtr);

    /**
     * delete that kernel!
     */
    void operator()(readdy::model::Kernel *);
};

/**
 * The KernelProvider is a singleton which can be accessed by getInstance()
 * and provides Kernels that can be added directly or loaded from directories.
 * If loaded from directories (with loadKernelsFromDirectory(string), the
 * specified directory will be scanned for shared libraries with the required
 * symbols, i.e., with an implementation of the Kernel class.
 */
class KernelProvider {
protected:
    /**
     * The constructor of KernelProvider. As it is a singleton, it is protected.
     */
    KernelProvider();

    /**
     * The destructor of KernelProvider.
     */
    virtual ~KernelProvider();

    /**
     * A protected method that determines if a path points to a shared library.
     *
     * @param path the path
     * @return True if the path points to a shared library, otherwise false.
     */
    bool isSharedLibrary(const std::string &path) const;

public:
    /**
     * pointer to a kernel instance equipped with an appropriate deleter
     */
    using kernel_ptr = std::unique_ptr<model::Kernel, KernelDeleter>;
    /**
     * raw pointer to a kernel instance
     */
    using raw_kernel_ptr = model::Kernel*;

    /**
     * no copying
     */
    KernelProvider(KernelProvider const &) = delete;

    /**
     * no copying
     */
    KernelProvider &operator=(KernelProvider const &) = delete;

    /**
     * no move
     */
    KernelProvider(KernelProvider &&) = delete;

    /**
     * no move
     */
    KernelProvider &operator=(KernelProvider &&) = delete;


    /**
     * Method that returns the singleton KernelProvider.
     *
     * @return The KernelProvider.
     */
    static KernelProvider &getInstance();

    /**
     * Method to load kernels (non-recursively) in shared libraries from a directory.
     *
     * @param directory the directory in which the shared libraries are located.
     */
    void loadKernelsFromDirectory(const std::string &directory);

    /**
     * Method that returns the currently available kernels. Output may change after a call
     * to loadKernelsFromDirectory(directory).
     *
     * @return the currently available kernels
     */
    std::vector<std::string> availableKernels() const;

    /**
     * Adds a kernel to the provider registry by providing a factory method and a name.
     * @param name the name
     * @param creator the factory method
     */
    void add(const std::string &name, const std::function<readdy::model::Kernel *()> &creator);

    /**
     * Method that allows to add a kernel to the KernelProvider by providing a path to a shared lib (containing an implementation of a kernel).
     *
     * @param sharedLib the path to the shared lib
     */
    void add(const std::string &sharedLib);

    /**
     * Method that gives the default kernel directory, i.e., where the kernel implementations are usually to be found.
     * First it is checked, if the environment variable 'READDY_PLUGIN_DIR' is set. In that case, the default kernel directory is the contents
     * of that environment variable.
     * Otherwise, the default kernel directory on unix systems is
     * \code{.unparsed}
     * /usr/local/readdy/lib/readdy_plugins
     * \endcode
     * and the default kernel directory on windows systems is
     * \code{.unparsed}
     * C:\\Program Files\ReaDDy\lib\readdy_plugins
     * \endcode
     *
     * @return the default kernel directory.
     */
    static const std::string getDefaultKernelDirectory();

    /**
     * Create a new kernel instance of specified name.
     * @param name the kernel name
     * @return a unique_ptr to the kernel instance, i.e., the caller has ownership
     */
    kernel_ptr create(std::string_view name) const;

private:
    /**
     * loads the name of a kernel that resides in a shared library
     * @param sharedLib the shared library
     * @return the kernel name
     */
    const std::string loadKernelName(const std::string &sharedLib);

    struct Impl;
    std::unique_ptr<Impl> pimpl;

};

NAMESPACE_END(plugin)
NAMESPACE_END(readdy)
