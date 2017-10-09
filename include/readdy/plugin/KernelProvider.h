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
    kernel_ptr create(const std::string &name) const;

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
