//
// Created by clonker on 07.03.16.
//

#include <iostream>
#include <readdy/plugin/Kernel.h>
#include <boost/range/iterator_range_core.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/dll.hpp>
#include <readdy/common/Utils.h>
#include <readdy/plugin/_internal/KernelPluginDecorator.h>
#include <boost/signals2/shared_connection_block.hpp>


namespace fs = boost::filesystem;
namespace utl = readdy::utils;
namespace plug = readdy::plugin;

plug::KernelProvider &plug::KernelProvider::getInstance() {
    static plug::KernelProvider instance;
    return instance;
}

plug::KernelProvider::KernelProvider(){
    fs::path path = fs::current_path();
    BOOST_LOG_TRIVIAL(debug) << "current path is " << path;
}

const std::string plug::KernelProvider::getDefaultKernelDirectory() {
    auto dir = getenv("READDY_PLUGIN_DIR");
    static std::string defaultDir;
    if (dir == NULL) {
        if (utl::isWindows()) {
            dir = getenv("PROGRAMFILES");
            if (dir == NULL) {
                defaultDir = "C:\\\\Program Files\\ReaDDy2\\lib\\readdy_plugins";
            } else {
                defaultDir = std::string(dir).append("\\ReaDDy2\\lib\\readdy_plugins");
            }
        } else {
            defaultDir = "/usr/local/readdy/lib/readdy_plugins";
        }
    } else {
        defaultDir = std::string(dir);
    }
    return defaultDir;
}

void plug::KernelProvider::loadKernelsFromDirectory(const std::string& directory) {
    const fs::path p(directory);
    if (fs::exists(p) && fs::is_directory(p)) {
        BOOST_LOG_TRIVIAL(debug) << "attempting to load plugins from directory " << p.string();
        //BOOST_LOG_TRIVIAL(debug) << "attempting to load plugins from directory " << p.string();
        BOOST_LOG_TRIVIAL(debug) << "current path: " << fs::current_path().string();
        for (auto dirEntry : boost::make_iterator_range(fs::directory_iterator(p), {})) {
            if (isSharedLibrary(dirEntry.path())) {
                add(dirEntry.path());
            } else {
                BOOST_LOG_TRIVIAL(debug) << "... skipping " << dirEntry.path().string() << " since it was no shared library.";
            }
        }
    } else {
        // TODO raise
        BOOST_LOG_TRIVIAL(debug) << "file [" << p.string() << "] did not exist or was a file.";
    }
    BOOST_LOG_TRIVIAL(debug) << "end of loadKernelsFromDirectory";
}

const std::string &plug::Kernel::getName() const {
    return this->name;
}

plug::Kernel::Kernel(const std::string name) : name(name) {
    BOOST_LOG_TRIVIAL(trace) << "creating kernel " << name;
}

readdy::plugin::Kernel::~Kernel() {
    BOOST_LOG_TRIVIAL(trace) << "destructing kernel \"" << name << "\"";
}

std::shared_ptr<readdy::plugin::Program> readdy::plugin::Kernel::createProgram(std::string name) {
    return nullptr;
}

std::vector<std::string> readdy::plugin::Kernel::getAvailablePrograms() {
    return std::vector<std::string>();
}

std::shared_ptr<readdy::model::KernelStateModel> readdy::plugin::Kernel::getKernelStateModel() {
    // todo
    throw std::runtime_error("todo");
}

std::shared_ptr<readdy::model::KernelContext> readdy::plugin::Kernel::getKernelContext() {
    // todo
    throw std::runtime_error("todo");
}

boost::signals2::connection readdy::plugin::Kernel::registerObservable(const std::shared_ptr<Observable> &observable) {
    // todo
    // todo replace copy of object by ptr
    boost::signals2::connection connection = signal.connect(std::bind(&Observable::evaluate, observable.get(), std::placeholders::_1, std::placeholders::_2));
    boost::signals2::shared_connection_block block {connection, false};
    return connection;
}

boost::signals2::connection readdy::plugin::Kernel::registerObservable(const ObservableType &observable, unsigned int stride) {
    // todo
    return signal.connect(observable);
}

std::vector<std::string> readdy::plugin::Kernel::getAvailableObservables() {
    return std::vector<std::string>();
}

std::shared_ptr<readdy::plugin::Observable> readdy::plugin::Kernel::createObservable(std::string name) {
    // todo
    throw std::runtime_error("todo");
}


bool readdy::plugin::KernelProvider::isSharedLibrary(const boost::filesystem::path &path) {
    const std::string s = path.string();
    return (s.find(".dll") != std::string::npos || s.find(".so") != std::string::npos || s.find(".dylib") != std::string::npos)
           && s.find(".lib") == std::string::npos
           && s.find(".exp") == std::string::npos
           && s.find(".pdb") == std::string::npos
           && s.find(".manifest") == std::string::npos
           && s.find(".rsp") == std::string::npos
           && s.find(".obj") == std::string::npos
           && s.find(".a") == std::string::npos;
}

void readdy::plugin::KernelProvider::add(const boost::filesystem::path &sharedLib) {
    using namespace readdy::plugin::_internal;
    auto shared = std::make_shared<KernelPluginDecorator>(sharedLib);
    plugins.emplace(std::make_pair(shared.get()->getName(), std::move(shared)));
}

void readdy::plugin::KernelProvider::add(const std::shared_ptr<Kernel> &&kernel) {
    const std::string name = kernel->getName();
    PluginProvider::add(name, std::move(kernel));
}


