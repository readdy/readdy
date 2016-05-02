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
#include <readdy/common/make_unique.h>
#include <readdy/plugin/Observable.h>
#include <readdy/plugin/_internal/ObservableFactory.h>

namespace fs = boost::filesystem;
namespace utl = readdy::utils;
namespace readdy {
    namespace plugin {
        struct Kernel::Impl {
            /**
             * The name of the kernel.
             */
            std::string name;
            /**
             * todo
             */
            std::unique_ptr<signal_t> signal = std::make_unique<signal_t>();
            /**
             * todo
             */
            std::unique_ptr<_internal::ObservableFactory> observableFactory;
            /**
             * todo
             */
            std::unordered_map<Observable*, boost::signals2::shared_connection_block> observableBlocks {};
        };

        KernelProvider &KernelProvider::getInstance() {
            static KernelProvider instance;
            return instance;
        }

        KernelProvider::KernelProvider(){
            fs::path path = fs::current_path();
            BOOST_LOG_TRIVIAL(debug) << "current path is " << path;
        }

        const std::string KernelProvider::getDefaultKernelDirectory() {
            auto&& dir = getenv("READDY_PLUGIN_DIR");
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

        void KernelProvider::loadKernelsFromDirectory(const std::string& directory) {
            const fs::path p(directory);
            if (fs::exists(p) && fs::is_directory(p)) {
                BOOST_LOG_TRIVIAL(debug) << "attempting to load plugins from directory " << p.string();
                BOOST_LOG_TRIVIAL(debug) << "current path: " << fs::current_path().string();
                for (auto&& dirEntry : boost::make_iterator_range(fs::directory_iterator(p), {})) {
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

        const std::string &Kernel::getName() const {
            return pimpl->name;
        }

        Kernel::Kernel(const std::string &name) : pimpl(std::make_unique<Kernel::Impl>()) {
            BOOST_LOG_TRIVIAL(trace) << "creating kernel " << name;
            pimpl->name = name;
            pimpl->observableFactory = std::make_unique<_internal::ObservableFactory>(this);
        }

        Kernel::~Kernel() {
            BOOST_LOG_TRIVIAL(trace) << "destructing kernel \"" << pimpl->name << "\"";
        }

        std::unique_ptr<Program> Kernel::createProgram(const std::string& name) const {
            throw std::runtime_error("todo, treat this properly (or better: make kernel abstract)");
        }

        std::vector<std::string> Kernel::getAvailablePrograms() const {
            return std::vector<std::string>();
        }

        readdy::model::KernelStateModel& Kernel::getKernelStateModel() const {
            // todo
            throw std::runtime_error("todo, treat this properly (or better: make kernel abstract)");
        }

        readdy::model::KernelContext& Kernel::getKernelContext() const {
            // todo
            throw std::runtime_error("todo");
        }

        boost::signals2::connection Kernel::registerObservable(Observable * const observable) {
            boost::signals2::connection connection = pimpl->signal->connect(std::bind(&Observable::evaluate, observable));
            boost::signals2::shared_connection_block block {connection, false};
            pimpl->observableBlocks[observable] = block;
            return connection;
        }

        boost::signals2::connection Kernel::registerObservable(const ObservableType &observable, unsigned int stride) {
            // todo wrap into observable class, set stride
            return pimpl->signal->connect(observable);
        }

        std::vector<std::string> Kernel::getAvailableObservables() {
            return pimpl->observableFactory->getRegisteredObservableNames();
        }

        std::unique_ptr<Observable> Kernel::createObservable(const std::string &name) {
            return pimpl->observableFactory->create(name);
        }

        std::tuple<std::unique_ptr<Observable>, boost::signals2::connection> Kernel::createAndRegisterObservable(const std::string &name, unsigned int stride) {
            // todo
            auto&& obs = createObservable(name);
            auto&& connection = registerObservable(obs.get());
            return std::make_tuple(std::move(obs), connection);
        }


        Kernel &Kernel::operator=(Kernel && rhs) = default;
        Kernel::Kernel(Kernel && rhs) = default;

        bool KernelProvider::isSharedLibrary(const boost::filesystem::path &path) const {
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

        void KernelProvider::add(const boost::filesystem::path &sharedLib) {
            using namespace _internal;
            const auto&& shared = std::make_shared<KernelPluginDecorator>(sharedLib);
            plugins.emplace(std::make_pair(shared.get()->getName(), std::move(shared)));
        }

        void KernelProvider::add(const std::shared_ptr<Kernel> &&kernel) {
            const std::string name = kernel->getName();
            PluginProvider::add(name, std::move(kernel));
        }
    }
}




