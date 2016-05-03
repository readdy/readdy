/**
 * << detailed description >>
 *
 * @file Kernel.cpp.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 02.05.16
 */

#include <readdy/common/make_unique.h>
#include <readdy/model/Kernel.h>
#include <readdy/model/_internal/ObservableFactory.h>
#include <atomic>

namespace readdy {
    namespace model {
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
            /**
             * todo
             */
            std::atomic<bool> hasModelTimeStepListener {false};
            std::function<void()> modelTimeStepListener = nullptr;
            boost::signals2::connection modelTimeStepListenerConnection;
        };

        const std::string &Kernel::getName() const {
            return pimpl->name;
        }

        Kernel::Kernel(const std::string &name) : pimpl(std::make_unique<Kernel::Impl>()) {
            BOOST_LOG_TRIVIAL(trace) << "creating kernel " << name;
            pimpl->name = name;
            pimpl->observableFactory = std::make_unique<_internal::ObservableFactory>(this);
            pimpl->modelTimeStepListener = [this] {
                const auto t = getKernelStateModel().getCurrentTimeStep();
                for(auto&& e : pimpl->observableBlocks) {
                    if(t % e.first->getStride() != 0) {
                        e.second.block();
                    } else {
                        e.second.unblock();
                    }
                }
                (*pimpl->signal)(t);
                return;
            };
        }

        Kernel::~Kernel() {
            if(pimpl->hasModelTimeStepListener)  {
                pimpl->modelTimeStepListenerConnection.disconnect();
            }
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
            if(!pimpl->hasModelTimeStepListener.exchange(true)) {

                pimpl->modelTimeStepListenerConnection = getKernelStateModel().addListener(pimpl->modelTimeStepListener);
            }
            boost::signals2::connection connection = pimpl->signal.get()->connect(std::bind(&Observable::evaluate, observable, std::placeholders::_1));
            boost::signals2::shared_connection_block block {connection, false};
            pimpl->observableBlocks[observable] = block;
            return connection;
        }

        void Kernel::evaluateObservablesAutomatically(bool evaluate) {

        }

        boost::signals2::connection Kernel::registerObservable(const ObservableType &observable, unsigned int stride) {
            // todo wrap into observable class, set stride
            return pimpl->signal->connect(observable);
        }

        std::vector<std::string> Kernel::getAvailableObservables() {
            return pimpl->observableFactory->getRegisteredObservableNames();
        }

        // todo: provide a templated version of this?
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
    }
}