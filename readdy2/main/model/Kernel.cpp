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
            std::unordered_map<ObservableBase *, boost::signals2::shared_connection_block> observableBlocks{};
            /**
             * todo
             */
            std::atomic<bool> hasModelTimeStepListener{false};
            std::function<void()> modelTimeStepListener = nullptr;
            boost::signals2::connection modelTimeStepListenerConnection;

            bool evaluateObservablesAutomatically = true;
        };

        const std::string &Kernel::getName() const {
            return pimpl->name;
        }

        Kernel::Kernel(const std::string &name) : pimpl(std::make_unique<Kernel::Impl>()) {
            BOOST_LOG_TRIVIAL(trace) << "creating kernel " << name;
            pimpl->name = name;
            pimpl->observableFactory = std::make_unique<_internal::ObservableFactory>(this);
            pimpl->modelTimeStepListener = [this] {
                if (pimpl->evaluateObservablesAutomatically) {
                    evaluateObservables();
                }
            };
        }

        Kernel::~Kernel() {
            if (pimpl->hasModelTimeStepListener) {
                pimpl->modelTimeStepListenerConnection.disconnect();
            }
        }

        std::unique_ptr<Program> Kernel::createProgram(const std::string &name) const {
            throw std::runtime_error("This method should not be called directly but overridden in a kernel implementation.");
        }

        std::vector<std::string> Kernel::getAvailablePrograms() const {
            return std::vector<std::string>();
        }

        readdy::model::KernelStateModel &Kernel::getKernelStateModel() const {
            throw std::runtime_error("This method should not be called directly but overridden in a kernel implementation.");
        }

        readdy::model::KernelContext &Kernel::getKernelContext() const {
            throw std::runtime_error("This method should not be called directly but overridden in a kernel implementation.");
        }

        boost::signals2::connection Kernel::registerObservable(ObservableBase *const observable) {
            if (!pimpl->hasModelTimeStepListener.exchange(true)) {
                pimpl->modelTimeStepListenerConnection = getKernelStateModel().addListener(pimpl->modelTimeStepListener);
            }
            boost::signals2::connection connection = pimpl->signal.get()->connect(std::bind(&ObservableBase::callback, observable, std::placeholders::_1));
            boost::signals2::shared_connection_block block{connection, false};
            pimpl->observableBlocks[observable] = block;
            return connection;
        }

        void Kernel::evaluateObservablesAutomatically(bool evaluate) {
            pimpl->evaluateObservablesAutomatically = evaluate;
        }

        std::tuple<std::unique_ptr<readdy::model::ObservableWrapper>, boost::signals2::connection> Kernel::registerObservable(const ObservableType &observable, unsigned int stride) {
            auto&& wrap = std::make_unique<ObservableWrapper>(this, observable, stride);
            auto&& connection = registerObservable(wrap.get());
            return std::make_tuple(std::move(wrap), connection);
        }

        std::vector<std::string> Kernel::getAvailableObservables() {
            return pimpl->observableFactory->getRegisteredObservableNames();
        }

        std::unique_ptr<ObservableBase> Kernel::createObservable(const std::string &name) {
            return pimpl->observableFactory->create(name);
        }

        std::tuple<std::unique_ptr<ObservableBase>, boost::signals2::connection> Kernel::createAndRegisterObservable(const std::string &name, unsigned int stride) {
            auto &&obs = createObservable(name);
            obs->setStride(stride);
            auto &&connection = registerObservable(obs.get());
            return std::make_tuple(std::move(obs), connection);
        }

        _internal::ObservableFactory &Kernel::getObservableFactory() const {
            return *pimpl->observableFactory;
        }

        void Kernel::evaluateObservables() {
            const auto t = getKernelStateModel().getCurrentTimeStep();
            for (auto &&e : pimpl->observableBlocks) {
                if (e.first->getStride() > 0 && t % e.first->getStride() != 0) {
                    e.second.block();
                } else {
                    e.second.unblock();
                }
            }
            (*pimpl->signal)(t);
        }

        void Kernel::evaluateAllObservables() {
            for (auto &&e : pimpl->observableBlocks) {
                e.second.unblock();
            }
            (*pimpl->signal)(getKernelStateModel().getCurrentTimeStep());
        }


        Kernel &Kernel::operator=(Kernel &&rhs) = default;

        Kernel::Kernel(Kernel &&rhs) = default;
    }
}



