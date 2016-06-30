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

            void ensureModelTimeStepListener(KernelStateModel * model) {
                if (!hasModelTimeStepListener.exchange(true) && modelTimeStepListener) {
                    modelTimeStepListenerConnection = model->addListener(modelTimeStepListener);
                }
            }
        };

        const std::string &Kernel::getName() const {
            return pimpl->name;
        }

        Kernel::Kernel(const std::string &name) : pimpl(std::make_unique<Kernel::Impl>()) {
            pimpl->name = name;
            pimpl->observableFactory = std::make_unique<_internal::ObservableFactory>(this);
            pimpl->modelTimeStepListener = [this] {
                if (pimpl->evaluateObservablesAutomatically) {
                    evaluateObservables();
                }
            };
            std::srand((unsigned int) std::time(0));
        }

        Kernel::~Kernel() {
            if (pimpl->hasModelTimeStepListener) {
                pimpl->modelTimeStepListenerConnection.disconnect();
            }
        }

        readdy::model::KernelStateModel &Kernel::getKernelStateModel() const {
            throw std::runtime_error("This method should not be called directly but overridden in a kernel implementation.");
        }

        readdy::model::KernelContext &Kernel::getKernelContext() const {
            throw std::runtime_error("This method should not be called directly but overridden in a kernel implementation.");
        }

        boost::signals2::scoped_connection Kernel::connectObservable(ObservableBase *const observable) {
            pimpl->ensureModelTimeStepListener(&getKernelStateModel());
            boost::signals2::scoped_connection connection (pimpl->signal.get()->connect(std::bind(&ObservableBase::callback, observable, std::placeholders::_1)));
            boost::signals2::shared_connection_block block{connection, false};
            pimpl->observableBlocks[observable] = block;
            return connection;
        }

        void Kernel::evaluateObservablesAutomatically(bool evaluate) {
            pimpl->evaluateObservablesAutomatically = evaluate;
        }

        std::tuple<std::unique_ptr<readdy::model::ObservableWrapper>, boost::signals2::scoped_connection> Kernel::registerObservable(const ObservableType &observable, unsigned int stride) {
            pimpl->ensureModelTimeStepListener(&getKernelStateModel());
            auto&& wrap = std::make_unique<ObservableWrapper>(this, observable, stride);
            auto&& connection = connectObservable(wrap.get());
            return std::make_tuple(std::move(wrap), std::move(connection));
        }

        readdy::model::_internal::ObservableFactory &Kernel::getObservableFactory() const {
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

        void Kernel::deconnectObservable(ObservableBase *const observable) {
            pimpl->observableBlocks.erase(observable);
        }

        std::vector<std::string> Kernel::getAvailablePotentials() const {
            return std::vector<std::string>();
        }

        std::unique_ptr<readdy::model::potentials::Potential> Kernel::createPotential(std::string &name) const {
            throw std::runtime_error("This method should not be called directly but overridden in a kernel implementation.");
        }

        potentials::PotentialFactory &Kernel::getPotentialFactory() const {
            throw std::runtime_error("This method should not be called directly but overridden in a kernel implementation.");
        }

        readdy::model::programs::ProgramFactory &Kernel::getProgramFactory() const {
            throw std::runtime_error("This method should not be called directly but overridden in a kernel implementation.");
        }

        readdy::model::reactions::ReactionFactory &Kernel::getReactionFactory() const {
            throw std::runtime_error("This method should not be called directly but overridden in a kernel implementation.");
        }

        void Kernel::addParticle(const std::string &type, const Vec3 &pos) {
            getKernelStateModel().addParticle({pos[0], pos[1], pos[2], getKernelContext().getParticleTypeID(type)});
        }


        Kernel &Kernel::operator=(Kernel &&rhs) = default;

        Kernel::Kernel(Kernel &&rhs) = default;
    }
}



