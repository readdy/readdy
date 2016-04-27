/**
 * << detailed description >>
 *
 * @file Observable.h
 * @brief << brief description >>
 * @author clonker
 * @date 22.04.16
 */

#ifndef READDY2_MAIN_OBSERVABLE_H
#define READDY2_MAIN_OBSERVABLE_H

#include <boost/signals2/signal.hpp>
#include "KernelContext.h"
#include "KernelStateModel.h"

namespace readdy {
    namespace model {

        class Observable {
        public:
            typedef boost::signals2::signal<void(const std::shared_ptr<KernelContext>, const std::shared_ptr<KernelStateModel>)> signal_t;

            Observable(unsigned int stride, signal_t &signal) : stride(stride) {
                connection = signal.connect(std::bind(&readdy::model::Observable::evaluate, this, std::placeholders::_1, std::placeholders::_2));
            };

            virtual ~Observable() {
                connection.disconnect();
            };

            void evaluate(const std::shared_ptr<KernelContext> &context, const std::shared_ptr<KernelStateModel> &model) {
                if(model->getCurrentTimeStep() % stride == 0) {
                    _evaluate(context, model);
                }
            };

        protected:
            unsigned int stride;
            virtual void _evaluate(const std::shared_ptr<KernelContext> &context, const std::shared_ptr<KernelStateModel> &model) = 0;
        private:
            boost::signals2::connection connection;
        };
    }
}
#endif //READDY2_MAIN_OBSERVABLE_H
