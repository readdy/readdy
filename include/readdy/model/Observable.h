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
#include <functional>

namespace readdy {
    namespace model {
        template<typename RetT, typename... ArgsT>
        class Observable<RetT(ArgsT...)> {

        protected:
            boost::signals2::signal<std::function<RetT(ArgsT...)>> signal;
            int stride;
        public:
            Observable(int stride) : stride(stride) {};
            virtual ~Observable() {};

            inline boost::signals2::connection connect(const signal::slot_type &observer) {
                return signal.connect(observer);
            }


        };
    }
}
#endif //READDY2_MAIN_OBSERVABLE_H
