/**
 * Contains various typedefs. Currently:
 *   - time_step_type: the type used for time steps
 *   - signal_t: the boost.signals2 type for observables
 *   - ObservableType: the std::function type used for observables
 *
 * @file Types.h
 * @brief Header file containing some basic typedefs.
 * @author clonker
 * @date 27.04.16
 */

#ifndef READDY_MAIN_TYPES_H
#define READDY_MAIN_TYPES_H

#include <boost/signals2/signal.hpp>

namespace readdy {
    namespace model {
        typedef unsigned long int time_step_type;

        typedef boost::signals2::signal<void(readdy::model::time_step_type)> signal_t;
        typedef signal_t::slot_type ObservableType;
    }
}
#endif //READDY_MAIN_TYPES_H
