/**
 * << detailed description >>
 *
 * @file Types.h
 * @brief << brief description >>
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
