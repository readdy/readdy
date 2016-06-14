/**
 * << detailed description >>
 *
 * @file py_observable.h
 * @brief << brief description >>
 * @author clonker
 * @date 14.06.16
 */

#ifndef READDY_MAIN_PY_OBSERVABLE_H
#define READDY_MAIN_PY_OBSERVABLE_H

#include <readdy/model/Kernel.h>
#include <boost/python/object.hpp>


namespace readdy {
    namespace py {
        class PyObservable : public readdy::model::ObservableBase{

        public:
            PyObservable(readdy::model::Kernel *const kernel, unsigned int stride, const boost::python::object &observableFun);
            PyObservable(const PyObservable&);
            PyObservable& operator=(const PyObservable&);
            virtual void evaluate() override;

        private:
            std::shared_ptr<boost::python::object> py_ptr;

        };
    }
}


#endif //READDY_MAIN_PY_OBSERVABLE_H
