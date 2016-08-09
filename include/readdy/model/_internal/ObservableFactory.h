/**
 * This header file contains the definition of the ObservableFactory. Its purpose is to create observable of different
 * types in the form of unique_ptrs. The actual implementation of an observable can be changed by specializing the
 * dispatcher for its type and invoking a virtual (and then: overridden) method within the factory.
 *
 * @file ObservableFactory.h
 * @brief This header file contains the definition of the ObservableFactory.
 * @author clonker
 * @date 29.04.16
 */

#ifndef READDY_MAIN_OBSERVABLEFACTORY_H
#define READDY_MAIN_OBSERVABLEFACTORY_H

#include <string>
#include <unordered_map>
#include <readdy/model/Observable.h>
#include <readdy/common/Utils.h>
#include <readdy/model/Observables.h>

namespace readdy {
    namespace model {
        class Kernel;
        namespace _internal {

            class ObservableFactory {
            public:
                ObservableFactory(Kernel *const kernel) : kernel(kernel) {};

                template<typename T, typename Obs1, typename Obs2>
                inline std::unique_ptr<T> create(Obs1 *obs1, Obs2 *obs2, unsigned int stride = 1) const {
                    return std::make_unique<T>(kernel, obs1, obs2, stride);
                };

                template<typename R, typename... Args>
                inline std::unique_ptr<R> create(unsigned int stride, Args... args) const {
                    return std::unique_ptr<R>(ObservableFactory::get_dispatcher<R, Args...>::impl(this, kernel, stride, std::forward<Args>(args)...));
                }

                virtual HistogramAlongAxisObservable* createAxisHistogramObservable(readdy::model::Kernel *const kernel, unsigned int stride,
                                                                                    std::vector<double> binBorders, std::vector<std::string> typesToCount, unsigned int axis) const {
                    // todo: provide default impl
                    throw std::runtime_error("Should be overridden (or todo: provide default impl)");
                }

                virtual NParticlesObservable* createNParticlesObservable(readdy::model::Kernel* const kernel, unsigned int stride, std::vector<std::string> typesToCount = {}) const {
                    throw std::runtime_error("should be overridden (or todo: provide default impl)");
                }

            protected:
                Kernel *const kernel;

                template<typename T, typename... Args> struct get_dispatcher;

                template<typename T, typename... Args> struct get_dispatcher {
                    static T *impl(const ObservableFactory * self, Kernel *const kernel, unsigned int stride, Args... args) {
                        // this only invokes the normal constructor
                        return new T(kernel, stride, std::forward<Args>(args)...);
                    };
                };

                template<typename... Args> struct get_dispatcher<readdy::model::HistogramAlongAxisObservable, Args...> {
                    static HistogramAlongAxisObservable *impl(const ObservableFactory * self, Kernel *const kernel, unsigned int stride, Args... args) {
                        return self->createAxisHistogramObservable(kernel, stride, std::forward<Args>(args)...);
                    }
                };

                template<typename... Args> struct get_dispatcher<readdy::model::NParticlesObservable, Args...> {
                    static NParticlesObservable* impl(const ObservableFactory* self, Kernel*const kernel, unsigned int stride, Args... args) {
                        return self->createNParticlesObservable(kernel, stride, std::forward<Args>(args)...);
                    }
                };
            };
        }
    }
}
#endif //READDY_MAIN_OBSERVABLEFACTORY_H
