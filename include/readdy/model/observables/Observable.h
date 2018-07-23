/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 ********************************************************************/


/**
 * The struct ObservableName is specialized for every observable type, giving access to a static string that corresponds
 * to an observable's name.
 *
 * ObservableBase is the base class of all observables. They have a stride and a current time step. Also, they provide
 * an abstract evaluate method, in which the subclasses should execute their calculations.
 *
 * Observable derives from ObservableBase and is templateized to a certain result type.
 *
 * Combiner takes two or more Observables and combines their results into a third result to avoid duplication of
 * work.
 *
 * The callback methods are responsible for triggering the evaluation and forwarding the results to the external callback
 * functions, which are optional.
 *
 * @file Observable.h
 * @brief Header file containing the definitions for ObservableName, ObservableBase, Observable and CombinerObservable.
 * @author clonker
 * @date 22.04.16
 */

#pragma once

#include <memory>

#include <h5rd/h5rd.h>

#include <readdy/common/common.h>
#include <readdy/common/signals.h>
#include <readdy/common/logging.h>
#include <readdy/common/tuple_utils.h>
#include <readdy/common/ReaDDyVec3.h>

NAMESPACE_BEGIN(readdy)

NAMESPACE_BEGIN(model)
class Kernel;

NAMESPACE_BEGIN(observables)
/**
 * The signal type for the observable-evaluation callback.
 */
using signal_type = readdy::signals::signal<void(time_step_type)>;
/**
 * Most general type of observable: A function pointer of type signal_type::slot_type.
 */
using observable_type = signal_type::slot_type;

/**
 * Base class for all observables, defining interfaces for serialization of results and striding of the evaluation.
 */
class ObservableBase {
public:
    /**
     * The type of the stride
     */
    using stride_type = readdy::stride_type;

    /**
     * Constructs an object of type ObservableBase. Needs a kernel and a stride, which is defaulted to 1, i.e.,
     * evaluation in every time step.
     * @param kernel the kernel
     * @param stride the stride
     */
    explicit ObservableBase(Kernel* kernel, stride_type stride = 1) : _stride(stride), kernel(kernel) {};

    /**
     * The stride at which the observable gets evaluated. Can be 0, which is equivalent to stride = 1.
     * @return the stride
     */
    stride_type stride() const {
        return _stride;
    }

    /**
     * The observable's current time step, i.e., the time step when it was last evaluated.
     * @return the observable's current time step
     */
    const time_step_type &currentTimeStep() const {
        return t_current;
    }

    ObservableBase(const ObservableBase&) = delete;
    ObservableBase& operator=(const ObservableBase&) = delete;
    ObservableBase(ObservableBase&&) = default;
    ObservableBase& operator=(ObservableBase&&) = delete;

    /**
     * The destructor.
     */
    virtual ~ObservableBase() = default;

    /**
     * Method that will trigger a callback with the current results if shouldExecuteCallback() is true.
     * @param t
     */
    virtual void callback(time_step_type t) {
        if (shouldExecuteCallback(t)) {
            firstCall = false;
            t_current = t;
            evaluate();
            if (writeToFile) append();
        }
    };

    /**
     * Method determining whether the observable should be evaluated at time step t.
     * @param t the time step
     * @return true, if we haven't been evaluated in the current time step yet and stride is 0 or a divisor of t
     */
    virtual bool shouldExecuteCallback(time_step_type t) const {
        return (t_current != t || firstCall) && (_stride == 0 || t % _stride == 0);
    }

    /**
     * Will be called automagically when shouldExecuteCallback() is true. Can also be called manually and will then
     * place the observable's results into the result member (in case of a readdy::model::Observable),
     * based on the current state of the system.
     */
    virtual void evaluate() = 0;

    /**
     * This should be called if the contents of the observable should be written into a file after evaluation.
     * @param file the file to write into
     * @param dataSetName the name of the data set, automatically placed under the group /readdy/observables
     * @param flushStride performance parameter, determining the hdf5-internal chunk size
     */
    void enableWriteToFile(File &file, const std::string &dataSetName, stride_type flushStride) {
        writeToFile = true;
        initializeDataSet(file, dataSetName, flushStride);
    }

    /**
     * Write now! Only has an effect if writing to file was enabled, see enableWriteToFile().
     */
    virtual void flush() = 0;

    virtual std::string type() const = 0;

protected:
    friend class readdy::model::Kernel;

    /**
     * Method that will be called upon registration on a kernel and that allows to
     * modify the simulation setup to the observable's needs.
     * @param kernel the kernel
     */
    virtual void initialize(Kernel * kernel) {};

    /**
     * Method that will be called once, if enableWriteToFile() is called and should create a readdy::io::DataSet that
     * the results are written into.
     * @param dataSetName the name of the data set to be created
     * @param flushStride the flush stride, more specifically the internal hdf5 chunk size
     */
    virtual void initializeDataSet(File &, const std::string &dataSetName, stride_type flushStride) = 0;

    /**
     * Called whenever result should be written into the file
     */
    virtual void append() = 0;

    /**
     * Stride at which the observable gets evaluated
     */
    stride_type _stride;
    /**
     * The kernel which created this observable
     */
    readdy::model::Kernel * kernel;
    /**
     * The current time step of the observable
     */
    time_step_type t_current = 0;
    /**
     * true if we should write to file, otherwise false
     */
    bool writeToFile = false;
    /**
     * this is only initially true and otherwise false
     */
    bool firstCall = true;
};

/**
 * Base class for builtin observables. Contains a result variable in which the current results will be stored, also
 * has a callback function which can be set externally and used to retrieve the data.
 *
 * @tparam Result Parameter describing the type of result in a time step. For instance,
 * positions would be of type std::vector<Vec3>.
 */
template<typename Result>
class Observable : public ObservableBase {
public:
    /**
     * type of the callback function
     */
    using callback_function = std::function<void(const Result &)>;
    /**
     * result type
     */
    typedef Result result_type;

    /**
     * Constructs an observable belonging to a kernel with a certain stride, default initializing the result.
     * @param kernel the kernel
     * @param stride the stride
     */
    Observable(Kernel * kernel, stride_type stride)
            : ObservableBase(kernel, stride), result() {
    }

    /**
     * Will return the observable's current result (i.e., the state of it's last evaluation).
     * @return the result
     */
    const result_type &getResult() {
        return result;
    }

    callback_function &callback() {
        return externalCallback;
    }
    
    const callback_function &callback() const {
        return externalCallback;
    }
    
    /**
     * Function that will evaluate the observable and trigger a callback if ObservableBase#shouldExecuteCallback()
     * is true.
     * @param t the time step
     */
    void callback(time_step_type t) override {
        if (shouldExecuteCallback(t)) {
            ObservableBase::callback(t);
            externalCallback(result);
        }
    }

protected:
    /**
     * the result variable, storing the current state
     */
    Result result;
    /**
     * the callback function
     */
    callback_function externalCallback = [](const Result /*unused*/) {};
};

/**
 * Combiner observable, that unifies a collection of observables and uses their results to compute something.
 * One should take care, that the parent observables strides fit with the stride of this observable or take into
 * account, that some values might be outdated.
 *
 * @tparam RESULT the result type, see readdy::model::observables::Observable.
 * @tparam PARENT_OBS the parent observables which will be evaluated before this observable gets evaluated
 */
template<typename RESULT, typename... PARENT_OBS>
class Combiner : public Observable<RESULT> {
public:

    /**
     * type of the stride
     */
    using stride_type = typename Observable<RESULT>::stride_type;

    /**
     * Constructs a combiner observable.
     * @param kernel the kernel it belongs to
     * @param stride a stride
     * @param parents the parent observables
     */
    Combiner(Kernel* kernel, stride_type stride, PARENT_OBS *... parents)
            : Observable<RESULT>(kernel, stride), parentObservables(std::forward<PARENT_OBS *>(parents)...) {}

    /**
     * If ObservableBase#shouldExecuteCallback() is true, this will call the Observable#callback() function
     * for each of its parents.
     * @param t the current time
     */
    void callback(time_step_type t) override {
        if (ObservableBase::shouldExecuteCallback(t)) {
            readdy::util::for_each_in_tuple(parentObservables, CallbackFunctor(ObservableBase::t_current));
            ObservableBase::callback(t);
        }
    }

    /**
     * Flushing is not supported, as file operations are not supported for combiner observables.
     */
    void flush() override {
        throw std::runtime_error("flush not supported for combiner observables");
    }

protected:
    /**
     * Not supported, see flush().
     * @param file the file
     * @param dataSetName data set name
     * @param flushStride flush stride
     */
    void initializeDataSet(File &file, const std::string &dataSetName, stride_type flushStride) override {
        throw std::runtime_error("not supported for combiner observables");
    }

    /**
     * Not supported, see flush().
     */
    void append() override {
        throw std::runtime_error("not supported for combiner observables");
    }

protected:
    /**
     * The parent observables, stored in a tuple
     */
    std::tuple<PARENT_OBS *...> parentObservables;
private:
    /**
     * Callback functor, which will call a parent observable's callback function.
     */
    struct CallbackFunctor {
        /**
         * its current time step
         */
        time_step_type currentTimeStep;

        /**
         * The current time step
         * @param currentTimeStep current time
         */
        explicit CallbackFunctor(time_step_type currentTimeStep) : currentTimeStep(currentTimeStep) {}

        /**
         * Calling the parent observable's callback function with #currentTimeStep.
         * @tparam T type of the parent observable
         * @param obs the parent observable
         */
        template<typename T>
        void operator()(T *const obs) {
            obs->callback(currentTimeStep);
        }
    };
};


NAMESPACE_END(observables)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
