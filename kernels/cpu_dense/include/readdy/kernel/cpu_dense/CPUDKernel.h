/**
 * << detailed description >>
 *
 * @file Kernel.h
 * @brief << brief description >>
 * @author clonker
 * @date 22.11.16
 */

#ifndef READDY_DENSE_KERNEL_H
#define READDY_DENSE_KERNEL_H

#include <readdy/model/Kernel.h>
#include <readdy/common/thread/Config.h>

#include "CPUDStateModel.h"

namespace readdy {
namespace kernel {
namespace cpu_dense {

class CPUDKernel : public readdy::model::Kernel {
public:
    static const std::string name;

    CPUDKernel();

    ~CPUDKernel();

    // factory method
    static readdy::model::Kernel* create();

    virtual readdy::model::programs::ProgramFactory &getProgramFactory() const override;

    virtual CPUDStateModel &getKernelStateModel() const override;

    virtual readdy::model::KernelContext &getKernelContext() const override;

    virtual readdy::model::potentials::PotentialFactory &getPotentialFactory() const override;

    virtual readdy::model::reactions::ReactionFactory &getReactionFactory() const override;

    virtual readdy::model::observables::ObservableFactory &getObservableFactory() const override;

    unsigned long getNThreads() const;

    void setNThreads(readdy::util::thread::Config::n_threads_t n);

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl;
};

}
}
}

extern "C" const char* name();

extern "C" readdy::model::Kernel* createKernel();


#endif //READDY_DENSE_KERNEL_H
