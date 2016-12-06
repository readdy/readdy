/**
 * << detailed description >>
 *
 * @file CPUKernel.h
 * @brief << brief description >>
 * @author clonker
 * @date 23.06.16
 */


#ifndef READDY_CPUKERNEL_CPUKERNEL_H
#define READDY_CPUKERNEL_CPUKERNEL_H

#include <readdy/model/Kernel.h>
#include <readdy/common/dll.h>
#include "CPUStateModel.h"

namespace readdy {
namespace kernel {
namespace cpu {


class CPUKernel : public readdy::model::Kernel {
public:
    static const std::string name;

    CPUKernel();

    ~CPUKernel();

    // factory method
    static readdy::model::Kernel* create();

    virtual readdy::model::programs::ProgramFactory &getProgramFactory() const override;

    virtual CPUStateModel &getKernelStateModel() const override;

    virtual readdy::model::KernelContext &getKernelContext() const override;

    virtual readdy::model::potentials::PotentialFactory &getPotentialFactory() const override;

    virtual readdy::model::reactions::ReactionFactory &getReactionFactory() const override;

    virtual readdy::model::_internal::ObservableFactory &getObservableFactory() const override;

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

#endif //READDY_CPUKERNEL_CPUKERNEL_H
