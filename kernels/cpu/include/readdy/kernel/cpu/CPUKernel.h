/**
 * << detailed description >>
 *
 * @file CPUKernel.h
 * @brief << brief description >>
 * @author clonker
 * @date 23.06.16
 */


#ifndef READDY_MAIN_CPUKERNEL_H
#define READDY_MAIN_CPUKERNEL_H

#define KERNEL_SINGLECPU_NO_EXPORT_ALIAS
#include <readdy/kernel/singlecpu/SingleCPUKernel.h>

// #define BOOST_DLL_FORCE_ALIAS_INSTANTIATION

namespace readdy {
    namespace kernel {
        namespace cpu {
            class CPUKernel : public singlecpu::SingleCPUKernel {
            public:
                static const std::string name;

                CPUKernel();

                // factory method
                static std::unique_ptr<CPUKernel> create();

                virtual readdy::model::programs::ProgramFactory &getProgramFactory() const override;

            private:
                struct Impl;
                std::unique_ptr<readdy::kernel::cpu::CPUKernel::Impl> pimpl;


            };


            BOOST_DLL_ALIAS(
                    readdy::kernel::cpu::CPUKernel::name,
                    name
            );
            // export factory method as "create_kernel"
            BOOST_DLL_ALIAS(
                    readdy::kernel::cpu::CPUKernel::create,
                    createKernel
            );
        }
    }
}

#endif //READDY_MAIN_CPUKERNEL_H
