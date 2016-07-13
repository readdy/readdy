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

#include <boost/dll.hpp>
#include <readdy/model/Kernel.h>
#include "CPUStateModel.h"

// #define BOOST_DLL_FORCE_ALIAS_INSTANTIATION

namespace readdy {
    namespace kernel {
        namespace cpu {
            class CPUKernel : public readdy::model::Kernel {
            public:
                static const std::string name;

                CPUKernel();

                unsigned int getNCores();

                // factory method
                static std::unique_ptr<CPUKernel> create();

                virtual readdy::model::programs::ProgramFactory &getProgramFactory() const override;

                virtual CPUStateModel &getKernelStateModel() const override;

                virtual readdy::model::KernelContext &getKernelContext() const override;

                virtual readdy::model::potentials::PotentialFactory &getPotentialFactory() const override;

                virtual readdy::model::reactions::ReactionFactory &getReactionFactory() const override;


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
