/**
 * todo
 *
 * @file KernelContext.h
 * @brief Container class for time independent information of the simulation.
 * @author clonker
 * @date 18.04.16
 * @todo write docs
 */

#ifndef READDY2_MAIN_KERNELCONTEXT_H
#define READDY2_MAIN_KERNELCONTEXT_H

namespace readdy {
    namespace model {
        class KernelContext {
        protected:
            double kBT;

        public:
            double getKBT() const;

            void setKBT(double kBT);
        };
    }
}
#endif //READDY2_MAIN_KERNELCONTEXT_H
