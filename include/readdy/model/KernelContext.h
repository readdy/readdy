/**
 * todo
 *
 * @file KernelContext.h
 * @brief Container class for time independent information of the KernelContext.
 * @author clonker
 * @date 18.04.16
 * @todo write docs, is kbt really time indep (or should be treated as such)?
 */

#ifndef READDY2_MAIN_KERNELCONTEXT_H
#define READDY2_MAIN_KERNELCONTEXT_H

#include <array>
#include <memory>

namespace readdy {
    namespace model {
        class KernelContext {
        public:
            double getKBT() const;

            void setKBT(double kBT);

            std::array<double, 3> getBoxSize() const;

            void setBoxSize(double dx, double dy, double dz);

            std::array<bool, 3> getPeriodicBoundary() const;

            void setPeriodicBoundary(bool pb_x, bool pb_y, bool pb_z);

            // ctor and dtor
            KernelContext();

            ~KernelContext();

            // move
            KernelContext(KernelContext &&rhs);

            KernelContext &operator=(KernelContext &&rhs);

            // copy
            KernelContext(const KernelContext &rhs);

            KernelContext &operator=(const KernelContext &rhs);

        private:
            struct Impl;
            std::unique_ptr<readdy::model::KernelContext::Impl> pimpl;
        };

    }
}
#endif //READDY2_MAIN_KERNELCONTEXT_H
