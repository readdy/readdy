/**
 * Header for the container of time dependent data in a kernel.
 *
 * @file KernelStateModel.h
 * @brief Defines KernelStateModel.
 * @author clonker
 * @date 18/04/16
 */

#ifndef READDY2_MAIN_KERNELSTATEMODEL_H
#define READDY2_MAIN_KERNELSTATEMODEL_H

namespace readdy {
    namespace model {
        class KernelStateModel {
        public:
            virtual ~KernelStateModel();

            virtual void updateModel(bool forces, bool distances) = 0;
        };
    }
}
#endif //READDY2_MAIN_KERNELSTATEMODEL_H
