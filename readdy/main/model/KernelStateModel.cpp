#include <readdy/model/KernelStateModel.h>
#include <readdy/common/make_unique.h>

/**
 * << detailed description >>
 *
 * @file KernelStateModel.h
 * @brief Implementation of KernelStateModel.
 * @author clonker
 * @date 18/04/16
 * @todo make proper reference to KernelStateModel.h?
 */

namespace readdy {
    namespace model {
        KernelStateModel::~KernelStateModel() = default;

        KernelStateModel::KernelStateModel() {
        }


        KernelStateModel &KernelStateModel::operator=(KernelStateModel &&rhs) = default;

        KernelStateModel::KernelStateModel(KernelStateModel &&rhs) = default;


    }
}

