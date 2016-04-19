/**
 * << detailed description >>
 *
 * @file KernelContext.cpp
 * @brief Implementation file of the KernelContext.
 * @author clonker
 * @date 18.04.16
 * @todo make proper reference to KernelContext.h, is kBT really indepdendent of t?
 */

#include <readdy/model/KernelContext.h>

using namespace readdy::model;

double KernelContext::getKBT() const  {
    return KernelContext::kBT;
}

void KernelContext::setKBT(double kBT) {
    KernelContext::kBT = kBT;
}



