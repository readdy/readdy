/**
 * << detailed description >>
 *
 * @file SingleCPUParticleDataImpl.h
 * @brief << brief description >>
 * @author clonker
 * @date 06.06.16
 */

#ifndef READDY_MAIN_SINGLECPUPARTICLEDATAIMPL_H
#define READDY_MAIN_SINGLECPUPARTICLEDATAIMPL_H

#include <algorithm>
#include <iostream>

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            namespace model {
                template<class T, class A>
                SingleCPUParticleData::skipping_iterator<T, A>::skipping_iterator(SingleCPUParticleData const *data, size_t pos, typename std::vector<T>::iterator it)
                        : it(it), data(data), pos(pos) {
                    skip();
                }

                template<class T, class A>
                SingleCPUParticleData::skipping_iterator<T, A>::skipping_iterator(const skipping_iterator <T, A> &rhs) {
                    it = rhs.it;
                    data = rhs.data;
                    pos = rhs.pos;
                }

                template<class T, class A>
                SingleCPUParticleData::skipping_iterator<T, A>::~skipping_iterator() {
                }

                template<class T, class A>
                SingleCPUParticleData::skipping_iterator<T, A> &SingleCPUParticleData::skipping_iterator<T, A>::operator=(const SingleCPUParticleData::skipping_iterator<T, A> &rhs) {
                    it = rhs.it;
                    return *this;
                }

                template<class T, class A>
                bool SingleCPUParticleData::skipping_iterator<T, A>::operator==(const SingleCPUParticleData::skipping_iterator<T, A> &rhs) {
                    return it == rhs.it;
                }

                template<class T, class A>
                bool SingleCPUParticleData::skipping_iterator<T, A>::operator!=(const skipping_iterator <T, A> &rhs) {
                    return it != rhs.it;
                }

                template<class T, class A>
                SingleCPUParticleData::skipping_iterator<T, A> &SingleCPUParticleData::skipping_iterator<T, A>::operator++() {
                    pos++;
                    it++;
                    skip();
                    return *this;
                }

                template<class T, class A>
                typename SingleCPUParticleData::skipping_iterator<T, A>::reference SingleCPUParticleData::skipping_iterator<T, A>::operator*() const {
                    return *it;
                }

                template<class T, class A>
                typename SingleCPUParticleData::skipping_iterator<T, A>::pointer SingleCPUParticleData::skipping_iterator<T, A>::operator->() const {
                    return &(*it);
                }

                template<class T, class A>
                void SingleCPUParticleData::skipping_iterator<T, A>::skip() {
                    while (std::find(data->deactivatedParticles->begin(), data->deactivatedParticles->end(), pos) != data->deactivatedParticles->end()) {
                        pos++;
                        it++;
                    }
                }

                template<class T, class A>
                SingleCPUParticleData::skipping_iterator<T, A> SingleCPUParticleData::skipping_iterator<T, A>::operator+(size_t offset) const {
                    auto copy = *this;
                    for (auto &&i = 0; i < offset; i++) {
                        copy.it++;
                    }
                    return copy;
                }

                template<class T, class A>
                size_t SingleCPUParticleData::skipping_iterator<T, A>::getInternalPosition() const {
                    return pos;
                }

            }
        }
    }
}

#endif //READDY_MAIN_SINGLECPUPARTICLEDATAIMPL_H
