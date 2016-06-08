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
                SingleCPUParticleData::const_skipping_iterator<T, A>::const_skipping_iterator(SingleCPUParticleData const *data, typename std::vector<T>::const_iterator begin,
                                                                                              typename std::vector<T>::const_iterator it)
                        : it(it), data(data), begin(begin) {
                    skip();
                }

                template<class T, class A>
                SingleCPUParticleData::const_skipping_iterator<T, A>::const_skipping_iterator(const const_skipping_iterator <T, A> &rhs) {
                    it = rhs.it;
                    data = rhs.data;
                    begin = rhs.begin;
                }

                template<class T, class A>
                SingleCPUParticleData::const_skipping_iterator<T, A>::~const_skipping_iterator() {
                }

                template<class T, class A>
                SingleCPUParticleData::const_skipping_iterator<T, A> &SingleCPUParticleData::const_skipping_iterator<T, A>::operator=(const const_skipping_iterator <T, A> &rhs) {
                    it = rhs.it;
                    return *this;
                }

                template<class T, class A>
                bool SingleCPUParticleData::const_skipping_iterator<T, A>::operator==(const const_skipping_iterator <T, A> &rhs) {
                    return it == rhs.it;
                }

                template<class T, class A>
                SingleCPUParticleData::const_skipping_iterator<T, A> &SingleCPUParticleData::const_skipping_iterator<T, A>::operator++() {
                    it++;
                    skip();
                    return *this;
                }

                template<class T, class A>
                bool SingleCPUParticleData::const_skipping_iterator<T, A>::operator!=(const const_skipping_iterator <T, A> &rhs) {
                    return it != rhs.it;
                }

                template<class T, class A>
                typename SingleCPUParticleData::const_skipping_iterator<T, A>::const_reference SingleCPUParticleData::const_skipping_iterator<T, A>::operator*() const {
                    return *it;
                }

                template<class T, class A>
                typename SingleCPUParticleData::const_skipping_iterator<T, A>::const_pointer SingleCPUParticleData::const_skipping_iterator<T, A>::operator->() const {
                    return &(*it);
                }

                template<class T, class A>
                void SingleCPUParticleData::const_skipping_iterator<T, A>::skip() {
                    auto &&idx = it - begin;
                    while (data->deactivatedParticles->find(idx) != data->deactivatedParticles->end()) {
                        ++idx;
                        ++it;
                    }
                }


                template<class T, class A>
                SingleCPUParticleData::skipping_iterator<T, A>::skipping_iterator(SingleCPUParticleData const *data, typename std::vector<T>::iterator begin, typename std::vector<T>::iterator it)
                        : it(it), data(data), begin(begin) {
                    skip();
                }

                template<class T, class A>
                SingleCPUParticleData::skipping_iterator<T, A>::skipping_iterator(const skipping_iterator <T, A> &rhs) {
                    it = rhs.it;
                    data = rhs.data;
                    begin = rhs.begin;
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
                    auto &&idx = it - begin;
                    while (data->deactivatedParticles->find(idx) != data->deactivatedParticles->end()) {
                        ++idx;
                        ++it;
                    }
                }

                template<class T, class A>
                void SingleCPUParticleData::skipping_iterator<T, A>::skip_reverse() {
                    auto &&idx = it - begin;
                    while (data->deactivatedParticles->find(idx) != data->deactivatedParticles->begin()) {
                        --idx;
                        --it;
                    }
                }

                template<class T, class A>
                SingleCPUParticleData::skipping_iterator<T, A> SingleCPUParticleData::skipping_iterator<T, A>::operator+(size_t offset) const {
                    auto copy = skipping_iterator<T, A>(*this);
                    copy += offset;
                    return copy;
                }

                template<class T, class A>
                size_t SingleCPUParticleData::skipping_iterator<T, A>::getInternalPosition() const {
                    return it - begin;
                }

                template<class T, class A>
                bool SingleCPUParticleData::skipping_iterator<T, A>::operator<(const skipping_iterator <T, A> &rhs) const {
                    return it < rhs.it;
                }

                template<class T, class A>
                bool SingleCPUParticleData::skipping_iterator<T, A>::operator>(const skipping_iterator <T, A> &rhs) const {
                    return it > rhs.it;
                }

                template<class T, class A>
                bool SingleCPUParticleData::skipping_iterator<T, A>::operator<=(const skipping_iterator &rhs) const {
                    return it <= rhs.it;
                }

                template<class T, class A>
                bool SingleCPUParticleData::skipping_iterator<T, A>::operator>=(const skipping_iterator &rhs) const {
                    return it >= rhs.it;
                }

                template<class T, class A>
                SingleCPUParticleData::skipping_iterator<T, A> SingleCPUParticleData::skipping_iterator<T, A>::operator++(int) {
                    auto &&copy = skipping_iterator<T, A>(*this);
                    operator++();
                    return copy;
                }

                template<class T, class A>
                SingleCPUParticleData::skipping_iterator<T, A> &SingleCPUParticleData::skipping_iterator<T, A>::operator--() {
                    --it;
                    skip_reverse();
                    return *this;
                }

                template<class T, class A>
                SingleCPUParticleData::skipping_iterator<T, A> SingleCPUParticleData::skipping_iterator<T, A>::operator--(int) {
                    auto &&copy = skipping_iterator<T, A>(*this);
                    operator--();
                    return copy;
                }

                template<class T, class A>
                SingleCPUParticleData::skipping_iterator<T, A> &SingleCPUParticleData::skipping_iterator<T, A>::operator+=(size_t offset) {
                    auto &&pos = getInternalPosition();
                    /*// how many indices to skip within offset
                    size_t skip = 0;
                    ++pos;
                    // how many indices fit so far
                    size_t fit = 0;
                    // next position to skip
                    auto &&lowerBound = data->deactivatedParticles->lower_bound(pos);
                    fit += *lowerBound - pos;
                    while (fit < offset && lowerBound != data->deactivatedParticles->end()) {
                        ++skip;
                        pos = *lowerBound + 1;
                        lowerBound = data->deactivatedParticles->lower_bound(pos);
                        fit += *lowerBound - pos;
                    }
                    std::advance(it, offset + skip);*/
                    auto&& a = data->deactivatedParticles->upper_bound(pos);
                    auto&& b = data->deactivatedParticles->lower_bound(pos+offset);
                    std::advance(it, offset + (b-a));
                    return *this;
                }

                template<class T, class A>
                typename SingleCPUParticleData::skipping_iterator<T, A>::difference_type SingleCPUParticleData::skipping_iterator<T,A>::operator-(skipping_iterator<T, A> rhs) const {
                    auto&& this_idx = it - begin;
                    auto&& that_idx = rhs.it - rhs.begin;
                    // wlog: this_idx <= that_idx
                    if(that_idx > this_idx) {
                        auto&& tmp = that_idx;
                        that_idx = this_idx;
                        this_idx = tmp;
                    }
                    // iterator pointing to the first ghost particle after this_deactivated_it
                    auto&& this_deactivated_it = data->deactivatedParticles->upper_bound(this_idx);
                    // iterator pointing to the first ghost particle after that_deactivated_it
                    auto&& that_deactivated_it = data->deactivatedParticles->upper_bound(that_idx);
                    return std::abs(that_idx - this_idx) - (that_deactivated_it - this_deactivated_it);
                }


                template<class T, class A>
                SingleCPUParticleData::skipping_iterator<T, A> operator+(size_t offs, const SingleCPUParticleData::skipping_iterator<T, A> &it) {
                    return it + offs;
                }

                template<class T, class A>
                SingleCPUParticleData::skipping_iterator<T, A> &SingleCPUParticleData::skipping_iterator<T, A>::operator-=(size_t offset) {
                    auto &&pos = getInternalPosition();
                    auto&& a = data->deactivatedParticles->lower_bound(pos-offset);
                    auto&& b = data->deactivatedParticles->upper_bound(pos) - 1;
                    std::advance(it, -offset-(b-a));
                    return *this;
                };

                template<class T, class A>
                SingleCPUParticleData::skipping_iterator<T, A> SingleCPUParticleData::skipping_iterator<T, A>::operator-(size_t offset) const {
                    auto copy = skipping_iterator<T, A>(*this);
                    copy -= offset;
                    return copy;
                }


                template<class T, class A>
                SingleCPUParticleData::skipping_iterator<T,A>::reference SingleCPUParticleData::skipping_iterator<T,A>::operator[](size_t idx) const {
                    // todo
                    return *(begin+idx);
                }


            }
        }
    }
}

#endif //READDY_MAIN_SINGLECPUPARTICLEDATAIMPL_H
