/**
 * << detailed description >>
 *
 * @file SingleCPUParticleData.h
 * @brief << brief description >>
 * @author clonker
 * @date 03.06.16
 */

#ifndef READDY_MAIN_SINGLECPUPARTICLEDATA_H
#define READDY_MAIN_SINGLECPUPARTICLEDATA_H

#include <memory>
#include <vector>
#include <readdy/model/Particle.h>
#include <set>
#include <boost/container/flat_set.hpp>

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            namespace model {

                class SingleCPUParticleData {
                public:

                    template<class T, class A = std::allocator<T>>
                    class skipping_iterator {
                    public:
                        typedef typename A::difference_type difference_type;
                        typedef typename A::value_type value_type;
                        typedef typename A::reference reference;
                        typedef typename A::pointer pointer;
                        typedef std::random_access_iterator_tag iterator_category;

                        skipping_iterator(SingleCPUParticleData const *, typename std::vector<T>::iterator, typename std::vector<T>::iterator);
                        skipping_iterator(const skipping_iterator &);
                        ~skipping_iterator();

                        skipping_iterator &operator=(const skipping_iterator &);

                        bool operator==(const skipping_iterator &);
                        bool operator!=(const skipping_iterator &);
                        bool operator<(const skipping_iterator&) const;
                        bool operator>(const skipping_iterator&) const;
                        bool operator<=(const skipping_iterator&) const;
                        bool operator>=(const skipping_iterator&) const;

                        skipping_iterator &operator++();
                        skipping_iterator operator+(size_t) const;
                        skipping_iterator operator++(int);
                        skipping_iterator& operator--();
                        skipping_iterator operator--(int);
                        skipping_iterator& operator+=(size_t);
                        skipping_iterator& operator-=(size_t);
                        skipping_iterator operator-(size_t) const;
                        difference_type operator-(skipping_iterator) const;

                        template<class Tf, class Af>
                        friend skipping_iterator<Tf,Af> operator+(size_t, const skipping_iterator<Tf,Af>&);

                        reference operator*() const;
                        pointer operator->() const;
                        reference operator[](size_t) const;

                        size_t getInternalPosition() const;

                    protected:
                        void skip();
                        void skip_reverse();

                        SingleCPUParticleData const *data;
                        typename std::vector<T>::iterator it;
                        typename std::vector<T>::iterator begin;
                    };

                    template<class T, class A = std::allocator<T>>
                    class const_skipping_iterator {
                    public:
                        typedef typename A::difference_type difference_type;
                        typedef typename A::value_type value_type;
                        typedef typename A::reference const_reference;
                        typedef typename A::pointer const_pointer;
                        typedef std::forward_iterator_tag iterator_category;

                        const_skipping_iterator(SingleCPUParticleData const *, typename std::vector<T>::const_iterator, typename std::vector<T>::const_iterator);
                        const_skipping_iterator(const const_skipping_iterator &);
                        ~const_skipping_iterator();

                        const_skipping_iterator &operator=(const const_skipping_iterator &);

                        bool operator==(const const_skipping_iterator &);
                        bool operator!=(const const_skipping_iterator &);

                        const_skipping_iterator& operator++();

                        const_reference operator*() const;

                        const_pointer operator->() const;

                    protected:
                        void skip();

                        SingleCPUParticleData const *data;
                        typename std::vector<T>::const_iterator it;
                        typename std::vector<T>::const_iterator begin;
                    };

                    // ctor / dtor
                    SingleCPUParticleData();

                    SingleCPUParticleData(unsigned int capacity);

                    ~SingleCPUParticleData();

                    // move
                    SingleCPUParticleData(SingleCPUParticleData &&rhs);

                    SingleCPUParticleData &operator=(SingleCPUParticleData &&rhs);

                    // copy
                    SingleCPUParticleData(const SingleCPUParticleData &rhs) = delete;

                    SingleCPUParticleData &operator=(const SingleCPUParticleData &rhs) = delete;

                    skipping_iterator<boost::uuids::uuid> begin_ids();
                    const_skipping_iterator<boost::uuids::uuid> begin_ids() const;
                    const_skipping_iterator<boost::uuids::uuid> cbegin_ids() const;
                    skipping_iterator<boost::uuids::uuid> end_ids();
                    const_skipping_iterator<boost::uuids::uuid> end_ids() const;
                    const_skipping_iterator<boost::uuids::uuid> cend_ids() const;

                    skipping_iterator<readdy::model::Vec3> begin_positions();
                    const_skipping_iterator<readdy::model::Vec3> begin_positions() const;
                    const_skipping_iterator<readdy::model::Vec3> cbegin_positions() const;
                    skipping_iterator<readdy::model::Vec3> end_positions();
                    const_skipping_iterator<readdy::model::Vec3> end_positions() const;
                    const_skipping_iterator<readdy::model::Vec3> cend_positions() const;

                    skipping_iterator<readdy::model::Vec3> begin_forces();
                    const_skipping_iterator<readdy::model::Vec3> begin_forces() const;
                    const_skipping_iterator<readdy::model::Vec3> cbegin_forces() const;
                    skipping_iterator<readdy::model::Vec3> end_forces();
                    const_skipping_iterator<readdy::model::Vec3> end_forces() const;
                    const_skipping_iterator<readdy::model::Vec3> cend_forces() const;

                    skipping_iterator<unsigned int> begin_types();
                    const_skipping_iterator<unsigned int> begin_types() const;
                    const_skipping_iterator<unsigned int> cbegin_types() const;
                    skipping_iterator<unsigned int> end_types();
                    const_skipping_iterator<unsigned int> end_types() const;
                    const_skipping_iterator<unsigned int> cend_types() const;

                    void swap(SingleCPUParticleData &rhs);

                    size_t size();

                    size_t max_size();

                    bool empty();

                    void addParticle(const readdy::model::Particle &particle);
                    void addParticles(const std::vector<readdy::model::Particle> &particles);

                    void removeParticle(const readdy::model::Particle &particle);

                    void removeParticle(const size_t index);

                    readdy::model::Particle operator[](const size_t index);

                    boost::container::flat_set<size_t> * getDeactivatedParticles() const;

                protected:
                    std::unique_ptr<std::vector<boost::uuids::uuid>> ids;
                    std::unique_ptr<std::vector<readdy::model::Vec3>> positions;
                    std::unique_ptr<std::vector<readdy::model::Vec3>> forces;
                    std::unique_ptr<std::vector<unsigned int>> type;
                    std::unique_ptr<boost::container::flat_set<size_t>> deactivatedParticles;
                    //size_t deactivated_index;
                };

            }
        }
    }
}

#include "impl/SingleCPUParticleDataImpl.h"

#endif //READDY_MAIN_SINGLECPUPARTICLEDATA_H
