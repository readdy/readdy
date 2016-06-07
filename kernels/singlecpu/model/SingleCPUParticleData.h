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
                        typedef std::forward_iterator_tag iterator_category;

                        skipping_iterator(SingleCPUParticleData const *data, size_t pos, typename std::vector<T>::iterator it);

                        skipping_iterator(const skipping_iterator &it);

                        ~skipping_iterator();

                        skipping_iterator &operator=(const skipping_iterator &rhs);

                        bool operator==(const skipping_iterator &rhs);

                        bool operator!=(const skipping_iterator &rhs);

                        skipping_iterator &operator++();

                        skipping_iterator operator+(size_t) const;

                        reference operator*() const;

                        pointer operator->() const;

                        size_t getInternalPosition() const;

                    protected:
                        void skip();

                        SingleCPUParticleData const *data;
                        typename std::vector<T>::iterator it;
                        size_t pos;
                    };

                    template<class T, class A = std::allocator<T>>
                    class const_skipping_iterator {
                    public:
                        typedef typename A::difference_type difference_type;
                        typedef typename A::value_type value_type;
                        typedef typename A::reference const_reference;
                        typedef typename A::pointer const_pointer;
                        typedef std::forward_iterator_tag iterator_category;

                        const_skipping_iterator(SingleCPUParticleData const *, size_t, typename std::vector<T>::const_iterator);
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
                        size_t pos;
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

                    void addParticles(const std::vector<readdy::model::Particle> &particles);

                    void removeParticle(const readdy::model::Particle &particle);

                    void removeParticle(const size_t index);

                    readdy::model::Particle operator[](const size_t index);

                    std::vector<size_t> * getDeactivatedParticles() const;

                protected:
                    std::unique_ptr<std::vector<boost::uuids::uuid>> ids;
                    std::unique_ptr<std::vector<readdy::model::Vec3>> positions;
                    std::unique_ptr<std::vector<readdy::model::Vec3>> forces;
                    std::unique_ptr<std::vector<unsigned int>> type;
                    std::unique_ptr<std::vector<size_t>> deactivatedParticles;
                };


            }
        }
    }
}

#include "impl/SingleCPUParticleDataImpl.h"

#endif //READDY_MAIN_SINGLECPUPARTICLEDATA_H
