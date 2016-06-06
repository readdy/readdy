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

                    skipping_iterator<boost::uuids::uuid> end_ids();

                    skipping_iterator<readdy::model::Vec3> begin_positions();

                    skipping_iterator<readdy::model::Vec3> end_positions();

                    skipping_iterator<readdy::model::Vec3> begin_forces();

                    skipping_iterator<readdy::model::Vec3> end_forces();

                    skipping_iterator<unsigned int> begin_types();

                    skipping_iterator<unsigned int> end_types();


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
