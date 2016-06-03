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

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            namespace model {

                class SingleCPUParticleData {
                public:


                    template<class T, class A = std::allocator<T>>
                    class particle_iterator {
                    public:

                        typedef A allocator_type;
                        typedef typename A::difference_type difference_type;
                        typedef typename A::value_type value_type;
                        typedef typename A::reference reference;
                        typedef typename A::pointer pointer;
                        typedef std::random_access_iterator_tag input_iterator; //or another tag

                        particle_iterator();
                        particle_iterator(const particle_iterator &rhs);
                        ~particle_iterator();

                        particle_iterator& operator=(const particle_iterator&);
                        bool operator==(const particle_iterator&) const;
                        bool operator!=(const particle_iterator&) const;

                        particle_iterator& operator++();
                        reference operator*() const;
                        pointer operator->() const;
                    };

                private:
                    struct Impl;
                    std::unique_ptr<Impl> pimpl;
                };
            }
        }
    }
}




#endif //READDY_MAIN_SINGLECPUPARTICLEDATA_H
