//
// Created by clonker on 08.03.16.
//

#ifndef READDY_MAIN_UTILS_H
#define READDY_MAIN_UTILS_H

#include <string>
#include <boost/smart_ptr/shared_ptr.hpp>
#include <ostream>
#include <iostream>
#include <memory>
#include <tuple>
#include <boost/functional/hash.hpp>

namespace readdy {
    namespace utils {
        std::string getOS();

        bool isWindows();

        namespace boost2std {
            template<typename T>
            boost::shared_ptr<T> make_shared_ptr(std::shared_ptr<T> &ptr) {
                return boost::shared_ptr<T>(ptr.get(), [ptr](T *) mutable { ptr.reset(); });
            }

            template<typename T>
            std::shared_ptr<T> make_shared_ptr(boost::shared_ptr<T> &ptr) {
                return std::shared_ptr<T>(ptr.get(), [ptr](T *) mutable { ptr.reset(); });
            }
        }

        namespace collections {
            template<typename MapType, typename KeyType = std::string>
            inline bool hasKey(const MapType &map, const KeyType &key) {
                return map.find(key) != map.end();
            }
        }

    }
}

namespace std {
    namespace {

        // Code from boost
        // Reciprocal of the golden ratio helps spread entropy
        //     and handles duplicates.
        // See Mike Seymour in magic-numbers-in-boosthash-combine:
        //     http://stackoverflow.com/questions/4948780

        template<class T>
        inline void hash_combine(std::size_t &seed, T const &v) {
            seed ^= hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }

        // Recursive template code derived from Matthieu M.
        template<class Tuple, size_t Index = std::tuple_size<Tuple>::value - 1>
        struct HashValueImpl {
            static void apply(size_t &seed, Tuple const &tuple) {
                HashValueImpl<Tuple, Index - 1>::apply(seed, tuple);
                hash_combine(seed, get<Index>(tuple));
            }
        };

        template<class Tuple>
        struct HashValueImpl<Tuple, 0> {
            static void apply(size_t &seed, Tuple const &tuple) {
                hash_combine(seed, get<0>(tuple));
            }
        };
    }

    template<typename ... TT>
    struct hash<std::tuple<TT...>> {
        size_t
        operator()(std::tuple<TT...> const &tt) const {
            size_t seed = 0;
            HashValueImpl<std::tuple<TT...> >::apply(seed, tt);
            return seed;
        }

    };
}

#endif //READDY_MAIN_UTILS_H
