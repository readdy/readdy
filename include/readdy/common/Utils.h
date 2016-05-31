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
            inline bool hasKey(const MapType &map, const KeyType& key) {
                return map.find(key) != map.end();
            }
        }
    }
}

#endif //READDY_MAIN_UTILS_H
