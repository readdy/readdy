/**
 * This header file contains definitions for various common utils. Currently:
 *   - getOS: returns a string corresponding to the executing operating system
 *   - isWindows: returns true if we are in windows
 *   - collections::hasKey: convenience method to check if a map contains a certain key
 *   - testing::getPluginsDirectory: Method that checks for some environment variables and then returns a potential
 *     directory in which the kernels are likely located.
 *
 * @file ObservableFactory.h
 * @brief Header file containing some common utils.
 * @author clonker
 * @date 08.03.16
 */

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

        namespace testing {
            std::string getPluginsDirectory();
        }

    }
}

#endif //READDY_MAIN_UTILS_H
