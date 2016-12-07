/**
 * << detailed description >>
 *
 * @file File.h
 * @brief << brief description >>
 * @author clonker
 * @date 31.08.16
 */

#ifndef READDY_MAIN_FILE_H
#define READDY_MAIN_FILE_H

#include <string>

namespace readdy {
    namespace io {
        class File {
            File(const std::string& path);
            virtual ~File();
            void flush();
            void close();
        };
    }
}
#endif //READDY_MAIN_FILE_H
