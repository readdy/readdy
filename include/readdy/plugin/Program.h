//
// Created by clonker on 08.04.16.
//

#ifndef READDY2_MAIN_PROGRAM_H
#define READDY2_MAIN_PROGRAM_H

#include <memory>

namespace readdy {
    namespace plugin {
        class Program {
        public:
            Program(const std::string name);

            const std::string getName() const;

        private:
            struct Impl;
            const std::unique_ptr<Impl> impl_ptr;
        };
    }
}


#endif //READDY2_MAIN_PROGRAM_H
