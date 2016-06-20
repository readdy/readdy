/**
 * << detailed description >>
 *
 * @file ProgramNames.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 20.06.16
 */
#include <readdy/model/Programs.h>
namespace readdy {
    namespace model {
        namespace _internal {
            const std::string ProgramName<TestProgram>::value = "TestProgram";
            const std::string ProgramName<AddParticleProgram>::value = "AddParticle";
            const std::string ProgramName<DiffuseProgram>::value = "Diffuse";
            const std::string ProgramName<UpdateStateModelProgram>::value = "UpdateStateModel";
        }
    }
}
