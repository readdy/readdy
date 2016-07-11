/**
 * << detailed description >>
 *
 * @file ProgramNames.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 20.06.16
 */
#include <readdy/model/programs/Programs.h>
namespace readdy {
    namespace model {
        namespace programs {
            namespace _internal {
                const std::string ProgramName<Test>::value = "Test";
                const std::string ProgramName<AddParticle>::value = "AddParticle";
                const std::string ProgramName<EulerBDIntegrator>::value = "Eulerian Brownian dynamics integrator";
                const std::string ProgramName<DefaultReactionProgram>::value = "PerformReactions";
                const std::string ProgramName<CalculateForces>::value = "Calculate forces";
                const std::string ProgramName<UpdateNeighborList>::value = "Update neighbor list";
            }
        }
    }
}
