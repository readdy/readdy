/**
 * This file contains mainly the simulation class. It can be
 * instantiated by using the default constructor and provides means to start a
 * simulation in a relatively simple manner.
 *
 * @file Kernel.h
 * @brief File containing the ReaDDy top level API.
 * @author clonker
 * @date 08.02.16
 */

#ifndef READDY_SIMULATION_H
#define READDY_SIMULATION_H

#include <readdy/model/Kernel.h>
#include <readdy/SimulationScheme.h>

namespace readdy {
    /**
     * Simulation is the focus of the high-level C++ API of ReaDDy.
     * This is where the system is set up and run for a certain number of
     * timesteps.
     * Things like temperature, boxsize, reactions and potentials belong to the
     * context and are given to the kernel when run() is called.
     */
    class Simulation {
    public:
        /**
         * The default constructor. Currently only instantiates the pimpl.
         */
        Simulation();

        /**
         * The destructor. Default behavior.
         */
        ~Simulation();

        /**
         * The move constructor: A simulation object is supposed to be moveable.
         */
        Simulation(Simulation &&rhs);

        /**
         * The move assign.
         */
        Simulation &operator=(Simulation &&rhs);

        /**
         * Method that returns the temperature the simulation is supposed to run at.
         * @return the temperature.
         */
        double getKBT() const;

        /**
         * Method to set the temperature of the system.
         * @param kBT the temperature.
         */
        void setKBT(double kBT);

        /**
         * Method to get the current box size.
         * @return the box size as Vec3 object.
         */
        readdy::model::Vec3 getBoxSize() const;

        /**
         * Method to set the box size.
         * @param dx length of the x-axis
         * @param dy length of the y-axis
         * @param dz length of the z-axis
         */
        void setBoxSize(double dx, double dy, double dz);
        /**
         * Method to set the box size.
         * @param boxSize a vector object with three components, indicating the three lengths of the axes.
         */
        void setBoxSize(const readdy::model::Vec3& boxSize);

        /**
         * Method that returns an array which indicates, if there are (partially) periodic boundaries.
         * @return an array of length 3, where each entry tells if the corresponding boundary is periodic.
         */
        std::array<bool, 3> getPeriodicBoundary() const;

        /**
         * Method to set which parts of the boundary are to be periodic.
         * @param periodic an array of length three with the corresponding entries.
         */
        void setPeriodicBoundary(std::array<bool, 3> periodic);

        /**
         * Registers a predefined observable with the kernel. A list of available observables can be obtained by
         * getAvailableObservables().
         * @param callbackFun a function which gets called upon evaluation of the observable
         * @param stride the stride argument which decides how often the observable gets called
         * @param args arguments for creation of the observable
         * @return a uuid with which the observable is associated
         */
        template<typename T, typename... Args>
        unsigned long registerObservable(const std::function<void(typename T::result_t)>& callbackFun,
                                              unsigned int stride, Args... args);

        /**
         * Registers an observable that implements the readdy::model::ObservableBase interface.
         * @param observable the observable
         * @return a uuid with which the observable is associated
         */
        unsigned long registerObservable(readdy::model::ObservableBase& observable);

        /**
         * Gives all available predefined observable names.
         * @return a vector containing all predefined observable names
         * @todo implement this (changed the factory type to constructor+dispatcher)
         */
        std::vector<std::string> getAvailableObservables();

        /**
         * Removes an observable by uuid.
         * @param uuid the uuid of the observable to be removed.
         */
        void deregisterObservable(const unsigned long uuid);

        /**
         * A method to access the particle positions of a certain type.
         * @param type the type
         * @return a vector containing the particle positions
         */
        std::vector<readdy::model::Vec3> getParticlePositions(std::string type);

        /**
         * A method to register a particle type
         * @param name the name (must be unique)
         * @param diffusionCoefficient the diffusion coefficient
         * @param radius the particle's radius, important for some potentials (like, e.g., harmonic repulsion)
         */
        void registerParticleType(const std::string &name, const double diffusionCoefficient, const double radius);

        /**
         * A method that allows to remove a certain potential type.
         * @param id the id of this potential
         * @todo test this thoroughly and see if the context can handle it with its internal maps
         */
        void deregisterPotential(const short id);

        //----------------------
        // Order 1 potentials
        //----------------------
        /**
         * Register a box potential.
         * @param the particle type for which the box potential should take effect
         * @param forceConstant the force constant of the box potential
         * @param origin the coordinate of the lower left corner of the box
         * @param extent the extent from the origin
         * @param considerParticleRadius a boolean that indicates if the particle radius should play a role
         *        when calculating force and energy or not
         * @return a uuid with which the potential can be removed again
         * @todo document this more thoroughly
         */
        const short registerBoxPotential(std::string particleType, double forceConstant,
                                                readdy::model::Vec3 origin, readdy::model::Vec3 extent,
                                                bool considerParticleRadius);
        /**
         * Register a potential of order 1, i.e., a potential that only depends on the position of one particle.
         * @param ptr a const pointer to a PotentialOrder1 instance that is to be evaluated
         * @param type the type for which the potential should be registered
         * @todo return uuid (?), this does only work for the single cpu kernel (-> descriptor language?)
         */
        void registerExternalPotentialOrder1(readdy::model::potentials::PotentialOrder1 *, const std::string &);

        //----------------------
        // Order 2 potentials
        //----------------------
        /**
         * Register a harmonic repulsion potential.
         * @param particleTypeA particle type A
         * @param particleTypeB particle type B
         * @param forceConstant the force constant
         * @return a uuid with which the potential can be removed again
         * @todo document this more thoroughly
         */
        const short registerHarmonicRepulsionPotential(std::string particleTypeA, std::string particleTypeB,
                                                              double forceConstant);

        /**
         * Register a weak interaction piecewise harmonic potential.
         * @param particleTypeA particle type A
         * @param particleTypeB particle type B
         * @param forceConstant the force constant
         * @param desiredParticleDistance the distance at which it is most favorable
         *        for the particles to be (w.r.t. this potential)
         * @param depth the depth of the energy well
         * @param noInteractionDistance the distance at which this potential has no effect anymore
         * @return a uuid with which the potential can be removed again
         * @todo document this more thoroughly, maybe make it available as a method of only two-three of the four: forceConstant, desiredDistance, depth, noInteractionDistance?
         */
        const short  registerWeakInteractionPiecewiseHarmonicPotential(
                std::string particleTypeA, std::string particleTypeB, double forceConstant,
                double desiredParticleDistance, double depth, double noInteractionDistance);
        /**
         * Register a potential of order 2, i.e., a potential that depends on the positions of a particle pair.
         * @param ptr a const pointer to a PotentialOrder2 instance that is to be evaluated
         * @param type1 one of the two types for which this potential should be registered
         * @param type2 one of the two types for which this potential should be registered
         * @todo return uuid (?), this does only work for the single cpu kernel (-> descriptor language?)
         */
        void registerPotentialOrder2(readdy::model::potentials::PotentialOrder2* ptr,
                                     const std::string &type1, const std::string &type2) {
            ensureKernelSelected();
            getSelectedKernel()->getKernelContext().registerExternalPotential(ptr, type1, type2);
        };

        //void registerReaction(const Reaction& reaction);
        //void registerReactionByDescriptor(const std::string descriptor);

        /**
         * A method to add a particle of a previously registered type to the system.
         * @param x the x coordinate
         * @param y the y coordinate
         * @param z the z coordinate
         * @param type the type of the particle
         */
        void addParticle(double x, double y, double z, const std::string& type);

        /**
         * Method that gives access to all the positions of all the particles in the system.
         * @return a vector containing all the positions
         */
        const std::vector<readdy::model::Vec3> getAllParticlePositions() const;

        /**
         * Method that allows to set a kernel for this simulation object. Most methods required that this method
         * is invoked beforehand.
         * @param kernel the kernel that is to be selected
         */
        void setKernel(const std::string& kernel);

        /**
         * Method that returns if a kernel is currently selected.
         * @return true if a kernel is selected, otherwise false.
         */
        bool isKernelSelected() const;

        /**
         * Method that returns the type of the currently selected kernel.
         * Raises if there is no kernel selected.
         * @return the kernel id string
         */
        const std::string& getSelectedKernelType() const;

        /**
         * Method to register a conversion reaction "A->B".
         * @param name the name of the reaction
         * @param from the type of A
         * @param to the type of B
         * @param rate the rate at which this reaction is to be performed
         * @return a uuid with which this reaction can be removed again
         * @todo implement removal of reactions
         */
        const short registerConversionReaction(const std::string &name, const std::string &from,
                                                             const std::string &to, const double rate);

        /**
         * Method to register an enzymatic reaction "A+C->B+C".
         * @param name the name of the reaction
         * @param catalyst the type of C
         * @param from the type of A
         * @param to the type of B
         * @param rate the rate at which this reaction is to be performed
         * @param eductDistance the distance at which B should be placed from C
         * @return a uuid with which this reaction can be removed again
         * @todo implement removal of reactions
         */
        const short registerEnzymaticReaction(const std::string &name, const std::string &catalyst,
                                                            const std::string &from, const std::string &to,
                                                            const double rate, const double eductDistance);

        /**
         * Method to register a fission reaction "A->B+C".
         * @param name the name of the reaction
         * @param from the type of A
         * @param to1 the type of B
         * @param to2 the type of C
         * @param rate the rate at which this reaction is to be performed
         * @param productDistance the distance at which the products are placed
         * @param weight1 the weight for particle B with respect to the product distance
         * @param weight2 the weight for particle C with respect to the product distance
         * @return a uuid with which this reaction can be removed again
         * @todo implement removal of reactions, explain the weights better
         */
        const short registerFissionReaction(const std::string &name, const std::string &from,
                                                          const std::string &to1, const std::string &to2,
                                                          const double rate, const double productDistance,
                                                          const double weight1 = 0.5, const double weight2 = 0.5);

        /**
         * Method to register a fusion reaction "A+B->C".
         * @param name the name of the reaction
         * @param from1 the type of A
         * @param from2 the type of B
         * @param to the type of C
         * @param rate the rate at which this reaction is to be performed
         * @param eductDistance the distance at which particles A and B become reactive
         * @param weight1 the weight of A with respect to the placement of C
         * @param weight2 the weight of B with respect to the placement of C
         * @return a uuid with which this reaction can be removed again
         * @todo implement removal of reactions, explain weights better
         */
        const short registerFusionReaction(const std::string &name, const std::string &from1,
                                                         const std::string &from2, const std::string &to,
                                                         const double rate, const double eductDistance,
                                                         const double weight1 = 0.5, const double weight2 = 0.5);

        /**
         * Method to register a decay reaction.
         * @param name the name of the reaction
         * @param particleType the type for which this decay should be performed
         * @param rate the rate
         * @return a uuid with which this reaction can be removed again
         * @todo implement removal of reactions
         */
        const short registerDecayReaction(const std::string &name, const std::string &particleType,
                                                        const double rate);

        void setTimeStep(const double);

        virtual void run(const readdy::model::observables::time_step_type steps, const double timeStep);

        double getRecommendedTimeStep(unsigned int N) const;

        template<typename SchemeType=readdy::api::ReaDDyScheme>
        readdy::api::SchemeConfigurator<SchemeType> runScheme(bool useDefaults = true);

    private:
        struct Impl;
        std::unique_ptr<readdy::Simulation::Impl> pimpl;
        void ensureKernelSelected() const;
        readdy::model::Kernel *const getSelectedKernel() const;
    };

    class NoKernelSelectedException : public std::runtime_error {
    public:
        NoKernelSelectedException(const std::string &__arg);
    };

struct Simulation::Impl {
    std::unordered_map<unsigned long, readdy::signals::scoped_connection> observableConnections {};
    std::unordered_map<unsigned long, std::unique_ptr<readdy::model::ObservableBase>> observables {};
    std::unique_ptr<readdy::model::Kernel> kernel;
    unsigned long counter = 0;
};


template<typename T, typename... Args>
unsigned long Simulation::registerObservable(const std::function<void(typename T::result_t)>& callbackFun, unsigned int stride, Args... args) {
    ensureKernelSelected();
    auto uuid = pimpl->counter++;
    auto obs = pimpl->kernel->createObservable<T>(stride, std::forward<Args>(args)...);
    obs->setCallback(callbackFun);
    auto connection = pimpl->kernel->connectObservable(obs.get());
    pimpl->observables.emplace(uuid, std::move(obs));
    pimpl->observableConnections.emplace(uuid, std::move(connection));
    return uuid;
}

template<typename SchemeType>
readdy::api::SchemeConfigurator<SchemeType> Simulation::runScheme(bool useDefaults) {
    ensureKernelSelected();
    return readdy::api::SchemeConfigurator<SchemeType>(getSelectedKernel(), useDefaults);
}

}

#endif //READDY_SIMULATION_H
