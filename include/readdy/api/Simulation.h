/********************************************************************
 * Copyright © 2016 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * This file is part of ReaDDy.                                     *
 *                                                                  *
 * ReaDDy is free software: you can redistribute it and/or modify   *
 * it under the terms of the GNU Lesser General Public License as   *
 * published by the Free Software Foundation, either version 3 of   *
 * the License, or (at your option) any later version.              *
 *                                                                  *
 * This program is distributed in the hope that it will be useful,  *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of   *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
 * GNU Lesser General Public License for more details.              *
 *                                                                  *
 * You should have received a copy of the GNU Lesser General        *
 * Public License along with this program. If not, see              *
 * <http://www.gnu.org/licenses/>.                                  *
 ********************************************************************/


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

#pragma once

#include <readdy/plugin/KernelProvider.h>
#include <readdy/model/Kernel.h>
#include <readdy/api/SimulationScheme.h>
#include <readdy/model/topologies/reactions/StructuralTopologyReaction.h>
#include "ObservableHandle.h"

NAMESPACE_BEGIN(readdy)
/**
 * Simulation is the focus of the high-level C++ API of ReaDDy.
 * This is where the system is set up and run for a certain number of
 * time steps.
 * Things like temperature, boxsize, reactions and potentials belong to the
 * context and are given to the kernel when run() is called.
 */
class Simulation {
    using particle = readdy::model::Particle;
public:
    using topology_reaction_mode = model::top::reactions::STRMode;

    /**
     * The default constructor. Instantiates the pimpl and decides if performance should be profiled.
     */
    explicit Simulation(bool profile = true);

    /**
     * The destructor. Default behavior.
     */
    ~Simulation();

    /**
     * The move constructor: A simulation object is supposed to be moveable.
     */
    Simulation(Simulation &&rhs) noexcept;

    /**
     * The move assign.
     */
    Simulation &operator=(Simulation &&rhs) noexcept;

    /**
     * Copy constructor, deleted.
     */
    Simulation(const Simulation &) = delete;

    /**
     * Copy assign, deleted
     */
    Simulation &operator=(const Simulation) = delete;

    /**
     * Method that returns the temperature the simulation is supposed to run at.
     * @return the temperature.
     */
    scalar getKBT() const;

    /**
     * Method to set the temperature of the system.
     * @param kBT the temperature.
     */
    void setKBT(scalar kBT);

    /**
     * Method to get the current box size.
     * @return the box size as Vec3 object.
     */
    readdy::model::Vec3 getBoxSize() const;

    topology_type_type registerTopologyType(const std::string& name,
                                            const std::vector<model::top::reactions::StructuralTopologyReaction> &reactions = {});

    readdy::model::TopologyParticle
    createTopologyParticle(const std::string &type, const readdy::model::Vec3 &pos) const;

    bool kernelSupportsTopologies() const;

    readdy::model::top::GraphTopology *addTopology(const std::string& type,
                                                   const std::vector<readdy::model::TopologyParticle> &particles);

    std::vector<readdy::model::top::GraphTopology *> currentTopologies();

    std::vector<model::Particle> getParticlesForTopology(const model::top::GraphTopology &topology) const;

    /**
     * Method to set the box size.
     * @param dx length of the x-axis
     * @param dy length of the y-axis
     * @param dz length of the z-axis
     */
    void setBoxSize(scalar dx, scalar dy, scalar dz);

    /**
     * Method to set the box size.
     * @param boxSize a vector object with three components, indicating the three lengths of the axes.
     */
    void setBoxSize(const readdy::model::Vec3 &boxSize);

    /**
     * Method that returns an array which indicates, if there are (partially) periodic boundaries.
     * @return an array of length 3, where each entry tells if the corresponding boundary is periodic.
     */
    const std::array<bool, 3>& getPeriodicBoundary() const;

    /**
     * Method to set which parts of the boundary are to be periodic.
     * @param periodic an array of length three with the corresponding entries.
     */
    void setPeriodicBoundary(const std::array<bool, 3> &periodic);

    /**
     * Allows to set an expected maximal number of particles in order to avoid reallocations.
     * @param n expected number of particles
     */
    void setExpectedMaxNParticles(std::size_t n);

    /**
     * Registers a predefined observable with the kernel. A list of available observables can be obtained by
     * getAvailableObservables().
     * @param stride the stride argument which decides how often the observable gets called
     * @param args arguments for creation of the observable
     * @return a uuid with which the observable is associated
     */
    template<typename T, typename... Args>
    ObservableHandle registerObservable(unsigned int stride, Args... args);

    /**
     * Registers a predefined observable with the kernel. A list of available observables can be obtained by
     * getAvailableObservables().
     * @param callbackFun a function which gets called upon evaluation of the observable
     * @param stride the stride argument which decides how often the observable gets called
     * @param args arguments for creation of the observable
     * @return a uuid with which the observable is associated
     */
    template<typename T, typename... Args>
    ObservableHandle registerObservable(const std::function<void(typename T::result_type)> &callbackFun,
                                        unsigned int stride, Args... args);

    /**
     * Registers an observable that implements the readdy::model::ObservableBase interface.
     * @param observable the observable
     * @return a uuid with which the observable is associated
     */
    ObservableHandle registerObservable(readdy::model::observables::ObservableBase &observable);

    /**
     * Gives all available predefined observable names.
     * @return a vector containing all predefined observable names
     * @todo implement this (changed the factory type to constructor+dispatcher)
     * @fixme
     */
    std::vector<std::string> getAvailableObservables();

    /**
     * Removes an observable by uuid.
     * @param uuid the uuid of the observable to be removed.
     */
    void deregisterObservable(ObservableHandle::observable_id uuid);

    void deregisterObservable(const ObservableHandle &uuid);

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
    particle::type_type
    registerParticleType(const std::string &name, scalar diffusionCoefficient, scalar radius,
                         readdy::model::particle_flavor flavor = readdy::model::particleflavor::NORMAL);

    /**
     * A method that allows to remove a certain potential type.
     * @param id the id of this potential
     * @todo test this thoroughly and see if the context can handle it with its internal maps
     */
    void deregisterPotential(short id);

    //----------------------
    // Order 1 potentials
    //----------------------
    /**
     * Register a box potential, which is used to confine particles to a cuboid volume. The energy function
     * increases quadratically with respect to the distance from the cuboid edges, resulting in a
     * harmonic repulsion.
     *
     * @param particleType the particle type for which the box potential should take effect
     * @param forceConstant the force constant determines the strength of repulsion
     * @param origin the coordinate of the lower left corner of the box
     * @param extent the extent from the origin
     * @param considerParticleRadius a boolean that indicates if the particle radius should play a role
     *        when calculating force and energy or not
     * @return a uuid with which the potential can be removed
     */
    const short registerBoxPotential(const std::string &particleType, scalar forceConstant,
                                     const readdy::model::Vec3 &origin, const readdy::model::Vec3 &extent,
                                     bool considerParticleRadius);

    /**
     * Register a sphere potential, which is used to confine particles inside a spherical volume. The energy function
     * increases quadratically with respect to the distance from the sphere edge, resulting in a harmonic repulsion.
     *
     * @param particleType the particle type for which the sphere potential should take effect
     * @param forceConstant the force constant determines the strength of repulsion
     * @param origin the center of the sphere
     * @param radius the extent of the sphere
     * @return a uuid with which the potential can be removed
     * @todo add a considerParticleRadius parameter
     */
    const short
    registerSphereInPotential(const std::string &particleType, scalar forceConstant, const readdy::model::Vec3 &origin,
                              scalar radius);

    /**
     * Register a sphere potential, which is used to confine particles outside a spherical volume. The energy function
     * increases quadratically with respect to the distance from the sphere edge, resulting in a harmonic repulsion.
     * 
     * @param particleType the particle type for which the potential should take effect
     * @param forceConstant the force constant determines the strength of interaction, like a spring constant
     * @param origin the center of the sphere
     * @param radius the extent of the sphere
     * @return a uuid with which the potential can be removed
     */
    const short
    registerSphereOutPotential(const std::string &particleType, scalar forceConstant, const readdy::model::Vec3 &origin,
                               scalar radius);

    /**
     * Register a spherical barrier potential. For positive height it represents a concentric barrier around the point origin
     * with a certain radius. The potential consists of multiple harmonic snippets.
     *
     * @param particleType the particle type for which the potential should take effect
     * @param origin the center of the sphere
     * @param radius the radius of the sphere
     * @param height the energetic height of the barrier, can be negative
     * @param width width of the barrier, behaves like full-width-at-half-maximum (FWHM)
     * @return a uuid with which the potential can be removed
     */
    const short
    registerSphericalBarrier(const std::string &particleType, const readdy::model::Vec3 &origin, scalar radius,
                             scalar height, scalar width);

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
    const short registerHarmonicRepulsionPotential(const std::string &particleTypeA, const std::string &particleTypeB,
                                                   scalar forceConstant);

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
    const short registerWeakInteractionPiecewiseHarmonicPotential(
            const std::string &particleTypeA, const std::string &particleTypeB, scalar forceConstant,
            scalar desiredParticleDistance, scalar depth, scalar noInteractionDistance);

    /**
    * Constructs a Lennard-Jones-type potential between two particle types A and B (where possibly A = B) of the form
    *
    * \f[ V_{\mbox{LJ}}(r) = k(\epsilon , n, m) \left[ \left(\frac{\sigma}{r}\right)^m - \left(\frac{\sigma}{r}\right)^n \right], \f]
    *
    * where n,m are exponent 1 and 2, respectively, with m > n.
    * If shift == true, it will be defined as
    *
    * \f[ V_{\mbox{LJ, shifted}}(r) = V_{\mbox{LJ}}(r) - V_{\mbox{LJ}}(r_{\mbox{cutoff}}) \f]
    *
    * for r <= cutoffDistance, which makes a difference in energy, but not in force.
    *
    * @param particleType1 particle type A
    * @param particleType2 particle type B
    * @param m first exponent
    * @param n second exponent
    * @param cutoffDistance the cutoff distance
    * @param shift if it should be shifted or not
    * @param epsilon the well depth
    * @param sigma the distance at which the inter-particle potential is zero
    * @return a uuid with which the potential can be removed again
    */
    const short
    registerLennardJonesPotential(const std::string &type1, const std::string &type2, unsigned int m, unsigned int n,
                                  scalar cutoff, bool shift, scalar epsilon, scalar sigma);

    /**
     * Constructs a potential that describes screened electrostatics with a hard-core repulsion between two
     * particle types A and B (where possibly A = B) of the form
     *
     * \f[ V(r) = C \frac{\exp(-\kappa r)}{r} + D\left(\frac{\sigma}{r}\right)^n, \f]
     *
     * where the first term is the electrostatic interaction, the constant C has the dimensions of an energy times distance. Its value
     * can be positive or negative and depends on the valencies of the particles (see Debye-Hueckel theory). \f$\kappa\f$ is
     * the inverse screening depth. The second term is a hard-core repulsion, that ensures
     * that the potential does not diverge to negative infinity.
     *
     * @param particleType1 particle type A
     * @param particleType2 particle type B
     * @param electrostaticStrength C
     * @param inverseScreeningDepth \f$\kappa\f$
     * @param repulsionStrength D
     * @param repulsionDistance \f$\sigma\f$
     * @param exponent n
     * @param cutoff the distance from which no energies and forces are calculated further
     * @return a uuid with which the potential can be removed again
     */
    const short
    registerScreenedElectrostaticsPotential(const std::string &particleType1, const std::string &particleType2,
                                            scalar electrostaticStrength,
                                            scalar inverseScreeningDepth, scalar repulsionStrength,
                                            scalar repulsionDistance, unsigned int exponent,
                                            scalar cutoff);

    void registerPotentialOrder1(readdy::model::potentials::PotentialOrder1 *ptr);

    /**
     * Register a potential of order 2, i.e., a potential that depends on the positions of a particle pair.
     * @param ptr a const pointer to a PotentialOrder2 instance that is to be evaluated
     * @param type1 one of the two types for which this potential should be registered
     * @param type2 one of the two types for which this potential should be registered
     * @todo return uuid (?), this does only work for the single cpu kernel (-> descriptor language?)
     */
    void registerPotentialOrder2(readdy::model::potentials::PotentialOrder2 *ptr);


    //void registerReactionByDescriptor(const std::string descriptor);

    /**
     * A method to add a particle of a previously registered type to the system.
     * @param x the x coordinate
     * @param y the y coordinate
     * @param z the z coordinate
     * @param type the type of the particle
     */
    void addParticle(const std::string &type, scalar x, scalar y, scalar z);

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
    void setKernel(const std::string &kernel);

    /**
     * Method that allows to set an already existing instance of a kernel for this simulation object.
     * @param kernel  the kernel
     */
    plugin::KernelProvider::raw_kernel_ptr setKernel(plugin::KernelProvider::kernel_ptr &&kernel);

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
    const std::string &getSelectedKernelType() const;

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
                                           const std::string &to, scalar rate);

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
                                          scalar rate, scalar eductDistance);

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
                                        scalar rate, scalar productDistance,
                                        scalar weight1 = 0.5, scalar weight2 = 0.5);

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
                                       scalar rate, scalar eductDistance,
                                       scalar weight1 = 0.5, scalar weight2 = 0.5);

    /**
     * Method to register a decay reaction.
     * @param name the name of the reaction
     * @param particleType the type for which this decay should be performed
     * @param rate the rate
     * @return a uuid with which this reaction can be removed again
     * @todo implement removal of reactions
     */
    const short registerDecayReaction(const std::string &name, const std::string &particleType,
                                      scalar rate);

    void registerSpatialTopologyReaction(const std::string &descriptor, scalar rate, scalar radius);

    void registerStructuralTopologyReaction(const std::string &topologyType,
                                            const model::top::reactions::StructuralTopologyReaction &reaction);

    const short
    registerCompartmentSphere(const std::unordered_map<std::string, std::string> &conversionsMap,
                              const std::string &name, const model::Vec3 &origin,
                              scalar radius, bool largerOrLess);

    const short registerCompartmentPlane(const std::unordered_map<std::string, std::string> &conversionsMap,
                                         const std::string &name,
                                         const model::Vec3 &normalCoefficients, scalar distanceFromPlane,
                                         bool largerOrLess);

    void
    configureTopologyBondPotential(const std::string &type1, const std::string &type2, scalar forceConstant,
                                   scalar length, api::BondType type = api::BondType::HARMONIC);

    void configureTopologyAnglePotential(const std::string &type1, const std::string &type2, const std::string &type3,
                                         scalar forceConstant, scalar equilibriumAngle,
                                         api::AngleType type = api::AngleType::HARMONIC);

    void configureTopologyTorsionPotential(const std::string &type1, const std::string &type2, const std::string &type3,
                                           const std::string &type4, scalar forceConstant, unsigned int multiplicity,
                                           scalar phi_0, api::TorsionType type = api::TorsionType::COS_DIHEDRAL);

    virtual void run(time_step_type steps, scalar timeStep);

    scalar getRecommendedTimeStep(unsigned int N) const;

    template<typename SchemeType=readdy::api::ReaDDyScheme>
    readdy::api::SchemeConfigurator<SchemeType> runScheme(bool useDefaults = true);

    bool singlePrecision() const;

    bool doublePrecision() const;

    /**
     * Access the root node of the performance measurement tree.
     * @return reference to root node of performance measurement
     */
    const util::PerformanceNode &performanceRoot();
private:
    struct Impl;
    std::unique_ptr<readdy::Simulation::Impl> pimpl;

    readdy::model::Kernel *const getSelectedKernel() const;

    void ensureKernelSelected() const;

};

class NoKernelSelectedException : public std::runtime_error {
public:
    explicit NoKernelSelectedException(const std::string &__arg);
};

NAMESPACE_END(readdy)

#include "bits/Simulation_misc.h"
