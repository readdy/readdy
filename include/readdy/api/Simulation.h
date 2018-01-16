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
 * simulation in a relatively simple manner. However, usage of the python top level API is recommended.
 *
 * @file Simulation.h
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
public:

    /**
     * Type of an observable callback corresponding to a certain observable
     */
    template<typename T>
    using observable_callback = typename std::function<void(typename T::result_type)>;

    /**
     * The default constructor. Instantiates the pimpl.
     */
    Simulation();

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
     * @see getKBT()
     */
    void setKBT(scalar kBT);

    /**
     * Method to get the current box size.
     * @return the box size as Vec3 object.
     */
    Vec3 getBoxSize() const;

    /**
     * Registers a new topology type which can be used to define structural and spatial topology reactions as well as
     * topology potential terms.
     * @param name the name of the type
     * @param reactions structual reactions of the type, can also be added with another call
     * @return internal id of the created type
     */
    topology_type_type registerTopologyType(const std::string& name,
                                            const std::vector<model::top::reactions::StructuralTopologyReaction> &reactions = {});

    /**
     * Creates a topology particle of a certain type at a position without adding it to the simulation box yet.
     * In order to instantiate it in the simulation it has to be used for creating a topology that owns this particular
     * particle instance.
     * @param type the type of topology flavored particle
     * @param pos the position
     * @return the created topology particle instance
     * @see addTopology()
     */
    readdy::model::TopologyParticle createTopologyParticle(const std::string &type, const Vec3 &pos) const;

    /**
     * Checks whether the selected kernel supports topologies
     * @return true if the kernel supports topologies
     */
    bool kernelSupportsTopologies() const;

    /**
     * Method for adding a new topology based on a topology type and topology particles. Node that initially the
     * particles are unconnected in the topology's connectivity graph, one has to make sure that there is a connection
     * through bonds between the particles upon simulation start.
     * @param type the type of topology to create
     * @param particles the topology's particles
     * @return a pointer to the instance of topology that was created
     * @see createTopologyParticle(type, pos)
     */
    readdy::model::top::GraphTopology *addTopology(const std::string& type,
                                                   const std::vector<readdy::model::TopologyParticle> &particles);

    /**
     * Method yielding a collection of pointers to the currently existing topologies in the simulation.
     * @return a vector of graph topology pointers
     */
    std::vector<readdy::model::top::GraphTopology *> currentTopologies();

    /**
     * Method yielding a collection of particles that corresponds to the particles currently contained in the given
     * topology instance.
     * @param topology the topology instance
     * @return a vector of particles
     */
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
    void setBoxSize(const Vec3 &boxSize);

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

    const model::observables::ObservableFactory &observe() const;

    /**
     * Registers a predefined observable with the kernel. A list of available observables can be obtained by
     * getAvailableObservables().
     * @tparam observable type
     * @param observable the observable
     * @return a uuid with which the observable is associated
     */
    template<typename T>
    ObservableHandle registerObservable(std::unique_ptr<T> observable, detail::is_observable_type<T>* = 0);

    /**
     * Registers a predefined observable with the kernel together with a callback.
     * A list of available observables can be obtained by getAvailableObservables().
     * @tparam T observable type
     * @param observable the observable instance
     * @param callback the callback
     * @return a observable handle that allows for post-hoc modification of the observable
     */
    template<typename T>
    ObservableHandle registerObservable(std::unique_ptr<T> observable, const observable_callback<T> &callback,
                                        detail::is_observable_type<T>* = 0);

    /**
     * Registers an observable that implements the readdy::model::ObservableBase interface.
     * @param observable the observable
     * @return a uuid with which the observable is associated
     */
    ObservableHandle registerObservable(readdy::model::observables::ObservableBase &observable);

    /**
     * Removes an observable by uuid.
     * @param uuid the uuid of the observable to be removed.
     */
    void deregisterObservable(ObservableHandle::observable_id uuid);

    /**
     * Removes an observable by its handle
     * @param handle the handle
     */
    void deregisterObservable(const ObservableHandle &handle);

    /**
     * A method to access the particle positions of a certain type.
     * @param type the type
     * @return a vector containing the particle positions
     */
    std::vector<Vec3> getParticlePositions(std::string type);

    /**
     * A method to register a particle type
     * @param name the name (must be unique)
     * @param diffusionCoefficient the diffusion coefficient
     * @param radius the particle's radius, important for some potentials (like, e.g., harmonic repulsion)
     * @param flavor the particle's flavor, one of NORMAL and TOPOLOGY; defaulting to NORMAL
     */
    particle_type_type
    registerParticleType(const std::string &name, scalar diffusionCoefficient,
                         readdy::model::particle_flavor flavor = readdy::model::particleflavor::NORMAL);

    /**
     * A method that allows to remove a certain potential type.
     * @param id the id of this potential
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
     * @return a uuid with which the potential can be removed
     */
    const short registerBoxPotential(const std::string &particleType, scalar forceConstant,
                                     const Vec3 &origin, const Vec3 &extent);

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
    registerSphereInPotential(const std::string &particleType, scalar forceConstant, const Vec3 &origin,
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
    registerSphereOutPotential(const std::string &particleType, scalar forceConstant, const Vec3 &origin,
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
    registerSphericalBarrier(const std::string &particleType, scalar height, scalar width, const Vec3 &origin,
                             scalar radius);

    //----------------------
    // Order 2 potentials
    //----------------------
    /**
     * Register a harmonic repulsion potential.
     * @param particleTypeA particle type A
     * @param particleTypeB particle type B
     * @param forceConstant the force constant
     * @param interactionDistance the interaction distance
     * @return a uuid with which the potential can be removed again
     * @todo document this more thoroughly
     */
    const short registerHarmonicRepulsionPotential(const std::string &particleTypeA, const std::string &particleTypeB,
                                                   scalar forceConstant, scalar interactionDistance);

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

    /**
     * Registers a pointer of order 1, i.e., a potential that acts as an external force on all particles of a
     * certain type.
     * @param ptr reference to the potential instance
     */
    void registerPotentialOrder1(readdy::model::potentials::PotentialOrder1 *ptr);

    /**
     * Register a potential of order 2, i.e., a potential that depends on the positions of a particle pair.
     * @param ptr a const pointer to a PotentialOrder2 instance that is to be evaluated
     */
    void registerPotentialOrder2(readdy::model::potentials::PotentialOrder2 *ptr);

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
    const std::vector<Vec3> getAllParticlePositions() const;

    /**
     * Method that allows to set a kernel for this simulation object. Most methods required that this method
     * is invoked beforehand.
     * @param kernel the kernel that is to be selected
     */
    void setKernel(const std::string &kernel);

    /**
     * Yields a modifyable reference to the current context of this simulation. Can be used to set time-independent
     * properties.
     * @return the context
     */
    model::Context &currentContext();

    /**
     * Yields a nonmodifyable reference to the current context of this simulation.
     * @return cref to the context
     */
    const model::Context &currentContext() const;

    /**
     * Method that allows to set an already existing instance of a kernel for this simulation object.
     * @param kernel the kernel
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
     */
    const short registerDecayReaction(const std::string &name, const std::string &particleType,
                                      scalar rate);

    /**
     * Method to register a spatial topology reaction. It is given by its descriptor, a rate, and a radius in which
     * it can happen.
     * @param descriptor the descriptor
     * @param rate the rate
     * @param radius the radius
     */
    void registerSpatialTopologyReaction(const std::string &descriptor, scalar rate, scalar radius);

    /**
     * Method to register a structural topology reaction. It is bound to a topology type.
     * @param topologyType the topology type
     * @param reaction the reaction, consisting out of a rate and a recipe-generating function
     */
    void registerStructuralTopologyReaction(const std::string &topologyType,
                                            const model::top::reactions::StructuralTopologyReaction &reaction);

    /**
     * Registers a sphere compartment.
     * @param conversionsMap the conversions
     * @param name the name
     * @param origin the origin
     * @param radius the radius
     * @param largerOrLess larger? less?
     * @return id
     * @todo more extensive documentation.
     */
    const short
    registerCompartmentSphere(const std::unordered_map<std::string, std::string> &conversionsMap,
                              const std::string &name, const Vec3 &origin,
                              scalar radius, bool largerOrLess);

    /**
     * Registers a plane compartment.
     * @param conversionsMap the conversions
     * @param name the name
     * @param normalCoefficients the diffusion coefficients
     * @param distanceFromPlane the distance from plane
     * @param largerOrLess larger? less?
     * @return id
     * @todo more extensive documentation
     */
    const short registerCompartmentPlane(const std::unordered_map<std::string, std::string> &conversionsMap,
                                         const std::string &name,
                                         const Vec3 &normalCoefficients, scalar distanceFromPlane,
                                         bool largerOrLess);

    /**
     * Configures a bond potential between a pair of particles with the given types. In order for the potential to come
     * into effect, the particles have to be connected in the connectivity graph of their topology. The types can be
     * reversed in order.
     * @param type1 first type of the tuple
     * @param type2 second type of the tuple
     * @param forceConstant the stiffness of the potential
     * @param length preferred distance between the particles with respect to this bond
     * @param type type of bond, defaults to harmonic bond
     */
    void configureTopologyBondPotential(const std::string &type1, const std::string &type2, scalar forceConstant,
                                        scalar length, api::BondType type = api::BondType::HARMONIC);

    /**
     * Configures an angle potential between a triple of particles with the given types. In order for the potential to
     * come into effect, the particles have to be connceted in the connectivity graph of their topology. The types
     * can be reversed in order.
     * @param type1 first type of the triple
     * @param type2 second type of the triple, type of the middle particle
     * @param type3 third type of the triple
     * @param forceConstant stiffness of the potential
     * @param equilibriumAngle the preferred angle between particles with respect to this potential
     * @param type the type of angle potential, defaults to harmonic angle
     */
    void configureTopologyAnglePotential(const std::string &type1, const std::string &type2, const std::string &type3,
                                         scalar forceConstant, scalar equilibriumAngle,
                                         api::AngleType type = api::AngleType::HARMONIC);

    /**
     * Configures a torsion potential between a quadruple of particles with the given types. In order for the potential
     * to come into effect, the particles have to be connected in the connectivity graph of their topology. The types
     * can be reversed in order.
     * @param type1 first type of the quadruple
     * @param type2 second type of the quadruple
     * @param type3 third type of the quadruple
     * @param type4 fourth type of the quadruple
     * @param forceConstant stiffness of the potential
     * @param multiplicity number of minima in the energy function
     * @param phi_0 reference torsion angle
     * @param type the type of torsion potential, defaults to cosine dihedral
     */
    void configureTopologyTorsionPotential(const std::string &type1, const std::string &type2, const std::string &type3,
                                           const std::string &type4, scalar forceConstant, unsigned int multiplicity,
                                           scalar phi_0, api::TorsionType type = api::TorsionType::COS_DIHEDRAL);

    /**
     * Runs the simulation using the default simulation loop for a given number of time steps and a given time step.
     * @param steps the number of steps
     * @param timeStep the time step
     * @see runScheme()
     */
    virtual void run(time_step_type steps, scalar timeStep);

    /**
     * This method yields a scheme configurator that can be used in fluent-api style to modify the simulation loop
     * with respect to the needs of the particular system.
     * @tparam SchemeType the type of simulation loop, can be one of ReaDDyScheme (default) and AdvancedScheme
     * @param useDefaults whether to initialize the configurator with the default values that are used in the simple run
     * @return the configurator object
     */
    template<typename SchemeType=readdy::api::ReaDDyScheme>
    readdy::api::SchemeConfigurator<SchemeType> runScheme(bool useDefaults = true);

    /**
     * Checks if the kernel is running on single precision.
     * @return true if it is running on single precision
     */
    bool singlePrecision() const;

    /**
     * Checks if the kernel is running on double precision.
     * @return true if it is running on double precision
     */
    bool doublePrecision() const;

    /**
     * Access the root node of the performance measurement tree.
     * @return reference to root node of performance measurement
     */
    const util::PerformanceNode &performanceRoot();

    /**
     * This method can be used to set kernel specific configuration. It takes a JSON string.
     * @param conf json
     */
    void setKernelConfiguration(const std::string &conf);
private:
    /**
     * the Impl struct
     */
    struct Impl;
    /**
     * reference to the impl
     */
    std::unique_ptr<readdy::Simulation::Impl> pimpl;

    /**
     * This method yields a pointer to the currently selected kernel, can be nullptr.
     * @return a pointer to the currently selected kernel
     */
    readdy::model::Kernel *const getSelectedKernel() const;

    /**
     * This method raises if currently no kernel is selected for this simulation instance
     */
    void ensureKernelSelected() const;

};

class NoKernelSelectedException : public std::runtime_error {
public:
    /**
     * Exception that gets raised if a method gets called that requires a kernel to be set on the simulation object.
     * @param __arg message
     */
    explicit NoKernelSelectedException(const std::string &__arg);
};

NAMESPACE_END(readdy)

#include "bits/Simulation_misc.h"
