---
title: Proposals
---

### General suggestions

- [ ] update API as below (and create top level python api):

```python
# define system here
system = ReactionDiffusionSystem(kernel="SingleCPU")

# create simulation (define algorithms to run the system)
simulation = system.simulation()  # use defaults. 
simulation.add_observable('rdf')  # actually part of the system, but configured through the simulation object
simulation.integrator = "EulerBDIntegrator"
simulation.compute_forces = False
simulation.reaction_scheduler = "UncontrolledApproximation"
simulation.evaluate_observables = False
simulation.run(10)  # does the configuration and runs the simulation

# continue simulation
simulation.run(10)

# second call RAISES exception, because you can only simulate a system once.
simulation = system.simulation()  # use defaults. 


# OR
simulation = system.simulation(integrator = "EulerBDIntegrator",
                               compute_forces = False,
...
)
simulation.run(10)
```

- [ ] suggestion: Allow geometry files as input for a box potential such that more complicated shapes can be realized with external tools
- [ ] implement CUDA kernel
    - meet up with Felix to discuss HALMD integration
- [ ] implement reactions with topologies
    - come up with convenient API to create / manipulate topologies
- [ ] improve reaction scheduler to gain more performance
    - filter particles out, that do not participate in any reaction
    - upon event creation, check if event is scheduled to happen in the current time interval
    - this introduces a bias on the probabilities of the remaining events (if there are intersections), try to balance that
- [ ] improve neighbor lists to gain more performance
    - verlet lists
- [ ] snapshotting
    - this point belongs together with the IO point
    - implement snapshotting using the observables framework
- [ ] implement IO (de-/serialization, dumping of trajectories into hdf5 files)
    - implement VMD plugin (see, e.g., lammps plugin on how to hide particles)
    - use and extend h5md?
    - use h5xx?
    - implement IO using the observables framework
- [ ] create benchmark (NCores x NParticles x rates)
    - maybe execute this benchmark automatically on some host
- [ ] domain decomposition (e.g., MPI)

### Topology reaction scheduling on GPUs

Let's assume we can do the following on the GPU
- Diffusion of normal particles and topologies
- Simple reactions, i.e. reactions between normal particles

Topology reactions can change the structure of topologies (e.g. polymerization, 
binding/unbinding to a complex). This cannot be done on the GPU. Instead those reactions
have to be performed on the CPU, which is in principle not a problem when those reactions
occur rarely. The actual problem is, that the __GPU cannot halt on its own__ when it find out that
a topology reaction should be performed. There are two ways of determining how long the GPU
should execute:
1. with a fixed time $\tau$
    - the GPU executes diffusion and normal reactions for a time $\tau$ which is much larger
    than the integration step and then returns
    - the CPU performs all possible topology reaction events based on its current state, 
    where reaction probabilities are $\mathrm{rate} * \tau$. This could be done with the fixed timestep
    version of our reaction schedulers
2. with a time $\tau$ sampled from a Gillespie algorithm
    - given one system state with a number of possible topology reactions events, 
    choose __one__ event and a corresponding $\tau
    - perform this reaction and then let the GPU run for $\tau$

### Aggregators and files

As of now observables describe only one point in time. Aggregators are observables that accumulate
data. Suggestion:
- Since aggregators accumulate the data, only they can clear the data. Thus the aggregator should also
be responsible for the file-writing process.
- Aggregators hold a file object, which is optional.
