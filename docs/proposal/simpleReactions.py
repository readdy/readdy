import readdy
import numpy.random
sim = readdy.Simulation()
sim.kBT = 2.
sim.boxsize = 5. # box will span from [-2.5,-2.5,-2.5] to [2.5,2.5,2.5]
sim.periodic_boundary = True

# Three species at once
sim.register_species(
	name=["A", "B", "C"],
	radius=[1., 1., 2.**(1./3.)],
	difffusion_coeff=[1., 1., 2.**(-1./3.)] )

# A and B react within distance of 2 with a rate of
# 420 s^-1 to form a C particle. The backwards reaction
# occurs with rate 0.3 s^-1.
sim.register_simple_reaction(
	name="binary",
	descr="A+{2}B{42e+1}<-->{3e-1}C")

def reaction(particle_A, particle_B, reaction_context):
	# perform reaction
	return Particle("type", ...)

sim.register_reaction(name="..", f=reaction, ...)


# place 666 A particles and 666 B particles at random positions
for i in range(666):
	sim.add_particle("A", 5.*numpy.random.random(3)-2.5)
	sim.add_particle("B", 5.*numpy.random.random(3)-2.5)

sim.register_observable_numbers(
	species=["A", "B", "C"],
	lagtime=100,
	save="results/particlenumbers.h5")

sim.configure_CUDA(...)
sim.run(maxtime=1e6, timestep=0.2)