import readdy # convenience objects and methods

sim = readdy.Simulation()

# set the temperature. 
sim.kBT = 1. 
#sim.boxsize = 10. # cubic box 
sim.boxsize = [10., 10., 200.] # non-cubic box
#sim.periodicboundary = True
# If one chooses to use non-periodic boundaries,
# the user has to ensure that particles do not leave
# the box e.g. by using repulsive walls.
sim.periodic_boundary = [True, True, False]

# Add two repulsive walls at the ends of the cuboid,
# both pointing inwards.
sim.register_repulsive_wall(
	normalvector=[0.,0.,+1.],
	point=[0.,0.,-100.] )
sim.register_repulsive_wall(
	normalvector=[0.,0.,-1.],
	point=[0.,0.,+100.] )

# First create a new particle species, then add one particle
# of this type to the simulation box.
# Returns (speciesID, speciesName)
sim.register_species(
	name="vesicle",
	radius=2.,
	diffusion_coeff=0.1) 

# Returns uniqueParticleID
sim.add_particle(
	species="vesicle", # optionally give speciesID
	position=[0.,0.,42.] ) 

# This observable will be computed and written to file on-the-fly
sim.register_observable_msd(
	species="vesicle", # optionally give speciesID
	#particles=0, # give uniqueID or list of IDs
	lagtime=10,
	save="results/meansquareddisplacement.h5")

# by default the CPU kernel will be used with default arguments:
# 1 node with 4 threads, or similar
sim.run(maxtime=10000, timestep=0.1)