---
layout: page
title: Development plan
permalink: /development-plan/
---

### Todos:
* update API as discussed (and create top level python api)
* implement CUDA kernel
* implement reactions with topologies
* improve reaction scheduler to gain more performance
* improve neighbor lists to gain more performance
* snapshotting
* implement IO (de-/serialization, dumping of trajectories into hdf5 files)
	* implement vmd plugin for that
* create benchmark (NCores x NParticles x rates)
* domain decomposition (e.g., MPI)