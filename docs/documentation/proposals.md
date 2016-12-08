---
layout: page
title: Proposals
---

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
1. with a fixed time `\tau`
    - the GPU executes diffusion and normal reactions for a time `\tau` which is much smaller
    than the integration step and then returns
    - the CPU performs all possible topology reaction events based on its current state, 
    where reaction probabilities are `rate * \tau`. This could be done with the fixed timestep
    version of our reaction schedulers
2. with a time `\tau` sampled from a Gillespie algorithm
    - given one system state with a number of possible topology reactions events, 
    choose __one__ event and a corresponding `\tau`
    - perform this reaction and then let the GPU run for `\tau`

### Aggregators and files

As of now observables describe only one point in time. Aggregators are observables that accumulate
data. Suggestion:
- Since aggregators accumulate the data, only they can clear the data. Thus the aggregator should also
be responsible for the file-writing process.
- Aggregators hold a file object, which is optional.
