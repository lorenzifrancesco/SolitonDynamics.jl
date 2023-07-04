# CondensateDynamics.jl dev notes
 We focus on obtaining
 - [x] ground states
[[dev_notes/ground_states]]
    - still to understand the strange slightly inexact results of the ground states with respect to analytics
    - [x] compute some final versions of the plots (mimick paper with [0.15, 0.4, 0.65] strengths)
    - NB: the 0.15 case does not correspond. Only for $\gamma = 0.65$ we have correspondence.
 - [x] tiles for the dynamics
[[dev_notes/tiles]]
- [x] Collapse 
[[dev_notes/collapse]]
 - Collapse for NPSE is consistent
 - Collapse for 3D GPE is measurable with density.
 - [ ] sigma 2
[[dev_notes/sigma2]]
## June 2023: Major refactoring
We changed radically the workflow, using scripts as standalone units in which we import CondensateDynamics
_NB: what equation are we actually simulating? We should [ ] write it down_

## July 2023: Finalization
### 30/06->2/07
- [ ] Obtain the plots with ground states and with the zooms
- [x] Fix the tiles / launch a prototype computation on the cloud
- [x] Launch the tiles simulations on the cloud
> results: unusual behaviour of 3D GPE for $\gamma=0.65$. We should check more carefully the reult tile by tile.
- [ ] collect all the parameters in a table (especially the barrier ones).
##### Problem: prepare_for_collision prepares in ground state but from time to time normalization check fail

